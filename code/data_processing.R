library(readxl)
library(dplyr)
library(lubridate)
library(stringr)
library(plyr)
library(imputeTS)

### 1) Precipitation (DAILY; remove Feb 29)
precip_raw <- read_excel("~file_path/BP_sikorski_precip.xlsx") %>%
  mutate(dtime = as.POSIXct(dtime * 86400, origin = "1970-01-01", tz = "UTC")) %>%
  select(dtime, precip_in = bpPrcp)

precip_raw$precip_in<-as.numeric(precip_raw$precip_in)

precip_daily <- precip_raw %>%
  mutate(date = as.Date(dtime)) %>%
  filter(format(date, "%m-%d") != "02-29") %>%
  select(date, precip_in)


### 2) Water level (HOURLY)
waterlevel_raw <- read_excel("~file_path/BPsealevel.xlsx") %>%
  dplyr::rename(
    dtime    = `Days since Jan1, 1970`,
    wl_obs   = `BP_NAVDsealevel (m)`,
    wl_tide  = `BP_NAVD_TIde (m)`,
    wl_resid = `BP_NonTide (m)`
  ) %>%
  mutate(
    time = as.POSIXct(dtime * 86400,
                      origin = "1970-01-01",
                      tz = "UTC")
  ) %>%
  select(time, wl_obs, wl_tide, wl_resid)


waterlevel_hourly_obs <- waterlevel_raw %>%
  mutate(time_raw = time,
         time = round_date(time_raw, "hour")) %>%
  arrange(time_raw) %>%
  group_by(time) %>%
  slice(1) %>%
  ungroup() %>%
  select(time, wl_obs, wl_tide, wl_resid) %>%
  filter(format(time, "%m-%d") != "02-29")

waterlevel_hourly_grid <- tibble(
  time = seq(min(waterlevel_hourly_obs$time), max(waterlevel_hourly_obs$time), by = "hour")
) %>%
  filter(format(time, "%m-%d") != "02-29")

waterlevel_hourly <- waterlevel_hourly_grid %>%
  left_join(waterlevel_hourly_obs, by = "time") %>%
  arrange(time)



### 3) Wind (HOURLY; pre-2000 from RDA + 2000+ from CSV)

# --- Pre-2000 wind from .Rda
env_wind <- new.env()
load("~file_path/Sikorsky-Bridgeport_merged1.Rda",
     envir = env_wind)
wind_pre2000_raw <- get(ls(env_wind), envir = env_wind)
rm(env_wind)

wind_pre2000_hourly <- wind_pre2000_raw %>%
  mutate(
    time_raw = as.POSIXct(
      sprintf("%04d-%02d-%02d %02d:%02d:%02d",
              year, month, day, hour, minute, second),
      tz = "UTC"
    ),
    time = round_date(time_raw, "hour")
  ) %>%
  arrange(time_raw) %>%
  group_by(time) %>%
  slice(1) %>%
  ungroup() %>%
  select(time, WSPD, WDIR) %>%
  filter(time <= as.POSIXct("1999-12-31 23:59:59", tz = "UTC")) %>%
  filter(format(time, "%m-%d") != "02-29")

# --- 2000+ wind from CSV
wind_2000plus_raw <- read.csv("~file_path/combined_data_filtered.csv",
                              stringsAsFactors = FALSE)
colnames(wind_2000plus_raw) <- c("Timestamp", "WDIR", "WSPD")

wind_2000plus_prep <- wind_2000plus_raw
wind_2000plus_prep$Timestamp <- trimws(as.character(wind_2000plus_prep$Timestamp))

wind_2000plus_prep$Timestamp <- ifelse(
  grepl("^\\d{4}-\\d{2}-\\d{2}$", wind_2000plus_prep$Timestamp),
  paste0(wind_2000plus_prep$Timestamp, " 00:00:00"),
  wind_2000plus_prep$Timestamp
)

wind_2000plus_prep$Timestamp <- trimws(wind_2000plus_prep$Timestamp)

wind_2000plus_prep$Timestamp <- as.POSIXct(
  wind_2000plus_prep$Timestamp,
  tz = "UTC",
  tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M")
)

wind_2000plus_hourly <- wind_2000plus_prep %>%
  mutate(
    time_raw = as.POSIXct(Timestamp, tz = "UTC"),
    time = round_date(time_raw, "hour")
  ) %>%
  arrange(time_raw) %>%
  group_by(time) %>%
  slice(1) %>%                 # first observation per rounded hour
  ungroup() %>%
  select(time, WSPD, WDIR) %>%
  filter(format(time, "%m-%d") != "02-29")

# --- Combine wind files
wind_hourly_obs <- bind_rows(wind_pre2000_hourly, wind_2000plus_hourly) %>%
  arrange(time) %>%
  group_by(time) %>%
  slice(1) %>%
  ungroup()

wind_hourly_grid <- tibble(
  time = seq(min(wind_hourly_obs$time), max(wind_hourly_obs$time), by = "hour")
) %>%
  filter(format(time, "%m-%d") != "02-29")

wind_hourly <- wind_hourly_grid %>%
  left_join(wind_hourly_obs, by = "time") %>%
  arrange(time)

#############################
# Restrict all final datasets to 1973–2023
#############################

# Daily precipitation
precip_daily <- precip_daily %>%
  filter(date >= as.Date("1973-01-01"),
         date <= as.Date("2023-12-31"))

# Hourly water level
waterlevel_hourly <- waterlevel_hourly %>%
  filter(time >= as.POSIXct("1973-01-01 00:00:00", tz = "UTC"),
         time <= as.POSIXct("2023-12-31 23:00:00", tz = "UTC"))

# Hourly wind
wind_hourly <- wind_hourly %>%
  filter(time >= as.POSIXct("1973-01-01 00:00:00", tz = "UTC"),
         time <= as.POSIXct("2023-12-31 23:00:00", tz = "UTC"))

precip_daily$precip_in<-as.numeric(precip_daily$precip_in)

# Add date parts 
wind_hourly <- wind_hourly %>%
  mutate(
    year   = year(time),
    month  = month(time),
    day    = day(time),
    hour   = hour(time),
    minute = minute(time),
    second = second(time)
  )

Data <- data.frame()

for (y1 in 1973:2023) {
  
  data_sub2 <- wind_hourly %>%
    filter(year == y1) %>%
    arrange(time)
  
  # -------------------------------------------------
  #  Clean WDIR before imputation
  # -------------------------------------------------
  data_sub2$WDIR[data_sub2$WDIR > 360] <- NA
  
  # Force paired missingness
  data_sub2$WSPD[!is.na(data_sub2$WSPD) & is.na(data_sub2$WDIR)] <- NA
  data_sub2$WDIR[ is.na(data_sub2$WSPD) & !is.na(data_sub2$WDIR)] <- NA
  
  # Identify NA runs
  rle1  <- rle(is.na(data_sub2$WSPD))
  rle11 <- rep(ifelse(rle1$values, rle1$lengths, 0), times = rle1$lengths)
  
  # Noise only for short gaps (<=5)
  Z1 <- rep(0, nrow(data_sub2))
  Z1[which(rle11 <= 5 & rle11 > 0)] <-
    rnorm(length(which(rle11 <= 5 & rle11 > 0)), 0, 1)
  
  sdWSPD <- sd(na.omit(data_sub2$WSPD))
  
  # -------------------------------------------------
  # Case 1: short gaps (<=5)
  # -------------------------------------------------
  c <- 1
  data_sub2$WSPD1 <- data_sub2$WSPD
  data_sub2$WSPD  <- na_interpolation(data_sub2$WSPD1, maxgap = 5) +
    (c * sdWSPD * Z1)
  
  data_sub2$WDIR1 <- data_sub2$WDIR
  data_sub2$WDIR  <- na_interpolation(data_sub2$WDIR1, maxgap = 5)
  
  # -------------------------------------------------
  # Case 2: long gaps (>5)
  # -------------------------------------------------
  rle1  <- rle(is.na(data_sub2$WSPD))
  rle11 <- rep(ifelse(rle1$values, rle1$lengths, 0), times = rle1$lengths)
  
  ind1 <- which(rle11 > 5)
  set.seed(5)
  sample1 <- sample(which(rle11 == 0), size = length(ind1), replace = TRUE)
  
  data_sub2$WSPD[ind1] <- data_sub2$WSPD[sample1]
  data_sub2$WSPD <- ifelse(data_sub2$WSPD <= 0.25, 0.25, data_sub2$WSPD)
  
  data_sub2$WDIR[ind1] <- data_sub2$WDIR[sample1]
  data_sub2$WDIR <- ifelse(data_sub2$WDIR <= 0, 0, data_sub2$WDIR)
  
  # -------------------------------------------------
  # Final WDIR handling
  # -------------------------------------------------
  data_sub2$WDIR <- data_sub2$WDIR %% 360      # force into [0,360)
  data_sub2$WDIR <- round_any(data_sub2$WDIR, 10)
  
  Data <- rbind(Data, data_sub2 %>% select(time, WSPD, WDIR))
}

wind_hourly_imputed <- Data %>%
  arrange(time)

wind_hourly_imputed$WDIR<-wind_hourly_imputed$WDIR%%360

waterlevel_hourly$wl_obs<-na_interpolation(waterlevel_hourly$wl_obs)
waterlevel_hourly$wl_resid<-na_interpolation(waterlevel_hourly$wl_resid)
waterlevel_hourly<-waterlevel_hourly[,c(1,2,4)]

## imputation of precipitation data
precip_daily$precip_in<-na_interpolation(precip_daily$precip_in)

### Keep only final objects
rm(list = setdiff(ls(), c("precip_daily", "waterlevel_hourly","wind_hourly_imputed")))

hourly_data <- full_join(
  waterlevel_hourly,
  wind_hourly_imputed,
  by = "time"
) %>%
  arrange(time)

#hourly_data <- hourly_data %>%
# mutate(
#    wdir_group = case_when(
#      is.na(WDIR)               ~ NA_real_,
#      WDIR >= 60  & WDIR <= 180 ~ 1,
#      WDIR >= 240 & WDIR <= 300 ~ 2,
#      TRUE                      ~ 3
#    )
#  )

hourly_data <- hourly_data %>%
  mutate(
    u_e = WSPD * cos(pi * (-13 + 630 - WDIR) / 180),
    u_n = WSPD * sin(pi * (-13 + 630 - WDIR) / 180)
  )

get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

daily_data <- hourly_data %>%
  dplyr::mutate(date = as.Date(time)) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    wdir_mode    = get_mode(wdir_group),
    #wl_obs_max   = max(wl_obs),
    d_t = max(wl_resid),
    neg_d_t = max(-wl_resid),
    S_t     = max(WSPD),
    u_t     = max(u_e),
    neg_u_t = max(-u_e),
    v_t     = max(u_n),
    .groups = "drop"
  ) %>%
  dplyr::arrange(date)

daily_data <- dplyr::left_join(
  daily_data,
  precip_daily,
  by = "date"
)

rm(list = setdiff(ls(), c("hourly_data", "daily_data")))

library(ggplot2)
library(patchwork)
library(lubridate)
library(colorspace)

# --- Prep
daily_data$date <- as.Date(daily_data$date)

# Keep ONLY 1973–2023 (safety; also prevents weird axis range)
daily_data <- daily_data %>%
  dplyr::filter(date >= as.Date("1973-01-01"),
                date <= as.Date("2023-12-31"))

# numeric time for smooth color mapping
daily_data$t_num <- as.numeric(daily_data$date)

# --- Print-friendly theme (big + readable)
th_print <- theme_bw(base_size = 18) +
  theme(
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# --- Rainbow color scale (includes yellow/orange; smooth like your example)
rainbow_like <- scale_color_gradientn(
  colors = colorspace::rainbow_hcl(12, c = 90, l = 55),
  guide  = "none"
)

# --- Fixed x-axis: exactly 1973–2023 (NO 2025)
xscale_fixed <- scale_x_date(
  limits      = as.Date(c("1973-01-01", "2023-12-31")),
  date_breaks = "5 years",
  date_labels = "%Y",
  expand      = c(0, 0)
)

# --- Panels
p_dt <- ggplot(daily_data, aes(x = date, y = d_t, color = t_num)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "Daily maximum sea-level residual",
    x = "Year",
    y = expression(bold(D[t])~"(m)")
  ) +
  xscale_fixed + rainbow_like + th_print

p_pr <- ggplot(daily_data, aes(x = date, y = precip_in, color = t_num)) +
  geom_line(linewidth = 0.6, alpha = 0.65) +  # slight alpha to avoid looking bolder
  labs(
    title = "Daily total precipitation",
    x = "Year",
    y = "Precipitation (in)"
  ) +
  xscale_fixed + rainbow_like + th_print

p_ut <- ggplot(daily_data, aes(x = date, y = u_t, color = t_num)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = expression("Daily maximum eastward wind component (" * u[t] * ")"),
    x = "Year",
    y = expression(bold(u[t])~"(m s"^{-1}*")")
  ) +
  xscale_fixed + rainbow_like + th_print

p_vt <- ggplot(daily_data, aes(x = date, y = v_t, color = t_num)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = expression("Daily maximum northward wind component (" * v[t] * ")"),
    x = "Year",
    y = expression(bold(v[t])~"(m s"^{-1}*")")
  ) +
  xscale_fixed + rainbow_like + th_print

# --- 2x2 grid
fig_2x2 <- (p_dt | p_pr) / (p_ut | p_vt) +
  plot_annotation(
    title = "Daily time series of sea level, precipitation, and wind components",
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5))
  )

# Display
fig_2x2

# Save (print-ready)
ggsave("daily_timeseries_2x2.png", fig_2x2, width = 14, height = 9, dpi = 300)
ggsave("daily_timeseries_2x2.pdf", fig_2x2, width = 14, height = 9)

library(ggplot2)

# main directory
out_dir <- "file_path"

# create subdirectory for scatterplots
scatter_dir <- file.path(out_dir, "Figure2_scatterplots")
dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)

# common y-axis range
ylims <- range(daily_data$d_t, na.rm = TRUE)

big_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  )

# (a) D_t vs P_t
fig_precip <- ggplot(daily_data, aes(x = precip_in, y = d_t)) +
  geom_point(shape = 19) +
  labs(
    title = "(a) Sea Level Residual vs Precipitation",
    x = "P_t (inches)",
    y = "D_t (m)"
  ) +
  coord_cartesian(ylim = ylims) +
  big_theme

ggsave(file.path(scatter_dir, "scatter_precip.png"),
       fig_precip, width = 7, height = 5, dpi = 300)

ggsave(file.path(scatter_dir, "scatter_precip.pdf"),
       fig_precip, width = 7, height = 5)


# (b) D_t vs -u_t
fig_westwind <- ggplot(daily_data, aes(x = neg_u_t, y = d_t)) +
  geom_point(shape = 19) +
  labs(
    title = "(b) Sea Level Residual vs Westward Wind",
    x = " -u_t (m s^-1)",
    y = "D_t (m)"
  ) +
  coord_cartesian(ylim = ylims) +
  big_theme

ggsave(file.path(scatter_dir, "scatter_westwind.png"),
       fig_westwind, width = 7, height = 5, dpi = 300)

ggsave(file.path(scatter_dir, "scatter_westwind.pdf"),
       fig_westwind, width = 7, height = 5)


# (c) D_t vs v_t
fig_northwind <- ggplot(daily_data, aes(x = v_t, y = d_t)) +
  geom_point(shape = 19) +
  labs(
    title = "(c) Sea Level Residual vs Northward Wind",
    x = "v_t (m s^-1)",
    y = "D_t (m)"
  ) +
  coord_cartesian(ylim = ylims) +
  big_theme

ggsave(file.path(scatter_dir, "scatter_northwind.png"),
       fig_northwind, width = 7, height = 5, dpi = 300)

ggsave(file.path(scatter_dir, "scatter_northwind.pdf"),
       fig_northwind, width = 7, height = 5)

fig 3

library(evd)

# main directory
out_dir <- "file_path"

# create subdirectory for return level curves
rl_dir <- file.path(out_dir, "Figure3_return_level_curves")
dir.create(rl_dir, showWarnings = FALSE, recursive = TRUE)


#-----------------------------
# Build bivariate datasets
#-----------------------------
Data1 <- data.frame(Pt = daily_data$precip_in, Dt = daily_data$d_t)
Data2 <- data.frame(neg_ut = daily_data$neg_u_t, Dt = daily_data$d_t)
Data3 <- data.frame(vt = daily_data$v_t, Dt = daily_data$d_t)


#-----------------------------
# (a) P_t vs D_t
#-----------------------------
k1 <- bvtcplot(Data1, pch=19, cex=0.1, cex.lab=2)$k0
thresh1 <- apply(Data1,2,sort,decreasing=TRUE)[(k1+5)/2,]

m1 <- fbvpot(Data1, thresh1, model="log")

png(file.path(rl_dir,"return_level_precip.png"), width=2000, height=1500, res=300)

plot(m1, which=3, p=c(0.9,0.95,0.99),
     col="orange", pch=19, cex=0.3,
     main="(a) Return Level Curves: P_t and D_t",
     xlab="P_t (inches)",
     ylab="D_t (m)",
     cex.lab=1.6, cex.axis=1.4, cex.main=1.4)

dev.off()

pdf(file.path(rl_dir,"return_level_precip.pdf"), width=7, height=5)

plot(m1, which=3, p=c(0.9,0.95,0.99),
     col="orange", pch=19, cex=0.3,
     main="(a) Return Level Curves: P_t and D_t",
     xlab="P_t (inches)",
     ylab="D_t (m)",
     cex.lab=1.6, cex.axis=1.4, cex.main=1.4)

dev.off()


#-----------------------------
# (b) -u_t vs D_t
#-----------------------------
k2 <- bvtcplot(Data2, pch=19, cex=0.1, cex.lab=2)$k0
thresh2 <- apply(Data2,2,sort,decreasing=TRUE)[(k2+5)/2,]

m2 <- fbvpot(Data2, thresh2, model="log")

png(file.path(rl_dir,"return_level_westwind.png"), width=2000, height=1500, res=300)

plot(m2, which=3, p=c(0.9,0.95,0.99),
     col="orange", pch=19, cex=0.3,
     main="(b) Return Level Curves: -u_t and D_t",
     xlab="-u_t (m s^-1)",
     ylab="D_t (m)",
     cex.lab=1.6, cex.axis=1.4, cex.main=1.4)

dev.off()

pdf(file.path(rl_dir,"return_level_westwind.pdf"), width=7, height=5)

plot(m2, which=3, p=c(0.9,0.95,0.99),
     col="orange", pch=19, cex=0.3,
     main="(b) Return Level Curves: -u_t and D_t",
     xlab="-u_t (m s^-1)",
     ylab="D_t (m)",
     cex.lab=1.6, cex.axis=1.4, cex.main=1.4)

dev.off()


#-----------------------------
# (c) v_t vs D_t
#-----------------------------
k3 <- bvtcplot(Data3, pch=19, cex=0.1, cex.lab=2)$k0
thresh3 <- apply(Data3,2,sort,decreasing=TRUE)[(k3+5)/2,]

m3 <- fbvpot(Data3, thresh3, model="log")

png(file.path(rl_dir,"return_level_northwind.png"), width=2000, height=1500, res=300)

plot(m3, which=3, p=c(0.9,0.95,0.99),
     col="orange", pch=19, cex=0.3,
     main="(c) Return Level Curves: v_t and D_t",
     xlab="v_t (m s^-1)",
     ylab="D_t (m)",
     cex.lab=1.6, cex.axis=1.4, cex.main=1.4)

dev.off()

pdf(file.path(rl_dir,"return_level_northwind.pdf"), width=7, height=5)

plot(m3, which=3, p=c(0.9,0.95,0.99),
     col="orange", pch=19, cex=0.3,
     main="(c) Return Level Curves: v_t and D_t",
     xlab="v_t (m s^-1)",
     ylab="D_t (m)",
     cex.lab=1.6, cex.axis=1.4, cex.main=1.4)

dev.off()

# Hourly data
write.csv(
  hourly_data,
  file = file.path(out_dir, "hourly_data.csv"),
  row.names = FALSE
)

# Daily data
write.csv(
  daily_data,
  file = file.path(out_dir, "daily_data.csv"),
  row.names = FALSE
)

