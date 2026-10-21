library(readxl)
library(dplyr)
library(lubridate)
library(stringr)
library(plyr)
library(imputeTS)

### 1) Precipitation (DAILY; remove Feb 29)
precip_raw <- read_excel("~/data/BP_sikorski_precip.xlsx") %>%
  mutate(dtime = as.POSIXct(dtime * 86400, origin = "1970-01-01", tz = "UTC")) %>%
  select(dtime, precip_in = bpPrcp)
precip_raw$precip_in <- as.numeric(precip_raw$precip_in)

precip_daily <- precip_raw %>%
  mutate(date = as.Date(dtime)) %>%
  filter(format(date, "%m-%d") != "02-29") %>%
  select(date, precip_in)

### 2) Water level (HOURLY; round to hour, 1 per hour, full grid, no Feb 29)
waterlevel_raw <- read_excel("~/data/BPsealevel.xlsx") %>%
  dplyr::rename(
    dtime    = `Days since Jan1, 1970`,
    wl_obs   = `BP_NAVDsealevel (m)`,
    wl_tide  = `BP_NAVD_TIde (m)`,
    wl_resid = `BP_NonTide (m)`
  ) %>%
  mutate(time = as.POSIXct(dtime * 86400, origin = "1970-01-01", tz = "UTC")) %>%
  select(time, wl_obs, wl_tide, wl_resid)

waterlevel_hourly_obs <- waterlevel_raw %>%
  mutate(time_raw = time, time = round_date(time_raw, "hour")) %>%
  arrange(time_raw) %>%
  group_by(time) %>% slice(1) %>% ungroup() %>%
  select(time, wl_obs, wl_tide, wl_resid) %>%
  filter(format(time, "%m-%d") != "02-29")

waterlevel_hourly_grid <- tibble(
  time = seq(min(waterlevel_hourly_obs$time), max(waterlevel_hourly_obs$time), by = "hour")
) %>% filter(format(time, "%m-%d") != "02-29")

waterlevel_hourly <- waterlevel_hourly_grid %>%
  left_join(waterlevel_hourly_obs, by = "time") %>%
  arrange(time)

### 3) Wind (HOURLY; pre-2000 from RDA + 2000+ from CSV)
env_wind <- new.env()
load("/Users/kapilchandak/Dropbox/Marine Sciences-2024/Code- Joint Probability of High Coastal Sea Level and High Precipitation/ImputationCode/Sikorsky-Bridgeport_merged1.Rda",
     envir = env_wind)
wind_pre2000_raw <- get(ls(env_wind), envir = env_wind)
rm(env_wind)

wind_pre2000_hourly <- wind_pre2000_raw %>%
  mutate(
    time_raw = as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:%02d",
                                  year, month, day, hour, minute, second), tz = "UTC"),
    time = round_date(time_raw, "hour")
  ) %>%
  arrange(time_raw) %>% group_by(time) %>% slice(1) %>% ungroup() %>%
  select(time, WSPD, WDIR) %>%
  filter(time <= as.POSIXct("1999-12-31 23:59:59", tz = "UTC")) %>%
  filter(format(time, "%m-%d") != "02-29")

wind_2000plus_raw <- read.csv("~/data/combined_data_filtered.csv", stringsAsFactors = FALSE)
colnames(wind_2000plus_raw) <- c("Timestamp", "WDIR", "WSPD")
wind_2000plus_prep <- wind_2000plus_raw
wind_2000plus_prep$Timestamp <- trimws(as.character(wind_2000plus_prep$Timestamp))
wind_2000plus_prep$Timestamp <- ifelse(
  grepl("^\\d{4}-\\d{2}-\\d{2}$", wind_2000plus_prep$Timestamp),
  paste0(wind_2000plus_prep$Timestamp, " 00:00:00"), wind_2000plus_prep$Timestamp)
wind_2000plus_prep$Timestamp <- trimws(wind_2000plus_prep$Timestamp)
wind_2000plus_prep$Timestamp <- as.POSIXct(wind_2000plus_prep$Timestamp, tz = "UTC",
                                           tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"))

wind_2000plus_hourly <- wind_2000plus_prep %>%
  mutate(time_raw = as.POSIXct(Timestamp, tz = "UTC"), time = round_date(time_raw, "hour")) %>%
  arrange(time_raw) %>% group_by(time) %>% slice(1) %>% ungroup() %>%
  select(time, WSPD, WDIR) %>%
  filter(format(time, "%m-%d") != "02-29")

wind_hourly_obs <- bind_rows(wind_pre2000_hourly, wind_2000plus_hourly) %>%
  arrange(time) %>% group_by(time) %>% slice(1) %>% ungroup()

wind_hourly_grid <- tibble(
  time = seq(min(wind_hourly_obs$time), max(wind_hourly_obs$time), by = "hour")
) %>% filter(format(time, "%m-%d") != "02-29")

wind_hourly <- wind_hourly_grid %>% left_join(wind_hourly_obs, by = "time") %>% arrange(time)

### Restrict all to 1973–2023
precip_daily <- precip_daily %>%
  filter(date >= as.Date("1973-01-01"), date <= as.Date("2023-12-31"))
waterlevel_hourly <- waterlevel_hourly %>%
  filter(time >= as.POSIXct("1973-01-01 00:00:00", tz = "UTC"),
         time <= as.POSIXct("2023-12-31 23:00:00", tz = "UTC"))
wind_hourly <- wind_hourly %>%
  filter(time >= as.POSIXct("1973-01-01 00:00:00", tz = "UTC"),
         time <= as.POSIXct("2023-12-31 23:00:00", tz = "UTC"))
precip_daily$precip_in <- as.numeric(precip_daily$precip_in)

wind_hourly <- wind_hourly %>%
  mutate(year = year(time), month = month(time), day = day(time),
         hour = hour(time), minute = minute(time), second = second(time))

Data <- data.frame()
set.seed(5)   # global seed, as in the paper

for (y1 in 1973:2023) {
  data_sub2 <- wind_hourly %>% filter(year == y1) %>% arrange(time)
  
  data_sub2$WDIR[data_sub2$WDIR > 360] <- NA
  data_sub2$WSPD[!is.na(data_sub2$WSPD) & is.na(data_sub2$WDIR)] <- NA
  data_sub2$WDIR[ is.na(data_sub2$WSPD) & !is.na(data_sub2$WDIR)] <- NA
  
  rle1  <- rle(is.na(data_sub2$WSPD))
  rle11 <- rep(ifelse(rle1$values, rle1$lengths, 0), times = rle1$lengths)
  
  ## [FLAG] record each hour's imputation type from the ORIGINAL gaps
  data_sub2$wind_imp <- ifelse(rle11 == 0, "obs",
                               ifelse(rle11 <= 5, "short", "long"))
  
  Z1 <- rep(0, nrow(data_sub2))
  Z1[which(rle11 <= 5 & rle11 > 0)] <- rnorm(length(which(rle11 <= 5 & rle11 > 0)), 0, 1)
  sdWSPD <- sd(na.omit(data_sub2$WSPD))
  
  # short gaps (<=5)
  cc <- 1
  data_sub2$WSPD1 <- data_sub2$WSPD
  data_sub2$WSPD  <- na_interpolation(data_sub2$WSPD1, maxgap = 5) + (cc * sdWSPD * Z1)
  data_sub2$WDIR1 <- data_sub2$WDIR
  data_sub2$WDIR  <- na_interpolation(data_sub2$WDIR1, maxgap = 5)
  
  # long gaps (>5)
  rle1  <- rle(is.na(data_sub2$WSPD))
  rle11 <- rep(ifelse(rle1$values, rle1$lengths, 0), times = rle1$lengths)
  ind1 <- which(rle11 > 5)
  sample1 <- sample(which(rle11 == 0), size = length(ind1), replace = TRUE)
  data_sub2$WSPD[ind1] <- data_sub2$WSPD[sample1]
  data_sub2$WSPD <- ifelse(data_sub2$WSPD <= 0.25, 0.25, data_sub2$WSPD)
  data_sub2$WDIR[ind1] <- data_sub2$WDIR[sample1]
  data_sub2$WDIR <- ifelse(data_sub2$WDIR <= 0, 0, data_sub2$WDIR)
  
  data_sub2$WDIR <- data_sub2$WDIR %% 360
  data_sub2$WDIR <- round_any(data_sub2$WDIR, 10)
  
  ## [FLAG] carry wind_imp through
  Data <- rbind(Data, data_sub2 %>% select(time, WSPD, WDIR, wind_imp))
}

wind_hourly_imputed <- Data %>% arrange(time)
wind_hourly_imputed$WDIR <- wind_hourly_imputed$WDIR %% 360

## [FLAG] mark which water-level residual hours were missing (before interpolation)
waterlevel_hourly$wl_resid_imp <- is.na(waterlevel_hourly$wl_resid)

waterlevel_hourly$wl_obs   <- na_interpolation(waterlevel_hourly$wl_obs)
waterlevel_hourly$wl_resid <- na_interpolation(waterlevel_hourly$wl_resid)
waterlevel_hourly <- waterlevel_hourly[, c("time","wl_obs","wl_resid","wl_resid_imp")]

## precip imputation (no missing -> no change)
precip_daily$precip_in <- na_interpolation(precip_daily$precip_in)

rm(list = setdiff(ls(), c("precip_daily", "waterlevel_hourly", "wind_hourly_imputed")))

hourly_data <- full_join(waterlevel_hourly, wind_hourly_imputed, by = "time") %>% arrange(time)

hourly_data <- hourly_data %>%
  mutate(wdir_group = case_when(
    is.na(WDIR)               ~ NA_real_,
    WDIR >= 60  & WDIR <= 180 ~ 1,
    WDIR >= 240 & WDIR <= 300 ~ 2,
    TRUE                      ~ 3),
    u_e = WSPD * cos(pi * (-13 + 630 - WDIR) / 180),
    u_n = WSPD * sin(pi * (-13 + 630 - WDIR) / 180))

get_mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

daily_data <- hourly_data %>%
  dplyr::mutate(date = as.Date(time)) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    wdir_mode = get_mode(wdir_group),
    d_t = max(wl_resid), neg_d_t = max(-wl_resid),
    S_t = max(WSPD), u_t = max(u_e), neg_u_t = max(-u_e), v_t = max(u_n),
    .groups = "drop") %>%
  dplyr::arrange(date)

daily_data <- dplyr::left_join(daily_data, precip_daily, by = "date")

## =====================================================================
## IMPUTATION IMPACT
## =====================================================================

## ---- (A) hourly-level imputation counts ----
hd <- hourly_data %>%
  mutate(date = as.Date(time)) %>%
  filter(date >= as.Date("1973-01-01"), date <= as.Date("2023-12-31"))
nh <- nrow(hd)
cat(sprintf("\nHourly records: %d\n", nh))
cat(sprintf("D_t (water level): %d hours interpolated (%.3f%%)\n",
            sum(hd$wl_resid_imp), 100 * mean(hd$wl_resid_imp)))
tw <- table(factor(hd$wind_imp, levels = c("obs","short","long")))
cat(sprintf("Wind: short-gap %d (%.3f%%), long-gap %d (%.3f%%), total %d (%.3f%%)\n",
            tw["short"], 100*tw["short"]/nh,
            tw["long"],  100*tw["long"]/nh,
            tw["short"]+tw["long"], 100*(tw["short"]+tw["long"])/nh))
## ---- (B) daily-maxima impact: did the max fall on an imputed hour? ----
imp <- hd %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    Dt_imp     = wl_resid_imp[which.max(wl_resid)],
    ut_type    = wind_imp[which.max(u_e)],
    negut_type = wind_imp[which.max(-u_e)],
    vt_type    = wind_imp[which.max(u_n)],
    .groups = "drop")
n <- nrow(imp)
cat(sprintf("\nDaily maxima: n = %d\n", n))
cat(sprintf("Precip : 0.000%% (no missing observations)\n"))
cat(sprintf("D_t    : %.3f%% of daily maxima fall on an interpolated hour\n", 100*mean(imp$Dt_imp)))
wind_pct <- function(col, lab) {
  tb <- table(factor(imp[[col]], levels = c("obs","short","long")))
  cat(sprintf("%-5s  : short = %.3f%%, long = %.3f%%, any = %.3f%%\n",
              lab, 100*tb["short"]/n, 100*tb["long"]/n, 100*(tb["short"]+tb["long"])/n))
}
wind_pct("ut_type",    "u_t")
wind_pct("negut_type", "-u_t")
wind_pct("vt_type",    "v_t")
