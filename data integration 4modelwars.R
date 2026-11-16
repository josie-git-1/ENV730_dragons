### Nov 15, 2026

## Data Integration: Building a Methane Flux Modeling Table

## from data sources to a model-ready table
  # target data product is a site-month modeling table

#set working directory

setwd("C:/Users/josie/Desktop/YSE/F26/ENV730/EDS730")

# Load libraries ----------------------------------------------------------
                    
library(sf)
library(terra)
library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(zoo)
library(lubridate)


  #writing help check functions for joins

count_missing <- function(data) {
  tibble(variable = names(data),
         missing_count = vapply(data, function(x) sum(is.na(x)), integer(1))) %>%
    arrange(desc(missing_count))
}

find_duplicate_keys <- function(data, ...) {
  data %>%
    count(..., name = "n") %>%
    filter(n > 1)
}


# Static site table -------------------------------------------------------

fluxnet_sites <- read.csv("data/FluxNet_Sites_2024.csv")
    
            #debug character display error
              #fluxnet_sites$SITE_NAME <- 
              #  iconv(fluxnet_sites$SITE_NAME, 
              #        from = "", 
              #        to = "UTF-8", 
              #        sub = "")

  #keep sites w CH4 data and valid coordinates      

fluxnet_ch4 <- fluxnet_sites %>%
  filter(FLUXNET.CH4 != "",
         !is.na(LOCATION_LONG),
         !is.na(LOCATION_LAT))

nrow(fluxnet_sites)

nrow(fluxnet_ch4)

  #convert methane site table to sf

fluxnet_ch4_sf <- st_as_sf(x = fluxnet_ch4,
                           coords = c("LOCATION_LONG", "LOCATION_LAT"),
                           crs = 4326,
                           remove = FALSE)

# Check the new attributes:

fluxnet_ch4_sf %>%
  st_drop_geometry() %>%
  select(SITE_ID, SITE_NAME) %>%
  head()

fluxnet_ch4_sf %>%
  st_drop_geometry() %>%
  select(SITE_ID, SITE_NAME) %>% summary


# Where are wetlands inundated? WAD2M -------------------------------------

wad2m_data <- rast("data/final_project/WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc")
range(names(wad2m_data))
#the layers are named as YYYY-MM

#building a monthly climatology, ignoring missing values in winter where snow and ice prevent inundation mapping
wad2m_month <- as.integer(substr(names(wad2m_all), 6, 7))


wad2m_climatology <- tapp(wad2m_all, index = wad2m_month, fun = mean, na.rm = TRUE)
names(wad2m_climatology) <- month.abb

#fill winter gaps with climatology ##and put it on the prediction grid

##to_ct_grid <- function(r) {
##  r %>%
##    resample(ct_grid, method = "bilinear") %>%
##    clamp(lower = 0, upper = 1, values = TRUE) %>%
##    mask(ct_grid)
##}

wad2m_2018 <- wad2m_all[[grepl("^2018", names(wad2m_all))]]

wad2m <- cover(wad2m_2018, wad2m_climatology) 
names(wad2m) <- month.abb

plot(wad2m[[c(1, 4, 7, 10)]])

#
tibble(month = month.abb,
       fraction_2018 = global(wad2m, "mean", na.rm = TRUE)[, 1],
       fraction_climatology = global(wad2m_climatology, "mean", na.rm = TRUE)[, 1])


#What to look for: WAD2M
#- Where: Higher fractions should follow river valleys (for example, the Connecticut River), coastal marshes, and low-lying terrain. Most upland cells are near zero.
#- When: Inundated area is largest in winter and spring, when evapotranspiration is low and snowmelt arrives, and smallest in late summer. That is the opposite of the flux-density cycle, so the two partly cancel in the budget.
#- 2018 vs. average: Compare fraction_2018 with fraction_climatology to see whether 2018 was a wet or dry year for wetlands.
#- Resolution: WAD2M is 0.25 degrees (about 25 km). Resampling to 5 km makes smoother maps but adds no new information about where individual wetlands are.



# Extract static raster predictors ----------------------------------------

  #extract soil data to each tower point

soil <- terra::rast("data/GlobalSoil_grids.tif")
names(soil)

    fluxnet_ch4_soil <- st_transform(fluxnet_ch4_sf, crs = crs(soil))

    soil_values <- terra::extract(soil, vect(fluxnet_ch4_soil)) %>%
    select(-ID) %>%
    rename(soil_bulk_density = BulkDensity,
           soil_ph = PH,
           soil_nitrogen = Nitrogen)
    

# Build the static site table ---------------------------------------------

    site_static <- fluxnet_ch4_sf %>%
      st_drop_geometry() %>%
      bind_cols(soil_values)
  
    site_static <- site_static %>%
      distinct(SITE_ID, .keep_all = TRUE)
    
    find_duplicate_keys(site_static, SITE_ID)
    
    site_static %>%
      select(SITE_ID, IGBP, soil_bulk_density, soil_ph, soil_nitrogen) %>%
      count_missing() %>%
      head(15)
    
    
### Dynamic monthly vegetation predictors
    modis_vi <- read.csv("data/ENV720-MOD13A3-061-results.csv")
    names(modis_vi)
    
    fluxnet_vi_monthly <- modis_vi %>%
      select(ID,
             Date,
             MOD13A3_061__1_km_monthly_EVI,
             MOD13A3_061__1_km_monthly_NDVI,
             MOD13A3_061__1_km_monthly_VI_Quality) %>%
      rename(SITE_ID = ID,
             EVI = MOD13A3_061__1_km_monthly_EVI,
             NDVI = MOD13A3_061__1_km_monthly_NDVI,
             QAQC = MOD13A3_061__1_km_monthly_VI_Quality) %>%
      mutate(date_formatted = as.Date(Date),
             YearMon = zoo::as.yearmon(date_formatted)) %>%
      filter(QAQC > 0 & QAQC != 65535 ) %>%
      distinct(SITE_ID, YearMon, .keep_all = TRUE)
    
    find_duplicate_keys(fluxnet_vi_monthly, SITE_ID, YearMon)
    
    range(fluxnet_vi_monthly$YearMon)
    
    
    
### Dynamic monthly inundation data
    
    inundation <- terra::rast("data/wad2m_data.tif")
    names(inundation)
    
    fluxnet_ch4_inundation <- st_transform(fluxnet_ch4_sf, crs = crs(inundation))
    
    inundation_values <- terra::extract(inundation, vect(fluxnet_ch4_inundation)) %>%
      select(-ID) %>%
      rename(soil_bulk_density = BulkDensity,
             soil_ph = PH,
             soil_nitrogen = Nitrogen)
    
### Monthly flux data
    
    
    load("data/FLUXNET_FLUXES.RDATA")
    
    fluxnet_flux <- FLUXNET.flux
    
    fluxnet_flux <- fluxnet_flux %>%
      mutate(SITE_ID = SITE)
    
    range(fluxnet_flux$YearMon)
    
    names(fluxnet_flux)
    
    find_duplicate_keys(fluxnet_flux, SITE_ID, YearMon)
    
### JOIN DYNAMIC PREDICTORS
    
    before_rows <- nrow(fluxnet_flux)
    
    flux_month_vi <- fluxnet_flux %>%
      left_join(fluxnet_vi_monthly, by = c("SITE_ID", "YearMon"))
    
    after_rows <- nrow(flux_month_vi)
    
    before_rows
    
    after_rows
    
    flux_month_vi %>%
      summarize(missing_ndvi = sum(is.na(NDVI)),
                missing_evi = sum(is.na(EVI)),
                total_rows = n())
    
    fluxnet_flux %>%
      anti_join(fluxnet_vi_monthly, by = c("SITE_ID", "YearMon")) %>%
      count(SITE_ID, sort = TRUE) %>%
      head()
    
    site_static_for_join <- site_static %>%
      select(SITE_ID, SITE_NAME, IGBP, LOCATION_LAT, LOCATION_LONG,
             soil_bulk_density, soil_ph, soil_nitrogen)
    
    before_rows <- nrow(flux_month_vi)
    
    flux_model_base <- flux_month_vi %>%
      left_join(site_static_for_join, by = "SITE_ID")
    
    after_rows <- nrow(flux_model_base)
    
    before_rows
    
    after_rows
    
    flux_model_base %>%
      filter(is.na(IGBP)) %>%
      count(SITE_ID, sort = TRUE) %>%
      head()
    
# Extract ERA5 and inundation values to tower sites
    
    ### ERA5
    
    era5_ppt_monthly <- rast("data/products/era5_ppt_monthly.tif")
    era5_tmean_monthly <- rast("data/products/era5_tmean_monthly.tif")
    
    nlyr(era5_ppt_monthly)
    
    nlyr(era5_tmean_monthly)
    
    names(era5_ppt_monthly) <- as.character(zoo::as.yearmon(names(era5_ppt_monthly)))
    names(era5_tmean_monthly) <- as.character(zoo::as.yearmon(names(era5_tmean_monthly)))
    
    fluxnet_ch4_era5 <- st_transform(fluxnet_ch4_with_context, crs = crs(era5_ppt_monthly))
    
    era5_ppt_long <- terra::extract(era5_ppt_monthly, vect(fluxnet_ch4_era5), ID = FALSE) %>%
      mutate(SITE_ID = fluxnet_ch4_era5$SITE_ID) %>%
      pivot_longer(-SITE_ID, names_to = "YearMon", values_to = "era5_p_f") %>%
      mutate(YearMon = zoo::as.yearmon(YearMon))
    
    era5_tmean_long <- terra::extract(era5_tmean_monthly, vect(fluxnet_ch4_era5), ID = FALSE) %>%
      mutate(SITE_ID = fluxnet_ch4_era5$SITE_ID) %>%
      pivot_longer(-SITE_ID, names_to = "YearMon", values_to = "era5_ta_f") %>%
      mutate(YearMon = zoo::as.yearmon(YearMon))
    
    era5_monthly <- era5_ppt_long %>%
      left_join(era5_tmean_long, by = c("SITE_ID", "YearMon"))
    
    find_duplicate_keys(era5_monthly, SITE_ID, YearMon)
    
    range(era5_monthly$YearMon)
    
    before_rows <- nrow(flux_model_base)
    
    flux_model_era5 <- flux_model_base %>%
      left_join(era5_monthly, by = c("SITE_ID", "YearMon"))
    
    after_rows <- nrow(flux_model_era5)
    
    before_rows
    
    after_rows
    
    #### inundation
    
        #lines 299-315 produced by Microsoft Copilot AI chat (prompt: https://copilot.cloud.microsoft/chat/share/eyJzaGFyZUlkIjoiYTJlZjlkNjItMTFlNC00NTAzLWFkOGEtNDg1YTk2YmNhNGI1IiwiY29udmVyc2F0aW9uSWQiOiI2MTM2MzE2Yy05NmNhLTQyNDYtODlkNS0yOTdiYzcyNGQ4ZmEifQ%3D%3D)
    
    flux_coords <- fluxnet_ch4_sf %>%
      st_coordinates() %>%
      as.data.frame()
    
    flux_lookup <- fluxnet_ch4_sf %>%
      st_drop_geometry() %>%
      select(SITE_ID) %>%
      bind_cols(flux_coords) %>%
      mutate(
        Long = round((X - 0.125) / 0.25) * 0.25 + 0.125,
        Lat  = round((Y - 0.125) / 0.25) * 0.25 + 0.125
      ) %>%
      select(SITE_ID, Long, Lat)
    
    lswt_site_monthly <- lswt_final %>%
      inner_join(flux_lookup, by = c("Long", "Lat"))

    nlyr(inundation_monthly)
    
    names(inundation_monthly) <- as.character(zoo::as.yearmon(names(inundation_monthly)))
    
    fluxnet_ch4_inundation <- st_transform(fluxnet_ch4_sf, crs = crs(inundation_monthly))
    
    inundation_long <- terra::extract(inundation_monthly, vect(fluxnet_ch4_inundation), ID = FALSE) %>%
      mutate(SITE_ID = fluxnet_ch4_inundation$SITE_ID) %>%
      pivot_longer(-SITE_ID, names_to = "YearMon", values_to = "inundation_p_f") %>%
      mutate(YearMon = zoo::as.yearmon(YearMon))
    
    inundation_monthly <- inundation_long %>%
      left_join(inundation_tmean_long, by = c("SITE_ID", "YearMon"))
    
    find_duplicate_keys(inundation_monthly, SITE_ID, YearMon)
    
    range(inundation_monthly$YearMon)
    
    before_rows <- nrow(flux_model_era5)
    
    flux_model_inundation <- flux_model_era5 %>%
      left_join(inundation_monthly, by = c("SITE_ID", "YearMon"))
    
    after_rows <- nrow(flux_model_inundation)
    
    before_rows
    
    after_rows
    
#build the model-ready table
    
    NE_methane_model_data <- flux_model_inundation %>%
      select(SITE_ID, SITE_NAME, YearMon,
             FCH4_F_gC,
             NDVI, EVI,
             IGBP,
             soil_ph, soil_nitrogen, soil_bulk_density,
             era5_p_f, era5_ta_f,
             Fw) %>%
      filter(!is.na(FCH4_F_gC))
    
    glimpse(NE_methane_model_data)
    
    nrow(NE_methane_model_data)
    
    n_distinct(NE_methane_model_data$SITE_ID)
    
    range(NE_methane_model_data$YearMon)

  #Data quality dashboard
    
    NE_methane_model_data %>%
      count_missing() %>%
      head(25)
    
    find_duplicate_keys(NE_methane_model_data, SITE_ID, YearMon)
    
    NE_methane_model_data %>%
      mutate(complete_predictors = complete.cases(P_F, TA_F, VPD_F,
                                                  NDVI, EVI,
                                                  soil_ph, soil_nitrogen,
                                                  climate_map, climate_mat,
                                                  elevation_m,
                                                  era5_p_f, era5_ta_f)) %>%
      ungroup() %>%
      reframe(.by = SITE_ID,
              total_rows = n(),
              complete_rows = sum(complete_predictors),
              complete_fraction = complete_rows / total_rows) %>%
      arrange(desc(complete_fraction), desc(complete_rows)) %>%
      head(15)
    
    ggplot(NE_methane_model_data, aes(x = FCH4_F_gC)) +
      geom_histogram(color = "black", fill = "white") +
      labs(title = "Distribution of monthly methane flux observations",
           x = "FCH4_F_gC",
           y = "Count")
    
    #Save work
    
    save(NE_methane_model_data, file = "data/products/NE_methane_model_data.RDATA")
    
    