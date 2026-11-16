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
              fluxnet_sites$SITE_NAME <- 
                iconv(fluxnet_sites$SITE_NAME, 
                      from = "", 
                      to = "UTF-8", 
                      sub = "")

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
    
    fluxnet_ch4_era5 <- st_transform(fluxnet_ch4_sf, crs = crs(era5_ppt_monthly))
    
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
    
    
            #### FIRST get data frame from wad2m
    
                    # Source - https://stackoverflow.com/a/71964103
                    # Posted by Robert Hijmans, modified by community. See post 'Timeline' for change history
                    # Retrieved 2026-11-15, License - CC BY-SA 4.0
                    
                    # source 2: https://towardsdatascience.com/how-to-crack-open-netcdf-files-in-r-and-extract-data-as-time-series-24107b70dcd/
                    
                    library(ncdf4)
                    library(tidyverse)
                    wad2m_nc_data <- nc_open("data/final_project/WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc")
                    print(wad2m_nc_data)
                    
                    attributes(wad2m_nc_data$var)
                    attributes(wad2m_nc_data$dim)
                    
                    # Get latitude and longitude with the ncvar_get function and store each into their own object:
                    
                    lat <- ncvar_get(wad2m_nc_data, "lat")
                    nlat <- dim(lat) #to check it matches the metadata: 720
                    
                    lon <- ncvar_get(wad2m_nc_data, "lon")
                    nlon <- dim(lon) #to check, should be 1440
                    
                    # Check your lat lon dimensions match the information in the metadata we explored before:
                    
                    print(c(nlon, nlat))
                    
                    # Get the time variable. Remember: our metadata said our time units are in seconds since 1981-01-01 00:00:00, so you will not see a recognizable date and time format, but a big number like "457185600". We will take care of this later
                    
                    time <- ncvar_get(wad2m_nc_data, "time")
                    head(time) # just to have a look at the numbers
                    tunits <- ncatt_get(wad2m_nc_data, "time", "units") #check units
                    nt <- dim(time) #should be 252
                    
                    #get the variable in "matrix slices"
                    lswt_array <- ncvar_get(wad2m_nc_data, "Fw") 
                    
                    fillvalue <- ncatt_get(wad2m_nc_data, "Fw", "_FillValue")
                    
                    dim(lswt_array) #to check; this should give you 1440 720 252
                    
                    #right away let's replace the nc FillValues with NAs
                    lswt_array[lswt_array==fillvalue$value] <- NA
                    lswt_array
                    
                    # human-readable time -----------------------------------------------------
                    
                    origin_date <- as.Date("1992-01-01")  # define origin date
                    
                    # Convert list to vector for easier processing
                    days_vec <- unlist(time)
                    
                    # Convert days to Date objects
                    dates <- origin_date + days_vec
                    
                    # Format into Year-Month (YYYY-MM)
                    year_month <- format(dates, "%Y-%m")
                    
                    # Output
                    print(year_month)
                    range(year_month)
                    
                    #Create 2D matrix of long, lat and time
                    
                    lonlattime <- as.matrix(expand.grid(lon,lat,year_month)) # this might take several seconds
                    
                    #reshape whole lswt_array
                    
                    lswt_vec_long <- as.vector(lswt_array)
                    length(lswt_vec_long)
                    
                    #Create data.frame
                    
                    lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
                    
                    head(lswt_obs)
                    
                    colnames(lswt_obs) <- c("Long","Lat","Date","Fw")
                    
                    head(lswt_obs)
                    
                    # Bold move time: remove all rows with NA in LSWT_final:
                    lswt_final <- na.omit(lswt_obs)
                    
                    dim(lswt_obs)
                    dim(lswt_final)
                    head(lswt_final)
                    
                    lswt_final <- lswt_final %>%
                      mutate(
                        Long = as.numeric(Long),
                        Lat = as.numeric(Lat)
                      )
                    
                    #lines 338-357 below are generated from Microsoft Copilot AI chat (prompt:https://copilot.cloud.microsoft/chat/share/eyJzaGFyZUlkIjoiYTJlZjlkNjItMTFlNC00NTAzLWFkOGEtNDg1YTk2YmNhNGI1IiwiY29udmVyc2F0aW9uSWQiOiI2MTM2MzE2Yy05NmNhLTQyNDYtODlkNS0yOTdiYzcyNGQ4ZmEifQ%3D%3D)
                    
                    
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
                      inner_join(
                        flux_lookup,
                        by = c("Long", "Lat"),
                        relationship = "many-to-many"
                      )
                    
          #now finish extracting to tower sites
          
          find_duplicate_keys(lswt_site_monthly, SITE_ID, Date)
          
          lswt_monthly <- lswt_site_monthly %>%
          mutate(
          YearMon = zoo::as.yearmon(Date)
          ) %>%
          select(SITE_ID, YearMon, Fw)
                    
          before_rows <- nrow(flux_model_era5)
                    
          flux_model_lswt <- flux_model_era5 %>%
          left_join(
          lswt_monthly,
          by = c("SITE_ID", "YearMon")
          )
                    
          after_rows <- nrow(flux_model_lswt)
                    
          before_rows
          after_rows
    
#build the model-ready table
    
    NE_methane_model_data <- flux_model_lswt %>%
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
      mutate(complete_predictors = complete.cases(NDVI, EVI,
                                                  soil_ph, soil_nitrogen,
                                                  era5_p_f, era5_ta_f,
                                                  Fw)) %>%
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
    
    