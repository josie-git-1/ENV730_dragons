
# To start random forest model using the NE_methane_model_data.RDATA, use this code

# -------------------------------------------------------------------------

#set working directory

setwd("C:/Users/josie/Desktop/YSE/F26/ENV730/EDS730")

# install packages and load libraries

install.packages("randomForest")
install.packages("tidyverse")
install.packages("GGally")

library(randomForest)
library(tidyverse)
library(GGally)

# load the modeling data

load("data/products/NE_methane_model_data.RDATA")
names(NE_methane_model_data)
glimpse(NE_methane_model_data)

#modeling checklist

  #before fitting any model, identify:
    #response variable: FCH4_F_gC;
    #climate predictors: era5_p_f, era5_ta_f;
    #vegetation predictors: NDVI, EVI;
    #soil predictors: soil_ph, soil_nitrogen, soil_bulk_density
    #ecosystem predictors: IGBP, Fw;
    #grouping variable for site-level assessment: SITE_NAME.

# (variables with _F mean data gaps were filled in the dataset to account for sampling errors)

#prepare the modeling table
model_data_org <- NE_methane_model_data %>%
  ungroup() %>%
  select(SITE_NAME, SITE_ID, YearMon, FCH4_F_gC,
         era5_p_f, era5_ta_f,
         soil_ph, soil_nitrogen, soil_bulk_density,
         NDVI, EVI,
         IGBP, Fw) %>%
  drop_na()

nrow(model_data_org)
n_distinct(model_data_org$SITE)
summary(model_data_org)

