#Author: Peter Hawman

#Purpose:
# 1. read in combined sensor QC'd data
# 2. fill gaps in NEE with established ML models for:
#     a. all directions, creek, interior


#Output:
# 1. complete timeseries of NEE for:
#     a. full footprint
#     b. creek to the south
#     c. interior to the north

#set wd
# setwd('Z:/FluxTower/GCELTER')
setwd('/Grad_Storage/FluxTower/GCELTER')

library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(tidyselect)
library(stringr)
library(ggplot2)
library(lubridate)
library(tidymodels) # machine learning tools and packages
library(doParallel) # for parallel processing
library(zoo)
library(xgboost)
library(multidplyr)

#set plotting theme
theme_set(theme_bw())

#read in data
fluxcombined <- read_csv(list.files('processed',pattern='c_GCE_Fluxes',full.names=T))
biophysical <- read_csv(list.files('processed',pattern='a_GCE_biophysical',full.names=T))
noaa_tides <- read_csv(list.files('processed',pattern='0_USGCE-noaa',full.names=T))

flux <- biophysical %>% rename(end_datetime=datetime) %>% left_join(fluxcombined,by='end_datetime') %>%
  left_join(noaa_tides %>% rename(end_datetime=datetime),by='end_datetime') %>%
  mutate(doy = yday(end_datetime),
         year = year(end_datetime),
         hour = hour(end_datetime)) %>%
  mutate(NEE_response = NEEcomb) %>%
  mutate(precip_12hr = rollapply(Total_Precip,width=23,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(PAR_3hr_cum = rollapply(Mean_PAR_Incident_Tower,width=7,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(Tair_3hr_cum = rollapply(Mean_Temp_Air,width=7,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(wind_cosine = cos(Mean_Wind_Direction),
         wind_sine = sin(Mean_Wind_Direction),
         doy_cosine = cos(doy),
         doy_sine = sin(doy),
         Solar_elv_sine = sin(Solar_elv),
         Solar_elv_cosine = cos(Solar_elv))

rm(fluxcombined,biophysical,noaa_tides)

#models
model_path <- 'models/gap-filling'


#### prep data for xgboost prediction ==============================================================================

#read in predictor names
predictor_list <- read_csv(paste(model_path,'1-model-predictors.csv',sep='/')) %>% pull(predictors)

model_data <- flux %>% select(NEEcomb,all_of(predictor_list)) %>%
  recipe(NEEcomb ~ .,data=.) %>%
  step_impute_bag(all_predictors()) %>%
  prep() %>% juice() 


### predict using xgboost models ===================================================================================

#convert to matrix for xgboost predictions
full_model_data <- model_data %>% select(all_of(predictor_list)) %>% as.matrix()

#list all models
model_list <- tibble(paths = list.files(model_path,recursive = F,full.names = T,pattern = '.rds'),
                     name = list.files(model_path,recursive = F,full.names = F,pattern = '.rds')) %>%
  mutate(model = sub('.*/','',name)) %>%
  mutate(num = as.numeric(gsub("\\D", "", model))) %>%
  arrange(num)

#predict using each model
for (i in 1:nrow(model_list)){
  #read model
  model <- xgb.load(model_list$paths[i])
  #predict
  flux$model <- predict(model,newdata = full_model_data)
  #update column name
  colnames(flux)[names(flux)=='model'] <- paste('model',i,sep='_')
  #remove model
  rm(model)
}

colnames(flux)

### average predictions and fill gaps ===================================================================================

cluster <- new_cluster(10)
cluster_library(cluster,'dplyr')

flux_avg <- flux %>%
  #get averages of model predictions
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    model_mean = mean(c_across(starts_with('model'))),
    model_median = median(c_across(starts_with('model')))) %>%
  collect() %>%
  ungroup()

flux_filled <- flux_avg %>%
  #fill gaps where missing NEE
  mutate(
    NEE_xgbGF = ifelse(is.na(NEEcomb),model_median,NEEcomb)) %>%
  # make separate columns for only "0" quality fluxes and gap-filling predictions
  mutate(
    NEE_xgbGF_HQ = ifelse(qc_co2_flux != 0, model_median,NEE_xgbGF)) %>%
  mutate(
    NEE_xgbGF_HQ = ifelse(is.na(NEE_xgbGF_HQ), model_median,NEE_xgbGF_HQ)) %>%
  #for high quality, use model predictions for NEE at night < 0
  mutate(
  NEE_xgbGF_HQ = case_when(Solar_elv < -10 & NEEcomb < 0 ~ model_median,
                                      TRUE ~ NEE_xgbGF_HQ)) %>%
  mutate(NEE_xgbGF_HQ = ifelse(between(NEE_xgbGF_HQ,-20,10),NEE_xgbGF_HQ,model_median)) %>%
  mutate(start_datetime = end_datetime - 1800) %>%
  select(start_datetime,end_datetime,PT_sensor_used:NEE_xgbGF_HQ) %>%
  arrange(start_datetime)

# 
# flux_filled %>%
#   pivot_longer(cols=c(ends_with('NEE_xgbGF_HQ'))) %>%
#   ggplot(aes(end_datetime,value))+
#   geom_line()+
#   facet_wrap(~name,ncol=1)



flux_filled %>% colnames()

summary(flux_filled)

# flux_filled %>%
#   mutate(roll_mean = rollapply(NEE_xgbGF_HQ,width=5,align = 'center',fill=NA,FUN=function(x) mean(x[-3], na.rm = TRUE),partial=T)) %>%
#   mutate(roll_sd = rollapply(NEE_xgbGF_HQ,width=5,align = 'center',fill=NA,FUN=function(x) sd(x[-3], na.rm = TRUE),partial=T)) %>%
#   ungroup() %>%
#   #remove outliers
#   mutate(lower_bound = roll_mean-roll_sd*1,
#          # upper_bound = roll_mean+roll_sd*1.5) %>%
#          upper_bound = roll_mean+roll_sd*1) %>%
#   mutate(NEE_xgbGF_HQ = ifelse(between(NEE_xgbGF_HQ,lower_bound,upper_bound),NEE_xgbGF_HQ,model_mean)) %>%
#   mutate(status = case_when(is.na(NEEcomb) ~ 'gap-filled',
#                             !is.na(NEEcomb) ~ 'measured')) %>%
#   filter(!year %in% c(2013,2018,2025)) %>%
#   # filter(between(doy,350,355)) %>%
#   mutate(minute = minute(end_datetime)) %>%
#   mutate(common_date = as.POSIXct(paste('3000',doy,hour,minute,sep='-'),format='%Y-%j-%H-%M')) %>%
#   filter(week(common_date)==40) %>%
#   # select(end_datetime,doy,hour,minute,common_date)
#   ggplot(aes(common_date,NEE_xgbGF_HQ,color=factor(year)))+
#   geom_line()+
#   geom_point(aes(shape=status),size=3)+
#   scale_shape_manual(values = c(1,18))


#save to file
write.csv(flux_filled,paste0('processed/d_GCE_Fluxes_qcd_comb_gapfilled_',
                                format(range(as.Date(flux_filled$end_datetime))[1],'%Y%m%d'),'_',
                                format(range(as.Date(flux_filled$end_datetime))[2],'%Y%m%d'),
                                '.csv'),row.names = F)


