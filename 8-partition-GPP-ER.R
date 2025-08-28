#Author: Peter Hawman

#Purpose:
# 1. read in gap-filled NEE
# 2. Partition NEE into ER and GPP:
#     a. all directions, creek, interior


#Output:
# 1. complete timeseries of GPP, ER, and NEE for:
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

#read in gap-filled data
flux <- read_csv(list.files('processed',pattern='d_GCE_Fluxes',full.names=T)) %>%
  mutate(wind_cosine = cos(Mean_Wind_Direction),
         wind_sine = sin(Mean_Wind_Direction),
         doy_cosine = cos(doy),
         doy_sine = sin(doy),
         hour_cosine = cos(hour),
         hour_sine = sin(hour),
         Solar_elv_sine = sin(Solar_elv),
         Solar_elv_cosine = cos(Solar_elv)) %>%
  mutate(date = as.Date(end_datetime-5*3600)) %>%
  mutate(NEE_response = NEEcomb)


#models
full_model_path <- 'models/partitioning'


#### predict ER for full footprint ==============================================================================

#nighttime average NEE from gap-filled data
night_NEE <- flux %>%
  filter(Solar_elv < -15) %>%
  group_by(date) %>%
  summarise(NEEnight = median(NEE_xgbGF_HQ)) %>%
  mutate(NEEnight_daylag = lag(NEEnight)) %>%
  mutate(NEEnight_weekavg = rollmean(NEEnight,k=7,align='right',fill=NA))

#daytime estiamtes of stuff
day_stuff <- flux %>%
  filter(Mean_PAR_Incident_Tower > 0) %>%
  group_by(date) %>%
  summarise(PAR_day_cum = sum(Mean_PAR_Incident_Tower),
            Temp_day_cum = sum(Mean_Temp_Air))


flux_together <- flux %>%
  left_join(night_NEE,by='date') %>%
  left_join(day_stuff,by='date')

#read in predictor names
predictor_list <- read_csv(paste(full_model_path,'1-model-predictors.csv',sep='/')) %>% pull(predictors)

model_data <- flux_together %>% select(NEE_response,all_of(predictor_list)) %>%
  recipe(NEE_response ~ .,data=.) %>%
  step_impute_bag(all_predictors()) %>%
  prep() %>% juice() 

summary(model_data)

#convert to matrix for xgboost predictions
full_model_data <- model_data %>% select(all_of(predictor_list)) %>% as.matrix()

#list all models
model_list <- tibble(paths = list.files(full_model_path,recursive = F,full.names = T,pattern = '.rds'),
                     name = list.files(full_model_path,recursive = F,full.names = F,pattern = '.rds')) %>%
  mutate(model = sub('.*/','',name)) %>%
  mutate(num = as.numeric(gsub("\\D", "", model))) %>%
  arrange(num)

#predict using each model
for (i in 1:nrow(model_list)){
  #read model
  model <- xgb.load(model_list$paths[i])
  #predict
  flux$ER_model <- predict(model,newdata = full_model_data)
  #update column name
  colnames(flux)[names(flux)=='ER_model'] <- paste('ER_model',i,sep='_')
  #remove model
  rm(model)
}

colnames(flux)


### average predictions and fill gaps ===================================================================================

cluster <- new_cluster(12)
cluster_library(cluster,'dplyr')

flux_er <- flux %>%
  #get averages of model predictions
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    ER_model_mean = mean(c_across(starts_with('ER_model'))),
    ER_model_median = median(c_across(starts_with('ER_model')))) %>%
  collect() %>%
  ungroup()

#partion NEE into GPP and ER for each sector
flux_partitioned <- flux_er %>%
  mutate(GPP = ER_model_mean-NEE_xgbGF_HQ) %>%
  #if dark, set GPP to Zero
  mutate(GPP_zero = ifelse(Solar_elv > -10, GPP,0)) %>%
  #if GPP negative, set to zero
  mutate(GPP_zero = ifelse(GPP_zero <0, 0,GPP_zero)) %>%
  #set respiration to model predicted during daytime, NEE during nighttme. 
  mutate(ER = ifelse(Solar_elv < -10, NEE_xgbGF_HQ,ER_model_mean)) %>%
  arrange(start_datetime)


flux_partitioned %>%
  pivot_longer(cols=c(starts_with('GPP_'))) %>%
  ggplot(aes(end_datetime,value))+
  geom_line()+
  facet_wrap(~name,ncol=1)

flux_partitioned %>%
  pivot_longer(cols=c(ER)) %>%
  ggplot(aes(end_datetime,value))+
  geom_line()+
  facet_wrap(~name,ncol=1)


#save to file
write.csv(flux_partitioned,paste0('processed/e_GCE_Fluxes_qcd_comb_gapfilled_partitioned_',
                             format(range(as.Date(flux_partitioned$end_datetime))[1],'%Y%m%d'),'_',
                             format(range(as.Date(flux_partitioned$end_datetime))[2],'%Y%m%d'),
                             '.csv'),row.names = F)


