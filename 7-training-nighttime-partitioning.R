library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(tidymodels) # machine learning tools and packages
library(doParallel) # for parallel processing
library(lubridate)
library(vip)
library(zoo)
library(xgboost)
# library(ranger)

setwd('/Grad_Storage/FluxTower/GCELTER/models/partitioning')

# setwd('z:/FluxTower/GCELTER/models/partitioning')

biophysical <- read_csv(list.files('../../processed',pattern='a_GCE',full.names = T))
noaa_tides <- read_csv(list.files('../../processed',pattern='0_USGCE',full.names=T))

#read in data
flux <- read_csv(list.files('../../processed',pattern='d_GCE',full.names=T))


  
flux_night <- flux %>%
  left_join(biophysical %>% rename(end_datetime=datetime)) %>%
  left_join(noaa_tides %>% rename(end_datetime=datetime)) %>%
  mutate(doy = yday(end_datetime),
         year = year(end_datetime),
         hour = hour(end_datetime)) %>%
  mutate(NEE_response = NEEcomb) %>%
  mutate(precip_12hr = rollapply(Total_Precip,width=23,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(Tair_3hr_cum = rollapply(Mean_Temp_Air,width=7,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(Tsoil_3hr_cum = rollapply(Mean_Temp_Soil_Tower,width=7,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  filter(qc_co2_flux == 0) %>%
  filter(Solar_elv < -10) %>%
  filter(Mean_PAR_Incident_Tower == 0 & PARTOA == 0) %>%
  filter(NEE_response >=0) %>%
  filter(NEE_response < 7.5) %>%
  mutate(wind_cosine = cos(Mean_Wind_Direction),
         wind_sine = sin(Mean_Wind_Direction),
         doy_cosine = cos(doy),
         doy_sine = sin(doy),
         Solar_elv_sine = sin(Solar_elv),
         Solar_elv_cosine = cos(Solar_elv))

# flux_night %>% ggplot(aes(end_datetime,NEE_response))+geom_point()

#nighttime average NEE from gap-filled data
night_NEE <- flux_night %>%
  filter(Solar_elv < -15) %>%
  mutate(end_datetime = end_datetime-5*3600) %>%
  group_by(date=as.Date(end_datetime)) %>%
  summarise(NEEnight = median(NEE_xgbGF)) %>%
  mutate(NEEnight_daylag = lag(NEEnight)) %>%
  mutate(NEEnight_weekavg = rollmean(NEEnight,k=7,align='right',fill=NA))
  # ggplot(aes(date,NEEnight))+geom_point()
  
#daytime estiamtes of stuff
day_stuff <- flux %>%
  filter(Mean_PAR_Incident_Tower > 0) %>%
  mutate(end_datetime = end_datetime-5*3600) %>%
  group_by(date=as.Date(end_datetime)) %>%
  summarise(PAR_day_cum = sum(Mean_PAR_Incident_Tower),
            Temp_day_cum = sum(Mean_Temp_Air))


flux_together <- flux_night %>% mutate(date = as.Date(end_datetime-5*3600) )%>% 
  left_join(night_NEE,by='date') %>%
  left_join(day_stuff,by='date')
# colnames(flux_together)

# flux_together %>% ggplot(aes(NEE_response))+geom_histogram()
# flux_together %>% ggplot(aes(end_datetime,NEE_response))+geom_point()

# summary(flux_together)

#### define predictors and response =======================================================================


# select variables for modeling
colnames(flux_together)

#define predictors and response
response <- 'NEE_response'
predictors <- c(
  #average night NEE
  'NEEnight',
  'NEEnight_daylag',
  'NEEnight_weekavg',
  #daytime things
  'PAR_day_cum',
  'Temp_day_cum',
  #air
  'wind_cosine',
  'wind_sine',
  'wind_speed',
  'Mean_Humidity',
  'VPD',
  'Mean_Temp_Air',
  'Tair_3hr_cum',
  #precip
  'Total_Precip',
  'precip_12hr',
  #soil conditions
  'Mean_Temp_Soil_Tower',
  'Mean_salinity',
  'water_level_combined',
  'hours_from_high_tide',
  'noaa_navd88',
  #seasonality and time
  'spring','summer','autumn','winter',
  'morning','afternoon','evening','night',
  'doy_cosine',
  'doy_sine')

#make tibble of selected variables
modeling_flux <- flux_together[c(response,predictors)] %>%
  filter(!is.na(NEE_response)) 


#### data splitting =====================================================
set.seed(42)
kmeans_dat <- modeling_flux %>% 
  select(spring:night)

Kclusters <- kmeans(kmeans_dat,10)

modeling_flux$kmeans <- as.factor(Kclusters$cluster)

# modeling_flux %>% ggplot(aes(datetime,NEE_response,color=k))+geom_point()


#### data imputation =====================================================
modeling_flux_imputed <- recipe(NEE_response ~ .,data=modeling_flux) %>%
  step_impute_bag(all_predictors()) %>%
  prep() %>% juice()

set.seed(42)
first_split <- initial_split(modeling_flux_imputed,prop=0.8,strata = kmeans)
#all training will be done on this data
training <- training(first_split)
#final model performance will be determined on this data
holdout <- testing(first_split) %>% select(-kmeans)

write.csv(holdout,'holdout-data.csv',row.names = F)
write.csv(training,'traning-data.csv',row.names=F)

#outer resampling number (split training data into n datasets) How many final models to make
n_outer <- 20


#### Model Tuning ================================================================================
#number of folds for cross validation
k <- 3
grid_size <- 15
further_iter <- 10

#create a list to save metrics from each resample
# holdout_perf_list <- list()


cl <- makeCluster(5)
registerDoParallel(cl)
for (i in 7:n_outer){
  
  set.seed(i)
  
  # random cv
  k_folds <- vfold_cv(training,v=k,
                      strata=kmeans)
  
  #set up model
  mod_spec <- boost_tree(trees = tune(), 
                         tree_depth = tune(), 
                         min_n = tune(), 
                         loss_reduction = tune(),
                         sample_size = tune(), 
                         mtry = tune(),
                         learn_rate = tune()) %>%
    set_engine('xgboost') %>%
    set_mode('regression')
  #recipe
  rec <-  recipe(NEE_response ~
                   #average night NEE
                   NEEnight+
                   NEEnight_daylag+
                   NEEnight_weekavg+
                   #daytime things
                   PAR_day_cum+
                   Temp_day_cum+
                   #air
                   wind_cosine+
                   wind_sine+
                   wind_speed+
                   Mean_Humidity+
                   VPD+
                   Mean_Temp_Air+
                   Tair_3hr_cum+
                   #precip
                   Total_Precip+
                   precip_12hr+
                   #soil conditions
                   Mean_Temp_Soil_Tower+
                   Mean_salinity+
                   water_level_combined+
                   hours_from_high_tide+
                   noaa_navd88+
                   #seasonality and time
                   spring+summer+autumn+winter+
                   morning+afternoon+evening+night+
                   doy_cosine+
                   doy_sine,
                 data=training)
  
  #create a workflow and add the recipe and model specifics
  mod_wf <- workflow() %>%
    add_model(mod_spec) %>%
    add_recipe(rec)
  
  #tune across the grid
  mod_res <- tune_bayes(mod_wf,
                        resamples=k_folds,
                        initial=grid_size,
                        iter=further_iter,
                        param_info = hardhat::extract_parameter_set_dials(mod_spec) %>%
                          update(mtry=finalize(mtry(),training)),
                        metrics = metric_set(rmse),
                        control = control_bayes(save_pred = F))
  
  #finalize model on training data
  fit_model <- finalize_workflow(mod_wf,select_best(mod_res,'rmse')) %>%
    fit(data=training)
  
  
  #save model
  ft <- fit_model %>% extract_fit_parsnip() %>% .$fit
  xgb.save(ft, paste(i,'model.rds',sep='-'))
  # saveRDS(ft, paste(i,'model.rds',sep='-'))
  
  
  #predict on testing data
  holdout <- holdout %>%
    mutate(pred = predict(fit_model,holdout)$.pred) 
  #save predictions for each model
  write.csv(holdout,paste(i,'model-predicted-holdout.csv',sep='-'),row.names=F)
  
  #model performance on holdout data
  holdout_regress <- holdout %>%
    nest(data=everything()) %>%
    mutate(fit = map(data,~lm(pred ~ NEE_response,data=.x)),
           glanced = map(fit,glance),
           tidied = map(fit,tidy),
           augmented = map(fit,augment))
  
  #put performace testing together
  holdout_perf <- bind_rows(holdout %>%
                              yardstick::rmse(NEE_response,pred) %>%
                              mutate(id = 'holdout') %>%
                              select(name=`.metric`,value=`.estimate`),
                            holdout_regress %>% 
                              unnest(tidied) %>%
                              mutate(term = ifelse(term == 'NEE_response','slope','intercept')) %>%
                              select(name=term,value=estimate),
                            
                            holdout_regress %>%
                              unnest(glanced) %>%
                              select(r.squared) %>%
                              pivot_longer(cols=`r.squared`),
                            
                            holdout_regress %>%
                              unnest(tidied) %>%
                              filter(term == 'NEE_response') %>%
                              unnest(augmented) %>%
                              summarise(nrmse = rmse_vec(NEE_response,pred)/(max(NEE_response,na.rm=T)-min(NEE_response,na.rm=T)),
                                        sspe = sum((NEE_response-pred)^2),
                                        Ubias = (n()*(mean(NEE_response)-mean(pred))^2)/sspe,
                                        Ub1 = (((mean(estimate)-1)^2)*sum((pred-mean(pred))^2))/sspe,
                                        Ue = sum((.fitted-NEE_response)^2)/sspe) %>%
                              pivot_longer(cols=c(nrmse:Ue))) %>%
    mutate(model = paste(i,'model',sep='-'))
  
  #save to list
  # holdout_perf_list[[i]] <- holdout_perf
  write.csv(holdout_perf,paste(i,'model-performance.csv',sep='-'))
  
  
  #save predictor list
  tibble(predictors = colnames(extract_mold(fit_model)$predictors)) %>%
    write.csv(paste(i,'model-predictors.csv',sep='-'),row.names=F)
  

  # variable importance
  # variable importance
  importance <- fit_model %>%
    extract_fit_parsnip() %>%
    vi(method='permute',
       target='NEE_response',train=rec %>% prep() %>% juice(),metric='rmse',
       pred_wrapper = function(object, newdata) predict(object, newdata = as.matrix(newdata)), #for xgboost, need to make data a matrix
       nsim = 20, parallel = T)
  write.csv(importance,paste(i,'model-feature-importance.csv',sep='-'),row.names=F)
}


#save performance table
# write.csv(bind_rows(holdout_perf_list),'holdout_perf.csv',row.names=F)

#other info including blocking variables and hyperparameter grid sizes
tibble(k=k,grid_size=grid_size,further_iter=further_iter) %>%
  write.csv('cv-grid-size.csv',row.names = F)
