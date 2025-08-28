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

setwd('/Grad_Storage/FluxTower/GCELTER/models/gap-filling')

#read in data

biophysical <- read_csv(list.files('../../processed',pattern='a_GCE',full.names = T))
noaa_tides <- read_csv(list.files('../../processed',pattern='0_USGCE',full.names=T))

flux <- read_csv(list.files('../../processed',pattern = 'c_GCE',full.names=T)) %>%
  left_join(biophysical %>% rename(end_datetime=datetime)) %>%
  left_join(noaa_tides %>% rename(end_datetime=datetime)) %>%
  mutate(doy = yday(end_datetime),
         year = year(end_datetime),
         hour = hour(end_datetime)) %>%
  mutate(NEE_response = NEEcomb) %>%
  mutate(precip_12hr = rollapply(Total_Precip,width=23,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(PAR_3hr_cum = rollapply(Mean_PAR_Incident_Tower,width=7,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  mutate(Tair_3hr_cum = rollapply(Mean_Temp_Air,width=7,align='right',partial=T,FUN=sum,na.rm=T)) %>%
  filter(qc_co2_flux == 0) %>%
  mutate(wind_cosine = cos(Mean_Wind_Direction),
         wind_sine = sin(Mean_Wind_Direction),
         doy_cosine = cos(doy),
         doy_sine = sin(doy),
         Solar_elv_sine = sin(Solar_elv),
         Solar_elv_cosine = cos(Solar_elv))
  

# summary(flux)

#### define predictors and response =======================================================================


# select variables for modeling
colnames(flux)

#define predictors and response
response <- 'NEE_response'
predictors <- c(
  #light
  'Mean_PAR_Incident_Tower',
  'PAR_3hr_cum',
  'CI',
  'Solar_elv_sine',
  'Solar_elv_cosine',
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
  'hour',
  'doy_cosine',
  'doy_sine',
  'year',
  #other
  'days_prior_maintenance',
  'days_since_maintenance')

#make tibble of selected variables
modeling_flux <- flux[c('end_datetime',response,predictors)] %>%
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
  step_rm(end_datetime) %>%
  step_impute_bag(all_predictors()) %>%
  prep() %>% juice()

modeling_flux_imputed$end_datetime <- modeling_flux$end_datetime

set.seed(42)
first_split <- initial_split(modeling_flux_imputed,prop=3/4,strata = kmeans)
#all training will be done on this data
training <- training(first_split) %>%
  select(-end_datetime)
#final model performance will be determined on this data
holdout <- testing(first_split) 

write.csv(holdout,'holdout-data.csv',row.names = F)

#outer resampling number (split training data into n datasets) How many final models to make
n_outer <- 20

#make outer split ids
training_nest <- training %>%
  #create random groups based on the number of outer folds
  group_by(kmeans) %>%
  mutate(outer_split = sample(1:n_outer,size = n(),replace=T)) %>%
  ungroup() 


#### Model Tuning ================================================================================
#number of folds for cross validation
k <- 5
grid_size <- 50
further_iter <- 10

#create a list to save metrics from each resample
# holdout_perf_list <- list()


cl <- makeCluster(20)
registerDoParallel(cl)
for (i in 1:20){
  
  set.seed(i)
  
  #filter to training data
  rs_training <- training_nest %>%
    filter(outer_split != i) %>%
    select(-outer_split)
  
  # random cv
  k_folds <- vfold_cv(rs_training,v=k,
                      strata=kmeans)
  
  #set up xgboost
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
  rec <-  recipe(NEE_response ~ #light
                   Mean_PAR_Incident_Tower+
                   PAR_3hr_cum+
                   CI+
                   Solar_elv_sine+
                   Solar_elv_cosine+
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
                   hour+
                   doy_cosine+
                   doy_sine+
                   year+
                   #other
                   days_prior_maintenance+
                   days_since_maintenance,
                 data=rs_training)
  
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
                          update(mtry=finalize(mtry(),rs_training)),
                        metrics = metric_set(rmse),
                        control = control_bayes(save_pred = F))
  
  #finalize model on training data
  fit_model <- finalize_workflow(mod_wf,select_best(mod_res,'rmse')) %>%
    fit(data=rs_training)
  
  #save fit model
  # saveRDS(fit_model,paste(i,'model.rds',sep='-'))
  
  #save model
  ft <- fit_model %>% extract_fit_parsnip() %>% .$fit
  xgb.save(ft, paste(i,'model.rds',sep='-'))
  
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
  importance <- fit_model %>%
    extract_fit_parsnip() %>%
    vi(method='permute',
       target='NEE_response',train=rec %>% prep() %>% juice(),metric='rmse',
       pred_wrapper = function(object, newdata) predict(object, newdata = as.matrix(newdata)), #for xgboost, need to make data a matrix
       nsim = 10, parallel = T)
  write.csv(importance,paste(i,'model-feature-importance.csv',sep='-'),row.names=F)
}


#save performance table
# write.csv(bind_rows(holdout_perf_list),'holdout_perf.csv',row.names=F)

#other info including blocking variables and hyperparameter grid sizes
tibble(k=k,grid_size=grid_size,further_iter=further_iter) %>%
  write.csv('cv-grid-size.csv',row.names = F)
