#Author: Peter Hawman

#Purpose:
# 1. read in gap-filled and partitioned fluxes
# 2. create daily fluxes
#     a. all directions, 
#     b. creek
#     c. interior

#Output:
# 1. daily fluxes and environment
#     a. full footprint
#     b. creek to the south
#     c. interior to the north

#set wd
setwd('Z:/FluxTower/GCELTER')

library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(tidyselect)
library(stringr)
library(ggplot2)
library(lubridate)
library(zoo)
library(VGAM) #for random draws from Laplace distribution

theme_set(ggthemes::theme_base())

#read in gap-filled and partitioned data
flux <- read_csv(list.files('processed',pattern='e_GCE',full.names = T))

colnames(flux)

#### daily NEE unceretainty ============================================================================

#set up dataframe
nee <- flux %>%
  mutate(NEE_meas = NEEcomb) %>%
  select(datetime=start_datetime,qc_co2_flux,NEE_meas,NEE_xgbGF_HQ,starts_with('model')) %>%
  #get the residuals between the meas NEE and averaged predicted value
  mutate(resid = model_median - NEE_meas)

nee %>% ggplot()+geom_histogram(aes(resid))


nee_gf_uncert <-  nee %>%
  filter(is.na(NEE_meas)) %>%
  pivot_longer(cols=c(contains('model'))) %>%
  filter(!str_detect(name,'me')) %>%
  group_by(datetime) %>%
  summarise(gapfill_variance = var(value)) %>%
  ungroup()

nee_gf_uncert %>% ggplot(aes(datetime,gapfill_variance))+geom_line()

#bin edges
# edges <- c(-25,seq(-14,0,5),seq(2,16,2))
nee %>% ggplot(aes(NEE_xgbGF_HQ))+geom_histogram()
edges <- c(seq(-20,15,5))
#number of bins
nbins <- length(edges)-1

#random uncertainty from random error for measured points 
out_list <- list()
for (i in 1:nbins) {
  
  #filter to bin i
  x <- nee %>% filter(NEE_meas >= edges[i] & NEE_meas < edges[i+1])
  
  n <- nrow(x) #number of meas
  bin_mean <- mean(x$NEE_meas,na.rm=T) #mean of bin
  
  std_resids <- sd(x$resid,na.rm=T)
  
  out_list[[i]] <- tibble(edge1=edges[i],
                          edge2=edges[i+1],
                          n=n,
                          NEE=bin_mean,
                          DEsigma=std_resids)
  
}

#initial variables 
init_vars <- bind_rows(out_list) %>%
  # mutate(DEsigma = ifelse(is.na(DEsigma),lag(DEsigma),DEsigma)) %>%
  # mutate(n = ifelse(n==0,lag(n),n)) %>%
  # mutate(DEsigma = ifelse(is.na(DEsigma),lead(DEsigma),DEsigma)) %>%
  # mutate(n = ifelse(is.na(n),lead(n),n)) %>%
  print()

#Random error estimate
set.seed(42)
draw_list <- list()
rand_list <- list()
for (i in 1:nrow(init_vars)) {
  
  #filter to bin i
  x <- nee %>% filter(NEE_xgbGF_HQ >= init_vars$edge1[i] & NEE_xgbGF_HQ < init_vars$edge2[i])
  
  #random draw from laplace distribution
  for (j in 1:500){

    draw <- rlaplace(n=nrow(x),location=0,scale=init_vars$DEsigma[i])

    y <- tibble(datetime=x$datetime,draw,j)

    draw_list[[j]] <- y
  }

  rand_list[[i]] <- bind_rows(draw_list)

}


#cumulative random uncertainty (variance) of measured values
nee_meas_uncert <- bind_rows(rand_list) %>%
  group_by(datetime) %>%
  summarise(measured_variance = var(draw)) %>%
  ungroup() %>%
  arrange(datetime)

nee %>% filter(!datetime %in% nee_meas_uncert$datetime)

nee_meas_uncert %>% ggplot()+geom_histogram(aes(measured_variance))
nee_meas_uncert %>% ggplot(aes(datetime,measured_variance))+geom_line()

#get the total uncertainty from measurements and gap-filling
NEE_uncert <- nee_meas_uncert %>%
  left_join(nee_gf_uncert,by='datetime') %>%
  replace(is.na(.),0) %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance)) %>%
  mutate(NEE_CI_95 = totalUncStd*2) %>% print()

nee %>%
  filter(week(datetime)==20) %>%
  left_join(NEE_uncert %>% select(datetime,NEE_CI_95),by='datetime') %>%
  ggplot(aes(datetime,NEE_xgbGF_HQ))+
  geom_errorbar(aes(ymin=NEE_xgbGF_HQ-NEE_CI_95,ymax=NEE_xgbGF_HQ+NEE_CI_95))+
  geom_point()+
  facet_wrap(~year(datetime),scales='free_x')

#### ER uncertainty ============================================================================

#Find random uncertainty using partition model residuals and add nee uncertainty
#UNCERTAINTY COMPOUNDS AS FLUXES ARE DERIVED

#set up dataframe
er <- flux %>%
  mutate(NEE_meas = NEEcomb) %>%
  mutate(day_night = ifelse(Solar_elv < -15, 'night','day')) %>%
  select(datetime=start_datetime,day_night,NEE_meas,ER,starts_with('ER_model')) %>%
  #get the residuals between the meas NEE and averaged predicted value
  mutate(resid = ER_model_median - NEE_meas)

#uncertainty from ER model predictions. The variance 20 predictions
er_model_uncert <- er %>%
  filter(day_night == 'day') %>%
  pivot_longer(cols=c(contains('model'))) %>%
  filter(!str_detect(name,'me')) %>%
  group_by(datetime) %>%
  summarise(resp_model_variance = var(value)) %>%
  ungroup() %>% print()

#combining uncertainties
ER_uncert <-nee_meas_uncert %>%
  left_join(nee_gf_uncert,by='datetime') %>%
  left_join(er_model_uncert,by='datetime') %>%
  replace(is.na(.),0) %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance + resp_model_variance)) %>%
  mutate(ER_CI_95 = totalUncStd*2) %>% print()

ER_uncert %>% summary()

er %>%
  filter(week(datetime)==20) %>%
  left_join(ER_uncert %>% select(datetime,ER_CI_95),by='datetime') %>%
  ggplot(aes(datetime,ER))+
  geom_hline(yintercept = 0, linetype='dotted')+
  geom_errorbar(aes(ymin=ER-ER_CI_95,ymax=ER+ER_CI_95))+
  geom_point()+
  facet_wrap(~year(datetime),scales='free_x')

#### daily GPP uncertainty ============================================================================
gpp <- flux %>%
  mutate(NEE_meas = NEEcomb) %>%
  mutate(day_night = ifelse(Solar_elv < -15, 'night','day')) %>%
  select(datetime=start_datetime,day_night,NEE_meas,NEE_xgbGF_HQ,ER,GPP,starts_with('ER_model'))

gpp_split_uncert <-  gpp %>%
  pivot_longer(cols=c(contains('ER_model'))) %>%
  filter(!str_detect(name,'me')) %>%
  mutate(pseudo_gpp = value-NEE_xgbGF_HQ) %>%
  group_by(datetime) %>%
  summarise(split_variance = var(pseudo_gpp)) %>%
  ungroup()

#combining uncertainties
GPP_uncert <- nee_meas_uncert %>%
  left_join(nee_gf_uncert,by='datetime') %>%
  left_join(er_model_uncert,by='datetime') %>%
  left_join(gpp_split_uncert,by='datetime') %>%
  replace(is.na(.),0) %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance + resp_model_variance + split_variance)) %>%
  mutate(GPP_CI_95 = totalUncStd*2) %>% print()

GPP_uncert %>% summary()

gpp %>%
  filter(week(datetime)==20) %>%
  left_join(GPP_uncert %>% select(datetime,GPP_CI_95),by='datetime') %>%
  ggplot(aes(datetime,GPP))+
  geom_hline(yintercept = 0, linetype='dotted')+
  geom_errorbar(aes(ymin=GPP-GPP_CI_95,ymax=GPP+GPP_CI_95))+
  geom_point()+
  facet_wrap(~year(datetime),scales='free_x')
#
#### PUT IT ALL TOGETHER =============================================================================


out <- flux %>%
  left_join(NEE_uncert %>% select(start_datetime=datetime,NEE_CI_95),by='start_datetime') %>%
  left_join(ER_uncert %>% select(start_datetime=datetime,ER_CI_95),by='start_datetime') %>%
  left_join(GPP_uncert %>% select(start_datetime=datetime,GPP_CI_95),by='start_datetime') %>%
  mutate(sensor = ifelse(is.na(NEEcomb),'gf',sensor))

out %>% count(sensor)

out %>% summary()

# #save to file
write.csv(out,paste0('processed/f_GCE_Fluxes_30mim_uncertainty_',
                                  format(range(as.Date(out$date))[1],'%Y%m%d'),'_',
                                  format(range(as.Date(out$date))[2],'%Y%m%d'),
                                  '.csv'),row.names = F)

