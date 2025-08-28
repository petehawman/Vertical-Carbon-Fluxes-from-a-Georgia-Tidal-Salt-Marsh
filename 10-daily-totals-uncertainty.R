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
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  mutate(NEE_meas = NEEcomb) %>%
  select(datetime,date,qc_co2_flux,NEE_meas,NEE_xgbGF_HQ,starts_with('model')) %>%
  #get the residuals between the meas NEE and averaged predicted value
  mutate(resid = model_median - NEE_meas)

nee %>% ggplot()+geom_histogram(aes(resid))


nee_gf_uncert <-  nee %>%
  filter(is.na(NEE_meas)) %>%
  pivot_longer(cols=c(contains('model'))) %>%
  filter(!str_detect(name,'me')) %>%
  group_by(date=as.Date(datetime),name) %>%
  summarise(gapfill_er = sum(value*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(gapfill_variance = var(gapfill_er)) %>%
  print()

nee_gf_uncert %>% ggplot(aes(date,gapfill_variance))+geom_line()+scale_x_date(date_breaks='1 year',date_labels='%Y')

#bin edges
# edges <- c(-25,seq(-14,0,5),seq(2,16,2))
edges <- seq(-20,15,5)
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
  group_by(date=as.Date(datetime),j) %>%
  summarise(measured_er = sum(draw*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(measured_variance = var(measured_er)) %>%
  print()

# nee_meas_uncert %>% ggplot()+geom_histogram(aes(measured_variance))
# nee_meas_uncert %>% ggplot(aes(date,measured_variance))+geom_line()

#get the total uncertainty from measurements and gap-filling
daily_NEE_uncert <- nee_meas_uncert %>%
  left_join(nee_gf_uncert,by='date') %>%
  replace(is.na(.),0) %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance)) %>%
  mutate(NEE_CI_95 = totalUncStd*2) %>% print()

flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(date) %>%
  summarise(
    #full footprint
    NEEgCmd = sum(NEE_xgbGF_HQ)*1800/1000000*44.01*0.272895) %>%
  left_join(daily_NEE_uncert %>% select(date,NEE_CI_95),by='date') %>%
  ggplot(aes(date,NEEgCmd))+
  geom_errorbar(aes(ymin=NEEgCmd-NEE_CI_95,ymax=NEEgCmd+NEE_CI_95))+
  geom_point()+
  facet_wrap(~year(date),scales='free_x')

#### daily ER uncertainty ============================================================================

#Find ER modoel uncertainty variance of model predictions
#UNCERTAINTY COMPOUNDS AS FLUXES ARE DERIVED

#set up dataframe
er <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  mutate(NEE_meas = NEEcomb) %>%
  mutate(day_night = ifelse(Solar_elv < -15, 'night','day')) %>%
  select(datetime,day_night,NEE_meas,ER,starts_with('ER_model')) %>%
  #get the residuals between the meas NEE and averaged predicted value
  mutate(resid = ER_model_median - NEE_meas)



#uncertainty from ER model predictions. The variance 20 predictions
er_model_uncert <- er %>%
  filter(day_night == 'day') %>%
  pivot_longer(cols=c(contains('model'))) %>%
  filter(!str_detect(name,'me')) %>%
  group_by(date=as.Date(datetime),name) %>%
  summarise(resp_model = sum(value*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(resp_model_variance = var(resp_model)) %>%
  ungroup() %>% print()

daily_ER_uncert <- nee_meas_uncert %>%
  left_join(nee_gf_uncert,by='date') %>%
  left_join(er_model_uncert,by='date') %>%
  replace(is.na(.),0) %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance + resp_model_variance)) %>%
  mutate(ER_CI_95 = totalUncStd*2) %>% print()


flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(date) %>%
  summarise(
    #full footprint
    ERgCmd = sum(ER)*1800/1000000*44.01*0.272895) %>%
  left_join(daily_ER_uncert %>% select(date,ER_CI_95),by='date') %>%
  ggplot(aes(date,ERgCmd))+
  geom_errorbar(aes(ymin=ERgCmd-ER_CI_95,ymax=ERgCmd+ER_CI_95))+
  geom_point()+
  facet_wrap(~year(date),scales='free_x')

#
#### daily GPP uncertainty ============================================================================
gpp <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  mutate(NEE_meas = NEEcomb) %>%
  mutate(day_night = ifelse(Solar_elv < -15, 'night','day')) %>%
  select(datetime,day_night,NEE_meas,NEE_xgbGF_HQ,ER,GPP,starts_with('ER_model'))

gpp_split_uncert <-  gpp %>%
  filter(day_night == 'day') %>%
  pivot_longer(cols=c(contains('ER_model'))) %>%
  filter(!str_detect(name,'me')) %>%
  mutate(pseudo_gpp = value-NEE_xgbGF_HQ) %>%
  group_by(date=as.Date(datetime),name) %>%
  summarise(split_er = sum(pseudo_gpp*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(split_variance = var(split_er)) %>%
  print()

#combining uncertainties
daily_total_uncert <-nee_meas_uncert %>%
  left_join(nee_gf_uncert,by='date') %>%
  left_join(er_model_uncert,by='date') %>%
  left_join(gpp_split_uncert,by='date') %>%
  replace(is.na(.),0) %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance + resp_model_variance + split_variance)) %>%
  mutate(GPP_CI_95 = totalUncStd*2) %>% print()

daily_total_uncert %>% summary()

flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(date) %>%
  summarise(GPPgCmd = sum(GPP)*1800/1000000*44.01*0.272895) %>%
  left_join(daily_total_uncert %>% select(date,GPP_CI_95),by='date') %>%
  ggplot(aes(date,GPPgCmd))+
  geom_errorbar(aes(ymin=GPPgCmd-GPP_CI_95,ymax=GPPgCmd+GPP_CI_95))+
  geom_point()+
  facet_wrap(~year(date),scales='free_x')

#
#### PUT IT ALL TOGETHER =============================================================================

#what proportion of each day is gap-filled?
prop_gf <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(date) %>%
  summarise(n_gf = sum(is.na(NEEcomb)),
            n = n()) %>%
  ungroup() %>% 
  mutate(prop_gf = n_gf/n)

#create a daily file
daily_C <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(date) %>%
  summarise(
    #full footprint
    avgNEEumolCO2s = mean(NEE_xgbGF_HQ),
    sdNEEumolCo2s = sd(NEE_xgbGF_HQ),
    avgGPPumolCO2s = mean(GPP),
    sdGPPumolCO2s = sd(GPP),
    avgERumolCO2s = mean(ER),
    sdERumolCO2s = sd(ER),
    NEEgCmd = sum(NEE_xgbGF_HQ)*1800/1000000*44.01*0.272895,
    GPPgCmd = sum(GPP)*1800/1000000*44.01*0.272895,
    ERgCmd = sum(ER)*1800/1000000*44.01*0.272895) %>%
  left_join(daily_total_uncert,by='date') %>%
  mutate(NEE_CI_95 = sqrt(measured_variance+gapfill_variance)*2,
         ER_CI_95 = sqrt(measured_variance+gapfill_variance+resp_model_variance)*2) %>%
  left_join(prop_gf %>% select(date,prop_gf),by='date') %>%
  select(date,prop_gf,starts_with('avg'),starts_with('sd'),starts_with('NEE'),starts_with('ER'),starts_with('GPP'))


#smoothed plot of fluxes
daily_C %>%
  select(-contains('umol'),-contains('CI_95')) %>%
  pivot_longer(cols=c(contains('gCmd')),names_to = 'meas',values_to = 'dailyC') %>%
  mutate(meas = str_sub(meas,end=-5)) %>%
  left_join(daily_C %>% select(-contains('umol'),-contains('gCmd')) %>% 
              pivot_longer(cols=c(contains('CI_95')),names_to='meas',values_to='uncert') %>%
              mutate(meas = str_sub(meas,end=-7)),
            by=c('date','meas')) %>% 
  arrange(date) %>%
  group_by(meas) %>%
  mutate(smoothed_C = rollapply(dailyC,width=14,align='center',FUN=mean,partial=T),
         smoothed_uncert = rollapply(uncert,width=14,align='center',FUN=mean,partial=T)) %>%
  ungroup() %>%
  ggplot(aes(date,dailyC,color=meas))+
  geom_hline(yintercept = 0, linetype = 'dotted')+
  geom_point(alpha = 0.15,show.legend=F)+
  geom_errorbar(aes(ymin=dailyC-uncert,ymax=dailyC+uncert),alpha=0.15)+
  geom_ribbon(aes(ymin = smoothed_C-smoothed_uncert, ymax=smoothed_C+smoothed_uncert,fill=meas),
              alpha=0.5,show.legend=F,color=NA)+
  geom_line(aes(y=smoothed_C,color = meas),linewidth=1)+
  annotate(geom='rect',xmin=as.Date('2018-01-01'),xmax=as.Date('2018-12-31'),
           ymin=-4,ymax=7,fill='white',alpha=0.65)+
  scale_color_manual(values = c('darkorange3','darkgreen','dodgerblue'))+
  scale_x_date(date_breaks = '1 year',date_labels = '%Y')+
  labs(x=NULL, y=bquote(g~C~m^-2~day^-1),color=NULL)+
       # caption='Lines and shading are 14-day rolling means ± 95% CI')+
  theme(legend.position=c(0.05,0.1),
        legend.background = element_blank(),
        plot.background = element_blank())

ggsave('processed/plots/flux-daily-timeseries.png',bg='white',device='png',
       units='cm',height=12, width=30)

# #save to file
write.csv(daily_C,paste0('processed/g_GCE_Fluxes_daily_carbon_',
                                  format(range(as.Date(daily_C$date))[1],'%Y%m%d'),'_',
                                  format(range(as.Date(daily_C$date))[2],'%Y%m%d'),
                                  '.csv'),row.names = F)


