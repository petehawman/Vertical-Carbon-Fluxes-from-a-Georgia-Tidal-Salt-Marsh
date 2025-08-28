#Author: Peter Hawman

#Purpose:
# 1. calculate interannual totals and uncertainty for
#    NEE, GPP, and ER
# 2. select to budegt:
#     a. full footprint, 
#     b. creek
#     c. interior

#Output:
# 1. interannual budgets

#set wd
setwd('Z:/FluxTower/GCELTER')

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(lubridate)
library(stringr)
library(ggplot2)
library(VGAM) #for random draws from Laplace distribution

theme_set(ggthemes::theme_base())


#read in data
flux <- read_csv(list.files('processed',pattern='e_GCE',full.names=T))

### annual NEE with uncertainty =================================================================================================


#set up dataframe
nee <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  mutate(NEE_meas = NEEcomb) %>%
  select(datetime,date,qc_co2_flux,NEE_meas,NEE_xgbGF_HQ,starts_with('model')) %>%
  #get the residuals between the meas NEE and averaged predicted value
  mutate(resid = model_median - NEE_meas)

#gap-filling uncertaintity
nee_gf_uncert <-  nee %>%
  filter(is.na(NEE_meas)) %>%
  pivot_longer(cols=c(contains('model'))) %>%
  filter(!str_detect(name,'me')) %>%
  group_by(year=year(datetime),name) %>%
  summarise(gapfill_er = sum(value*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(gapfill_variance = var(gapfill_er)) %>%
  ungroup()

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
  group_by(year=year(datetime),j) %>%
  summarise(meas_er = sum(draw*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>% 
  group_by(year) %>%
  summarise(measured_variance = var(meas_er)) %>%
  print()



#### year ER uncertainty ============================================================================

#Find ER modoel uncertainty variance of model predictions
#UNCERTAINTY COMPOUNDS AS FLUXES ARE DERIVED

#set up dataframe
er <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
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
  group_by(year=year(datetime),name) %>%
  summarise(resp_model = sum(value*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(resp_model_variance = var(resp_model)) %>%
  ungroup()


#
#### year GPP uncertainty ============================================================================
gpp <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  mutate(date = as.Date(datetime)) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  mutate(NEE_meas = NEEcomb) %>%
  mutate(day_night = ifelse(Solar_elv < -15, 'night','day')) %>%
  select(datetime,day_night,NEE_meas,NEE_xgbGF_HQ,ER,GPP,starts_with('ER_model'))

gpp_split_uncert <-  gpp %>%
  pivot_longer(cols=c(contains('ER_model'))) %>%
  filter(!str_detect(name,'me')) %>%
  mutate(pseudo_gpp = value-NEE_xgbGF_HQ) %>%
  group_by(year=year(datetime),name) %>%
  summarise(split_er = sum(pseudo_gpp*1800/1000000*44.01*0.272895)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(split_variance = var(split_er))

#combining uncertainties
year_total_uncert <- nee_gf_uncert %>%
  left_join(nee_meas_uncert,by='year') %>%
  left_join(er_model_uncert,by='year') %>%
  left_join(gpp_split_uncert,by='year') %>%
  mutate(totalUncStd = sqrt(gapfill_variance  + measured_variance + resp_model_variance + split_variance)) %>%
  mutate(GPP_CI_95 = totalUncStd*2) %>% print()


year_total_uncert %>%
  filter(year != 2018) %>%
  select(year,measured_variance, gapfill_variance,resp_model_variance,split_variance) %>%
  pivot_longer(cols=c(-year)) %>%
  group_by(year) %>%
  mutate(cumulative = cumsum(value)) %>%
  mutate(normalized = cumulative/max(cumulative)) %>%
  ungroup() %>%
  mutate(name = factor(name,levels=c('measured_variance', 'gapfill_variance','resp_model_variance','split_variance'),
                       labels=c('random measurment\nuncertainty', 'gap-filling model\nuncertainty',
                                'ER model\nuncertainty','partitioning\nuncertainty'))) %>%
  group_by(name) %>%
  summarise(mean = mean(normalized),
            sd = sd(normalized)) %>%
  # filter(year == 2022) %>%
  ggplot(aes(name,mean,group=1))+
  geom_ribbon(aes(ymin=0,ymax=mean),fill='grey70')+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),fill=NA,color='black',linetype='dashed')+
  geom_line(linewidth=1)+
  scale_y_continuous(labels=scales::percent)+
  labs(x=NULL,y='cumulative uncertainty')


year_total_uncert %>% summary()

write.csv(year_total_uncert,'processed/i_GCE_annual_uncertanty-sources.csv',row.names=F)

#
#### PUT IT ALL TOGETHER =============================================================================

#what proportion of each day is gap-filled?
prop_gf <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(year=year(datetime)) %>%
  summarise(n_gf = sum(is.na(NEEcomb)),
            n = n()) %>%
  ungroup() %>% 
  mutate(prop_gf = n_gf/n) %>%
  print()

year_C <- flux %>%
  mutate(datetime = start_datetime-5*3600) %>%
  filter(year(datetime) %in% c(2014:2024)) %>%
  group_by(year=year(datetime)) %>%
  summarise(
    NEEgCmd = sum(NEE_xgbGF_HQ)*1800/1000000*44.01*0.272895,
    GPPgCmd = sum(GPP)*1800/1000000*44.01*0.272895,
    ERgCmd = sum(ER)*1800/1000000*44.01*0.272895) %>%
  left_join(year_total_uncert,by='year') %>% 
  mutate(NEE_CI_95 = sqrt(measured_variance+gapfill_variance)*2,
         ER_CI_95 = sqrt(measured_variance+gapfill_variance+resp_model_variance)*2) %>%
  left_join(prop_gf %>% select(year,prop_gf),by='year') %>%
  select(year,prop_gf,starts_with('NEE'),starts_with('ER'),starts_with('GPP')) %>%
  print()

### plot and save to file ===============================================================================================



plot_budgets <- year_C %>%
  select(-contains('CI_95')) %>%
  pivot_longer(cols=c(contains('gCmd')),names_to = 'meas',values_to = 'budget') %>%
  mutate(meas = str_sub(meas,end=-5)) %>%
  left_join(year_C %>% select(-contains('gCmd')) %>% 
              pivot_longer(cols=c(contains('CI_95')),names_to='meas',values_to='CI_95') %>%
              mutate(meas = str_sub(meas,end=-7)),
            by=c('year','meas'))

plot_budgets %>%
  ggplot(aes(year,budget,fill=meas,color=meas))+
  geom_bar(stat='identity',position = 'dodge',alpha=0.5)+
  geom_errorbar(aes(ymin=budget-CI_95,ymax=budget+CI_95),width=0.5,show.legend=F,position=position_dodge(0.95))+
  geom_text(data = plot_budgets %>% filter(meas == 'NEE'),
            aes(label=round(budget)),position=position_dodge(width=0.9),show.legend=F,vjust=3,hjust=-0.3,fontface='bold')+
  geom_text(data = plot_budgets %>% filter(meas == 'ER'),
            aes(label=round(budget)),position='dodge',show.legend=F,vjust=-2,hjust=1.7,fontface='bold')+
  geom_text(data = plot_budgets %>% filter(meas == 'GPP'),
            aes(label=round(budget)),position='dodge',show.legend=F,vjust=-2,hjust=0.5,fontface='bold')+
  scale_fill_manual(values = c('darkorange3','darkgreen','dodgerblue'))+
  scale_color_manual(values = c('darkorange3','darkgreen','dodgerblue'))+
  scale_x_continuous(breaks = seq(2014,2024,1))+
  annotate(geom='rect',xmin=2017.5,xmax=2018.5,
           ymin=-Inf,ymax=Inf,fill='white',alpha=0.65)+
  geom_hline(yintercept = 0,linetype='dotted')+
  ylim(c(-400,1200))+
  labs(x=NULL,y=bquote(g~C~m^-2~y^-1),color=NULL,fill=NULL)+#,caption='Error bars are ± 95% CI')+
  theme(legend.position = c(0.05,0.9),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor.x = element_line(color='grey30',linewidth = 0.25,linetype='dashed'))

ggsave('processed/plots/flux-annual.png',bg='white',device='png',
       units='cm',height=12, width=30)


#save to file
write.csv(year_C,paste0('processed/h_GCE_annual_carbon_flux_',
                         min(year_C$year),'_',
                         max(year_C$year),
                         '.csv'),row.names = F)
