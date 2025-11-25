#Author: Peter Hawman

#Purpose:
# 1. read in processed 30min data fluxes (eddypro full output files)
# 2. join biophysical data
# 3. set QC values and date removal
# 4. conduct QC for fluxes
#   a. remove days prior and post maintenance (2 days)
#   b. EddyPro QC falgs
#   c. remove extreme values
#   d. u* filtering

#Output:
# 1. flux dataset ready for combination with flux1

#set wd
setwd('Z:/FluxTower/GCELTER')

#load packages
library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(tidyselect)
library(stringr)
library(ggplot2)
library(lubridate) # for making dates easier
library(REddyProc) # for estmating u* threshold
library(zoo)       # for rolling functions
library(fuzzyjoin) # for joining between date ranges (ustar filtering)

#set plotting theme
theme_set(theme_bw())

#### Read in data ====================================================================================================


#read in biophysical
biophysical <- read_csv(list.files('processed',pattern='a_GCE_biophysical',full.names=T),
                        col_types = cols(.default=col_double(),datetime=col_datetime(),
                                         tide_stage=col_character(),PT_sensor_used=col_character()))


#select which sensor to process
flux_sensor <- 'flux2'


#read in fluxes
#list files
(flux_list <- list.files(paste0('processed/results/',flux_sensor),pattern = 'full_output',full.names=T,recursive = T))
#function to read in flux full output files and select desired variables
flux_read <- function(x) {
  read_csv(x,skip=3,col_names=as.character(read.csv(x,nrows=1,skip=1,header=F,stringsAsFactors = F))) %>%
    mutate(datetime = as.POSIXct(paste(as.Date(date,'%m/%d/%Y'),time),tz='GMT')) %>%
    select(datetime,
           co2_mole_fraction,co2_mixing_ratio,
           h2o_mole_fraction,h2o_mixing_ratio,
           H,qc_H,rand_err_H, # sensible heat
           LE,qc_LE,rand_err_LE, # latent heat
           co2_flux,qc_co2_flux,rand_err_co2_flux, # co2 flux
           h2o_flux,qc_h2o_flux,rand_err_h2o_flux, # water vapor flux
           H_strg,LE_strg,co2_strg,h2o_strg, # flux storage 
           ET, # evapotranspiration 
           wind_speed,wind_dir,`u*`,L,`(z-d)/L`,v_var, # extra variables for footprint analysis
           x_peak:`x_90%`, #EddyPro footprint predictions
           RH,VPD_flux = VPD)
}

#open each flux full ouput and put in a list
flux_data_list <- lapply(flux_list,flux_read)
#bind all listed dataframes into one
flux <- bind_rows(flux_data_list) %>%
  arrange(datetime)

#replace with NA
flux[flux==-9999] <- NA

#create NEE variable as co2 flux + co2 storage and LE and H
flux <- flux %>%
  mutate(NEE = co2_flux + co2_strg,
         LE = LE + LE_strg,
         H = H + H_strg)

# flux %>%
#   pivot_longer(cols=c(NEE,LE,H)) %>%
#   ggplot(aes(datetime,value))+
#   geom_line()+
#   facet_wrap(~name,scales='free_y')

#create a vector of years for this analysis. This will be used to filter the data
years <- unique(year(flux$datetime))

#### join biophysical data =====================================================================================

#join biophysical data and fluxes to specific date range
full_data <- tibble(datetime = seq.POSIXt(from = range(biophysical$datetime)[1],
                                          to = range(biophysical$datetime)[2],
                                          by = '30 mins')) %>%
  left_join(biophysical,by='datetime') %>%
  left_join(flux,by='datetime') %>%
  #filter to years specified above
  filter(year(datetime) %in% years)

length(seq.POSIXt(from=full_data$datetime[1],to=full_data$datetime[nrow(full_data)],by='30 min'))
nrow(full_data)

# full_data %>%
#   ggplot(aes(datetime,NEE))+
#   geom_point()+
#   ylim(c(-100,100))


#### Data removal near maintenance ===============================================================================

#make a table of maintenance data based on data file splits
maintenance <- tibble(file = flux_list) %>%
  mutate(start_date = as.Date(str_split_i(file,'-',i=2),format='%Y%m%d'),
         end_date = as.Date(str_split_i(file,'-',i=3),format='%Y%m%d'))

maint_list <- list()
for (i in 1:length(maintenance$start_date)){
  
  x <- tibble(date = seq.Date(from = maintenance$start_date[i], to = maintenance$end_date[i],by = '1 day'))
  
  y <- x %>%
    mutate(days_prior_maintenance = nrow(x)-row_number(),
           days_since_maintenance = row_number()-1)
  
  maint_list[[i]] <- y
}

maint_dates <- bind_rows(maint_list)

#join to flux data and set data 2 days prior to maintenance to NA
full_data_maint <- full_data %>%
  mutate(date = as.Date(datetime)) %>%
  left_join(maint_dates, by='date',relationship = 'many-to-many') %>%
  select(-date) %>% 
  mutate(days_prior_maintenance = ifelse(is.na(days_prior_maintenance),0,days_prior_maintenance),
         days_since_maintenance = ifelse(is.na(days_since_maintenance),0,days_since_maintenance)) %>%
  mutate(NEE_maint = ifelse(days_prior_maintenance < 2, NA, NEE),
         LE_maint = ifelse(days_prior_maintenance < 2, NA, LE),
         H_maint = ifelse(days_prior_maintenance < 2, NA, H))

#### QC Flags ========================================================================================================

qc_value <- c(0,1)

#plot data by QC value

# full_data_maint %>% 
#   ggplot(aes(x=datetime))+
#   geom_point(aes(y=NEE_maint,color=factor(qc_co2_flux)))+
#   scale_x_datetime(date_breaks='1 year',date_labels = '%Y')+
#   ylim(c(-20,20))

#if CO2_qc is not acceptable, set to NA
full_data_qc <- full_data_maint %>%
  mutate(NEE_maint_qc = ifelse(qc_co2_flux %in% qc_value, NEE_maint, NA))  %>% # set any NEE with QC != qc_value to NA
  #latent and sensible heat
  mutate(LE_maint_qc = ifelse(qc_LE %in% qc_value, LE_maint, NA),
         H_maint_qc = ifelse(qc_H %in% qc_value, H_maint, NA)) 
  #remove any extreme values

# full_data_qc %>% 
#   ggplot(aes(x=datetime))+
#   geom_point(aes(y=NEE_maint),color='red')+
#   geom_point(aes(y=NEE_maint_qc))+
#   scale_x_datetime(date_breaks='1 year',date_labels = '%Y')+
#   ylim(c(-20,20))

#### Data removal (visual/co2 mixing ratio) =====================================================================================

#select table containing dates to remove based on comparison and/or CO2 mixing ratio
removal_dates <- read_csv('metadata/by_date_removal_co2mixing-comparison-looksbad.csv',
                          col_types = cols(start_date=col_date(format = '%m/%d/%Y'),
                                           end_date=col_date(format = '%m/%d/%Y'))) %>%
  filter(tower == flux_sensor) # filter to current flux sensor

#read in sections of data that need to be removed based on mixing ratio and comparison by David

removal_dates

# full_data_qc %>%
#   # filter(year(datetime) %in% c(2017)) %>%
#   ggplot()+
#   geom_rect(data=removal_dates,
#             aes(xmin=as.POSIXct(paste(start_date,'00:00',format='%Y-%m-%d %H:%M')),
#                 xmax = as.POSIXct(paste(end_date,'00:00',format='%Y-%m-%d %H:%M')),
#                 ymin=-Inf,ymax=Inf),
#             fill='red',alpha=0.3)+
#   geom_point(aes(x=datetime,y=NEE_maint_qc))+
#   # lims(y=c(-20,10))+
#   scale_x_datetime(date_breaks = '1 month',date_label='%Y-%m')+
#   theme(axis.text.x = element_text(angle=45,hjust=1))
# 
# 


#join date ranges with bad data and set to NA
full_data_qc_removal <- full_data_qc %>%
  mutate(date = as.Date(datetime)) %>%
  fuzzy_left_join(removal_dates,
                by=c('date'='start_date','date'='end_date'),
                match_fun=list(`>=`,`<=`)) %>% 
  mutate(NEE_maint_qc_removal = if_else(!is.na(start_date),NA,NEE_maint_qc)) %>%
  mutate(LE_maint_qc_removal = if_else(!is.na(start_date),NA,LE_maint_qc)) %>%
  mutate(H_maint_qc_removal = if_else(!is.na(start_date),NA,H_maint_qc))



#plot
# full_data_qc_removal %>%
#   # filter(as.Date(datetime)==removal_dates$start_date[2]+1) %>% pull(NEE_maint_qc_removal)
#   filter(year(datetime)==2024) %>%
#   ggplot(aes(datetime,NEE_maint_qc))+
#   geom_point(color='red')+
#   geom_point(aes(y=NEE_maint_qc_removal))+
#   scale_x_datetime(date_breaks = '1 year',date_labels='%Y')


#### extreme value filtering ============================================================================================

#remove extreme values based on NEE
full_data_qc_removal_ext <- full_data_qc_removal %>%
  mutate(NEE_maint_qc_removal_ext = ifelse(between(NEE_maint_qc_removal,-20,15),NEE_maint_qc_removal,NA),
         LE_maint_qc_removal_ext = ifelse(between(NEE_maint_qc_removal,-20,15),LE_maint_qc_removal,NA),
         H_maint_qc_removal_ext = ifelse(between(NEE_maint_qc_removal,-20,15),H_maint_qc_removal,NA))

full_data_qc_removal_ext %>%
  ggplot(aes(x=datetime))+
  geom_point(aes(y=NEE_maint_qc_removal),color='red',size=0.25)+
  geom_point(aes(y=NEE_maint_qc_removal_ext),color='black',size=0.25)+
  geom_hline(yintercept = -20,linetype='dotted')+
  geom_hline(yintercept = 15,linetype='dotted')+
  scale_x_datetime(date_breaks='1 year',date_labels = '%Y')


#### u* filtering ========================================================================================================
# u* threshold determination and filtering
# vignette("uStarCases")
#get data ready for REddyProc. Needs DateTime, NEE, Rg, Tair, VPD, Ustar in df
fluxes_reddy_prior_2018 <- full_data_qc_removal_ext %>%
  mutate(Mean_Shortwave_Net = Mean_PAR_Incident_Tower/2) %>%
  select(DateTime=datetime,NEE=NEE_maint_qc_removal_ext,Rg=Mean_Shortwave_Net,Tair=Mean_Temp_Air,VPD,Ustar=`u*`) %>%
  filter(year(DateTime)<2018) %>%
  as.data.frame()

fluxes_reddy_prior_2018 %>% mutate(gap = DateTime - lag(DateTime)) %>% as_tibble() %>% arrange(gap) %>% filter(gap == 0)

#convert to whatever object REddyProc likes
fluxes_reddyproc_prior_2018 <- sEddyProc$new('gce',fluxes_reddy_prior_2018,c('NEE','Rg','Tair','VPD','Ustar'))

#ustar threshold determination with bootstrapping
ustarthres_prior2018 <- fluxes_reddyproc_prior_2018$sEstUstarThold()


#do same for post 2018
fluxes_reddy_post_2018 <- full_data_qc_removal_ext %>%
  mutate(Mean_Shortwave_Net = Mean_PAR_Incident_Tower/2) %>%
  select(DateTime=datetime,NEE=NEE_maint_qc_removal_ext,Rg=Mean_Shortwave_Net,Tair=Mean_Temp_Air,VPD,Ustar=`u*`) %>%
  filter(year(DateTime)>2018) %>%
  as.data.frame()

#convert to whatever object REddyProc likes
fluxes_reddyproc_post_2018 <- sEddyProc$new('gce',fluxes_reddy_post_2018,c('NEE','Rg','Tair','VPD','Ustar'))

#ustar threshold determination with bootstrapping
ustarthres_post2018 <- fluxes_reddyproc_post_2018$sEstUstarThold()

#combine two sets of years around 2018
ustar_seasonal_values <- bind_rows(ustarthres_prior2018 %>%
                                     as_tibble() %>%
                                     drop_na(),
                                   ustarthres_post2018 %>%
                                     as_tibble() %>%
                                     drop_na()) %>%

  mutate(year = str_sub(season,end=4),
         month = str_sub(season,start=6)) %>%
  mutate(start_date = as.Date(paste(year,month,'01'),format='%Y %m %d')) %>%
  mutate(end_date = lead(start_date)) %>%
  select(start_date,end_date,uStar)


#remove data if less than u* threshold
full_data_qc_removal_ext_ustar <- full_data_qc_removal_ext %>%
  mutate(date = as.Date(datetime)) %>%
  fuzzy_left_join(ustar_seasonal_values,
                  by=c('date'='start_date','date'='end_date'),
                  match_fun=list(`>=`,`<`)) %>%
  mutate(NEE_maint_qc_removal_ext_ustar = ifelse(`u*`<uStar,NA,NEE_maint_qc_removal_ext),
         LE_maint_qc_removal_ext_ustar = ifelse(`u*`<uStar,NA,LE_maint_qc_removal_ext),
         H_maint_qc_removal_ext_ustar = ifelse(`u*`<uStar,NA,H_maint_qc_removal_ext)) 


full_data_qc_removal_ext_ustar %>%
  ggplot(aes(x=datetime))+
  geom_point(aes(y=NEE_maint_qc_removal_ext),color='red')+
  geom_point(aes(y=NEE_maint_qc_removal_ext_ustar),color='black')+
  scale_x_datetime(date_breaks='1 year',date_labels = '%Y')

#### export file =====================================================================================

#check column names
colnames(full_data_qc_removal_ext_ustar)

#rename NEE to org and remove unwanted columns
out <- full_data_qc_removal_ext_ustar %>%
  rename(NEE_org = NEE,
         end_datetime=datetime) %>%
  mutate(start_datetime = end_datetime-1800) %>%
  #reorder a bit and remove biophysical data
  select(start_datetime,end_datetime,
         ET:days_since_maintenance,
         co2_mole_fraction:h2o_mixing_ratio,
         co2_flux,qc_co2_flux,rand_err_co2_flux,co2_strg,starts_with('NEE'),
         h2o_flux,qc_h2o_flux,rand_err_h2o_flux,h2o_strg,
         LE,qc_LE,rand_err_LE,LE_strg,starts_with('LE_'),
         H,qc_H,rand_err_H,H_strg,starts_with('H_'))

colnames(out)

#final plot before writing
out %>%
  pivot_longer(cols=c(contains('_ext'))) %>%
  ggplot(aes(end_datetime,value))+
  geom_point()+
  facet_wrap(~name,scales='free_y',ncol=1)

out %>% count(end_datetime,sort=T)

write.csv(out,paste0('processed/b_NEW_',
                     flux_sensor,'_GCE_Fluxes_qcd_',
                     format(range(as.Date(out$end_datetime))[1],'%Y%m%d'),'_',
                     format(range(as.Date(out$end_datetime))[2],'%Y%m%d'),
                     '.csv'),row.names = F)


