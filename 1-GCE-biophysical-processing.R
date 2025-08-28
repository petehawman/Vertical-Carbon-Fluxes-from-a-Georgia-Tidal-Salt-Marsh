#Author: Peter Hawman

#Purpose:
# 1. reads in 5-min combined variables GCE dataset
# 2. QC data based on flags
# 3. Additional eyeball QC
# 4. Create 30-min dataset
# 5. Adds solar elevation and top-of-atmosphere PAR
# 6. Fixes calibration issues with PAR and removes some bad data
# 7. Fills data gaps in air temp, humidity, and PAR with Marsh Landing data
# 8. Corrects high humidity readings and calcs VPD
# 9. Adds salinity data from Lower Duplin
# 10. Create fuzzy variables for time of day and season
# 11. Exports 30-min biophysical data file


#libraries needed
library(tidyverse) #for the best data
library(lubridate) #for date and times!
library(zoo) #for rolling means
library(suntools) #for calculating Solar Elevation
library(photobiology) # for find_peaks function
library(multidplyr) #multi-core processing

#set clusters for multidplyr
cluster <- new_cluster(7)

#### Read in 5min combined dataset =====================================================================

#this gets the column names
cnames <- as.character(read.csv('biophysical-data/combined-5min/Combined_5min_20130709-20210809.csv',skip=2,nrows=1,header=F,stringsAsFactors=F))
#read in data
biophysical_1 <- read_csv('biophysical-data/combined-5min/Combined_5min_20130709-20210809.csv',
              #assign column names
              col_names = cnames,
              #skip header info
              skip=5,
              #set column types
              col_types = do.call(paste0,as.list(ifelse(cnames=='Date','T',ifelse(cnames %in% str_subset(cnames,'Flag'),'c','d')))))

biophysical_2 <- read_csv('biophysical-data/combined-5min/Combined_FluxTower_5min_20200401-20211231.csv',
                          #assign column names
                          col_names = cnames,
                          #skip header info
                          skip=5,
                          #set column types
                          col_types = do.call(paste0,as.list(ifelse(cnames=='Date','T',ifelse(cnames %in% str_subset(cnames,'Flag'),'c','d')))))

biophysical_3 <- read_csv('biophysical-data/combined-5min/Combined_FluxTower_5min_20220124-20250225.csv',
                          #assign column names
                          col_names = cnames,
                          #skip header info
                          skip=5,
                          #set column types
                          col_types = do.call(paste0,as.list(ifelse(cnames=='Date','T',ifelse(cnames %in% str_subset(cnames,'Flag'),'c','d')))))

#make a seqence of 5 min timestamps through range of data so we have continous timestamps
date_df <- tibble(Date = seq.POSIXt(range(biophysical_1$Date)[1],range(biophysical_3$Date)[2],by='5 min'))

#bind and join to all timestamps
biophysical <- date_df %>%
  left_join(bind_rows(biophysical_1,
          biophysical_2 %>% filter(!Date %in% biophysical_1$Date),
          biophysical_3),by='Date')

biophysical %>%
  count(Date) %>% arrange(desc(n))

range(biophysical$Date)

tibble(Date = seq.POSIXt(from=range(biophysical$Date)[1],to=range(biophysical$Date)[2],by='5 min')) %>%
  # filter(!Date %in% biophysical$Date) %>%
  nrow()

nrow(biophysical)

#### Use Flags to QC data ==============================================================================

#get flags in long format
flags <- biophysical %>%
  pivot_longer(cols = starts_with('Flag'),names_to='flags',values_to = 'flag_values') %>%
  select(Date,flags,flag_values) 
  
#get measurements in long format
measures <- biophysical %>%
  pivot_longer(cols= c(-starts_with('Flag'),-Date), names_to='meas',values_to = 'meas_values') %>%
  select(Date,meas,meas_values)

#bind columns of flags and measurements
biophsyical_flagged <- bind_cols(flags,measures) %>% 
  #fix duplicate date column
  rename(Date=Date...1) %>%
  select(-Date...4) %>%
  #if there is a flag (E or Q) set measure to NA
  mutate(meas_values = ifelse(is.na(flag_values),meas_values,NA)) %>%
  #remove flag columns
  select(-flags,-flag_values) %>%
  #spread measures back out
  pivot_wider(names_from=meas,values_from=meas_values)

summary(biophsyical_flagged)


##################### Remove bad data not captured from QC ##############################################

#variable names
biophsyical_flagged %>% colnames()

#plot soil temperatures
biophsyical_flagged %>%
  filter(year(Date) %in% c(2021:2023)) %>%
  # filter(month(Date)==4) %>%
  select(-contains('Delta')) %>%
  pivot_longer(cols=c(contains('Soil'))) %>%
  ggplot(aes(Date,value,color=name,group=name))+
  geom_line()+
  geom_vline(xintercept = as.POSIXct('2021-05-10 00:00',tz='GMT'))


#remove odd soil Temp
biophsyical_flagged_ST <- biophsyical_flagged %>%
  mutate(Mean_Temp_Soil_Marsh = ifelse(between(Date,as.POSIXct('2020-04-17 12:00',tz='GMT'),as.POSIXct('2020-04-19 20:00',tz='GMT')),
                                               NA, Mean_Temp_Soil_Marsh)) %>%
  mutate(Mean_Temp_Soil_Marsh = ifelse(between(Date,as.POSIXct('2020-07-24 20:00',tz='GMT'),as.POSIXct('2020-07-25 01:00',tz='GMT')),
                                       NA, Mean_Temp_Soil_Marsh)) %>%
  mutate(Mean_Temp_Soil_Marsh = ifelse(between(Date,as.POSIXct('2020-08-14 17:00',tz='GMT'),as.POSIXct('2020-08-14 22:00',tz='GMT')),
                                       NA, Mean_Temp_Soil_Marsh)) %>%
  mutate(Mean_Temp_Soil_Marsh = ifelse(between(Date,as.POSIXct('2020-11-15 12:00',tz='GMT'),as.POSIXct('2020-11-15 17:00',tz='GMT')),
                                       NA, Mean_Temp_Soil_Marsh)) %>%
  mutate(Mean_Temp_Soil_Tower = ifelse(between(Date,as.POSIXct('2021-05-10 00:00',tz='GMT'),as.POSIXct('2022-02-22 00:00',tz='GMT')),
                                       NA, Mean_Temp_Soil_Tower)) %>%
  mutate(Mean_Temp_Soil_Marsh = ifelse(between(Date,as.POSIXct('2021-05-10 00:00',tz='GMT'),as.POSIXct('2022-01-01',tz='GMT')), 
                                       NA, Mean_Temp_Soil_Marsh)) 

  # filter(year(Date) %in% c(2021)) %>%
  # filter(month(Date)==c(11)) %>%
  # select(-contains('Delta')) %>%
  # pivot_longer(cols=c(contains('Soil'))) %>%
  # ggplot(aes(Date,value,color=name,group=name))+
  # geom_line()
  # geom_vline(xintercept = as.POSIXct('2022-02-22 00:00',tz='GMT'))
  # geom_vline(xintercept = as.POSIXct('2020-11-15 17:00',tz='GMT'))


#plot air temperatures
# biophsyical_flagged_ST %>%
#   pivot_longer(cols=c(contains('Air'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()+
#   facet_wrap(~name)

#plot humidity
# biophsyical_flagged_ST %>%
#   pivot_longer(cols=c(contains('Humidity'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()

  
#### tides and water table ==============================================================================================

biophsyical_flagged_ST %>%
  # filter(as.Date(Date)==as.Date('2019-05-14')) %>%
  # filter(month(Date)==05 & year(Date)==2019) %>%
  filter(year(Date) %in% c(2017:2018)) %>%
  pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh)) %>%
  drop_na(value) %>%
  ggplot(aes(Date,value,color=name))+
  geom_line()

#remove bad sections of data
biophsyical_flagged_ST_water <- biophsyical_flagged_ST %>%
  
  #marsh water table
  mutate(Mean_Water_Level_Marsh = ifelse(Date < as.POSIXct('2017-09-12 00:00',tz='GMT'),
                                         Mean_Water_Level_Marsh+0.0331,Mean_Water_Level_Marsh)) %>%
  mutate(Mean_Water_Level_Marsh = ifelse(between(Date, left = as.POSIXct('2017-09-11 00:00',tz='GMT'),right = as.POSIXct('2017-09-12 00:00',tz='GMT')),
                                         NA,Mean_Water_Level_Marsh)) %>%
  
  
  mutate(Mean_Water_Level_Marsh = ifelse(between(Date, left = as.POSIXct('2017-09-12 00:00',tz='GMT'),right = as.POSIXct('2017-12-20 19:10',tz='GMT')),
                                         Mean_Water_Level_Marsh+0.0825,Mean_Water_Level_Marsh)) %>%
  mutate(Mean_Water_Level_Marsh = ifelse(between(Date, left = as.POSIXct('2019-08-04 17:30',tz='GMT'),right = as.POSIXct('2019-12-20 23:30',tz='GMT')),
                                         NA,Mean_Water_Level_Marsh)) %>%
  mutate(Mean_Water_Level_Marsh = ifelse(between(as.Date(Date),left=as.Date('2019-06-24'),right = as.Date('2019-06-25')),
                                         NA,Mean_Water_Level_Marsh)) %>%
  mutate(Mean_Water_Level_Marsh = ifelse(between(as.Date(Date),left=as.Date('2021-07-01'),right = as.Date('2022-01-23')),
                                         NA,Mean_Water_Level_Marsh)) %>%
  
  
  #creek water table
  mutate(Mean_Water_Level_Creek = ifelse(as.Date(Date) < as.Date('2013-11-01'),
                                         NA,Mean_Water_Level_Creek)) %>%
  mutate(Mean_Water_Level_Creek = ifelse(between(as.Date(Date), left = as.Date('2017-07-01'),right = as.Date('2017-12-20')),
                                         NA,Mean_Water_Level_Creek)) %>%
  mutate(Mean_Water_Level_Creek = ifelse(year(Date)==2017,NA,Mean_Water_Level_Creek)) %>%
  mutate(Mean_Water_Level_Creek = ifelse(as.Date(Date) >= as.Date('2021-05-09'),NA,Mean_Water_Level_Creek)) %>%
  mutate(Mean_Water_Level_Creek = ifelse(between(Date,as.POSIXct('2018-03-20 19:00',tz='GMT'),as.POSIXct('2018-03-21 00:20',tz='GMT')),
                                         NA,Mean_Water_Level_Creek)) %>%
  mutate(Mean_Water_Level_Creek = ifelse(between(Date,as.POSIXct('2018-07-30 21:00',tz='GMT'),as.POSIXct('2018-07-31 00:00',tz='GMT')),
                                         NA,Mean_Water_Level_Creek)) %>%
  mutate(Mean_Water_Level_Creek = ifelse(between(Date,as.POSIXct('2018-07-30 17:40',tz='GMT'),as.POSIXct('2018-07-30 19:00',tz='GMT')),
                                         NA,Mean_Water_Level_Creek)) 
  

# crap_dates <- c(as.Date('2018-03-20'), as.Date('2018-07-30'),as.Date('2019-12-30'),as.Date('2021-05-06'))

biophsyical_flagged_ST_water %>%
  # filter(as.Date(Date) %in% crap_dates) %>%
  filter(as.Date(Date) %in% c(as.Date('2018-07-30'),as.Date('2018-07-31'))) %>%
  # filter(year(Date) %in% c(2021,2022)) %>%
  pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh)) %>%
  ggplot(aes(Date,value,color=name))+
  geom_line()+
  geom_point()+
  scale_x_datetime(date_breaks = '1 hour',date_labels='%H')


biophsyical_flagged_ST_water %>%
  # filter(as.Date(Date)==as.Date('2019-05-14')) %>%
  # filter(month(Date) %in% c(9)) %>%
  filter(year(Date) %in% c(2017:2021)) %>%
  # mutate(Mean_Water_Level_Marsh = ifelse(between(Date, left = as.POSIXct('2017-09-12 00:00',tz='GMT'),right = as.POSIXct('2017-12-20 19:10',tz='GMT')),
  #                                        Mean_Water_Level_Marsh+0.0825,Mean_Water_Level_Marsh)) %>%
  # filter(day(Date) %in% c(10:12)) %>% #select(Date,Mean_Water_Level_Marsh) %>% view()
  pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh)) %>%
  ggplot(aes(Date,value,color=name))+
  geom_line()+
  geom_vline(xintercept = as.POSIXct('2017-09-12 00:00',tz='GMT'))+
  geom_vline(xintercept = as.POSIXct('2017-09-11 00:00',tz='GMT'))

#### create 30min averages (totals for precipitation) ===================================================

#cluster zoo lirbary for rollapply
cluster_library(cluster, 'zoo')

#make 30min averages
biophsyical_30min_means <- biophsyical_flagged_ST_water %>%
  pivot_longer(c(-Date,-contains(c('Total','Wind_Direction'))),names_to = 'variables',values_to = 'values') %>%
  group_by(variables) %>%
  partition(cluster) %>%
  mutate(rollingmeans = rollapply(values,width = 7, FUN=mean, na.rm=T, align='right',fill=NA, partial=T)) %>%
  collect() %>%
  ungroup() %>%
  pivot_wider(id_cols = Date,names_from=variables,values_from=rollingmeans) %>%
  filter(minute(Date) == 00 | minute(Date) == 30)

#make 30min totals for precip
biophysical_30min_totals <- biophsyical_flagged_ST_water %>%
  select(Date,Total_Precip) %>%
  mutate(Total_Precip = rollapply(Total_Precip,width = 6, FUN=sum, align='right',fill=NA, partial=T)) %>%
  filter(minute(Date) == 00 | minute(Date) == 30)

#average wind direction. Only doing vector averaging for direction, will use scalar average from above for speed
mean_wind <- biophsyical_flagged_ST_water %>%
  select(Date,Wind_Direction,Mean_Wind_Speed) %>%
  #calculate vector components
  mutate(u = Mean_Wind_Speed*sin(pi*Wind_Direction/180),
         v = Mean_Wind_Speed*cos(pi*Wind_Direction/180)) %>%
  #30min means of vector compoonents
  mutate(u_mean = rollapply(u,width=7,FUN=mean,na.rm=T,align='right',fill=NA,partial=T),
         v_mean = rollapply(v,width=7,FUN=mean,na.rm=T,align='right',fill=NA,partial=T)) %>%
  filter(minute(Date)== 00 | minute(Date) == 30) %>%
  #average direction from 30min vector means
  mutate(Mean_Wind_Direction = (atan2(u_mean,v_mean)*180/pi)) %>% 
  #if degrees are negative, add 360
  mutate(Mean_Wind_Direction = ifelse(Mean_Wind_Direction < 0, Mean_Wind_Direction+360,Mean_Wind_Direction)) %>%
  select(Date,Mean_Wind_Direction)

# biophsyical_flagged_ST_water %>%
#   left_join(mean_wind %>% select(Date,Mean_Wind_Direction),by='Date') %>%
#   filter(year(Date)==2017 & yday(Date) %in% 43:60) %>%
#   ggplot(aes(Date,Wind_Direction))+
#   geom_point()+
#   geom_point(aes(y=Mean_Wind_Direction),color='red')

#put them together
biophsyical_30min <- biophsyical_30min_means %>%
  left_join(biophysical_30min_totals,by='Date') %>%
  left_join(mean_wind,by='Date')

# biophsyical_30min_means %>%
#   filter(as.Date(Date) %in% c(as.Date('2019-05-13'),as.Date('2019-05-14'),as.Date('2019-05-15'))) %>%
#   # filter(month(Date)==05 & year(Date)==2019) %>%
#   pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh)) %>%
#   ggplot(aes(Date,value,color=name))+
#   geom_line()

#### Plot 30min data ====================================================================================

#plot soil temperatures
# biophsyical_30min %>%
#   select(-contains('Delta')) %>%
#   pivot_longer(cols=c(contains('Soil'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()

#plot air temperatures
# biophsyical_30min %>%
#   pivot_longer(cols=c(contains('Air'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()+
#   facet_wrap(~name)

#plot humidity
# biophsyical_30min %>%
#   pivot_longer(cols=c(contains('Humidity'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()+
#   facet_wrap(~name)

#plot PAR
# biophsyical_30min %>%
#   select(-contains('Total')) %>%
#   pivot_longer(cols=c(contains('PAR'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()

#plot shortwave radiation
# biophsyical_30min %>%
#   select(-contains('Total')) %>%
#   pivot_longer(cols=c(contains('Shortwave'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()+
#   facet_wrap(~name)

#plot water tables
# biophsyical_30min %>%
#   pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh)) %>%
#   ggplot(aes(Date,value,color=name))+
#   geom_line()+
#   facet_wrap(~name)

#### Add solar elevation info and TOA PAR ========================================================================

#Use maptools package to calculate Solar Elevation angle using lat/long of flux tower
#matrix containing lat and long
geolocation <- matrix(c(-81.28233, 31.44401),nrow=1,ncol=2,dimnames = list(c('1'),c('x','y')))

#get solar positions. The first column is solar azimuth, the second is solar elevation. 
#Adjust time by 15 min to match with GCE right aligned averaging
solar.postion <- as.data.frame(solarpos(crds=geolocation,dateTime = biophsyical_30min$Date-.25*60*60,POSIXct.out=TRUE)) %>% #
  rename(azimuth=V1,elevation=V2) #rename columns

#add Solar elevation angle to dataset
biophsyical_30min$Solar_elv <- solar.postion$elevation

biophsyical_30min_Solar <- biophsyical_30min %>%
  mutate(DOY = yday(Date),
         Itoa = 1370*(1+0.033*cos(360*DOY/365))*sin(Solar_elv*pi/180),
         PARTOA = 0.48*Itoa*4.57) %>%
  #set PARTOA at night to 0
  mutate(PARTOA=ifelse(Solar_elv < 0,0,PARTOA)) %>%
  mutate(Mean_PAR_Incident_Tower = ifelse(Solar_elv < -5,0,Mean_PAR_Incident_Tower)) %>%
  mutate(Mean_PAR_Incident_Tower = ifelse(Mean_PAR_Incident_Tower < 0, 0, Mean_PAR_Incident_Tower)) %>%
  select(-Itoa)

# biophsyical_30min_Solar %>%
#   select(-contains('Total'),-PARTOA) %>%
#   pivot_longer(cols=c(contains('PAR'))) %>%
#   ggplot(aes(Date,value,color=name,group=name))+
#   geom_line()

#### fix PAR incident issues in Oct & Nov 2014 and drift =========================================================

october_days <- seq.Date(from=as.Date('2014-10-19'),to=as.Date('2014-10-28'),by='1 day')
november_days <- as.Date(c('2014-11-01','2014-11-02','2014-11-04','2014-11-05','2014-11-07','2014-11-08','2014-11-09','2014-11-10'))

#use Mean_Shortwave_Net_Tower to fill in the gaps and replace PAR measures for these days

#create a linear model of PAR ~ Shortwave
summary(lm(Mean_PAR_Incident_Tower~Mean_Shortwave_Net_Tower,data=biophsyical_30min_Solar))
radiation_fit <- lm(Mean_PAR_Incident_Tower~Mean_Shortwave_Net_Tower,data=biophsyical_30min_Solar)

biophsyical_30min_PARcorrected <- biophsyical_30min_Solar %>%
  #create PAR predictions
  mutate(PAR_predicted = predict(radiation_fit,newdata=.)) %>% 
  #for october days, replace with PAR predictions
  mutate(Mean_PAR_Incident_Tower = ifelse(as.Date(Date) %in% october_days,PAR_predicted,Mean_PAR_Incident_Tower)) %>% 
  #for november days, replace with PAR predictions
  mutate(Mean_PAR_Incident_Tower = ifelse(as.Date(Date) %in% november_days,PAR_predicted,Mean_PAR_Incident_Tower)) %>%
  #where PAR is NA, replace with PAR predictions
  mutate(Mean_PAR_Incident_Tower = ifelse(is.na(Mean_PAR_Incident_Tower),PAR_predicted,Mean_PAR_Incident_Tower)) %>%
  select(-PAR_predicted) %>%
  #make sure PAR is 0 when the sun is down
  mutate(Mean_PAR_Incident_Tower = ifelse(Solar_elv < -5, 0, Mean_PAR_Incident_Tower)) %>%
  mutate(Mean_PAR_Incident_Tower = ifelse(Mean_PAR_Incident_Tower < 0, 0, Mean_PAR_Incident_Tower))

#correct uncalibrated PAR. (See PAR-2019-calibration.R for details)
#for PAR from 2017 and beyond, correct with linear model
biophsyical_30min_PARcorrected <- biophsyical_30min_PARcorrected %>%
  mutate(Mean_PAR_Incident_Tower = ifelse(year(Date) >= 2017,
                                    ifelse(Mean_PAR_Incident_Tower > 0,
                                           Mean_PAR_Incident_Tower*1.060160 + 36.454558,
                                           Mean_PAR_Incident_Tower),Mean_PAR_Incident_Tower)) 

#before
# biophsyical_30min_PARcorrected %>%
#   ggplot(aes(Date,Mean_PAR_Incident_Tower))+
# geom_line()

#after
# biophsyical_30min_PARcorrected %>%
#   ggplot(aes(Date,Mean_PAR_Incident_Tower))+
#   geom_line()


#
#### bring in Marsh landing for filling gaps ================================================================

#custom function to read in data from Marsh Landing
marsh_reader <- function(x){
  cnames <- as.character(read.csv(x,
                                  skip=2,
                                  nrows=1,
                                  header=F,
                                  stringsAsFactors = F))
  
  #read in data and apply column names. Making sure columns are classed appropriately
  dat <- read_csv(x,
                  skip=5,
                  col_names = cnames,
                  col_types = do.call(paste0,as.list(ifelse(cnames=='Date','T',ifelse(cnames %in% str_subset(cnames,'Flag'),'c','d'))))) %>%
    select(-NESDIS_ID,-NWSLID,-Year:-Minute) %>%
    rename(datetime=Date)
}

#list files
marshlanding_ls <- list.files('biophysical-data/Marsh-Landing/',pattern = '.csv',full.names = T)

#read them in using custom funciton
marshlanding_tibs <- lapply(marshlanding_ls,marsh_reader)

#bidn together and get rid of some uneeded columns
marshlanding <- bind_rows(marshlanding_tibs) %>% 
  select(datetime,Temp_Air,Flag_Temp_Air,Humidity,Flag_Humidity,Mean_PAR,Flag_Mean_PAR,Precipitation,Flag_Precipitation,
         Wind_Direction, Flag_Wind_Direction,Wind_Speed,Flag_Wind_Speed)

#set values that have flags to NA
flags <- marshlanding %>%
  pivot_longer(cols = starts_with('Flag'),names_to='flags',values_to = 'flag_values') %>%
  select(datetime,flags,flag_values)

meas <- marshlanding %>%
  pivot_longer(cols= c(-starts_with('Flag'),-datetime), names_to='meas',values_to = 'meas_values') %>%
  select(datetime,meas,meas_values)

marshinglanding_qcd <- bind_cols(flags,meas) %>%
  rename(Date=datetime...1) %>%
  select(-datetime...4) %>%
  #if there is a flag set measure to NA
  mutate(meas_values = ifelse(is.na(flag_values),meas_values,NA)) %>%
  #remove flag columns
  select(-flags,-flag_values) %>%
  #spread measures back out
  pivot_wider(names_from=meas,values_from=meas_values)


#update names of columns to distinguish them
colnames(marshinglanding_qcd) <- paste0('ML_',colnames(marshinglanding_qcd))
marshinglanding_qcd <- rename(marshinglanding_qcd,Date=ML_Date)

#turn to 30-min mean from 15
cluster_library(cluster,'zoo')

marshlanding_30min<- marshinglanding_qcd %>%
  select(-contains(c('Precip','Total','Direction'))) %>% #remove precipiation and total radiation measures
  pivot_longer(c(-Date)) %>%
  group_by(name) %>%
  partition(cluster) %>%
  mutate(rollingmeans = rollapply(value,width = 3, FUN=mean,na.rm=T, align='right',fill=NA, partial=T)) %>%
  collect() %>%
  ungroup() %>%
  pivot_wider(id_cols = Date,names_from=name,values_from=rollingmeans) %>%
  filter(minute(Date) == 00 | minute(Date) == 30)

#for precip, get 30min totals
marshlanding_30min_precip <- marshinglanding_qcd %>%
  select(Date,ML_Precipitation) %>%
  mutate(rollingtotal = rollapply(ML_Precipitation,width = 2, FUN=sum,na.rm=T, align='right',fill=NA, partial=T)) %>% 
  filter(minute(Date) == 00 | minute(Date) == 30)

#average wind direction. Only doing vector averaging for direction, will use scalar average from above for speed
ml_mean_wind <- marshinglanding_qcd %>%
  select(Date,ML_Wind_Direction,ML_Wind_Speed) %>%
  #calculate vector components
  mutate(u = ML_Wind_Speed*sin(pi*ML_Wind_Direction/180),
         v = ML_Wind_Speed*cos(pi*ML_Wind_Direction/180)) %>%
  #30min means of vector compoonents
  mutate(u_mean = rollapply(u,width=7,FUN=mean,na.rm=T,align='right',fill=NA,partial=T),
         v_mean = rollapply(v,width=7,FUN=mean,na.rm=T,align='right',fill=NA,partial=T)) %>%
  filter(minute(Date)== 00 | minute(Date) == 30) %>%
  #average direction from 30min vector means
  mutate(ML_Wind_Direction = (atan2(u_mean,v_mean)*180/pi)) %>% 
  #if degrees are negative, add 360
  mutate(ML_Wind_Direction = ifelse(ML_Wind_Direction < 0, ML_Wind_Direction+360,ML_Wind_Direction)) %>%
  select(Date,ML_Wind_Direction)

#add precip back to 30min means
marshlanding_30min$ML_Precipitation <- marshlanding_30min_precip$rollingtotal
marshlanding_30min$ML_Mean_Wind_Direction <- ml_mean_wind$ML_Wind_Direction

marshlanding_30min %>%
  pivot_longer(cols=c(-Date)) %>%
  ggplot(aes(Date,value,color=name))+
  geom_line()+
  facet_wrap(~name,scales='free_y')

#fill in gaps for air temp, humidity, PAR
biophsyical_30min_PARcorrected_ML <- biophsyical_30min_PARcorrected %>%
  left_join(marshlanding_30min,by='Date') %>%
  mutate(Mean_Temp_Air = ifelse(is.na(Mean_Temp_Air),ML_Temp_Air,Mean_Temp_Air),
         Mean_Humidity = ifelse(is.na(Mean_Humidity),ML_Humidity,Mean_Humidity),
         Mean_PAR_Incident_Tower = ifelse(is.na(Mean_PAR_Incident_Tower),ML_Mean_PAR,Mean_PAR_Incident_Tower),
         Total_Precip = ifelse(is.na(Total_Precip),ML_Precipitation,Total_Precip),
         Mean_Wind_Direction = ifelse(is.na(Mean_Wind_Direction),ML_Mean_Wind_Direction,Mean_Wind_Direction),
         Mean_Wind_Speed = ifelse(is.na(Mean_Wind_Speed),ML_Wind_Speed,Mean_Wind_Speed)) %>%
  select(-contains('ML_'))



#### Humidity correction and VPD Calc ===================================================================
#calculate VPD
biophsyical_30min_PARcorrected_ML_VPD <- biophsyical_30min_PARcorrected_ML %>%
  #using the GCE_flux humidity sensor (which goes over 100)
  #if >00, set to 100
  mutate(Mean_Humidity = ifelse(Mean_Humidity > 100, 100, Mean_Humidity)) %>%
  #calculate VPD
  mutate(SVP = 6.112*exp(17.67*Mean_Temp_Air/(Mean_Temp_Air+243.5)),
         AVP = ifelse(Mean_Humidity > 0,Mean_Humidity*SVP/100,NA),
         VPD = (SVP - AVP)/10) %>%
  select(-SVP,-AVP) %>%
  mutate(VPD = ifelse(VPD < 0, 0, VPD))


#### Create cloudiness index ===============================================================================

biophsyical_30min_PARcorrected_ML_VPD_CI <- biophsyical_30min_PARcorrected_ML_VPD %>%
  mutate(CI = ifelse(Mean_PAR_Incident_Tower > PARTOA,0,1-Mean_PAR_Incident_Tower/PARTOA)) %>% # add cloudiness index
  mutate(CI = ifelse(Mean_PAR_Incident_Tower == 0, 0, CI))  


#### Add salinity data form Lower Duplin ==================================================================
#Lower Duplin data for 5-min salinity measurements
duplin_list <- list.files('biophysical-data/Lower-Duplin/',full.names = T)

read_dat <- function(x){
  y <- read_csv(x,skip=5,col_names=as.character(read.csv(x,nrows=1,skip=2,header=F,stringsAsFactors = F))) %>%
    select(datetime=Date,Salinity,Flag_Salinity)
}

#create list of dfs
df_list <- lapply(duplin_list,read_dat)
#combine into one df
duplin_15 <- do.call(rbind,df_list)

#set any data with flags to NA
duplin_15_selected <- duplin_15 %>%
  mutate(salinity = ifelse(is.na(Flag_Salinity)==TRUE,Salinity,NA)) %>%
  select(Date=datetime,salinity)

#create 30min averages from 15min
duplin_30min <- duplin_15_selected %>%
  mutate(Mean_salinity = rollapply(salinity,width = 3, FUN=mean,na.rm=T, align='right',fill=NA, partial=T)) %>%
  filter(minute(Date) == 00 | minute(Date) == 30) %>%
  select(Date,Mean_salinity)

#join to full dataset
biophsyical_30min_PARcorrected_ML_VPD_CI_salinity <- biophsyical_30min_PARcorrected_ML_VPD_CI %>%
  left_join(duplin_30min,by='Date')

#### create fuzzy variables ==============================================================================================

biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity %>%
  mutate(month = as.integer(month(Date)), # add month variable 
         hour = hour(Date) + minute(Date)/60) %>% # add hour variable
  #spring: starts in january, peaks in april, ends in july
  mutate(spring = lfl::triangle(month,1,4,7), 
         #summer: starts in april, peaks in july, ends in october
         summer = lfl::triangle(month,4,7,10), 
         #autumn: starts in july, peaks in october, ends in january
         autumn = lfl::triangle(ifelse(month == 1, month+12, month),7,10,13), 
         #winter: starts in october, peaks in january, ends in april
         winter = lfl::triangle(ifelse(month <= 4, month+12,month),10,13,16)) %>% 
  mutate(morning = lfl::triangle(hour,8,14,20),
         afternoon = lfl::triangle(ifelse(hour <=2, hour+24, hour),14,20,26),
         evening = lfl::triangle(ifelse(hour <= 8, hour+24, hour),20,26,32),
         night = lfl::triangle(ifelse(hour <= 14, hour + 24, hour),25,32,38))

#### make a combined water level between marsh and creek PTs =========================================================

biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy %>%
  mutate(marsh_water = ifelse(Mean_Water_Level_Marsh < 0.84,NA,Mean_Water_Level_Marsh)) %>%
  #create a variable to identify where the measure is coming from
  mutate(PT_sensor_used = ifelse(is.na(marsh_water),'creek', 'marsh')) %>%
  #if marsh water NA, add in the Creek data
  mutate(water_level_combined = ifelse(is.na(marsh_water),Mean_Water_Level_Creek,marsh_water)) %>%
  mutate(water_level_combined = ifelse(is.na(Mean_Water_Level_Creek),Mean_Water_Level_Marsh,water_level_combined)) %>%
  #Last step, make sure that idenifer for sensor used is NA if no data
  mutate(PT_sensor_used = ifelse(is.na(water_level_combined),NA,PT_sensor_used)) %>%
  select(-marsh_water) %>%
  
  mutate(water_level_combined = ifelse(between(Date,left=as.POSIXct('2018-12-01 21:30',tz='GMT'),
                                               right = as.POSIXct('2018-12-04 06:30',tz='GMT')),NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(between(Date,left=as.POSIXct('2018-07-19 15:30',tz='GMT'),
                                               right = as.POSIXct('2018-07-20 20:00',tz='GMT')),NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(between(Date,left=as.POSIXct('2019-06-24 17:00',tz='GMT'),
                                               right = as.POSIXct('2019-06-24 18:30',tz='GMT')),NA,water_level_combined)) %>%
  
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2018-06-02 21:30',tz='GMT'),
                                       NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2018-06-02 21:30',tz='GMT'),
                                       na.approx(water_level_combined),water_level_combined)) %>%
  
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2018-03-20 23:00',tz='GMT'),
                                       NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2018-03-20 23:00',tz='GMT'),
                                       na.approx(water_level_combined),water_level_combined)) %>%
  
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2020-10-29 20:00',tz='GMT'),
                                       NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2020-10-29 20:00',tz='GMT'),
                                       na.approx(water_level_combined),water_level_combined)) %>%
  
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2018-09-30 12:00',tz='GMT'),
                                       NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2020-09-30 12:00',tz='GMT'),
                                       na.approx(water_level_combined),water_level_combined)) %>%
  
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2020-11-10 16:00',tz='GMT'),
                                       NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(Date == as.POSIXct('2020-11-10 16:00',tz='GMT'),
                                       na.approx(water_level_combined),water_level_combined)) %>%
  # mutate(water_level_combined = ifelse(as.Date(Date) >= as.Date('2021-05-08'),NA,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(year(Date) == 2017,Mean_Water_Level_Marsh,water_level_combined)) %>%
  mutate(water_level_combined = ifelse(between(Date,as.POSIXct('2021-05-06 15:00',tz='GMT'),as.POSIXct('2021-05-06 19:00',tz='GMT')),Mean_Water_Level_Creek,water_level_combined))

# crap_dates <- c(as.Date('2018-03-20'), as.Date('2018-07-30'),as.Date('2019-12-30'),as.Date('2021-05-06'))

# biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
#   filter(as.Date(Date) %in% crap_dates) %>%
#   # filter(as.Date(Date) %in% c(as.Date('2018-03-20'),as.Date('2018-03-21'))) %>%
#   # filter(year(Date) %in% c(2021,2022)) %>%
#   pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh,water_level_combined)) %>%
#   ggplot(aes(Date,value,color=name))+
#   geom_line()+
#   scale_x_datetime(date_breaks = '1 hour',date_labels='%H')+
#   facet_grid(name~as.Date(Date),scales='free')

biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
  # filter(as.Date(Date) %in% c(as.Date('2019-05-13'),as.Date('2019-05-14'),as.Date('2019-05-15'))) %>%
  # filter(as.Date(Date) %in% c(as.Date('2019-05-14'),as.Date('2019-05-15'))) %>%
  # filter(month(Date)==05 & year(Date)==2019) %>%
  # filter(year(Date) %in% c(2017:2018)) %>%
  pivot_longer(cols=c(Mean_Water_Level_Creek,Mean_Water_Level_Marsh,water_level_combined)) %>% #filter(is.na(value))
  ggplot(aes(Date,value,color=name))+
  geom_line()+
  geom_hline(yintercept = 0.84)+facet_wrap(~name,ncol=1)

# biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
#   # filter(as.Date(Date) %in% c(as.Date('2019-05-13'),as.Date('2019-05-14'),as.Date('2019-05-15'))) %>%
#   filter(as.Date(Date) %in% c(as.Date('2019-05-14'),as.Date('2019-05-15'))) %>%
#   # filter(month(Date)==05 & year(Date)==2019) %>%
#   ggplot(aes(Date,water_level_combined,color=PT_sensor_used))+
#   geom_point()+
#   geom_hline(yintercept = 0.84)

# 
# biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
#   # filter(year(Date)==2021 & month(Date) == 5) %>%
#   ggplot(aes(Date,water_level_combined))+
#   geom_line()+
#   geom_vline(xintercept = as.POSIXct('2021-05-08 00:00',tz='GMT'),color='red')
  # scale_x_datetime(date_breaks = '1 month')+
  # theme(axis.text.x = element_text(angle=45,hjust=1))+
  # facet_wrap(~year(Date),scales='free_x')
# 
# biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
#   # filter(year(Date)==2021) %>%
#   ggplot(aes(Date,Mean_Water_Level_Creek))+
#   geom_line()
# 
# biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
#   # filter(year(Date)==2021) %>%
#   ggplot(aes(Date,Mean_Water_Level_Marsh))+
#   geom_line()


#### Try to fill some NA's using spline or approx ==============================================================

#cluster zoo lirbary for rollapply
cluster_library(cluster, 'zoo')
# 
biophsyical_30min_PARcorrected_filled <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
  pivot_longer(cols=c(-Date,-PT_sensor_used)) %>%
  group_by(name) %>%
  partition(cluster) %>%
  mutate(approxed = na.approx(value,maxgap = 12,na.rm=F)) %>%
  mutate(value = ifelse(is.na(value),approxed,value)) %>%
  collect() %>%
  ungroup() %>%
  select(-approxed) %>%
  pivot_wider(names_from=name,values_from=value)

# summary(biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT)
# summary(biophsyical_30min_PARcorrected_filled)
# 
# biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT %>%
#   pivot_longer(cols=c(-Date,-PT_sensor_used)) %>%
#   count(is.na(value))
# 
# biophsyical_30min_PARcorrected_filled %>%
#   pivot_longer(cols=c(-Date,-PT_sensor_used)) %>%
#   count(is.na(value))

#### Label data as high, low, rising, or ebbing tide stages ===================================================================================

#get the high tide times
highs <- biophsyical_30min_PARcorrected_filled %>% 
  mutate(peaks = find_peaks(water_level_combined,ignore_threshold = 0.2,span=13)) %>%
  mutate(tide_stage = ifelse(peaks == TRUE, 'high',NA)) %>%
  filter(tide_stage == 'high') %>%
  select(Date,water_level_combined,tide_stage)


#empty list to store each row
low_list <- list()
#for each high tide filter to data between and find the minimum water level, this is low tide
for (i in 1:nrow(highs)){
  
  out <- biophsyical_30min_PARcorrected_filled %>%
    filter(between(Date,highs$Date[i],highs$Date[i+1])) %>%
    filter(water_level_combined == min(water_level_combined,na.rm=T)) %>%
    mutate(tide_stage = 'low') %>%
    select(Date,water_level_combined,tide_stage) 
  
  
  
  low_list[[i]] <- out[1,] #if there are two equal mins, just grab the first
}

#combine low tide data
lows <- bind_rows(low_list)

#combine low and high tides
highs_lows <- bind_rows(highs,lows) %>% arrange(Date) %>% drop_na()

#check no two highs or lows are in a row (checking for missing highs or lows)
highs_lows %>%
  mutate(lead = lead(tide_stage)) %>%
  filter(tide_stage == lead)

#count highs for each day
highs_lows %>%
  filter(tide_stage=='high') %>%
  group_by(as.Date(Date)) %>%
  count() %>% arrange(desc(n))

#save high and low table
write.csv(highs_lows,'data/processed/tide-highs-lows.csv',row.names = F)



#assign rising or ebbing between low and high tides
rising_ebbing_list <- list()
for (i in 1:nrow(highs_lows)){
  
  if (highs_lows$tide_stage[i] == 'high'){
    
    out <- biophsyical_30min_PARcorrected_filled %>%
      select(Date,water_level_combined) %>%
      filter(Date > highs_lows$Date[i] & Date < highs_lows$Date[i+1]) %>%
      mutate(tide_stage = 'ebbing')
    
  } else if (highs_lows$tide_stage[i] == 'low'){
    
    out <- biophsyical_30min_PARcorrected_filled %>%
      select(Date,water_level_combined) %>%
      filter(Date > highs_lows$Date[i] & Date < highs_lows$Date[i+1]) %>%
      mutate(tide_stage = 'rising')
    
  }
  
  rising_ebbing_list[[i]] <- out
}

#bind together
stages <- bind_rows(rising_ebbing_list) %>%
  bind_rows(highs_lows) %>%
  arrange(Date) 

#join stages back to full data
biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides <- biophsyical_30min_PARcorrected_filled %>%
  left_join(stages %>% select(-water_level_combined),by='Date')

#### flood tide characteristics =========================================================================================

tide_ids_list <- list()
for (i in 1:(nrow(lows))){
  
  tide <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides %>%
    filter(Date >= lows$Date[i] & Date < lows$Date[i+1]) %>%
    select(Date,water_level_combined,tide_stage) %>% 
    mutate(tide_id = i) 
  
  high_tide <- tide %>%
    filter(tide_stage=='high')
  
  out <- tide %>%
    mutate(high_tide_level = high_tide$water_level_combined,
           high_tide_time = high_tide$Date) %>%
    mutate(hours_from_high_tide = Date-high_tide_time) %>%
    select(-water_level_combined,-tide_stage,-high_tide_time)
  
  #convert units of time difference to seconds to they are all the same
  units(out$hours_from_high_tide) <- 'secs'
  #change to numeric and convert to decimal hours
  out$hours_from_high_tide <- as.numeric(out$hours_from_high_tide)/3600
  
  tide_ids_list[[i]] <- out
}

biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides %>%
  left_join(bind_rows(tide_ids_list),by='Date') 


#if hours from high are out of range -7.5 - 9, set to NA
biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char_ranged <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char %>%
  mutate(high_tide_level = ifelse(!between(hours_from_high_tide,-7.5,9),
                                  NA,high_tide_level)) %>%
  mutate(hours_from_high_tide = ifelse(!between(hours_from_high_tide,-7.5,9),
                                       NA,hours_from_high_tide)) %>%
  mutate(tide_stage = ifelse(!between(hours_from_high_tide,-7.5,9),
                             NA,tide_stage)) %>%
  mutate(water_level_combined = ifelse(!between(hours_from_high_tide,-7.5,9),
                                       NA,water_level_combined)) 

biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char_ranged %>%
  ggplot(aes(Date,water_level_combined))+
  geom_line(color='gray50')+
  geom_point(data=biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char_ranged %>%
               filter(tide_stage %in% c('high','low')),
             aes(Date,water_level_combined,color=tide_stage))


#### Set up file for export ================================================================================
#reorder variables and write to csv

colnames(biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char_ranged)


out <- biophsyical_30min_PARcorrected_ML_VPD_CI_salinity_fuzzy_WT_tides_char_ranged %>%
  rename(datetime=Date)
  
out %>%
  # filter(year(datetime)%in%c(2018:2021)) %>%
  pivot_longer(cols=c(Mean_PAR_Incident_Tower,Mean_Humidity,Mean_Temp_Air,VPD,CI,water_level_combined,Total_Precip,
                      Mean_Water_Level_Creek)) %>%
  ggplot(aes(datetime,value,color=name))+
  geom_line()+
  facet_wrap(~name,scales='free_y')


#write to csv with date in name
write.csv(out,
          paste0('processed/a_GCE_biophysical_variables_',
                 format(range(as.Date(out$datetime))[1],'%Y%m%d'),'_',
                 format(range(as.Date(out$datetime))[2],'%Y%m%d'),
                 '.csv'),
          row.names = F)


