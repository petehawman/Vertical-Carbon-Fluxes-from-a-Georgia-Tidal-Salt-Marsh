#Author: Peter Hawman

#Purpose:
# 1. read in processed and QC'd 30min data fluxes from flux1 and flux2
# 2. join when missing
# 3. average when coincident


#Output:
# 1. a combined file with fluxes from both and averaged when coincident

#set wd
setwd('Z:/FluxTower/GCELTER')


library(tidyverse)
library(foreach)
library(doParallel)

registerDoParallel(10)

theme_set(theme_bw())

#read in flux 1 and 2
flux1 <- read_csv(list.files('processed',pattern='b_flux1',full.names = T))

flux2 <- read_csv(list.files('processed',pattern='b_flux2',full.names = T))

#isolate nee for easier working
f1 <- flux1 %>% select(start_datetime,
                       wind_speed:v_var,RH:NEE_org,
                       co2_mole_fraction:co2_strg,
                       days_prior_maintenance,days_since_maintenance,
                       NEE_maint_qc_removal_night_ext_ustar) %>%
  mutate(sensor = 'flux1')

f2 <- flux2 %>% select(start_datetime,
                       wind_speed:v_var,RH:NEE_org,
                       co2_mole_fraction:co2_strg,
                       days_prior_maintenance,days_since_maintenance,
                       NEE_maint_qc_removal_night_ext_ustar) %>%
  mutate(sensor = 'flux2') %>%
  filter(start_datetime %in% f1$start_datetime)

results <- foreach(i=1:nrow(f1)) %dopar% {
  
  library(tidyverse)
  
  x <- f1[i,]
  y <- f2[i,]
  
  #if both have a measurement average everyting
  if (!is.na(x$NEE_maint_qc_removal_night_ext_ustar) & !is.na(y$NEE_maint_qc_removal_night_ext_ustar)){
    
    nee <-  bind_rows(x,y) %>%
        group_by(start_datetime) %>%
        summarise(sd = sd(NEE_maint_qc_removal_night_ext_ustar),
               mean = mean(NEE_maint_qc_removal_night_ext_ustar)) %>%
        ungroup() %>%
        mutate(removal = ifelse(sd < 5,'keep','throwout'))
    
    other_avg <- bind_rows(x,y) %>%
      select(-wind_speed,-wind_dir,-sensor,-NEE_maint_qc_removal_night_ext_ustar) %>%
      pivot_longer(cols=c(-start_datetime)) %>%
      group_by(name,start_datetime) %>%
      summarise(mean = mean(value)) %>%
      ungroup() %>%
      pivot_wider(id_cols=start_datetime,names_from=name,values_from=mean)
    
    wind_avg <- bind_rows(x,y) %>% 
      mutate(u = wind_speed*sin(pi*wind_dir/180),
             v = wind_speed*cos(pi*wind_dir/180)) %>%
      group_by(start_datetime) %>%
      summarise(wind_speed_a = mean(wind_speed),
                u_mean = mean(u),
                v_mean = mean(v)) %>%
      mutate(wind_dir_a = (atan2(u_mean,v_mean)*180/pi)) %>% 
      #if degrees are negative, add 360
      mutate(wind_dir_a = ifelse(wind_dir_a < 0, wind_dir_a+360,wind_dir_a)) %>%
      select(-u_mean,-v_mean) %>%
      rename(wind_dir=wind_dir_a,wind_speed=wind_speed_a)
    
    out <- nee %>% 
      left_join(other_avg,by='start_datetime') %>%
      left_join(wind_avg,by='start_datetime') %>%
      mutate(sensor = 'mean')
    
    #if flux1 has meas but flux2 doesn't, use flux 1
  } else if (!is.na(x$NEE_maint_qc_removal_night_ext_ustar) & is.na(y$NEE_maint_qc_removal_night_ext_ustar)){
    
    out <- x
    
    #if flux1 doesn't have meas but flux2 does, use flux 2
  } else if (is.na(x$NEE_maint_qc_removal_night_ext_ustar) & is.na(!y$NEE_maint_qc_removal_night_ext_ustar)) {
    
    out <- y 
    
    #if both are missing choose the one with the least amount of NA's
  } else{

    out <- tibble(which=c('x','y'),na_count = c(sum(is.na(x)),sum(is.na(y))),
                  data = c(nest(x),nest(y))) %>% 
      filter(na_count == min(na_count)) %>%
      select(data) %>%
      unnest(data) %>% unnest(data)
    
  }
  
  return(out)
  
}

stopImplicitCluster()

#bind rows. this is the merged fluxes
combined <- bind_rows(results) %>% 
  mutate(NEEcomb = ifelse(!is.na(mean),mean,NEE_maint_qc_removal_night_ext_ustar)) %>%
  mutate(NEEcomb = case_when(removal == 'throwout' ~ NA,
                             removal != 'throwout' ~ NEEcomb,
                             is.na(removal) ~ NEEcomb)) 
  
combined %>% 
  ggplot(aes(start_datetime,NEEcomb,color=sensor))+
  geom_point()

complete <- combined %>% 
  bind_rows(flux2 %>% select(start_datetime,
                             wind_speed:v_var,RH:NEE_org,
                             co2_mole_fraction:co2_strg,
                             days_prior_maintenance,days_since_maintenance,
                             NEEcomb=NEE_maint_qc_removal_night_ext_ustar) %>%
              mutate(sensor = 'flux2') %>%
              filter(!start_datetime %in% combined$start_datetime)) %>%
  mutate(end_datetime = start_datetime + 1800) %>%
  select(start_datetime,end_datetime,wind_speed:v_var,RH:NEE_org,
         co2_mole_fraction:co2_strg,
         days_prior_maintenance,days_since_maintenance,
         NEEcomb,sensor)

colnames(complete)

complete %>%
  ggplot(aes(start_datetime,NEEcomb,color=sensor))+
  geom_point()

write.csv(complete,paste0('processed/c_GCE_Fluxes_qcd_comb_',
                     format(range(as.Date(complete$end_datetime))[1],'%Y%m%d'),'_',
                     format(range(as.Date(complete$end_datetime))[2],'%Y%m%d'),
                     '.csv'),row.names = F)

