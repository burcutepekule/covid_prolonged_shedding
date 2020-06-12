# Setup
rm(list=ls())
library(rstan)
library(zoo)
library(Rcpp)
library(lubridate)
library(cowplot)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("/net/cephfs/home/btepek/COVID_IC/setup.R")

args <- commandArgs(TRUE)
print(args)

r_c_in <- as.integer(args[1])
m_relax_in <- as.integer(args[2])
end_index_in <- as.integer(args[3])

# country_name    = "Italy"; # end_index_in    = 6;
# country_name    = "United States of America"; # end_index_in    = 6;
# country_name    = "Spain"; # end_index_in    = 10;
country_name    = "France"; # end_index_in    = 9;
# country_name    = "China"; # end_index_in    = 14;
# country_name    = "Greece"; # end_index_in    = 7;

data_all        = (read.csv("/net/cephfs/home/btepek/COVID_IC/data/DATA_WHO.csv") %>% tbl_df())
data_country    = data_all %>% filter(ADM0_NAME==country_name)
dates_country   = data_country$date_epicrv;
allDates        = as.Date(unlist(dates_country), format = "%Y-%m-%dT00:00:00.000Z");
cases_data_all  = data_country$NewCase;
deaths_data_all = data_country$NewDeath;

# get rid of the last day 
allDates        = allDates[-length(allDates)];
cases_data_all  = cases_data_all[-length(cases_data_all)];
deaths_data_all = deaths_data_all[-length(deaths_data_all)];

date_data      = as.Date(allDates[1]-1);
date_end       = as.Date(tail(allDates, n=1));
data_start     = ymd(date_data);
data_end       = ymd(date_end);

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

if((country_name=="Spain")==TRUE){
  print("Calculating for Spain")
  popsize_data_all = 46752096;
  data_control = ymd("2020-03-14")
  date_control = ymd("2020-03-14")
  date_relax   = ymd("2020-05-09")
}
if((country_name=="Italy")==TRUE){
  print("Calculating for Italy")
  popsize_data_all = 60475018;
  data_control = ymd("2020-03-08")
  date_control = ymd("2020-03-08")
  date_relax   = ymd("2020-04-14")
}
if((country_name=="France")==TRUE){
  print("Calculating for France")
  popsize_data_all = 65252509;
  data_control = ymd("2020-03-17")
  date_control = ymd("2020-03-17")
  date_relax   = ymd("2020-05-11")
}

if((country_name=="China")==TRUE){
  print("Calculating for China")
  popsize_data_all = 1438523180;
  data_control = ymd("2020-01-23")
  date_control = ymd("2020-01-23")
  date_relax   = ymd("2020-04-08")
}

if((country_name=="United States of America")==TRUE){
  print("Calculating for United States of America")
  popsize_data_all = 330722412;
  data_control = ymd("2020-03-20")
  date_control = ymd("2020-03-30")
  date_relax   = ymd("2020-05-15")
  deaths_data_all[106] = -1*deaths_data_all[106]; #I THINK THERE IS A MISTAKE IN THE TABLE
}


# after lockdown before relaxation time
vecIndex_begin  = which(abs(allDates-date_control) == min(abs(allDates - date_control)));
vecIndex_end    = which(abs(allDates-date_relax) == min(abs(allDates - date_relax)));
vecIndex_diff   = vecIndex_end-vecIndex_begin+1;

vecIndex        = vecIndex_begin+end_index_in*5;
use_data_end    = allDates[vecIndex];
cases_data_use  = cases_data_all[1:vecIndex];
deaths_data_use = deaths_data_all[1:vecIndex];
allDates_use    = allDates[1:vecIndex];

cases_data_all = cases_data_all[1:vecIndex_end];
deaths_data_all= deaths_data_all[1:vecIndex_end];
allDates       = allDates[1:vecIndex_end];
date_end       = date_relax; #use relaxation as the end of data

# 1 day before the observation of the first case
daily_cases_data_use   = append(cases_data_use,0,after=0);
daily_deaths_data_use  = append(deaths_data_use,0,after=0)
daily_cases_data_all   = append(cases_data_all,0,after=0)
daily_deaths_data_all  = append(deaths_data_all,0,after=0)
allDates_agg           = append(allDates,data_start,after=0);


agg_data_all     = tibble(allDates_agg,daily_cases_data_all,daily_deaths_data_all);
allDates_agg_use = append(allDates_use,data_start,after=0);
agg_data_use     = tibble(allDates_agg_use,daily_cases_data_use,daily_deaths_data_use);
observed_data    = c(daily_cases_data_use,daily_deaths_data_use);
observed_dates   = c(allDates_agg_use,allDates_agg_use);


write.table(agg_data_use, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/agg_data_use_",end_index_in,".csv"))
write.table(agg_data_all, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/agg_data_all_",end_index_in,".csv"))

days2add     = vecIndex_end-vecIndex+7;
date_simul   = use_data_end + days2add;


print(c(r_c_in,m_relax_in,end_index_in))

data_list = list(
  pop_t=popsize_data_all,
  tswitch=as.numeric(date_control-date_data+1),
  trelax=as.numeric(date_relax-date_data+1),
  r_end=1,
  m_relax=m_relax_in/100,
  gamma_c=1/14,
  D=as.numeric(use_data_end-date_data+1),
  k_daily_cases  = daily_cases_data_use,
  k_daily_deaths = daily_deaths_data_use,
  
  p_pi        = c(1,999),
  p_R0        = c(2.5,0.5),
  p_tau       = 1/2.5,
  p_gamma_s   = 1/2.5,
  p_gamma_H   = 1/12,
  p_gamma_ICU = 1/12,
  p_eps_H     = c(0.08,0.02),
  p_eps_H2ICU = c(0.4,0.08),
  p_eps_x_ICU = c(0.4,0.08),
  p_r_d_s     = c(0.2,0.03),
  p_r_d_c     = c(0.075,0.015),
  p_r_lock    = c(1,1),
  p_phi       = 1/100,
  
  t0=0,
  t_data=1,
  S=as.numeric(date_simul-date_data+1),
  E=days2add, 
  ts=1:as.numeric(use_data_end-date_data+1),
  ts_pred=1:as.numeric(date_simul-date_data+1),
  
  r_c=r_c_in
)

# M_model_BT     = stan_model("/net/cephfs/home/btepek/COVID_IC/models/model_DENMARK.stan")
M_model_BT     = readRDS("/net/cephfs/home/btepek/COVID_IC/models/model_DENMARK.rds")
T_modelBT      = sampling(M_model_BT,data = data_list,warmup=500,iter=1000,chains=8,init=0.5)

# save(T_modelBT, file =paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/RDATA/T_modelBT_rc_",r_c_in,"_mrelax_",m_relax_in,"_endidx_",end_index_in,".RData"))

pp = c("fitted_k_daily_cases","fitted_k_daily_deaths")
model_output_fitted = summary(T_modelBT,pp)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$D,2),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp,each=data_list$D)) %>%
  mutate(populations=factor(type,levels=pp,
                            labels=c("Number of daily cases","Number of daily deaths"))) 

agg_data_all = tibble(all_dates=observed_dates,all_data=observed_data,populations=model_output_fitted$populations);
write.table(model_output_fitted, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/CSVS/output_fit_rc_",r_c_in,"_mrelax_",m_relax_in,"_endidx_",end_index_in,".csv"))


ggplot() +
  geom_ribbon(data=model_output_fitted,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=populations),alpha=.5) +
  geom_line(data=model_output_fitted,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ populations ,scales="free",nrow=2) +
  geom_vline(xintercept=date_control,linetype=2) +
  geom_vline(xintercept=date_end,linetype=2) +
  scale_colour_manual(values=c("grey20","black"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  scale_x_date(date_breaks="2 weeks",date_labels="%b %d") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_point(data=agg_data_all,aes(x=all_dates,y=all_data,fill=populations),shape=21,size=2,colour = "black", fill = "white") +
  labs(x="Days",y="N")+
  scale_y_continuous(trans = 'log10')
ggsave(paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/FIGS/figure_fit_rc_",r_c_in,"_mrelax_",m_relax_in,".png"),width = 20, height = 10)


parameterSummary = summary(T_modelBT, c("R0","tau","gamma_s","gamma_H","gamma_ICU","eps_H","eps_H2ICU","eps_x_ICU","r_d_s","r_d_c","r_lock","shift_lock","m_lock","phi"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
write.table(parameterSummary$summary, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/CSVS/parameters_fit_rc_",r_c_in,"_mrelax_",m_relax_in,"_endidx_",end_index_in,".csv"))


pp3 = c("predicted_daily_cases","predicted_daily_deaths")
model_output_all = summary(T_modelBT,pp3)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$S,2),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp3,each=data_list$S)) %>%
  mutate(populations=factor(type,levels=pp3,
                            labels=c("Number of daily cases","Number of daily deaths")))
write.table(model_output_all, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/CSVS/output_all_rc_",r_c_in,"_mrelax_",m_relax_in,"_endidx_",end_index_in,".csv"))


ggplot() +
  geom_ribbon(data=model_output_all,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=populations),alpha=.5) +
  geom_line(data=model_output_all,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ populations ,scales="free",nrow=2) +
  geom_vline(xintercept=date_control,linetype=2) +
  geom_vline(xintercept=date_end,linetype=2) +
  scale_colour_manual(values=c("grey20","black"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  scale_x_date(date_breaks="2 weeks",date_labels="%b %d") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_point(data=agg_data_all,aes(x=all_dates,y=all_data,fill=populations),shape=21,size=2,colour = "black", fill = "white") +
  labs(x="Days",y="N")+
  scale_y_continuous(trans = 'log10')
ggsave(paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_FRANCE/FIGS/figure_all_rc_",r_c_in,"_mrelax_",m_relax_in,"_endidx_",end_index_in,".png"),width = 20, height = 10)

