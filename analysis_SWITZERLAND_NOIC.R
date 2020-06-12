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

path2save = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_SWITZERLAND");
dir.create(path2save)
dir.create(paste0(path2save,"/CSVS/"))
dir.create(paste0(path2save,"/FIGS/"))
dir.create(paste0(path2save,"/RDATA/"))

args <- commandArgs(TRUE)
print(args)

r_c_in <- as.integer(args[1])
m_relax_in <- as.integer(args[2])

source("/net/cephfs/home/btepek/COVID_IC/prepare_data.R")
source("/net/cephfs/home/btepek/COVID_IC/setup.R")

# no issue with CH in terms of NA values, but might need to use na.locf at one point for other cantons
canton_name = "CH";
cases_data  = cases_data_all %>% select(canton_name);
deaths_data = deaths_data_all %>% select(canton_name);
icu_data    = icu_data_all %>% select(canton_name);
hosp_data   = hosp_data_all %>% select(canton_name);
pop_size    = as.numeric(popsize_data_all[popsize_data_all$name==canton_name,2]) 

daily_cases_data = append(diff(as.numeric(unlist(cases_data)),lag=1),as.numeric(cases_data[1,1]),after=0);
daily_deaths_data= append(diff(as.numeric(unlist(deaths_data)),lag=1),as.numeric(deaths_data[1,1]),after=0);
icu_data         = as.numeric(unlist(icu_data));
hosp_data        = as.numeric(unlist(hosp_data));
daily_cases_data = append(daily_cases_data,0,after=0);
daily_deaths_data= append(daily_deaths_data,0,after=0);
icu_data         = append(icu_data,0,after=0);
hosp_data        = append(hosp_data,0,after=0);
allDates_agg     = append(allDates,data_start,after=0);
agg_data         = tibble(allDates_agg,daily_cases_data,daily_deaths_data,icu_data,hosp_data);
observed_data    = c(daily_cases_data,hosp_data,icu_data,daily_deaths_data);
observed_dates   = c(allDates_agg,allDates_agg,allDates_agg,allDates_agg);

days2add     = 26*7;
date_simul   = date_end + days2add;

print(c(r_c_in,m_relax_in))

data_list = list(
  pop_t=pop_size,
  tswitch=as.numeric(date_control-date_data+1),
  trelax=as.numeric(date_relax-date_data+1),
  r_end=1,
  m_relax=m_relax_in/100,
  gamma_c=1/14,
  D=as.numeric(date_end-date_data+1),
  k_daily_cases  = daily_cases_data,
  k_hospit       = hosp_data,
  k_icu          = icu_data,
  k_daily_deaths = daily_deaths_data,
  
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
  ts=1:as.numeric(date_end-date_data+1),
  ts_pred=1:as.numeric(date_simul-date_data+1),
  
  r_c=r_c_in
)

# M_model_BT   = stan_model("/net/cephfs/home/btepek/COVID_IC/models/model_BT.stan")
M_model_BT     = readRDS("/net/cephfs/home/btepek/COVID_IC/models/model_BT.rds")
T_modelBT      = sampling(M_model_BT,data = data_list,warmup=500,iter=1000,chains=8,init=0.5)

# save(T_modelBT, file =paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_SWITZERLAND/RDATA/T_modelBT_rc_",r_c_in,"_mrelax_",m_relax_in,".RData"))

parameterSummary = summary(T_modelBT, c("R0","tau","gamma_s","gamma_H","gamma_ICU","eps_H","eps_H2ICU","eps_x_ICU","r_d_s","r_d_c","r_lock","shift_lock","m_lock","phi"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
write.table(parameterSummary$summary, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_SWITZERLAND/CSVS/parameters_fit_rc_",r_c_in,"_mrelax_",m_relax_in,".csv"))

pp3 = c("predicted_daily_cases","predicted_current_hospit","predicted_current_icu","predicted_daily_deaths")
model_output_all = summary(T_modelBT,pp3)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$S,4),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp3,each=data_list$S)) %>%
  mutate(populations=factor(type,levels=pp3,
                            labels=c("Number of daily cases","Hospital admissions","ICU admissions","Number of daily deaths")))
write.table(model_output_all, file = paste0("/net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_SWITZERLAND/CSVS/output_all_rc_",r_c_in,"_mrelax_",m_relax_in,".csv"))

