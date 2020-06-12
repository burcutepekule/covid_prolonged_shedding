# Setup
source("/net/cephfs/home/btepek/COVID_IC/setup.R")
library(readxl)
library(xtable)

# load the data 
cases_data_all   = (read.csv("/net/cephfs/home/btepek/COVID_IC/data/covid19-cases-switzerland/covid19_cases_switzerland_openzh.csv") %>% tbl_df())
deaths_data_all  = (read.csv("/net/cephfs/home/btepek/COVID_IC/data/covid19-cases-switzerland/covid19_fatalities_switzerland_openzh.csv") %>% tbl_df())
hosp_data_all    = (read.csv("/net/cephfs/home/btepek/COVID_IC/data/covid19-cases-switzerland/covid19_hospitalized_switzerland_openzh.csv") %>% tbl_df())
icu_data_all     = (read.csv("/net/cephfs/home/btepek/COVID_IC/data/covid19-cases-switzerland/covid19_icu_switzerland_openzh.csv") %>%tbl_df())
popsize_data_all = read_xlsx("/net/cephfs/home/btepek/COVID_IC/data/covid19-cases-switzerland/cantonPopSizes.xlsx") %>%tbl_df()

allDates     = as.Date(cases_data_all$Date);
date_data    = as.Date("2020-02-24");
date_end     = as.Date(tail(allDates, n=1));
data_start   = ymd(date_data)
data_end     = ymd(date_end)
data_control = ymd("2020-03-17")
date_control = ymd("2020-03-17")
date_relax   = ymd("2020-04-27")



