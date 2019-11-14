library(tidyverse)
library(skimr)
library(survival)

full_data <- read_csv('data/imputed_data.csv', col_types = cols(
    .default = col_double(),
    fac = col_factor(),
    province = col_factor(),
    fac_type = col_factor(),
    educ = col_factor(),
    educ_imp = col_logical(),
    enroll_HIV_stage_imp = col_logical(),
    hh_inc = col_factor(),
    hh_inc_imp = col_logical(),
    male_imp = col_logical(),
    marital = col_factor(),
    marital_imp = col_logical(),
    mpr_days_imp = col_logical(),
    mpr_imp = col_logical()
))

txt_return_times <- full_data %>% filter(traced == 1) %>% pull(T_return)
txt_events <- full_data %>% filter(traced == 1) %>% pull(event_return)

txt_surv <- survival::Surv(time = txt_return_times, event = txt_events, type = 'right')
txt_km <- survival::survfit(txt_surv ~ 1, type = 'kaplan-meier', conf.int = .95)
plot(txt_km)

control_return_times <- full_data %>% filter(traced == 0) %>% pull(T_return)
control_events <- full_data %>% filter(traced == 0) %>% pull(event_return)

control_surv <- survival::Surv(time = control_return_times, event = control_events, type = 'right')
control_km <- survival::survfit(control_surv ~ 1, type = 'kaplan-meier', conf.int = .95)
plot(control_km)

