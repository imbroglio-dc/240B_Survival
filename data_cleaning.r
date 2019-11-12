library(tidyverse); library(lubridate)
set.seed(0)
# briefly uses packages haven, mice

setwd("~/projects/Zambia/")
# Load & Format Data ---------------------------------------------------------------
ltfu <- haven::read_dta("inbound/upd_analysis_art.dta") %>% 
    filter(lost == TRUE) %>% 
    left_join(., haven::read_dta("inbound/mpr.dta"), "ARTPatID") %>% 
    left_join(., haven::read_dta("inbound/prob_tracer.dta"), "ARTPatID") %>% 
    transmute(pat_id = as.numeric(as_factor(ARTPatID)),
              province = as_factor(Province), 
              fac_type = as_factor(Facilitytype), 
              fac_size = as.numeric(fac_size),
              fac = as.character(haven::as_factor(CurrentClinic)), 
              txt_date = DBClose, 
              apt_date = AptDate, 
              visit_date = visitdate, 
              last_visit = lastvisit,
              age = time_length(txt_date - dob, "year"), 
              yrs_ART = time_length(txt_date - ARTStartD, "year"), 
              yrs_enrolled = time_length(txt_date - DateOfEnroll, "year"),
              yrs_enroll_cd4 = time_length(txt_date - enrollcd4Dt, "year"),
              yrs_init_cd4 = time_length(txt_date - initiationCD4dt, "year"),
              yrs_last_cd4 = time_length(txt_date - lastcd4dt, "year"),
              educ = factor(case_when(Educ == 0 ~ "none",
                                      Educ %in% 1:6 ~ "1-6",
                                      Educ %in% 7:12 ~ "7-12",
                                      Educ == 13 ~ "college/univ"), 
                            levels = c("none", "1-6", "7-12", "college/univ"),
                            ordered = T), 
              enroll_cd4 = as.numeric(enrollcd4),
              init_cd4 = as.numeric(initiationCD4), 
              last_cd4 = as.numeric(lastcd4),
              enroll_HIV_stage = factor(HIVEnrollStage, ordered = T),
              hh_inc = factor(HHInc, 1:5, 
                              c("<50K", "50-99K", "100-199K", "200-499K", ">500K"), 
                              ordered = T), 
              male = as.numeric(Male), 
              marital = haven::as_factor(Marital), 
              mpr = as.numeric(mpr), 
              mpr_days = time_length(txt_date - mprd, "day"), 
              traced = case_when(sampled == 1 ~ 1, T ~ 0), # traced iff sampled == 1
              trace_result = haven::as_factor(tracing_result), 
              chart_rev_delay = time_length(a2_form_date - txt_date, "day"), 
              tracing_delay = time_length(b13_interview_date - txt_date, "day"),
              transfer_days = time_length(transferdate - txt_date, "day"),
              p1 = as.numeric(p1),
              p2 = as.numeric(p2),
              tracer = b01_staff_id) %>% 
    arrange(pat_id, visit_date, apt_date)



# facility level ----------------------------------------------------------
facilities <- ltfu %>% 
    dplyr::select(pat_id, fac, province, fac_type, fac_size, txt_date, p1, p2) %>% 
    distinct()

# Account for facility splitting done during sampling ---------------------
facilities <- facilities %>% 
    mutate(fac = case_when(fac_size == 1306 ~ "Kafue District Hospital_A", 
                           fac_size == 1319 ~ "Kafue District Hospital_B", 
                           T ~ fac), 
           fac = as_factor(fac))

# end of study ------------------------------------------------------------
facilities <- facilities %>% 
    mutate(txt_date = time_length(max(ltfu$visit_date, na.rm = T) - txt_date, "day")) %>% 
    rename(study_end = txt_date)

# probability of sampling -------------------------------------------------




# time relalated ----------------------------------------------------------
dates <- ltfu %>% dplyr::select(pat_id, txt_date, apt_date, visit_date, last_visit)

# LTFU history variables --------------------------------------------------
dates <- dates %>% 
    mutate_at(vars(contains("date")), lubridate::date) %>%
    distinct() %>% 
    transmute(pat_id = pat_id,
              visit = time_length(visit_date - txt_date, "day"),
              apt = time_length(apt_date - txt_date, "day"),
              last = time_length(last_visit - txt_date, "day"),
              last_visit = time_length(last_visit - txt_date, "day")) %>% 
    gather(key = "tmp", value = "visit", visit, last_visit) %>% 
    select(-tmp) %>% 
    gather(key = "type", value = 'days', -pat_id, -last) %>%
    distinct() %>% arrange(pat_id, days, type)

# periods of ltfu in the 2 years prior to t_0
ltfu_hx <- dates %>% 
    filter(days <= last) %>%
    distinct() %>%
    arrange(pat_id, days, type) %>% 
    mutate(pat_id = pat_id,
           diff = c(diff(days, lag=1), NaN),
           type = c(paste0(head(type, -1), " -> ", tail(type, -1)), "NA"),
           next_pat = c(diff(pat_id), 1)) %>% 
    filter(str_ends(type, "visit") | next_pat == 1) %>% 
    mutate(diff = case_when(next_pat == 1 ~ -last, 
                            T ~ diff)) %>% 
    filter((str_starts(type, "apt") & diff > 90 ) | 
               (str_starts(type, "visit") & diff > 180) | 
               next_pat == 1) %>% 
    group_by(pat_id) %>% 
    summarise(days_ltfu = sum(diff), lost_eps = n(), 
              days_last = -unique(last) - 90) %>% 
    mutate(days_ltfu = days_ltfu - 90*lost_eps)

# outcome = visits after txt assigment
outcomes <- dates %>% 
    filter(days >= 0) %>% 
    transmute(pat_id = pat_id,
              last = last,
              type = c("NA", paste0(head(type, -1), " -> ", tail(type, -1))),
              days = days,
              diff = c(NaN, diff(days, lag=1)),
              next_pat = c(1, diff(pat_id)) > 0) %>% 
    filter(!(str_ends(type, "apt") & next_pat == T)) %>% 
    mutate(next_pat = as.numeric(c(1, diff(pat_id)) > 0),
           diff = case_when(next_pat == T ~ days - last, 
                            T ~ diff)) %>% 
    filter(next_pat == 1 | 
               (type == "apt -> visit" & diff > 90) | 
               (type == "visit -> visit" & diff > 180)) %>% 
    arrange(pat_id, days) %>% 
    dplyr::select(pat_id, days) %>% 
    group_by(pat_id) %>% 
    summarise(T_return = head(days, 1), 
              T_ltfu = case_when(length(head(days, 2)) == 1 ~ Inf, 
                                 T ~ tail(head(days, 2), 1))) %>% 
    left_join(distinct(dplyr::select(dates, pat_id)), .) %>% 
    mutate(T_return = case_when(is.na(T_return) ~ Inf, 
                                T ~ T_return),
           T_ltfu = case_when(is.na(T_ltfu) ~ Inf, 
                              T ~ T_ltfu))



# individual level --------------------------------------------------------
individuals <- ltfu %>% 
    dplyr::select(pat_id, age, yrs_ART, yrs_enrolled, yrs_ART, educ,
                  enroll_HIV_stage, hh_inc, male, marital, mpr, mpr_days, 
                  traced, trace_result, chart_rev_delay, tracing_delay, 
                  transfer_days, tracer) %>% 
    distinct()

# negative chart review delay ---------------------------------------------
# 1 case, -192 days before txt date
# chart review delay ~ half a year to a full year!
individuals$chart_rev_delay[individuals$chart_rev_delay < 0] <- NA

# deal with cd4s ----------------------------------------------------------
# concerns about time ordering of cd4 measurements relative to t_0 and each other
# will remove for now - may circle back if subject matter experts think it's important
individuals <- individuals %>% dplyr::select_at(vars(-contains("cd4")))

# filter unrealistic ages -------------------------------------------------
# ages between 15 and 85 (semi-arbitrary cut-offs) at time of txt assignment
# deleted 33 (5 traced) but could leave for imputation
individuals <-  filter(individuals, age < 85 & age > 15)


# fill in tracer ----------------------------------------------------------
individuals <- individuals %>% 
    mutate(tracer = factor(case_when(is.na(tracer) ~ 0, 
                              T ~ tracer)))



# consolidate -------------------------------------------------------------
data <- inner_join(inner_join(outcomes, facilities), inner_join(ltfu_hx, individuals)) %>% 
    mutate(T_return = case_when(is.infinite(T_return) ~ study_end, 
                                T ~ T_return), 
           T_ltfu = case_when(is.infinite(T_ltfu) ~ study_end, 
                                T ~ T_ltfu))

rm(list = ls(pattern="[^data]"))
setwd("~/projects/240B_Survival/")
write_csv(data, path = "data/cleaned_data.csv")




# imputation --------------------------------------------------------------
library(mice)

imputation <- data %>% 
    dplyr::select(-pat_id, -T_return, -T_ltfu, -fac, -study_end, -p1,
                  -p2, -trace_result, -chart_rev_delay, -tracing_delay, -transfer_days,
                  -tracer, -traced) %>% 
    mice(., m = 5, maxit = 25, visitSequence = "monotone")

imputed_data <- mice::complete(imputation, action = 1) %>% 
    full_join(data, ., c('province', 'fac_type', 'fac_size', 'days_ltfu', 'lost_eps',
                         'days_last', 'age', 'yrs_ART', 'yrs_enrolled')) %>% 
    mutate_at(vars(ends_with(".x")), is.na) %>% 
    rename_at(vars(ends_with(".x")), ~str_replace(., "\\.x", "_imp")) %>% 
    rename_at(vars(ends_with(".y")), ~str_remove(., "\\.y"))

imputed_data <- imputed_data %>% 
    mutate(event_return = as.numeric(T_return != study_end),
           event_ltfu = as.numeric(T_ltfu != study_end)) %>% 
    dplyr::select(T_return, event_return, T_ltfu, event_ltfu, traced, everything()) %>% 
    dplyr::select(-pat_id, -study_end, -p1, -p2, -days_ltfu, -trace_result,
                  -chart_rev_delay, -tracing_delay, -transfer_days)

imputed_data <- imputed_data[, c(1:9, 9+order(colnames(imputed_data)[-c(1:9)]))]

write_csv(imputed_data, "./data/imputed_data.csv")




# Figuring things out -----------------------------------------------------
tmp <- data %>% select(fac, province, fac_type, p1, fac_size) %>% distinct %>% 
    mutate_at(c("p1"), ~case_when(is.na(.) ~ 0, 
                                        T ~ .)) %>% 
    arrange(fac) %>% 
    group_by(fac, province, fac_type, fac_size) %>% 
    summarise_all(~max(., na.rm = T)) %>% 
    ungroup %>% group_by(province, fac_type) %>% 
    arrange(province, fac_type)

tmp <- haven::read_dta("../Zambia/inbound/lostbyclinic.dta") %>%
    rename(fac = clinic) %>% full_join(tmp, .)

tmp <- dplyr::select(data, fac) %>% group_by_all() %>% summarise(n_ltfu = n()) %>% 
    full_join(tmp, .)

tmp %>% group_by(province, fac_type) %>% 
    mutate_at(c("fac_size", "totallost", "n_ltfu"), ~1/(sum(.) / (2*.))) %>% 
    head(20)

