library(tidyverse)
library(skimr)
library(MOSS)

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

data <- full_data %>% sample_n(300)

sl_lib_g <- c("SL.mean", "SL.glm")
sl_lib_censor <- c("SL.mean", "SL.glm")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.step.forward")

initial_sl_fit <- MOSS::initial_sl_fit(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    W = data %>% select(
        province, fac_type, fac_size, 
        age, educ, educ_imp, enroll_HIV_stage, enroll_HIV_stage_imp, 
        hh_inc, hh_inc_imp, male, male_imp,
        marital, marital_imp, mpr, mpr_days, mpr_imp, tracer, yrs_ART, yrs_enrolled
    ),
    t_max = data %>% select(T_return) %>% max(),
    sl_treatment = sl_lib_g,
    sl_censoring = sl_lib_censor,
    sl_failure = sl_lib_failure
    )

initial_sl_fit$density_failure_1$hazard_to_survival()
k_grid <- 1:max(data$T_return)
initial_sl_fit$density_failure_1$t <- k_grid
initial_sl_fit$density_failure_0$t <- k_grid

moss_hazard_fit <- MOSS_hazard$new(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    density_failure = initial_sl_fit$density_failure_1,
    density_censor = initial_sl_fit$density_censor_1,
    g1W = initial_sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
)

psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(
    epsilon = 1, max_num_interation = 2, verbose = TRUE, method = 'l2'
)

moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)
moss_hazard_fit_1$display(type = 'survival')

survival_curve_estimate <- as.vector(moss_hazard_fit_1$survival)

eic_fit <- eic$new(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    density_failure = moss_hazard_fit_1$density_failure,
    density_censor = moss_hazard_fit_1$density_censor,
    g1W = moss_hazard_fit_1$g1W,
    psi = survival_curve_estimate,
    A_intervene = 1
)

eic_matrix <- eic_fit$all_t(k_grid = k_grid)
std_err <- compute_simultaneous_ci(eic_matrix)
upper_bound <- survival_curve_estimate + 1.96 * std_err
lower_bound <- survival_curve_estimate - 1.96 * std_err
print(survival_curve_estimate)
print(upper_bound)
print(lower_bound)
