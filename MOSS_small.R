library(tidyverse)
library(MOSS)
source(file = "my_init_sl_fit.R")
set.seed(0)

full_data <- read_rds("data/imputed_data.RDS")
data <- dplyr::select(full_data, 
                      province, fac_type, fac_size, 
                      enroll_HIV_stage, enroll_HIV_stage_imp, mpr, mpr_imp, 
                      yrs_ART, yrs_enrolled, lost_eps, 
                      traced, T_return, event_return) %>% 
    dplyr::filter(T_return > 0) %>% 
    sample_n(1e3) %>% 
    mutate(event_return = case_when(T_return > 730 ~ 0, 
                                    TRUE ~ event_return),
           T_return = case_when(T_return > 730 ~ 730 + T_return %% 180, 
                                TRUE ~ T_return), 
           T_return = T_return %/% 14)

sl_lib_g <- c("SL.glm")
sl_lib_censor <- c("SL.mean", "SL.glm", "sl_xgboost")
sl_lib_failure <- c("SL.mean", "SL.glm", "sl_xgboost", "SL.step.forward")

sl_fit <- my_init_sl_fit(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    W = dplyr::select(data,
                      province, fac_type,
                      enroll_HIV_stage, enroll_HIV_stage_imp, mpr, mpr_imp, 
                      yrs_ART, yrs_enrolled, lost_eps),
    t_max = max(data$T_return),
    sl_failure = sl_lib_failure,
    sl_censoring = sl_lib_censor,
    sl_treatment = sl_lib_g,
    cvControl = list(V = 5, stratifyCV = FALSE,
                     shuffle = TRUE, validRows = NULL)
)


sl_fit$density_failure_1$hazard_to_survival()
saveRDS(sl_fit, file = "data/sl_fit_small.RDS")

sl_fit <- read_rds("data/sl_fit_small.RDS")
k_grid <- 1:max(data$T_return)
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid


# Treated -----------------------------------------------------------------

moss_hazard_fit <- MOSS_hazard$new(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
)


psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(
    epsilon = 1e-2, max_num_interation = 50, verbose = TRUE, method = 'glm'
)


moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)
moss_hazard_fit_1$display(type = 'survival')

survival_curve_estimate <- as.vector(moss_hazard_fit_1$survival)
eic_fit <- eic$new(
    A = data$traced,
    T_tilde = data$T_return,
    Delta = data$event_return,
    density_failure = moss_hazard_fit$density_failure,
    density_censor = moss_hazard_fit$density_censor,
    g1W = moss_hazard_fit$g1W,
    psi = survival_curve_estimate,
    A_intervene = 1
)

eic_matrix <- eic_fit$all_t(k_grid = k_grid)
std_err <- compute_simultaneous_ci(eic_matrix)
upper_bound <- survival_curve_estimate + 1.96 * std_err
lower_bound <- survival_curve_estimate - 1.96 * std_err

out_df <- data.frame(days = 14*(1:length(upper_bound)),
           lost = survival_curve_estimate, 
           u = upper_bound, 
           l = lower_bound, 
           type = rep("A = 1", length(upper_bound)))
treated_plot <- out_df %>% 
    ggplot(aes(x = days, y = lost)) +
    # Add a ribbon with the confidence band
    geom_smooth(aes(ymin = l, ymax = u), stat = "identity") +
    xlab("Days") +
    ylab("Proportion Lost")



# Control -----------------------------------------------------------------

moss_hazard_fit <- MOSS_hazard$new(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 0,
    k_grid = k_grid
)


psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(
    epsilon = 1e-2, max_num_interation = 50, verbose = TRUE, method = 'glm'
)


moss_hazard_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)
moss_hazard_fit_0$display(type = 'survival')

survival_curve_estimate <- as.vector(moss_hazard_fit_0$survival)
eic_fit <- eic$new(
    A = data$traced,
    T_tilde = data$T_return,
    Delta = data$event_return,
    density_failure = moss_hazard_fit$density_failure,
    density_censor = moss_hazard_fit$density_censor,
    g1W = moss_hazard_fit$g1W,
    psi = survival_curve_estimate,
    A_intervene = 0
)

eic_matrix <- eic_fit$all_t(k_grid = k_grid)
std_err <- compute_simultaneous_ci(eic_matrix)
upper_bound <- survival_curve_estimate + 1.96 * std_err
lower_bound <- survival_curve_estimate - 1.96 * std_err


out_df_0 <- data.frame(days = 14*(1:length(upper_bound)),
                       lost = survival_curve_estimate, 
                       u = upper_bound, 
                       l = lower_bound, 
                       type = rep("A = 0", length(upper_bound))) %>% 
    rbind(out_df, .)


# Combined Survival Plot --------------------------------------------------

combined_plot <- out_df_0 %>% 
    ggplot(aes(x = days, y = lost)) +
    # Add a ribbon with the confidence band
    geom_smooth(aes(ymin = l, ymax = u, fill = type, colour = type), 
                stat = "identity") +
    xlab("Days") +
    ylab("Proportion Lost")

