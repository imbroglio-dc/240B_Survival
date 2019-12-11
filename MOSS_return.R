library(tidyverse)
library(MOSS)
source(file = "my_init_sl_fit.R")
set.seed(0)

options(na.action = na.omit)

full_data <- read_rds("data/imputed_data.RDS")
data <- dplyr::select(full_data, 
                      province, fac_type, fac_size,
                      mpr, mpr_imp, yrs_ART, lost_eps, 
                      traced, T_return, event_return) %>% 
    sample_n(15e3) %>% 
    mutate(T_return = pmin(74, T_return %/% 14 + 1))

sl_lib_g <- c("SL.glm")
sl_lib_censor <- c("SL.mean", "SL.glm", "sl_xgboost")
sl_lib_failure <- c("SL.mean", "SL.glm", "sl_xgboost", "SL.step.forward")

sl_fit <- my_init_sl_fit(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    W = dplyr::select(data,
        province, fac_type, fac_size, 
        mpr, mpr_imp, yrs_ART, lost_eps
    ),
    t_max = max(data$T_return),
    sl_failure = sl_lib_failure, 
    sl_censoring = sl_lib_censor,
    sl_treatment = sl_lib_g,
    cvControl = list(V = 10, stratifyCV = FALSE, 
                     shuffle = TRUE, validRows = NULL)
)

k_grid <- sl_fit$density_failure_1$t
# sl_fit$density_failure_1$t <- k_grid
# sl_fit$density_failure_0$t <- k_grid

sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()

saveRDS(sl_fit, file = "data/sl_fit.RDS")
gc()



# A = 1 -----------------------------------------------------------------
hazard_fit_1 <- MOSS_hazard$new(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
)

hazard_1 <- hazard_fit_1$iterate_onestep(
    epsilon = 5e-2, max_num_interation = 15, verbose = TRUE, method = 'glm'
)
saveRDS(hazard_1, file = "data/hazard_1.RDS")

surv_est <- as.vector(survival_curve$new(t = k_grid, survival = hazard_1)$survival)

eic_fit <- eic$new(
    A = data$traced,
    T_tilde = data$T_return,
    Delta = data$event_return,
    density_failure = hazard_fit_1$density_failure,
    density_censor = hazard_fit_1$density_censor,
    g1W = hazard_fit_1$g1W,
    psi = surv_est,
    A_intervene = 1
)

eic_matrix <- eic_fit$all_t(k_grid = k_grid)
std_err <- compute_simultaneous_ci(eic_matrix)
upper_bound <- pmin(surv_est + 1.96 * std_err, 1)
lower_bound <- surv_est - 1.96 * std_err

out_df <- data.frame(days = 14 * (k_grid - 1),
                     lost = surv_est, 
                     u = upper_bound, 
                     l = lower_bound, 
                     type = "A = 1")



# A = 0 -----------------------------------------------------------------
hazard_fit_0 <- MOSS_hazard$new(
    T_tilde = data$T_return,
    Delta = data$event_return,
    A = data$traced,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0,
    k_grid = k_grid
)

hazard_0 <- hazard_fit_0$iterate_onestep(
    epsilon = 5e-2, max_num_interation = 15, verbose = TRUE, method = 'glm'
)
saveRDS(hazard_0, file = "data/hazard_0.RDS")

surv_est <- as.vector(survival_curve$new(t = k_grid, survival = hazard_0)$survival)

eic_fit <- eic$new(
    A = data$traced,
    T_tilde = data$T_return,
    Delta = data$event_return,
    density_failure = hazard_fit_0$density_failure,
    density_censor = hazard_fit_0$density_censor,
    g1W = hazard_fit_0$g1W,
    psi = surv_est,
    A_intervene = 0
)

eic_matrix <- eic_fit$all_t(k_grid = k_grid)
std_err <- compute_simultaneous_ci(eic_matrix)
upper_bound <- pmin(surv_est + 1.96 * std_err, 1)
lower_bound <- surv_est - 1.96 * std_err

out_df <- data.frame(days = 14 * (k_grid - 1),
                     lost = surv_est, 
                     u = upper_bound, 
                     l = lower_bound, 
                     type = "A = 0") %>% 
    rbind(out_df, .)


# Combined Survival Plot --------------------------------------------------

combined_plot <- out_df %>% 
    filter(days < 700) %>% 
    ggplot(aes(x = days, y = lost)) +
    # Add a ribbon with the confidence band
    geom_smooth(aes(ymin = l, ymax = u, fill = type, colour = type), 
                stat = "identity") +
    xlab("Days") +
    ylab("Proportion Lost") + 
    theme_minimal()

save.image("data/MOSS_return.RData")
