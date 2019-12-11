my_init_sl_fit <- function (T_tilde, Delta, A, W, t_max, sl_failure = c("SL.glm"), 
                            sl_censoring = c("SL.glm"), sl_treatment = c("SL.glm"), 
                            cvControl = list(), gtol = 0.001) 
{
    ftime <- T_tilde
    ftype <- Delta
    trt <- A
    adjustVars <- W
    t_0 <- t_max
    trtOfInterest <- 0:1
    SL.ftime <- sl_failure
    SL.ctime <- sl_censoring
    SL.trt <- sl_treatment
    cv.Control <- cvControl
    adjustVars <- data.frame(adjustVars)
    ftypeOfInterest <- unique(ftype)
    n <- length(ftime)
    id <- seq_len(n)
    dat <- data.frame(id = id, ftime = ftime, ftype = ftype, 
                      trt = trt)
    if (!is.null(adjustVars)) 
        dat <- cbind(dat, adjustVars)
    nJ <- length(ftypeOfInterest)
    allJ <- sort(unique(ftype[ftype != 0]))
    ofInterestJ <- sort(ftypeOfInterest)
    ntrt <- length(trtOfInterest)
    uniqtrt <- sort(trtOfInterest)
    message("estimating treatment")
    trtOut <- survtmle::estimateTreatment(dat = dat, ntrt = ntrt, cvControl = cv.Control, 
                                          uniqtrt = uniqtrt, adjustVars = adjustVars, SL.trt = SL.trt, 
                                          verbose = TRUE, returnModels = TRUE, gtol = gtol)
    dat <- trtOut$dat
    trtMod <- trtOut$trtMod
    dataList <- survtmle::makeDataList(dat = dat, J = allJ, ntrt = ntrt, 
                                       uniqtrt = uniqtrt, t0 = t_0, bounds = NULL)
    message("estimating censoring")
    censOut <- tryCatch({
        survtmle::estimateCensoring(dataList = dataList, ntrt = ntrt, cvControl = cv.Control,  
                                    uniqtrt = uniqtrt, t0 = t_0, verbose = TRUE, adjustVars = adjustVars, 
                                    SL.ctime = SL.ctime, glm.family = "binomial", returnModels = TRUE, 
                                    gtol = gtol)
    }, error = function(cond) {
        message("censoring sl error")
        NULL
    })
    if (is.null(censOut)) {
        censOut <- list()
        censOut$dataList <- dataList
        censOut$dataList$obs[, "G_dC"] <- 1
        censOut$dataList$"0"[, "G_dC"] <- 1
        censOut$dataList$"1"[, "G_dC"] <- 1
        is_sl_censoring_converge <- FALSE
        dataList <- censOut$dataList
    }
    else {
        dataList <- censOut$dataList
        ctimeMod <- censOut$ctimeMod
        is_sl_censoring_converge <- TRUE
    }
    message("estimating hazards")
    estOut <- survtmle::estimateHazards(dataList = dataList, cvControl = cv.Control, 
                                        J = allJ, verbose = TRUE, bounds = NULL, adjustVars = adjustVars, 
                                        SL.ftime = SL.ftime, glm.family = "binomial", returnModels = TRUE)
    dataList <- estOut$dataList
    ftimeMod <- estOut$ftimeMod
    suppressWarnings(if (all(dataList[[1]] == "convergence failure")) {
        return("estimation convergence failure")
    })
    g_1 <- dat$g_1
    g_0 <- dat$g_0
    d1 <- dataList$`1`
    d0 <- dataList$`0`
    haz1 <- d1[, c("id", "t", "Q1Haz")]
    haz1 <- tidyr::spread(haz1, t, Q1Haz)
    haz1$id <- NULL
    haz0 <- d0[, c("id", "t", "Q1Haz")]
    haz0 <- tidyr::spread(haz0, t, Q1Haz)
    haz0$id <- NULL
    S_Ac_1 <- d1[, c("id", "t", "G_dC")]
    S_Ac_1 <- tidyr::spread(S_Ac_1, t, G_dC)
    S_Ac_1 <- S_Ac_1[, -1]
    S_Ac_0 <- d0[, c("id", "t", "G_dC")]
    S_Ac_0 <- tidyr::spread(S_Ac_0, t, G_dC)
    S_Ac_0 <- S_Ac_0[, -1]
    density_failure_1 <- survival_curve$new(t = seq(range(ftime)[1], 
                                                    range(ftime)[2]), hazard = haz1)
    density_failure_0 <- survival_curve$new(t = seq(range(ftime)[1], 
                                                    range(ftime)[2]), hazard = haz0)
    density_censor_1 <- survival_curve$new(t = seq(range(ftime)[1], 
                                                   range(ftime)[2]), survival = S_Ac_1)
    density_censor_0 <- survival_curve$new(t = seq(range(ftime)[1], 
                                                   range(ftime)[2]), survival = S_Ac_0)
    return(list(density_failure_1 = density_failure_1, density_failure_0 = density_failure_0, 
                density_censor_1 = density_censor_1, density_censor_0 = density_censor_0, 
                g1W = g_1[, 1]))
}



# sl_xgboost --------------------------------------------------------------

sl_xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 800, 
                        max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                        nthread = 6, verbose = FALSE, save_period = NULL, ...) 
{
    SuperLearner:::.SL.require("xgboost")
    if (packageVersion("xgboost") < 0.6) 
        stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
    if (!is.matrix(X)) {
        X = model.matrix(~. - 1, X)
    }
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    if (family$family == "gaussian") {
        model = xgboost::xgboost(data = xgmat, objective = "reg:linear", 
                                 nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                 eta = shrinkage, verbose = verbose, nthread = nthread, 
                                 params = params, save_period = save_period)
    }
    if (family$family == "binomial") {
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                                 nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                 eta = shrinkage, verbose = verbose, nthread = nthread, 
                                 params = params, save_period = save_period)
    }
    if (family$family == "multinomial") {
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                                 nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                 eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                                 nthread = nthread, params = params, save_period = save_period)
    }
    if (!is.matrix(newX)) {
        newX = model.matrix(~. - 1, newX)
    }
    pred = predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
}

