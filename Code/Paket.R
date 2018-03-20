# Funktion zum Berechnen von Wahrscheinlichkeit 
calculate_p <- function(X, y, beta) {
    
    eta <- X %*% beta
    exp_eta <- exp(eta)
    
    p <- exp_eta / (1 + exp_eta)
    
    return(as.vector(p))
}


# Funktion zum Berechnen des Deviance Residuals
calculate_Deviance_Residuals <- function(X, y, beta) {
    
    eta <- X %*% beta #
    s <- y
    s[s == 0] = -1
    
    d = s * sqrt(-2*((y * eta) - (log(1 + exp(eta)))))
    
    return(as.numeric(d))
}


# Funktion zum Berechnung des Maximalen Likelihoods
maxLikeEst <- function(X, y) {
    
    beta <- rep(0, times = ncol(X))
    M <- diag(nrow = nrow(X))
    tolerance <- exp(-20)
    diff <- 10 * abs(tolerance)
    maxIteration <- 5000
    i <- 0
    
    while (diff > tolerance || i < maxIteration) {
        
        p <- calculate_p(X, y, beta)
        M <- diag(p * (1 - p))
        beta_change <- solve(t(X) %*% M %*% X) %*% t(X) %*% (y - p)
        beta <- beta + beta_change
        diff <- sum(abs(beta_change))
        i <- i + 1
    }
    
    # degrees of freedom = observations - params
    dfRes <- nrow(X) - ncol(X)   
    dfNull <- nrow(X) - 1 
    
    # cov matrix
    vcov <- solve(t(X) %*% M %*% X) 
    
    # Deviance Residual
    devianceResidual <- calculate_Deviance_Residuals(X, y, beta)
    
    # Maximumwert der Log Likelihood Funktion
    maxLogLikeValue <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))
    
    result <- list(coefficients = beta,
                   vcov = vcov,
                   devianceResidual = devianceResidual,
                   dfRes = dfRes,
                   dfNull = dfNull,
                   maxLogLikeValue = maxLogLikeValue) 
    
    return(result)
    
}


# Funktion zum Erstellen des Logit-Modells mit der Klasse LogitMod
logitMod <- function(formula, data) {
    
    modelFrame <- model.frame(formula, data)
    X <- model.matrix(formula, modelFrame)
    y <- model.response(modelFrame)
    
    result <- maxLikeEst(X, y)
    
    # Null Modell
    nullModell <- maxLikeEst(X = matrix(rep(1, times = nrow(X)), ncol = 1), y = y)
    
    result$nullModell <- nullModell
    result$formula <- formula
    result$call <- match.call()
    result$X <- X
    result$y <- y
    
    class(result) <- "logitMod"
    
    return(result)
    
}
logitModell <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)


# Funktion für die Print-Methode der Klasse logitMod
print.logitMod <- function(x, ...){
    
    cat("Call: ", paste0(deparse(x$call)), fill = TRUE)
    
    cat("\n\nCoefficients:\n")
    
    print.default(format(coef(x)[,1], digits = 4L),
                  print.gap = 1L, quote = FALSE, right = TRUE)
    
    cat("\nDegrees of Freedom: ", x$dfNull, " Total (i.e. Null); ",
        x$dfRes, " Residual")
    
    # Berechnung von null deviance, residual deviance & aic
    nullDeviance <- sum(x$nullModell$devianceResidual^2)
    x$nullDeviance <- nullDeviance
    devianceResidual <- sum(x$devianceResidual^2)
    x$devianceResidual <- devianceResidual
    x_AIC <- (-2*x$maxLogLikeValue + 2*ncol(x$X))
    x$AIC <- x_AIC
    
    cat("\nNull Deviance:\t", round(nullDeviance,1))
    cat("\nResidual Deviance:", round(devianceResidual,1), "\t", 
        "AIC: ", round(x_AIC,1))
    
    # invisibly return linMod object
    invisible(x)
    
}
logitModell <- logitMod(admit ~ gre + gpa + rank, testData)
print(logitModell)


# Funktion für die Summary-Methode der Klasse logitMod
summary.logitMod <- function(x, ...) {
    
    # Koeffizienten Standardfehler
    betaStandardError <- as.matrix(sqrt(diag(x$vcov)))
    x$betaStandardError <- betaStandardError
    
    # z-Statistik
    zStat <- x$coefficients / betaStandardError
    x$zStat <- zStat
    
    # p-Werte
    pValue <- 2 * pnorm(-abs(zStat))
    x$pValue <- pValue
    
    # sigCode <- pValue
    # sigCode[0 <= sigCode && sigCode < 0.001] <- "***"
    # sigCode[0.001 <= sigCode && sigCode < 0.01] <- "**"
    # sigCode[0.01 <= sigCode && sigCode < 0.5] <- "*"
    # sigCode[0.05 <= sigCode && sigCode < 0.1] <- "."
    # sigCode[0.1 <= sigCode && sigCode <= 1] <- ""
    # x$sigCode <- sigCode
    
    # Zusammenfassung der Werte für die Koeffizienten
    x$coefficients <- cbind("Estimate" = x$coefficients[,],
                            "Std. error" = x$betaStandardError[,],
                            "z value" = x$zStat[,],
                            "Pr(>|z|)" = x$pValue[,])
    #" " = x$sigCode[,])
    
    # Berechnung von nullDeviance, residualDeviance & aic
    nullDeviance <- sum(x$nullModell$devianceResidual^2)
    x$nullDeviance <- nullDeviance
    residualDeviance <- sum(x$devianceResidual^2)
    x$residualDeviance <- residualDeviance
    x_AIC <- (-2*x$maxLogLikeValue + 2*ncol(x$X))
    x$AIC <- x_AIC
    
    class(x) <- "summary.logitMod"
    
    return(x)
    
}
print.summary.logitMod <- function(x, ...) {
    
    cat("Call: ", deparse(x$call), fill = TRUE)
    
    cat("\nDeviance Residuals:\n")
    print.summaryDefault(summary(x$devianceResidual)[-4], digits = 4L)
    
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, signif.legend = TRUE, digits = 4L)
    
    cat("\n(Dispersion parameter for binomial family taken to be 1)\n")
    
    cat("\n    Null deviance: ", 
        round(x$nullDeviance,2), " on ", x$dfNull, " degrees of freedom\n")
    cat("Residual deviance: ", 
        round(x$residualDeviance,2), " on ", x$dfRes, " degrees of freedom\n")
    cat("AIC: ", round(x$AIC,2))
    
    cat("\n\nNumber of Fisher Scoring iterations: ", "\n")
    
    # invisibly return summary 
    invisible(x)
}


# Funktion für die Plot-Methode der Klasse logitMod
