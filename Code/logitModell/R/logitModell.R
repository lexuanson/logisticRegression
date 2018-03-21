#' @title Maximum Likelihood
#' @description This function calculates the maximum likelihood for binary logistic regression
#' @param y a matrix/vector containing the dependent variable
#' @param X a matrix containing the intercept and all independent variables
#' @return a list with maximum likelihood estimation results
#' @examples  
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' maxLikeEst(X = model.matrix(testModell, testModelFrame), y = model.response(testModelFrame))
#' @export
maxLikeEst <- function(y, X) {
    
    # initialisiere beta
    beta <- rep(0, times = ncol(X))
    
    # initialisiere ...
    M <- diag(nrow = nrow(X))
    
    # setze Abbruchskriterien
    tolerance <- exp(-10)
    diff <- 10 * abs(tolerance)
    maxIteration <- 1000
    i <- 0
    
    # Solange Abbruchskriterien nicht erreicht
    while (diff > tolerance || i < maxIteration) {
        
        # berechne die Wahrscheinlichkeit
        eta <- X %*% beta
        exp_eta <- exp(eta)
        p <- as.vector(exp_eta / (1 + exp_eta))
        
        # aktualisiere ...
        M <- diag(p * (1 - p))
        
        # berechne die Änderung von Beta
        beta_change <- solve(t(X) %*% M %*% X) %*% t(X) %*% (y - p)
        
        # aktualisiere Beta
        beta <- beta + beta_change
        
        # berechne die totale Änderung von Beta
        diff <- sum(abs(beta_change))
        
        # nächste Iteration
        i <- i + 1
    }
    
    # Freiheitsgrad = Anzahl an Beobachtungen - Anzahl an Parameter
    dfRes <- nrow(X) - ncol(X)  # ResiduensFreiheitsgrad 
    dfNull <- nrow(X) - 1       # Nullmodell Freiheitsgrad
    
    # Kovarianzmatrix
    vcov <- solve(t(X) %*% M %*% X) 
    
    # Devianz Residual
    s <- y
    s[s == 0] = -1
    devianceResidual = as.numeric(s * sqrt(-2*((y * eta) - (log(1 + exp(eta))))))
    
    # Maximumwert der Log Likelihood Funktion
    maxLogLikeValue <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))
    
    # Liste der zurückgegebenen Werte
    result <- list(coefficients = beta,
                   vcov = vcov,
                   devianceResidual = devianceResidual,
                   dfRes = dfRes,
                   dfNull = dfNull,
                   maxLogLikeValue = maxLogLikeValue) 
    
    return(result)
    
}


#' @title Interface for an alternative logistic regression implementation
#' @description This function computes coefficients of a binary logistic regression
#' and contructs an object of "logitMod" class
#' @param formula a formula object. Matrix or array are not accepted
#' @param data a data frame contains all variables
#' @return a list of class "logitMod" containing the coefficients, a null model, formula, call,
#' the dependent variable and all independent variables.
#' @examples 
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' @export
logitMod <- function(formula, data) {
    
    # initialisiere y (Zielvariable) und X (Matrix enthält alle erklärenden Variablen)
    modelFrame <- model.frame(formula, data)
    X <- model.matrix(formula, modelFrame)
    y <- model.response(modelFrame)
    
    # erstelle eine Liste als Ergebnis der maxLikeEst-Funktion
    result <- maxLikeEst(y, X)
    
    # erstelle das Null Modell und speichere es in die Ergebnisliste
    nullModell <- maxLikeEst(X = matrix(rep(1, times = nrow(X)), ncol = 1), y = y)
    result$nullModell <- nullModell
    
    # speichere notwendige Parameter in die Ergebnisliste
    result$formula <- formula
    result$call <- match.call()
    result$X <- X
    result$y <- y
    
    # ordne die Ergebnisliste der Klasse "logitMod" zu
    class(result) <- "logitMod"
    
    return(result)
    
}



#' @title Printing method for logitMod estimations
#' @description Printing method for class "logitMod"
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return a standard print output equivalent to the built in binary logistic regression
#' @examples 
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' print(logm)
#' @export
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



#' @title Summary method for "logitMod" estimations
#' @description Summary method for class "logitMod"
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return a list of all necessary values equivalent to the summary output of the
#' built in binary logistic regression
#' @examples 
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' summary(logm)
#' @export
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



#' @title Printing method for the summary of "logitMod" estimations
#' @description Printing method for summary of class "logitMod"
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return an equivalent output to the summary of the built in binary logistic regression
#' @examples 
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' summary(logm)
#' @export
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



#' @title Plotting method for "logitMod" objects
#' @description Plotting method for objects of class "logitMod"
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return equivalent plots to those of the built in binary logistic regression
#' @examples 
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' plot(logm)
#' @export
plot.logitMod <- function(x, ...) {
    
    #1
    plot(y = x$devianceResidual, x = (x$X %*% x$coefficients),
         main = "Residuals vs Fitted",
         ylab = "Residuals",
         xlab = paste("Predicted Values\n",
                      deparse(x$call)))
    abline(a = 0, b = 0, lty = 3)
    
    #2
    qqnorm(x$devianceResidual, 
           main = "Normal Q-Q", 
           ylab = "Std. deviance resid.",
           xlab = paste("Theoretical Quantiles\n", deparse(x$call))) 
    qqline(x$devianceResidual, lty = 3)
    
    #3
    plot(y = sqrt(abs(x$devianceResidual)), x = (x$X %*% x$coefficients),
         main = "Scale Location",
         ylab = expression(sqrt("|Std. deviance resid.|")),
         xlab = paste("Predicted Values\n", deparse(x$call)))
    
}
