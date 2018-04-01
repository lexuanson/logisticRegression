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
    
    # falls y nicht als 0-1/Variable eingegeben wird
    if (!(0 %in% y && 1 %in% y)) {
        y <- factor(y, labels = c(0,1))
    }
    y <- as.numeric(as.character(y))
    
    # initialisiere beta
    beta <- rep(0, times = ncol(X))
    
    # initialisiere ...
    W <- diag(nrow = nrow(X))
    
    # setze Abbruchskriterien
    tolerance <- exp(-6)
    diff <- 10 * abs(tolerance)
    maxIteration <- 100
    i <- 0

    # Solange Abbruchskriterien nicht erreicht
    while (diff > tolerance) {
        
        # berechne die Wahrscheinlichkeit
        eta <- X %*% beta
        exp_eta <- exp(eta)
        p <- as.vector(exp_eta / (1 + exp_eta))
        
        # aktualisiere ...
        W <- diag(p * (1 - p))
        
        # berechne die Änderung von Beta
        betaChange <- solve(t(X) %*% W %*% X) %*% t(X) %*% (y - p)
        
        # aktualisiere Beta
        beta <- beta + betaChange
        
        # berechne die totale Änderung von Beta
        diff <- abs(sum(betaChange))
        
        # nächste Iteration
        i <- i + 1
        
        if (i > maxIteration) {
            stop("Der Algorithmus konvergiert nicht!")
        }
    }
    
    # Freiheitsgrad = Anzahl an Beobachtungen - Anzahl an Parameter
    dfRes <- nrow(X) - ncol(X)  # ResiduensFreiheitsgrad 
    dfNull <- nrow(X) - 1       # Nullmodell Freiheitsgrad
    
    # Devianz Residual
    s <- y
    s[s == 0] = -1
    devianceResiduals = as.numeric(s * sqrt(-2*((y * eta) - (log(1 + exp(eta))))))
    #devianceResiduals = -2 * as.numeric(crossprod(y,eta) - sum(log(1 + exp(eta))))    
    
    # Kovarianzmatrix
    vcov <- solve(t(X) %*% W %*% X) 
    
    # Maximumwert der Log Likelihood Funktion
    maxLogLikeValue <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))
    
    # Liste der zurückgegebenen Werte
    result <- list(coefficients = beta,
                   vcov = vcov,
                   devianceResiduals = devianceResiduals,
                   dfRes = dfRes,
                   dfNull = dfNull,
                   maxLogLikeValue = maxLogLikeValue,
                   anzahlIteration = i,
                   p = p,
                   W = W) 
    
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
#' @importFrom stats model.frame model.matrix model.response 
#' @export
logitMod <- function(formula, data) {
    
    # initialisiere y (Zielvariable) und X (Matrix enthält alle erklärenden Variablen)
    modelFrame <- model.frame(formula, data)
    X <- model.matrix(formula, modelFrame)
    y <- model.response(modelFrame)
    
    # falls y nicht als 0-1/Variable eingegeben wird
    if (!(0 %in% y && 1 %in% y)) {
        y <- factor(y, labels = c(0,1))
    }
    y <- as.numeric(as.character(y))
    
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
    
    # Berechnung von nullDeviance, residualDeviance & aic
    nullDeviance <- -2 * result$nullModell$WaxLogLikeValue
    result$nullDeviance <- nullDeviance
    residualDeviance <- -2 * result$WaxLogLikeValue
    result$residualDeviance <- residualDeviance
    x_AIC <- (-2*result$WaxLogLikeValue + 2*ncol(result$X))
    result$AIC <- x_AIC
    
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
#' @importFrom stats coef
#' @export
print.logitMod <- function(x, ...){
    
    cat("Call: ", paste0(deparse(x$call)), fill = TRUE)
    
    cat("\nCoefficients:\n")
    
    print.default(format(coef(x)[,1], digits = 4L),
                  print.gap = 1L, quote = FALSE, right = TRUE)
    
    cat("\nDegrees of Freedom: ", x$dfNull, " Total (i.e. Null); ",
        x$dfRes, " Residual")
    
    cat("\nNull Deviance:\t", round(x$nullDeviance,1))
    cat("\nResidual Deviance:", round(x$residualDeviance,1), "\t", 
        "AIC: ", round(x$AIC,1), "\n") 
    
    # invisibly return linMod object
    invisible(x)
    
}


#' @title Summary method for "logitMod" estimations
#' @description Summary method for class "logitMod"
#' @param object an object of class "logitMod"
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
#' @importFrom stats pnorm
#' @export
summary.logitMod <- function(object, ...) {
    
    # Koeffizienten Standardfehler
    object$betaStandardError  <- as.matrix(sqrt(diag(object$vcov)))
    
    # z-Statistik
    zStat <- object$coefficients / object$betaStandardError
    object$zStat <- zStat
    
    # p-Werte
    pValue <- 2 * pnorm(-abs(zStat))
    object$pValue <- pValue
    
    # Zusammenfassung der Werte für die Koeffizienten
    object$coefficients <- cbind("Estimate" = object$coefficients[,],
                            "Std. error" = object$betaStandardError[,],
                            "z value" = object$zStat[,],
                            "Pr(>|z|)" = object$pValue[,])
    
    class(object) <- "summary.logitMod"
    
    return(object)
    
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
#' @importFrom stats printCoefmat 
#' @export
print.summary.logitMod <- function(x, ...) {
    
    cat("Call: ", deparse(x$call), fill = TRUE)
    
    cat("\nDeviance Residuals:\n")
    print.summaryDefault(summary(x$devianceResiduals)[-4], digits = 4L)
    
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, signif.legend = TRUE, digits = 4L)
    
    cat("\n(Dispersion parameter for binomial family taken to be 1)\n")
    
    cat("\n    Null deviance: ", 
        round(x$nullDeviance,2), " on ", x$dfNull, " degrees of freedom\n")
    cat("Residual deviance: ", 
        round(x$residualDeviance,2), " on ", x$dfRes, " degrees of freedom\n")
    cat("AIC: ", round(x$AIC,2))
    
    cat("\n\nNumber of Fisher Scoring iterations: ", x$anzahlIteration, "\n")
    
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
#' @importFrom graphics plot abline 
#' @importFrom stats qqnorm qqline
#' @export
plot.logitMod <- function(x, ...) {
    
    #1
    plot(y = x$devianceResiduals, x = (x$X %*% x$coefficients),
         main = "Residuals vs Fitted",
         ylab = "Residuals",
         xlab = paste("Predicted Values\n",
                      deparse(x$call)))
    abline(a = 0, b = 0, lty = 3)
    
    #2
    qqnorm(x$devianceResiduals,
       main = "Normal Q-Q",
       ylab = "Std. deviance resid.",
       xlab = paste("Theoretical Quantiles\n", deparse(x$call)))
    qqline(x$devianceResiduals, lty = 3)

    #3
    plot(y = sqrt(abs(x$devianceResiduals)), x = (x$X %*% x$coefficients),
     main = "Scale Location",
     ylab = expression(sqrt("|Std. deviance resid.|")),
     xlab = paste("Predicted Values\n", deparse(x$call)))
    
    
    #4
    pearsonResidual <- (x$y - x$p)/sqrt(x$p*(1 - x$p))
    leverage <- diag(sqrt(x$W) %*% x$X %*% (solve(t(x$X) %*% x$W %*% x$X)) %*% t(x$X) %*% sqrt(x$W))
    plot(y = pearsonResidual, 
         x = leverage,
         main = "Residual vs Leverage",
         ylab = "Std. Pearson resid.",
         xlab = paste("Leverage\n", deparse(x$call)))

}

