## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----Preparing data------------------------------------------------------
#install.packages("logitModell")
library(logitModell)

testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
testData <- testData[1:100,]
testData$rank <- factor(testData$rank)
testFormula <- as.formula("admit ~ gre + gpa + rank")
testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)

# Alternatively one can use the following code to test the result on another data set. This data set is default built in R.

# data(cats, package = "MASS")
# testData <- cats[1:100,]
# #testData$Sex <- factor(testData$Sex, labels = c(0,1))
# testFormula <- as.formula("Sex ~ Bwt + Hwt")
# testModelFrame <- model.frame(Sex ~ Bwt + Hwt, testData)

## ----Fit logit model-----------------------------------------------------
fittedLogModell <- logitMod(formula = testFormula, data = testData)
standardModell <- glm(formula = testFormula, data = testData, family = "binomial")

## ----echo = FALSE--------------------------------------------------------
all(all.equal(standardModell$coefficients,
              as.numeric(fittedLogModell$coefficients),
              check.attributes = FALSE),
    all.equal(vcov(standardModell), 
              fittedLogModell$vcov,
              check.attributes = FALSE, tolerance = exp(-10)),
    all.equal(residuals(standardModell),
              as.numeric(fittedLogModell$devianceResidual),
              check.attributes = FALSE, tolerance = exp(-10)),
    all.equal(standardModell$df.residual, fittedLogModell$dfRes),
    all.equal(standardModell$df.null, fittedLogModell$dfNull))

## ----Print method--------------------------------------------------------
print(fittedLogModell)

## ----Summary method------------------------------------------------------
summary(fittedLogModell)

## ----Plot method, fig.width = 4, fig.height = 3, fig.align = "center"----
plot(fittedLogModell)

