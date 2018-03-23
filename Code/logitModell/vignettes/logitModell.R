## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----Preparing data------------------------------------------------------
#install.packages("logitModell")
library(logitModell)

testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
testData$rank <- factor(testData$rank)
testModell <- as.formula("admit ~ gre + gpa + rank")
testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)

# Alternatively one can use the following code to test the result on 
# another data set. This daat set is default in R.

# data(cats, package = "MASS")
# testData <- cats[1:100,]
# testData$Sex <- factor(testData$Sex, labels = c(0,1))
# #testData$Sex <- as.numeric(as.character(testData$Sex))
# testModell <- as.formula("Sex ~ Bwt + Hwt")
# testModelFrame <- model.frame(Sex ~ Bwt + Hwt, testData)

## ----Fit logit model-----------------------------------------------------
fittedLogModell <- logitMod(formula = testModell, data = testData)

## ----Print method--------------------------------------------------------
fittedLogModell

## ----Summary method------------------------------------------------------
summary(fittedLogModell)

## ----Plot method---------------------------------------------------------
plot(fittedLogModell)

