library(logitModell)

testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
testData$rank <- factor(testData$rank)
testModell <- as.formula("admit ~ gre + gpa + rank")
testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#maxLikeEst(y = model.response(testModelFrame), X = model.matrix(testModell, testModelFrame))
