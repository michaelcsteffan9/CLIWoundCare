library(randomForest)
require(caTools)
library(dplyr)
library(tidyr)

#Data
data <- read.csv("/Users/michael/Desktop/analyticData.csv", header=T)
data_rf<-data[c(6:36,40)]
data_rf <- lapply(data_rf, factor)

#Imputation
set.seed(0408)
data.imp<-rfImpute(ampFreeSurv2yr~ .,data_rf)


#Omitting
data_rf<-data.frame(data_rf)
data.omit<-drop_na(data_rf)


#Training Set(s)
sample.imp = sample.split(data.imp$ampFreeSurv2yr, SplitRatio = .8)
train.imp = subset(data.imp, sample.imp == TRUE)
test.imp  = subset(data.imp, sample.imp == FALSE)

sample.omit = sample.split(data.omit$ampFreeSurv2yr, SplitRatio = .8)
train.omit = subset(data.omit, sample.omit == TRUE)
test.omit  = subset(data.omit, sample.omit == FALSE)


#Random Forest(s)
rf.imp <- randomForest(ampFreeSurv2yr~ .,data=train.imp)
rf.omit <- randomForest(ampFreeSurv2yr~ .,data=train.omit, na.action = na.pass)

importance(rf.imp)
importance(rf.omit)


#Predictions
pred.imp = predict(rf.imp, test.imp)
pred.omit = predict(rf.omit, test.omit)


#Confusion Matrix (True Pos and True Neg)
cm.imp = table(test.imp[,1], pred.imp)
cm.omit = table(test.omit[,1], pred.imp)


