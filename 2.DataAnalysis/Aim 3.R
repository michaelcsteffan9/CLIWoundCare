library(randomForest)
require(caTools)
library(dplyr)
library(tidyr)
library(caret)

#Data
data <- read.csv("/Users/michael/Desktop/analyticData.csv", header=T)
data_rf<-data[c(3:7,10:40,44)]
data_rf <- lapply(data_rf, factor)
data_rf$ageAtIndexProcedure=as.integer(data_rf$ageAtIndexProcedure)
data_rf$maxRutherfordClass=factor(data_rf$maxRutherfordClass,ordered=T)
data_rf$ischemia=factor(data_rf$ischemia,ordered=T)
data_rf$woundClass=factor(data_rf$woundClass,ordered=T) 


#Omitting
data_rf<-data.frame(data_rf)
data.omit<-drop_na(data_rf)


#Training Set(s)
set.seed(0412)
sample.omit = sample.split(data.omit$ampFreeSurv2yr, SplitRatio = .8)
train.omit = subset(data.omit, sample.omit == TRUE)
test.omit  = subset(data.omit, sample.omit == FALSE)


#randomforest Tuning
X <- train.omit[,c(-37)]
Y <- train.omit$ampFreeSurv2yr
rf_tune<-tuneRF(X, Y, ntreeTry=500, stepFactor=1.25, improve=0.05, trace=T, plot=T,doBest=T)
print(rf_tune)
plot(rf_tune)

#CARET Tuning
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random", p=1)
rf_caret <- train(ampFreeSurv2yr~ .,data=train.omit, method="rf", metric='Accuracy', 
                   tuneLength=10, trControl=control)
print(rf_caret)
plot(rf_caret)


#Random Forest(s) and Importance
rf.omit <- randomForest(ampFreeSurv2yr~ .,data=train.omit, na.action = na.pass, mtry=6)
importance(rf.omit)


#Predictions
pred.omit = predict(rf.omit, test.omit)


#Confusion Matrix (True Pos and True Neg)
cm.omit = table(test.omit[,37], pred.omit)
