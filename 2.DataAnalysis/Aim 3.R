library(randomForest)
require(caTools)
library(dplyr)
library(tidyr)
library(caret)

#Data
data <- read.csv("/Users/michael/Desktop/analyticData.csv", header = T)
data.rf <- data[c(3:7,10:40,44)] 


#Grouping Hypertension, Diabetes and Renal Disease
data.rf$anyHtn <- ifelse(data.rf$HYPERTENSION.COMPLICATED == 1
                          |data.rf$HYPERTENSION.UNCOMPLICATED == 1, 1, 0)

data.rf$anyDM <- ifelse(data.rf$DIABETES.COMPLICATED == 1|data.rf$DIABETES == 1
                       |data.rf$DIABETES.UNCOMPLICATED == 1, 1, 0)

data.rf$anyRD <- ifelse(data.rf$RENAL.DISEASE.CKD == 1|data.rf$RENAL.DISEASE.COMPLICATED == 1
                       |data.rf$RENAL.DISEASE.ESRD == 1|data.rf$RENAL_DISEASE == 1, 1, 0)


#Binary maxRutherfordClass to 5-6 and 4
data.rf$maxRutherfordClass_bin <- ifelse(data.rf$maxRutherfordClass >= 5, 1, 0)


#Grouping ischemia to 0-2 and 3
data.rf$ischemia_bin <- ifelse(data.rf$ischemia > 2, 1, 0)


#Binary for PATIENT_RACE_1 
data.rf$PATIENT_RACE_bin <- ifelse(data.rf$PATIENT_RACE_1 == 'WHITE OR CAUCASIAN', 1, 0)


#Binarize woundClass
data.rf$woundClass_0 <- ifelse(data.rf$woundClass == 0, 1, 0)
data.rf$woundClass_1 <- ifelse(data.rf$woundClass == 1, 1, 0)
data.rf$woundClass_2 <- ifelse(data.rf$woundClass == 2, 1, 0)
data.rf$woundClass_3 <- ifelse(data.rf$woundClass == 3, 1, 0)


#Adjusting Variables
names <- c(6:9,13:47)
data.rf[,names] <- lapply(data.rf[,names], factor)
data.rf$ischemia <- factor(data.rf$ischemia, ordered = T)
data.rf$woundClass <- factor(data.rf$woundClass, ordered = T) 
data.rf$maxRutherfordClass <- factor(data.rf$maxRutherfordClass, ordered = T)
str(data.rf)


#Separating into Two Datasets Based on Primary Group
data.WM <- data.rf[data.rf[,5] == 'WM',]
data.ReV <- data.rf[data.rf[,5] == 'Revasc',]







#Training Set(s)
set.seed(75)
sample.WM <- sample.split(data.WM$ampFreeSurv2yr, SplitRatio = 0.9)
train.WM <- subset(data.WM, sample.WM == TRUE)
test.WM <- subset(data.WM, sample.WM == FALSE)

sample.ReV <- sample.split(data.ReV$ampFreeSurv2yr, SplitRatio = 0.9)
train.ReV <- subset(data.ReV, sample.ReV == TRUE)
test.ReV <- subset(data.ReV, sample.ReV == FALSE)


#Random Forest(s)
WM.1 <-train.WM[c(1:4,6:37)]    #With Original Variables
ReV.1 <-train.ReV[c(1:4,6:37)]
WM.2 <-train.WM[c(1,4,6:9,13:16,23:24,29:47)]    #With Adjustments
ReV.2 <-train.ReV[c(1,4,6:9,13:16,23:24,29:47)]

rf.WM.1 <- randomForest(ampFreeSurv2yr~ ., data = WM.1, na.action = na.pass, mtry = 6)
rf.ReV.1 <- randomForest(ampFreeSurv2yr~ ., data = ReV.1, na.action = na.pass, mtry = 6)
rf.WM.2 <- randomForest(ampFreeSurv2yr~ ., data = WM.2, na.action = na.pass, mtry = 6)
rf.ReV.2 <- randomForest(ampFreeSurv2yr~ ., data = ReV.2, na.action = na.pass, mtry = 6)


#Importance
varImpPlot(rf.WM.1, sort=T)
varImpPlot(rf.WM.2, sort=T)
varImpPlot(rf.ReV.1, sort=T)
varImpPlot(rf.ReV.2, sort=T)


#Predictions
pred.WM.1 <- predict(rf.WM.1, test.WM)
pred.WM.2 <- predict(rf.WM.2, test.WM)
pred.ReV.1 <- predict(rf.ReV.1, test.ReV)
pred.ReV.2 <- predict(rf.ReV.2, test.ReV)


#Confusion Matrix (True Pos and True Neg)
cm.WM.1 <- table(test.WM[,37], pred.WM.1)
cm.WM.2 <- table(test.WM[,37], pred.WM.2)
cm.ReV.1 <- table(test.ReV[,37], pred.ReV.1)
cm.ReV.2 <- table(test.ReV[,37], pred.ReV.2)







####################
#randomforest Tuning
####################
#Comorbidities As Originally Defined

#WM Group
X1.WM <- train.WM[c(1:36)]
Y1.WM <- train.WM$ampFreeSurv2yr
rf.tune1.WM <- tuneRF(X1.WM, Y1.WM, ntreeTry = 500, stepFactor = 1.5, improve = 0.05, trace = T, 
                      plot = T, doBest = T)
print(rf.tune1.WM)
plot(rf.tune1.WM)

#ReVasc Group
X1.ReV <- train.ReV[c(1:36)]
Y1.ReV <- train.ReV$ampFreeSurv2yr
rf.tune1.ReV <- tuneRF(X1.ReV, Y1.ReV, ntreeTry = 500, stepFactor = 1.5, improve = 0.05, trace = T, 
                       plot = T, doBest = T)
print(rf.tune1.ReV)
plot(rf.tune1.ReV)


#With Grouping and Adjustment (ETHNICITY and DEMENTIA.NA (<10% prev) dropped)

#WM Group
X2.WM <- train.WM[c(1,4:9,13:16,23,24,29:36,38:46)]
Y2.WM <- train.WM$ampFreeSurv2yr
rf.tune2.WM <- tuneRF(X2.WM, Y2.WM, ntreeTry = 500, stepFactor = 1.5, improve = 0.05, trace = T, 
                      plot = T, doBest = T)
print(rf.tune2.WM)
plot(rf.tune2.WM)


#ReVasc Group
X2.ReV <- train.ReV[c(1,4:9,13:16,23,24,29:36,38:46)]
Y2.ReV <- train.ReV$ampFreeSurv2yr
rf.tune2.ReV <- tuneRF(X2.ReV, Y2.ReV, ntreeTry = 500, stepFactor = 1.5, improve = 0.05, trace = T, 
                       plot = T, doBest = T)
print(rf.tune2.ReV)
plot(rf.tune2.ReV)

##############
#CARET Tuning
##############
#Comorbidities As Originally Defined

#WM Group
control.WM.1 <- trainControl(method = "repeatedcv", number = 5, repeats = 3, search = "random", p = 1)

rf_caret <- train(ampFreeSurv2yr~ .,data = WM.1, method = "rf", metric = 'Accuracy', 
                  tuneLength = 10, trControl = control.WM.1)
print(rf_caret)
plot(rf_caret)

#ReVasc Group
control.ReV.1 <- trainControl(method = "repeatedcv", number = 5, repeats = 3, search = "random", p = 1)
rf_caret <- train(ampFreeSurv2yr~ .,data = ReV.1, method = "rf", metric = 'Accuracy', 
                  tuneLength = 10, trControl = control.ReV.1)
print(rf_caret)
plot(rf_caret)


#With Grouping and Adjustment (ETHNICITY and DEMENTIA.NA (<10% prev) dropped)

#WM Group
control.WM.2 <- trainControl(method = "repeatedcv", number = 5, repeats = 3, search = "random", p = 1)
rf_caret <- train(ampFreeSurv2yr~ .,data = WM.2, method = "rf", metric = 'Accuracy', 
                  tuneLength = 10, trControl = control.WM.2)
print(rf_caret)
plot(rf_caret)


#ReVasc Group
control.ReV.2 <- trainControl(method = "repeatedcv", number = 5, repeats = 3, search = "random", p = 1)
ReV.2 <-train.ReV[c(1:46)]
rf_caret <- train(ampFreeSurv2yr~ .,data = ReV.2, method = "rf", metric = 'Accuracy', 
                  tuneLength = 10, trControl = control.ReV.2)
print(rf_caret)
plot(rf_caret)




