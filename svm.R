library(ggplot2)
library(e1071)
library(stringr)
library(rpart)
library(rattle)
library(caret)
library(plyr)

#
#PREPROCESAMIENTO DE DATOS
#

#cargo set de datos
datos <- read.csv("nci60_binary_class_training_set.csv")
datos_testset <- read.csv("nci60_non_coding_test_set.csv")
#View(datos)

str(datos) # busco columnas que puedan crear conflicto(que no son numero) con rowMeans()

datos$LC.NCI.H23[is.na(datos$LC.NCI.H23)] <- 0 # lleno de 0's columna conflictiva


datos$Labels=as.factor(datos$Labels)

datos$ensembl_gene_id <- NULL
datos$gene_biotype <- NULL
#datos$chromosome_name <- NULL

datos_testset$ensembl_gene_id <- NULL
datos_testset$gene_biotype <- NULL
#datos_testset$chromosome_name <- NULL
#datos_testset$hgnc_id <- NULL
datos_testset$LC.NCI.H23[is.na(datos_testset$LC.NCI.H23)] <- 0 # lleno de 0's columna conflictiva

datos <- na.omit(datos) # se omiten filas con NA
str(datos)


# setting the seed for reproducible results
set.seed(420) 
# separate the partitions between training and test set
trainIndex <- createDataPartition(datos$Labels, 
                                  p = .8, 
                                  list = FALSE, 
                                  times = 1)

training_set <- datos[ trainIndex, ]
test_set <- datos[ -trainIndex, ]

#
# FIN DE PREPROCESAMIENTO
#


#CASO 1 SVM del profe
# create svm model

svm_model <- svm(Labels~., data=training_set, type='C-classification',kernel='radial')

summary(svm_model)
# prediccion sobre training set
pred_train <-predict(svm_model,training_set)

mean(pred_train==training_set$Labels)

conf.mat.svm.train <- confusionMatrix(pred_train,training_set$Labels)

conf.mat.svm.train

# prediccion sobre test set
pred_test <-predict(svm_model,test_set)
summary(pred_test)
conf.mat.svm.test <- confusionMatrix(pred_test, test_set$Labels)

conf.mat.svm.test

# percentage of testset predicted correctly by svm
mean(pred_test==test_set$Labels)
#

#TUNE(grid search), dura bastante

svm.tune.test <- tune.svm(Labels~.,data = datos, gamma=c(0.001, 0.01, 0.1, 1, 5, 10, 20, 50, 100, 150, 200), cost = (1:2))

summary(svm.tune.test)
plot(svm.tune.test)
svm.tune.test$best.parameters





# CASO 2
#svm de internet, https://dataaspirant.com/2017/01/19/support-vector-machine-classifier-implementation-r-caret-package/
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

set.seed(3233)
svm_Radial <- train(Labels ~., data = training_set, method = "svmRadial",
                      trControl=trctrl,
                      preProcess = c("center", "scale"),
                      tuneLength = 10)
svm_Radial
plot(svm_Radial)

test_pred_Radial <- predict(svm_Radial, newdata = test_set)
confusionMatrix(test_pred_Radial, test_set$Labels)



# grid search, dura bastante
grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                           C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                                   1, 1.5, 2,5))
set.seed(3233)

svm_Radial_Grid <- train(Labels ~., data = training_set, method = "svmRadial",
                           trControl=trctrl,
                           preProcess = c("center", "scale"),
                           tuneGrid = grid_radial,
                           tuneLength = 10)

svm_Radial_Grid

test_pred_Radial_Grid <- predict(svm_Radial_Grid, newdata = test_set)
 
confusionMatrix(test_pred_Radial_Grid, test_set$Labels )



#SVM, classification, probando modelo con datos sin clasificar

#Predict Output 
predict(svm_Radial_Grid , newdata = datos_testset, interval = 'prediction')
