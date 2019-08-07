library(ggplot2)
library(ggtheme)
library(e1071)
library(stringr)
library(rpart)
library(rattle)
library(caret)
library(plyr)
library(Boruta)
library(UBL)
library(Metrics)
library(dplyr)
library(mlbench)
library(DMwR)
#============================================
#PREPROCESAMIENTO DE DATOS
#============================================
#cargo set de datos
datos <- read.csv("C:/Users/oscar/Documents/Aprendizaje Automatico CI2600/ultima entrega/nci60_binary_class_training_set.csv")
#View(datos)
str(datos) # busco columnas que puedan crear conflicto(que no son numero) con rowMeans()
datos$Labels=as.factor(datos$Labels)
# se quitan variables que sobran
datos$ensembl_gene_id <- NULL
datos$gene_biotype <- NULL
datos$X <- NULL
datos$entrezgene <- row.names(datos$entrezgene)
datos$entrezgene <- NULL
datos$chromosome_name <- NULL
datos$LC.NCI.H23 <-NULL
levels(datos$Labels) <- c("non_driver", "driver")

# se hace un shuffle de los elementos
set.seed(156) 
rand <- sample(nrow(datos))
rand
datos<-datos[rand, ]
datos <- na.omit(datos) # se omiten filas con NA
str(datos)
#log2?
datos[, 1:59] <- log(datos[1:59], 2)

#============================================
# PARTICIONAMIENTO
#============================================
# setting the seed for reproducible results
set.seed(420) 
# separate the partitions between training and test set
trainIndex <- createDataPartition(datos$Labels, 
                                  p = .8, 
                                  list = FALSE, 
                                  times = 1)

training_set <- datos[ trainIndex, ]
test_set <- datos[ -trainIndex, ]
#BALANCEO
#============================================
# Balanceo SMOTE, OVER 100%
#============================================
smoted_data <- SMOTE(Labels~., training_set, perc.over=100)

training_set <- smoted_data

#FEATURE SELECTION
#============================================
#RFE (~30 minutos) https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
#============================================ 
# ensure the results are repeatable
set.seed(7)

# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(datos[,1:59], datos[,60], sizes=c(1:59), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))
#=============================================
# Feature selection lvq
#=============================================
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Labels~., data=datos, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)
#============================================
#highly correlated https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
#============================================ 
# ensure the results are repeatable
set.seed(7)
# calculate correlation matrix
correlationMatrix <- cor(datos[,1:59])
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)











#modelos:



#svm de internet, https://dataaspirant.com/2017/01/19/support-vector-machine-classifier-implementation-r-caret-package/
#=========================================
#SVM RADIAL
#=========================================
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

set.seed(3233)
svm_Radial <- train(Labels ~ LC.NCI.H522 + BR.MCF7 + LE.HL.60.TB. + CO.KM12 + RE.TK.10+
                      LE.CCRF.CEM + OV.SK.OV.3 + CNS.SF.268 + RE.CAKI.1 + LC.NCI.H460 + OV.OVCAR.4
                    + OV.IGROV1 + CO.HCT.15 + ME.UACC.62 + CO.SW.620 + RE.SN12C + RE.786.0 + OV.OVCAR.3+
                      LE.K.562 + BR.BT.549 + RE.UO.31 + ME.M14 + OV.NCI.ADR.RES + ME.SK.MEL.2 + ME.MALME.3M
                    + ME.SK.MEL.28 + BR.MDA.MB.231 + ME.MDA.N + LC.HOP.92 + BR.HS.578T + RE.ACHN + OV.OVCAR.8+
                      RE.A498 + CNS.SNB.19 + CO.HT29 + ME.SK.MEL.5, 
                    data = training_set, method = "svmRadial",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)
svm_Radial
plot(svm_Radial)

# RESULTADO SET DE ENTRENAMIENTO
test_pred_Radial1 <- predict(svm_Radial, newdata = training_set)
confusionMatrix(test_pred_Radial1, training_set$Labels)
# RESULTADO SET DE PRUEBA
test_pred_Radial2 <- predict(svm_Radial, newdata = test_set)
confusionMatrix(test_pred_Radial2, test_set$Labels)



# grid search, dura bastante
grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                           C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                                 1, 1.5, 2,5))
set.seed(3233)

svm_Radial_Grid <- train(Labels ~LC.NCI.H522 + BR.MCF7 + LE.HL.60.TB. + CO.KM12 + RE.TK.10+
                           LE.CCRF.CEM + OV.SK.OV.3 + CNS.SF.268 + RE.CAKI.1 + LC.NCI.H460 + OV.OVCAR.4
                         + OV.IGROV1 + CO.HCT.15 + ME.UACC.62 + CO.SW.620 + RE.SN12C + RE.786.0 + OV.OVCAR.3+
                           LE.K.562 + BR.BT.549 + RE.UO.31 + ME.M14 + OV.NCI.ADR.RES + ME.SK.MEL.2 + ME.MALME.3M
                         + ME.SK.MEL.28 + BR.MDA.MB.231 + ME.MDA.N + LC.HOP.92 + BR.HS.578T + RE.ACHN + OV.OVCAR.8+
                           RE.A498 + CNS.SNB.19 + CO.HT29 + ME.SK.MEL.5,
                         data = training_set, method = "svmRadial",
                         trControl=trctrl,
                         preProcess = c("center", "scale"),
                         tuneGrid = grid_radial,
                         metric = "Kappa",
                         tuneLength = 10)

svm_Radial_Grid



# RESULTADO SET DE ENTRENAMIENTO
test_pred_Radial_Grid2 <- predict(svm_Radial_Grid, newdata = training_set)
confusionMatrix(test_pred_Radial_Grid2, training_set$Labels )
# RESULTADO SET DE PRUEBA
test_pred_Radial_Grid1 <- predict(svm_Radial_Grid, newdata = test_set)
confusionMatrix(test_pred_Radial_Grid1, test_set$Labels )




#=========================================
#SVM LINEAR
#=========================================
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

set.seed(3233)
svm_Radial <- train(Labels ~ LC.NCI.H522 + BR.MCF7 + LE.HL.60.TB. + CO.KM12 + RE.TK.10+
                      LE.CCRF.CEM + OV.SK.OV.3 + CNS.SF.268 + RE.CAKI.1 + LC.NCI.H460 + OV.OVCAR.4
                    + OV.IGROV1 + CO.HCT.15 + ME.UACC.62 + CO.SW.620 + RE.SN12C + RE.786.0 + OV.OVCAR.3+
                      LE.K.562 + BR.BT.549 + RE.UO.31 + ME.M14 + OV.NCI.ADR.RES + ME.SK.MEL.2 + ME.MALME.3M
                    + ME.SK.MEL.28 + BR.MDA.MB.231 + ME.MDA.N + LC.HOP.92 + BR.HS.578T + RE.ACHN + OV.OVCAR.8+
                      RE.A498 + CNS.SNB.19 + CO.HT29 + ME.SK.MEL.5, 
                    data = training_set, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)
svm_Radial

# RESULTADO SET DE ENTRENAMIENTO
test_pred_Radial1 <- predict(svm_Radial, newdata = training_set)
confusionMatrix(test_pred_Radial1, training_set$Labels)
# RESULTADO SET DE PRUEBA
test_pred_Radial2 <- predict(svm_Radial, newdata = test_set)
confusionMatrix(test_pred_Radial2, test_set$Labels)




#=========================================

#SVM, classification, probando modelo con datos sin clasificar

#Predict Output 
predict(svm_Radial_Grid , newdata = datos_testset, interval = 'prediction')
