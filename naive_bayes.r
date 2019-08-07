#Getting started with Naive Bayes
#Install the package
#install.packages(“e1071”)
#Loading the library
library(caret)
library(e1071)
library(DMwR)

#Se cargan datos con la última columna modifcada (1 -> YES, -1 -> NO)
datos = read.csv("nci60_binary_class_training_set.csv")
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

Naive_Bayes_Model=naiveBayes(Labels ~ LC.NCI.H522 + BR.MCF7 + LE.HL.60.TB. + CO.KM12 + RE.TK.10+
                               LE.CCRF.CEM + OV.SK.OV.3 + CNS.SF.268 + RE.CAKI.1 + LC.NCI.H460 + OV.OVCAR.4
                             + OV.IGROV1 + CO.HCT.15 + ME.UACC.62 + CO.SW.620 + RE.SN12C + RE.786.0 + OV.OVCAR.3+
                               LE.K.562 + BR.BT.549 + RE.UO.31 + ME.M14 + OV.NCI.ADR.RES + ME.SK.MEL.2 + ME.MALME.3M
                             + ME.SK.MEL.28 + BR.MDA.MB.231 + ME.MDA.N + LC.HOP.92 + BR.HS.578T + RE.ACHN + OV.OVCAR.8+
                               RE.A498 + CNS.SNB.19 + CO.HT29 + ME.SK.MEL.5,
                             data=training_set)

summary(Naive_Bayes_Model)
NB_Predictions=predict(Naive_Bayes_Model,training_set)
confusionMatrix(NB_Predictions,training_set$Labels)

NB_Predictions_test=predict(Naive_Bayes_Model,test_set)
confusionMatrix(NB_Predictions_test,test_set$Labels)



#NB, probando modelo con datos_testset sin clasificar
#Predict Output 
datos_testset = read.csv("nci60_non_coding_test_set.csv")
datos_testset$ensembl_gene_id <- NULL
datos_testset$gene_biotype <- NULL
datos_testset$X <- NULL
datos_testset$entrezgene <- row.names(datos_testset$entrezgene)
datos_testset$entrezgene <- NULL
datos_testset$chromosome_name <- NULL
datos_testset$LC.NCI.H23 <-NULL
datos_testset$hgnc_id<-NULL
datos_testset <- na.omit(datos_testset)
datos_testset <- log(datos_testset, 2)
summary(predict(Naive_Bayes_Model , newdata = datos_testset, interval = 'prediction'))

