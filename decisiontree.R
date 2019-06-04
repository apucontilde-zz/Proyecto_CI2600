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
datos <- read.csv("C:/Users/oscar/Documents/Aprendizaje Automatico CI2600/proyecto iteracion 2/nci60_binary_class_training_set.csv")
datos_testset <- read.csv("C:/Users/oscar/Documents/Aprendizaje Automatico CI2600/proyecto iteracion 2/nci60_non_coding_test_set.csv")
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
datos_testset$hgnc_id <- NULL
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

str(datos_testset)








#CART, classification
set.seed(101)
model <- rpart(Labels ~ .,data= datos, method='class')

plotcp(model)
summary(model) 

plot(model, uniform=TRUE,
     main="Classification Tree for Labels")
text(model, use.n=TRUE, all=TRUE, cex=.8)

#prediccion sobre particion de prueba

model.pred <- predict(model,type='class')

mean(model.pred == datos$Labels)

cart.conf.mat <- confusionMatrix(model.pred, datos$Labels)

cart.conf.mat


#CART, classification, probando modelo con datos sin clasificar

fit <- rpart(Labels ~ ., data = datos,method="class")  #grow tree 
summary(fit)
#Predict Output 
predict(fit , newdata = datos_testset, interval = 'prediction')

