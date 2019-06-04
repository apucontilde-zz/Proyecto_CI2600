#Getting started with Naive Bayes
#Install the package
#install.packages(“e1071”)
#Loading the library
library(e1071)

#Se cargan datos con la última columna modifcada (1 -> YES, -1 -> NO)
BreastCancer_df = read.csv("nci60_binary_class_training_set.csv")

final_test = read.csv("nci60_non_coding_test_set.csv")

BreastCancer_df$LC.NCI.H23[is.na(BreastCancer_df$LC.NCI.H23)] <- 0 # lleno de 0's columna conflictiva

BreastCancer_df$Labels=as.factor(BreastCancer_df$Labels)

BreastCancer_df$ensembl_gene_id <- NULL
BreastCancer_df$gene_biotype <- NULL
BreastCancer_df$chromosome_name <- NULL



# setting the seed for reproducible results
set.seed(420) 
# separate the partitions between training and test set
trainIndex <- createDataPartition(datos$Labels, 
                                  p = .8, 
                                  list = FALSE, 
                                  times = 1)

training_set <- BreastCancer_df[ trainIndex, ]
test_set <- BreastCancer_df[ -trainIndex, ]

Naive_Bayes_Model=naiveBayes(Labels~.,data=training_set)

summary(Naive_Bayes_Model)
NB_Predictions=predict(Naive_Bayes_Model,test_set)
summary(NB_Predictions)
#Confusion matrix to 
table(NB_Predictions,test_set$Labels)
NB_Predictions_final=predict(Naive_Bayes_Model,final_test)
summary(NB_Predictions_final)

