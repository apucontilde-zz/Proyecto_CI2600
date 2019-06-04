#Getting started with Naive Bayes
#Install the package
#install.packages(“e1071”)
#Loading the library
library(e1071)

#Se cargan datos con la última columna modifcada (1 -> YES, -1 -> NO)
BreastCancer_df = read.csv("nci60_binary_class_training_set_yn.csv")

Naive_Bayes_Model=naiveBayes(Labels~.,data=BreastCancer_df)

summary(Naive_Bayes_Model)
NB_Predictions=predict(Naive_Bayes_Model,BreastCancer_df)
summary(NB_Predictions)
#Confusion matrix to 
table(NB_Predictions,BreastCancer_df$Labels)