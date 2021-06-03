# The function splits ExpressionSet object (or matrix) to test and train set and create partition to apply SNFS

train_test <- function (Dataforanalyze, Class_info, train_ratio){
  
  if(class(Dataforanalyze)=="ExpressionSet") {
    
    for(i in 1:length(varLabels(Dataforanalyze))){
      if (varLabels(Dataforanalyze)[i]==Class_info){
        Class_index=i
      }
    }
    
    varLabels(Dataforanalyze)[Class_index]="type"
    
    
    
    realratio<-summary(Dataforanalyze$type)[1]/summary(Dataforanalyze$type)[2]
    if (realratio>=1){
      accepted_error=0.5
    }
    if ((realratio<1)&(realratio>0.1)){
      accepted_error=0.1
    }
    
    ## setting training size according to the given train_ratio
    smp_size <- floor(train_ratio * as.integer(ncol(Dataforanalyze)))
    
    ## select training samples from whole data
    train_ind <- sample(seq_len(as.integer(ncol(Dataforanalyze))), replace=FALSE, size = smp_size)
    train <- Dataforanalyze[,train_ind ]
    
    ## creating test data
    test <- Dataforanalyze[,-train_ind ]
    
    ## compare the train/test class ratio to the real data's and reproducing if it is not acceptable
    ttratio<-summary(train$type)[1]/summary(train$type)[2]
    if ((ttratio>realratio+accepted_error) | (ttratio+accepted_error<realratio)){
      train_ind <- sample(seq_len(as.integer(ncol(Dataforanalyze))), replace=FALSE, size = smp_size)
      train <- Dataforanalyze[,train_ind ]
      test <- Dataforanalyze[,-train_ind ]
      ttratio<-summary(train$type)[1]/summary(train$type)[2]
    }
    
    
    trainandtestlist<-list(train_ind, train, test, Dataforanalyze)
    names(trainandtestlist)<-c("Indices of training samples", "Training dataset", "Test dataset", "Data")
    
    return(trainandtestlist)
    
  }
  
  
  else {
    
    
    for(i in 1:length(rownames(Dataforanalyze))){
      if (rownames(Dataforanalyze)[i]==Class_info){
        Class_index=i
      }
    }
    
    rownames(Dataforanalyze)[Class_index]="type"
    
    Dataforanalyze=list(as.factor(Dataforanalyze[Class_index,]), Dataforanalyze[-Class_index,])
    names(Dataforanalyze)=c("type", "expression")
    
    
    realratio<-summary(Dataforanalyze$type)[[1]]/summary(Dataforanalyze$type)[[2]]
    if (realratio>=1){
      accepted_error=0.5
    }
    if ((realratio<1)&(realratio>0.1)){
      accepted_error=0.1
    }
    
    ## setting training size according to the given train_ratio
    smp_size <- floor(train_ratio * as.integer(ncol(Dataforanalyze$expression)))
    
    ## select training samples from whole data
    train_ind <- sample(seq_len(as.integer(ncol(Dataforanalyze$expression))), replace=FALSE, size = smp_size)
    train <- Dataforanalyze$expression[,train_ind]
    train_class<-Dataforanalyze$type[train_ind]
    
    ## creating test data
    if(train_ratio<1){
      test <- Dataforanalyze$expression[,-train_ind]
      test_class<-Dataforanalyze$type[-train_ind]
    }
    
    ## compare the train/test class ratio to the real data's and reproducing if it is not acceptable
    ttratio<-summary(train_class)[[1]]/summary(train_class)[[2]]
    
    
    if ((ttratio>realratio+accepted_error) | (ttratio+accepted_error<realratio)){
      train_ind <- sample(seq_len(as.integer(ncol(Dataforanalyze$expression))), replace=FALSE, size = smp_size)
      train <- Dataforanalyze$expression[,train_ind]
      train_class<-Dataforanalyze$type[train_ind]
      if(train_ratio<1){
        test <- Dataforanalyze$expression[,-train_ind]
        test_class<-Dataforanalyze$type[-train_ind]
      }
      ttratio<-summary(train_class)[[1]]/summary(train_class)[[2]]
    }
    
    if(train_ratio<1){
      trainandtestlist<-list(train_ind, train, train_class, test, test_class, Dataforanalyze)
      names(trainandtestlist)<-c("Indices of training samples", "Training dataset", "Trainin dataset class info", "Test dataset", "Test dataset class info", "Data")
    }
    else {
      trainandtestlist<-list(train_ind, train, train_class, character(0), character(0), Dataforanalyze)
      names(trainandtestlist)<-c("Indices of training samples", "Training dataset", "Trainin dataset class info", "-", "-", "Data")
      
    }
    
    return(trainandtestlist)
    
  }
  
  
}