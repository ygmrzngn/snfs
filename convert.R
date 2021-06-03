# Before SNFS the data must be convert to ExpressionSet class in order to apply the analysis. 
# Therefore, the function converts the exprSet class object to ExpressionSet object

convert.exprSet.to.ExpressionSet <- function(exprSet, sample_index){
  
  require(plyr)
  
  if (class(exprSet)=="exprSet"){
    
    Assay_Data=slot(exprSet, "exprs")
    colnames(Assay_Data)=sample_index
    Assay_Data=Assay_Data[,order(as.integer(colnames(Assay_Data)))] 
    
    Pheno_Data=slot(slot(exprSet, "phenoData"),"pData")
    Pheno_Data=arrange(Pheno_Data,sample_index)
    rownames(Pheno_Data)=colnames(Assay_Data)
    
    Meta_Data=data.frame(row.names=colnames(Pheno_Data))
    
    Annotated_Data_Frame=new("AnnotatedDataFrame", data=Pheno_Data, varMetadata=Meta_Data)
    Annotated_Data_Frame@varMetadata$labelDescription=(slot(slot(exprSet, "phenoData"),"varLabels"))
    
    Experiment_Data=slot(exprSet, "description")
    
    Data_Annotation=slot(exprSet, "annotation")
    
    DataclassExpressionSet=new("ExpressionSet", exprs=Assay_Data, phenoData=Annotated_Data_Frame, experimentData=Experiment_Data, annotation=Data_Annotation)
    
    return(DataclassExpressionSet)
    
  }
  
  else {
    print("*exprSet* sinifindan verinin *ExpressionSet* sinifindan veriye donusturulme islemi icin sinifi *exprSet* olan bir data kullaniniz")
  }
};convert<-function(Dataforanalyze, sample_index){
  
  if(class(Dataforanalyze)[1]=="ExpressionSet"){
    Dataforanalyze<-Dataforanalyze
  }
  if(class(Dataforanalyze)[1]=="exprSet"){
    Dataforanalyze <-convert.exprSet.to.ExpressionSet(Dataforanalyze, sample_index)
  }
  else 
    print("function requires exprSet or ExpressionSet class object")
}