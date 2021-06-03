SNFS<- function(partition, train_ratio, selectiontype1, selectiontype2, selection_ratio, detect_com, metric_for, biomarker_ratio, skip) {
  
  #This function works for data sets contain only two classes. Hence, in this particular code chunk evaluates the number of classes
  if (length(summary(partition[[4]][,partition[[1]]]$type))>=3){
    stop("The data set contains more than two classes. Therefore, the execution of the function stopped.")
  }
  
  #data prep for classification with all features
  data_class=data.frame(partition[[4]]$type); colnames(data_class)=c("type")
  combined_data_for_filtering=cbind(data_class, (t(exprs(partition[[4]]))))
  
  #splitting train and test
  train_set<- combined_data_for_filtering[ partition[[1]],]
  test_set <- combined_data_for_filtering[-partition[[1]],]
  
  attach(train_set, warn.conflicts = FALSE)
  #classification with SVM
  x <- train_set[,-1]
  y <- train_set$type
  
  
  #proposed method step 1: ranking
  
  ###############################
  ###       Ranking Data      ###
  ###############################
  require("genefilter")
  require("FSelector")
  require("e1071")
  
  
  if(((selectiontype1=="ig")||(selectiontype2 =="ig"))){
    #ranking with information gain filter
    weights_filter_ig<-information.gain(type~., train_set)
    ranked_list_filter_ig=(weights_filter_ig[order(-weights_filter_ig$attr_importance), , drop = FALSE])
    rank_filter_ig=matrix(0,dim(x)[2],1)
    for (i in 1:dim(x)[2]){
      rank_filter_ig[i]=i
    }
    list_filter_ig=(cbind(ranked_list_filter_ig,rank_filter_ig))
    colnames(list_filter_ig)=c("importance", "rank from ig filter")
  }
  
  
  
  if(((selectiontype1=="cs")||(selectiontype2 =="cs"))){
    #ranking with chi-square filter
    weights_filter_cs<-chi.squared(type~., train_set)
    
    ranked_list_filter_cs=(weights_filter_cs[order(-weights_filter_cs$attr_importance), , drop = FALSE])
    rank_filter_cs=matrix(0,dim(x)[2],1)
    for (i in 1:dim(x)[2]){
      rank_filter_cs[i]=i
    }
    list_filter_cs=(cbind(ranked_list_filter_cs,rank_filter_cs))
    colnames(list_filter_cs)=c("importance", "rank from chi-square filter")
    
  }
  
  
  
  if(((selectiontype1=="svmrfe")||(selectiontype2 =="svmrfe"))){
    #ranking with SVM-RFE (with RFE or SQRT-RFE)
    ranked_list_svmrfe<-svmrfeFeatureRanking(x,y)
    #ranked_list_svmrfe<-svm.sqrt.rfe(x,y)
    
    rank_svmrfe=matrix(0,dim(x)[2],1)
    
    for (i in 1:dim(x)[2]){
      rank_svmrfe[i]=i
    }
    
    list_svmrfe=(cbind(ranked_list_svmrfe,rank_svmrfe))
    colnames(list_svmrfe)=c("gene number", "rank from svm-rfe")
    list_svmrfe=list_svmrfe[order(ranked_list_svmrfe),,drop=FALSE]
    rownames(list_svmrfe)=rownames(t(x))
  }
  
  
  
  if(((selectiontype1=="rf")||(selectiontype2 =="rf"))){
    #ranking with random-forest
    weights_rf<-random.forest.importance(type~., train_set , importance.type = 1)
    
    ranked_list_rf=(weights_rf[order(-weights_rf$attr_importance), , drop = FALSE])
    rank_rf=matrix(0,dim(x)[2],1)
    for (i in 1:dim(x)[2]){
      rank_rf[i]=i
    }
    list_rf=(cbind(ranked_list_rf,rank_rf))
    colnames(list_rf)=c("importance", "rank from random-forest")
    
  }
  
  
  
  if(((selectiontype1=="rfrfe")||(selectiontype2 =="rfrfe"))){
    #ranking with RF-RFE
    ranked_list_rfrfe<-rf.sqrt.rfe(x,y)
    
    rank_rfrfe=matrix(0,dim(x)[2],1)
    
    for (i in 1:dim(x)[2]){
      rank_rfrfe[i]=i
    }
    
    list_rfrfe=(cbind(ranked_list_rfrfe,rank_rfrfe))
    colnames(list_rfrfe)=c("gene number", "rank from rf-rfe")
    list_rfrfe=list_rfrfe[order(ranked_list_rfrfe),,drop=FALSE]
    rownames(list_rfrfe)=rownames(t(x))
  }
  
  
  
  
  #proposed method step 2: merging the rankings and selecting top % genes
  
  if(((selectiontype1=="ig")|| (selectiontype2 =="ig"))&& ((selectiontype1=="svmrfe")|| (selectiontype2 =="svmrfe"))){
    ranked_data=(cbind((list_svmrfe[rownames(list_filter_ig),]), list_filter_ig))
    final_rank_sum=ranked_data$`rank from svm-rfe`+ranked_data$`rank from ig filter`
    
  }
  
  if(((selectiontype1=="ig")|| (selectiontype2 =="ig"))&& ((selectiontype1=="cs")|| (selectiontype2 =="cs"))){
    ranked_data=(cbind((list_filter_cs[rownames(list_filter_ig),]), list_filter_ig))
    final_rank_sum=ranked_data$`rank from chi-square filter`+ranked_data$`rank from ig filter`
    
  }
  
  if(((selectiontype1=="ig")|| (selectiontype2 =="ig"))&& ((selectiontype1=="rf")|| (selectiontype2 =="rf"))){
    ranked_data=(cbind((list_rf[rownames(list_filter_ig),]), list_filter_ig))
    final_rank_sum=ranked_data$`rank from random-forest`+ranked_data$`rank from ig filter`
    
  }
  
  if(((selectiontype1=="cs")|| (selectiontype2 =="cs"))&& ((selectiontype1=="svmrfe")|| (selectiontype2 =="svmrfe"))){
    ranked_data=(cbind((list_svmrfe[rownames(list_filter_cs),]), list_filter_cs))
    final_rank_sum=ranked_data$`rank from svm-rfe`+ranked_data$`rank from chi-square filter`
    
  }
  
  if(((selectiontype1=="cs")|| (selectiontype2 =="cs"))&& ((selectiontype1=="rf")|| (selectiontype2 =="rf"))){
    ranked_data=(cbind((list_rf[rownames(list_filter_cs),]), list_filter_cs))
    final_rank_sum=ranked_data$`rank from random-forest`+ranked_data$`rank from chi-square filter`
    
  }
  
  if(((selectiontype1=="rf")|| (selectiontype2 =="rf"))&& ((selectiontype1=="svmrfe")|| (selectiontype2 =="svmrfe"))){
    ranked_data=(cbind((list_svmrfe[rownames(list_rf),]), list_rf))
    final_rank_sum=ranked_data$`rank from svm-rfe`+ranked_data$`rank from random-forest`
    
  }
  
  if(((selectiontype1=="rf")|| (selectiontype2 =="rf"))&& ((selectiontype1=="rfrfe")|| (selectiontype2 =="rfrfe"))){
    ranked_data=(cbind((list_rfrfe[rownames(list_rf),]), list_rf))
    final_rank_sum=ranked_data$`rank from rf-rfe`+ranked_data$`rank from random-forest`
    
  }
  
  if(((selectiontype1=="cs")|| (selectiontype2 =="cs"))&& ((selectiontype1=="rfrfe")|| (selectiontype2 =="rfrfe"))){
    ranked_data=(cbind((list_rfrfe[rownames(list_filter_cs),]), list_filter_cs))
    final_rank_sum=ranked_data$`rank from rf-rfe`+ranked_data$`rank from chi-square filter`
    
  }
  
  if(((selectiontype1=="ig")|| (selectiontype2 =="ig"))&& ((selectiontype1=="rfrfe")|| (selectiontype2 =="rfrfe"))){
    ranked_data=(cbind((list_rfrfe[rownames(list_filter_ig),]), list_filter_ig))
    final_rank_sum=ranked_data$`rank from rf-rfe`+ranked_data$`rank from ig filter`
    
  }
  
  
  
  final_ranked_data=(cbind(ranked_data, final_rank_sum))
  head(final_ranked_data)
  
  #ordering genes according to ranks
  ordered_genes=(final_ranked_data[order(final_ranked_data$final_rank_sum), , drop = FALSE])
  
  #expression data according to the ranks
  exprsmerged=exprs(partition[[4]][,partition[[1]]])[rownames(ordered_genes),]
  
  
  
  #selecting %selection_ratio percent of the total feature
  selection=(ceiling(dim(partition[[4]])[[1]]*selection_ratio))
  top_percent=exprsmerged[1:selection,]
  top_percent=top_percent[,order(partition[[4]][,partition[[1]]]$type)]
  top_percent_class=partition[[4]][,partition[[1]]]$type[order(partition[[4]][,partition[[1]]]$type)]
  
  
  if(skip==FALSE){
    ###############################
    ###      Aggregate Data     ###
    ###############################
    
    #seperate groups (binary class problem) 
    top_percent1=top_percent[,1:summary(partition[[4]][,partition[[1]]]$type)[1]]
    top_percent2=top_percent[,summary(partition[[4]][,partition[[1]]]$type)[1]+1:summary(partition[[4]][,partition[[1]]]$type)[2]]
    
    
    aggregateddatatop1=matrix(0,selection,1)
    for (i in 1:selection){
      aggregateddatatop1[i,1]= mean(top_percent1[i,])
    }
    aggregateddatatop2=matrix(0,selection,1)
    for (i in 1:selection){
      aggregateddatatop2[i,1]= mean(top_percent2[i,])
    }
    
    
    top_percenta=cbind(aggregateddatatop1,aggregateddatatop2)
    rownames(top_percenta)=rownames(exprsmerged[1:selection,])
    colnames(top_percenta)=c("ALL","AML")
    
    
    
    
    
    ###############################
    ###     Bi-Clustering       ###
    ###############################
    
    adjacency_matrix=matrix(0,selection,selection)
    rownames(adjacency_matrix)=rownames(top_percenta)
    colnames(adjacency_matrix)=rownames(top_percenta)
    
    #If it is desired an iterative process can be added to clustering for ex. this procedure can be repeated 5 times
    
    cluster1=kmeans(top_percenta, 3)
    
    cluster_member1=matrix(cluster1$cluster,dim(adjacency_matrix)[1],1)
    rownames(cluster_member1)=names(cluster1$cluster)
    colnames(cluster_member1)=c("Cluster Membership")
    cluster_member1
    for (i in 1:dim(cluster_member1)[1]) {
      which_cluster=cluster_member1[i,1]
      for (j in 1:dim(cluster_member1)[1]) {
        if (cluster_member1[j,1]==which_cluster) {
          adjacency_matrix[i,j]=adjacency_matrix[i,j]+1
        }
      }
    }
    
    
    cluster2=kmeans(top_percenta, 4)
    
    cluster_member2=matrix(cluster2$cluster,dim(adjacency_matrix)[1],1)
    rownames(cluster_member2)=names(cluster2$cluster)
    colnames(cluster_member2)=c("Cluster Membership")
    cluster_member2
    for (i in 1:dim(cluster_member2)[1]) {
      which_cluster=cluster_member2[i,1]
      for (j in 1:dim(cluster_member2)[1]) {
        if (cluster_member2[j,1]==which_cluster) {
          adjacency_matrix[i,j]=adjacency_matrix[i,j]+1
        }
      }
    }
    
    
    cluster3=kmeans(top_percenta, 5)
    
    cluster_member3=matrix(cluster3$cluster,dim(adjacency_matrix)[1],1)
    rownames(cluster_member3)=names(cluster3$cluster)
    colnames(cluster_member3)=c("Cluster Membership")
    cluster_member3
    for (i in 1:dim(cluster_member3)[1]) {
      which_cluster=cluster_member3[i,1]
      for (j in 1:dim(cluster_member3)[1]) {
        if (cluster_member3[j,1]==which_cluster) {
          adjacency_matrix[i,j]=adjacency_matrix[i,j]+1
        }
      }
    }
    
    
    diag(adjacency_matrix)=0
    
    
    
    unweighted_adjacency_matrix=adjacency_matrix
    for (i in 1:dim(adjacency_matrix)[1]){
      for(j in 1:dim(adjacency_matrix)[2]) {
        
        if(adjacency_matrix[i,j]==0){
          unweighted_adjacency_matrix[i,j]=0
        }
        else {
          unweighted_adjacency_matrix[i,j]=1
        } 
        
      }
    }
    
    
    require(igraph)
    SNAmat=as.matrix(adjacency_matrix)                                  # coerces the data set as a matrix
    SNAgraph=graph.adjacency(SNAmat,mode="undirected", weighted="TRUE") # this will create an 'igraph object'
    #plot.igraph(SNAgraph, vertex.label=V(SNAgraph)$name)
    
  }
  
  
  if(skip==TRUE){
    
    #Calculate Similarity Measurements
    
    #1.Correlation as Similarity Measurement
    cormatrix<- cor(t(top_percent), method="pearson")     
    
    #absolute value for negative correlations *Referance: Steve Horvath, Weighted Network Analysis: Applications in Genomics and System Biology Chapter 5 pg.91-121, Springer
    abscormatrix<-abs(cormatrix)                                       #also mentioned in chapter 7 as "Unsigned correlation network adjacency matrix"
    
    
    #or transform to the range [0,1] 
    #for (i in 1:dim(cormatrix)[1]){
    #    for(j in 1:dim(cormatrix)[2]){
    #      
    #        cormatrix[i,j]=(0.5+0.5*cormatrix[i,j])                   #also mentioned in chapter 7 as "Signed correlation network adjacency matrix"
    #    }
    #  }
    
    
    #zero for negative correlations (biased and unsuitable for genetic expression data *Gates et al.)
    #for (i in 1:dim(cormatrix)[1]){
    #  for(j in 1:dim(cormatrix)[2]){
    #    if(cormatrix[i,j]<0){
    #      cormatrix[i,j]=0
    #    }
    #  }
    #}
    
    
    #Thresholding
    
    #Calculate Hard Thresholding Cut-off: Retaining only top %1  of the relationships method (Easy)
    thresholdmatrix<- abscormatrix; diag(thresholdmatrix)<-NA; thresholdmatrix<-thresholdmatrix[!is.na(thresholdmatrix)]
    cutoff_t=as.numeric(quantile(thresholdmatrix, probs=c(0.99), na.rm = TRUE))
    
    #Hard thresholding: Bonferroni p-value method?
    #Hard thresholding: Spectral?
    
    #Hard thresholding: Arbitrary?
    #cutoff_t=0.30
    
    
    
    
    #NOTE: Filtering invalidates the Scale-free topology assumption. Therefore we can not use the method in Zhang, Bin and Horvath, (2005) "A general framework for weighted gene co-expression network analysis", Statistical Applications in Genetics and Molecular Biology: 4:1,1-43.
    
    #Soft thresholding: Default depending on the sample size
    #if(length(top_percent_class)<20){
    #  power_t=18
    #}
    #if((length(top_percent_class)<=30)&&(length(top_percent_class)>=20)){
    #  power_t=16
    #}
    #if((length(top_percent_class)<=40)&&(length(top_percent_class)>30)){
    #  power_t=14
    #}
    #if(length(top_percent_class)>40){
    #  power_t=12
    #}
    
    #Soft thresholding: Arbitrary
    power_t=1
    
    
    #Coverting object to weighted network (with soft thresholding)
    adjacency_matrix<-(abscormatrix)^(power_t)
    
    
    #Coverting object to network (Hard thresholding)
    unweighted_adjacency_matrix=matrix(0,dim(adjacency_matrix)[1],dim(adjacency_matrix)[2])
    rownames(unweighted_adjacency_matrix)=rownames(adjacency_matrix)
    colnames(unweighted_adjacency_matrix)=rownames(adjacency_matrix)
    
    for (i in 1:dim(adjacency_matrix)[1]){
      for(j in 1:dim(adjacency_matrix)[2]){
        if(adjacency_matrix[i,j]>cutoff_t){
          unweighted_adjacency_matrix[i,j]=1
        }
        else {
          unweighted_adjacency_matrix[i,j]=0
        }
      }
    }
    
    require("igraph")
    SNAmat=as.matrix(adjacency_matrix)    
    SNAgraph<-graph.adjacency(SNAmat, weighted=TRUE, mode="undirected", diag=FALSE)
    
    
  }
  
  
  #NOTE: following code requires Force Atlas 2 function to draw sociogram with force atlas layout.
  
  
  #proposed method step 4: community detection, biomarker selection
  
  if(detect_com=="walktrap"){
    imc_walktrap <-cluster_walktrap(SNAgraph)
    community_list=communities(imc_walktrap)
    community_membership=membership(imc_walktrap)
    community_index=community_list
    #modularity(imc_walktrap)
  }
  
  if(detect_com=="infomap"){
    imc_infomap  <-cluster_infomap(SNAgraph)
    community_list=communities(imc_infomap)
    community_membership=membership(imc_infomap)
    community_index=community_list
    #modularity(imc_infomap)
  }
  
  if(detect_com=="louvain"){
    imc_louvain   <-cluster_louvain(SNAgraph)
    community_list=communities(imc_louvain)
    community_membership=membership(imc_louvain)
    community_index=community_list
    #modularity(imc_louvain)
  }
  
  if(detect_com=="fast"){
    imc_fastgreedy<- cluster_fast_greedy(SNAgraph)
    community_list=communities(imc_fastgreedy)
    community_membership=membership(imc_fastgreedy)
    community_index=community_list
    #modularity(imc_fastgreedy)
  }
  
  
  ##############################
  ### Plotting Communities
  ##############################
  
  community_number=length(community_list); community_number
  if(community_number>8){
    stop("There are more than eight communities. Therefore the execution of the function stops.")
  }
  
  community_list
  #require(sna)
  layout_identical<-layout.forceatlas2(SNAgraph, iterations = 500, plotstep=10)
  #sna_graph<-plot.igraph(SNAgraph, layout=layout_nicely, mark.groups = community_list, mark.shape=1, mark.col = rainbow(length(community_list), alpha = 0.3), mark.border=rainbow(length(community_list),alpha=1), mark.expand= 25, vertex.label=V(SNAgraph)$name, vertex.size=5, vertex.color="red", vertex.frame.color="red", vertex.shape="sphere", vertex.label.cex=0.7, vertex.label.color="black", vertex.label.font=2, edge.color="grey", edge.width=0.1, margin=0.2, frame="TRUE")
  #sna_graph<-plot.igraph(SNAgraph, layout=layout_identical, mark.groups = community_list, mark.shape=1, mark.col = rainbow(length(community_list), alpha = 0.3), mark.border=rainbow(length(community_list),alpha=1), mark.expand= 25, vertex.label=V(SNAgraph)$name, vertex.size=5, vertex.color="red", vertex.frame.color="red", vertex.shape="sphere", vertex.label.cex=0.7, vertex.label.color="black", vertex.label.font=2, edge.color="grey", edge.width=0.1, margin=0.2, frame="TRUE")
  
  
  #coords <- layout_(SNAgraph, as_star())
  #pdf(file = "communities.pdf", title="Communities")
  #plot.igraph(SNAgraph, layout=layout_identical, mark.groups = community_list, mark.shape=1, mark.col = rainbow(length(community_list), alpha = 0.3), mark.border=rainbow(length(community_list),alpha=1), mark.expand= 25, vertex.label=V(SNAgraph)$name, vertex.size=5, vertex.color="red", vertex.frame.color="red", vertex.shape="sphere", vertex.label.cex=0.7, vertex.label.color="black", vertex.label.font=2, edge.color="grey", edge.width=0.1, margin=0.2, frame="TRUE")
  #dev.off()
  
  
  
  
  ##################################
  #calculating network measurements#
  ##################################
  
  #unweighted degree of centrality
  unweighted_degree=matrix(0,dim(unweighted_adjacency_matrix)[1],1)
  rownames(unweighted_degree)=rownames(unweighted_adjacency_matrix)
  colnames(unweighted_degree)=c("unweighted degree centrality")
  for (i in 1:dim(unweighted_adjacency_matrix)[1]){
    unweighted_degree[i,1]=sum(unweighted_adjacency_matrix[i,])
  }
  
  
  #weighted degree centrality with degree/node strength
  weighted_degree=matrix(0,dim(adjacency_matrix)[1],3)
  rownames(weighted_degree)=rownames(adjacency_matrix)
  colnames(weighted_degree)=c("node strength or degree", "membership", "weighted degree centrality")
  
  #library(tnet)
  #degree_w(SNA_matrix)                                       #the function gives exactly the same results with our calculation of node strength
  
  for (i in 1:dim(adjacency_matrix)[1]){
    weighted_degree[i,1]=sum(adjacency_matrix[i,])-1          #Sum of the total weights minus 1. Because weights of the diagonal elements are omitted.
    weighted_degree[i,2]=as.integer(community_membership)[i]
    
    if(weighted_degree[i,2]==1){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"1")
    }
    
    if(weighted_degree[i,2]==2){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"2")
    }
    
    if(weighted_degree[i,2]==3){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"3")
    }
    
    if(weighted_degree[i,2]==4){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"4")
    }
    
    if(weighted_degree[i,2]==5){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"5")
    }
    
    if(weighted_degree[i,2]==6){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"6")
    }
    
    if(weighted_degree[i,2]==7){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"7")
    }
    
    if(weighted_degree[i,2]==8){
      weighted_degree[i,3]=weighted_degree[i,1]/length(community_index$"8")
    }
    
  }
  
  
  
  #intra-degree of centrality
  intra=matrix(0,dim(weighted_degree)[1],dim(weighted_degree)[1])
  rownames(intra)=rownames(adjacency_matrix)
  colnames(intra)=colnames(adjacency_matrix)
  for (i in 1:dim(intra)[1]) {
    
    for (j in 1:dim(intra)[2]) {
      
      if ((weighted_degree[j,2]==weighted_degree[i,2])&&(unweighted_adjacency_matrix[i,j]>0)){
        intra[i,j]=1
      }
      
      else {
        intra[i,j]=0
      }
      
    }
    
  }
  intra_degree=matrix(rowSums(intra),dim(weighted_degree)[1],1)
  rownames(intra_degree)=rownames(adjacency_matrix)
  
  #out-of-degree centrality
  out_of=matrix(0,dim(adjacency_matrix)[1],dim(adjacency_matrix)[2])
  rownames(out_of)=rownames(adjacency_matrix)
  colnames(out_of)=colnames(adjacency_matrix)
  for (i in 1:dim(adjacency_matrix)[1]) {
    for (j in 1:dim(adjacency_matrix)[2]) {
      if ((weighted_degree[j,2]!=weighted_degree[i,2])&&(adjacency_matrix[i,j]>0)){
        out_of[i,j]=1
      }
      
      else {
        out_of[i,j]=0
      }
      
    }
    
  }
  outof_degree=matrix(rowSums(out_of),dim(out_of)[1],1)
  rownames(outof_degree)=rownames(adjacency_matrix)
  
  #degree of domesticity
  epsilon=0.0001
  degree_domesticity=(intra_degree+epsilon)/(outof_degree+epsilon)
  
  
  #coverage
  coverage=matrix(0,dim(adjacency_matrix)[1],1)
  rownames(coverage)=rownames(adjacency_matrix)
  colnames(coverage)=c("coverage")
  for (i in 1:dim(adjacency_matrix)[1]){
    if(as.integer(community_membership)[i]==1){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"1")
    }
    if(as.integer(community_membership)[i]==2){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"2")
    }
    if(as.integer(community_membership)[i]==3){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"3")
    }
    
    if(as.integer(community_membership)[i]==4){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"4")
    }
    
    if(as.integer(community_membership)[i]==5){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"5")
    }
    if(as.integer(community_membership)[i]==6){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"6")
    }
    if(as.integer(community_membership)[i]==7){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"7")
    }
    if(as.integer(community_membership)[i]==8){
      coverage[i,1]=intra_degree[i,1]/length(community_index$"8")
    }
    
    
  }
  
  
  #identity
  identity=matrix(ordered_genes[1:selection,5])
  rownames(identity)=rownames(ordered_genes)[1:selection]
  
  
  
  #determining the biomarkers based on the network metrics  
  metrics<-cbind(weighted_degree[rownames(coverage),2], unweighted_degree[rownames(coverage),], weighted_degree[rownames(coverage),3], intra_degree[rownames(coverage),], outof_degree[rownames(coverage),], degree_domesticity[rownames(coverage),],identity[rownames(coverage),],coverage)
  colnames(metrics)=c("membership","unweighted degree centrality", "weighted degree centrality", "intra-community UDC", "out-of-community UDC", "degree of domesticity", "identity", "coverage")
  
  
  if(metric_for==1){s_metric="unweighted degree centrality"};if(metric_for==2){s_metric="weighted degree centrality"};if(metric_for==3){s_metric="intra-community UDC"};if(metric_for==4){s_metric="out-of-community UDC"};if(metric_for==5){s_metric="degree of domesticity"}
  
  
  
  if (s_metric=="unweighted degree centrality"){
    index_of_metric=2
  }
  
  if (s_metric=="weighted degree centrality"){
    index_of_metric=3
  }
  
  if (s_metric=="intra-community UDC"){
    index_of_metric=4
  }
  
  if (s_metric=="out-of-community UDC"){
    index_of_metric=5
  }
  
  if (s_metric=="degree of domesticity"){
    index_of_metric=6
  }
  
  
  
  selectedbio=metrics[order(metrics[,1],-metrics[,index_of_metric]),]
  
  if(community_number==2){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
  }
  if(community_number==3){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
    com3=split(selectedbio[,index_of_metric], selectedbio[,1])$"3"
  }
  if(community_number==4){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
    com3=split(selectedbio[,index_of_metric], selectedbio[,1])$"3"
    com4=split(selectedbio[,index_of_metric], selectedbio[,1])$"4"
  }
  if(community_number==5){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
    com3=split(selectedbio[,index_of_metric], selectedbio[,1])$"3"
    com4=split(selectedbio[,index_of_metric], selectedbio[,1])$"4"
    com5=split(selectedbio[,index_of_metric], selectedbio[,1])$"5"
  }
  if(community_number==6){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
    com3=split(selectedbio[,index_of_metric], selectedbio[,1])$"3"
    com4=split(selectedbio[,index_of_metric], selectedbio[,1])$"4"
    com5=split(selectedbio[,index_of_metric], selectedbio[,1])$"5"
    com6=split(selectedbio[,index_of_metric], selectedbio[,1])$"6"
  }
  if(community_number==7){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
    com3=split(selectedbio[,index_of_metric], selectedbio[,1])$"3"
    com4=split(selectedbio[,index_of_metric], selectedbio[,1])$"4"
    com5=split(selectedbio[,index_of_metric], selectedbio[,1])$"5"
    com6=split(selectedbio[,index_of_metric], selectedbio[,1])$"6"
    com7=split(selectedbio[,index_of_metric], selectedbio[,1])$"7"
  }
  if(community_number==8){
    com1=split(selectedbio[,index_of_metric], selectedbio[,1])$"1"
    com2=split(selectedbio[,index_of_metric], selectedbio[,1])$"2"
    com3=split(selectedbio[,index_of_metric], selectedbio[,1])$"3"
    com4=split(selectedbio[,index_of_metric], selectedbio[,1])$"4"
    com5=split(selectedbio[,index_of_metric], selectedbio[,1])$"5"
    com6=split(selectedbio[,index_of_metric], selectedbio[,1])$"6"
    com7=split(selectedbio[,index_of_metric], selectedbio[,1])$"7"
    com8=split(selectedbio[,index_of_metric], selectedbio[,1])$"8"
  }
  
  
  if(community_number==2){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
  }
  if(community_number==3){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
    select3=attributes(com3[1:ceiling(length(com3)*biomarker_ratio)])$names
  }
  if(community_number==4){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
    select3=attributes(com3[1:ceiling(length(com3)*biomarker_ratio)])$names
    select4=attributes(com4[1:ceiling(length(com4)*biomarker_ratio)])$names
  }
  if(community_number==5){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
    select3=attributes(com3[1:ceiling(length(com3)*biomarker_ratio)])$names
    select4=attributes(com4[1:ceiling(length(com4)*biomarker_ratio)])$names
    select5=attributes(com5[1:ceiling(length(com5)*biomarker_ratio)])$names
  }
  if(community_number==6){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
    select3=attributes(com3[1:ceiling(length(com3)*biomarker_ratio)])$names
    select4=attributes(com4[1:ceiling(length(com4)*biomarker_ratio)])$names
    select5=attributes(com5[1:ceiling(length(com5)*biomarker_ratio)])$names
    select6=attributes(com6[1:ceiling(length(com5)*biomarker_ratio)])$names
  }
  if(community_number==7){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
    select3=attributes(com3[1:ceiling(length(com3)*biomarker_ratio)])$names
    select4=attributes(com4[1:ceiling(length(com4)*biomarker_ratio)])$names
    select5=attributes(com5[1:ceiling(length(com5)*biomarker_ratio)])$names
    select6=attributes(com6[1:ceiling(length(com5)*biomarker_ratio)])$names
    select7=attributes(com7[1:ceiling(length(com5)*biomarker_ratio)])$names
  }
  if(community_number==8){
    select1=attributes(com1[1:ceiling(length(com1)*biomarker_ratio)])$names
    select2=attributes(com2[1:ceiling(length(com2)*biomarker_ratio)])$names
    select3=attributes(com3[1:ceiling(length(com3)*biomarker_ratio)])$names
    select4=attributes(com4[1:ceiling(length(com4)*biomarker_ratio)])$names
    select5=attributes(com5[1:ceiling(length(com5)*biomarker_ratio)])$names
    select6=attributes(com6[1:ceiling(length(com5)*biomarker_ratio)])$names
    select7=attributes(com7[1:ceiling(length(com5)*biomarker_ratio)])$names
    select8=attributes(com8[1:ceiling(length(com5)*biomarker_ratio)])$names
    
  }
  
  
  
  if(community_number==2){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,]); rownames(selected_metrics)=c(select1,select2)
  }
  if(community_number==3){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,], selectedbio[select3,]); rownames(selected_metrics)=c(select1,select2,select3)
  }
  if(community_number==4){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,], selectedbio[select3,], selectedbio[select4,]); rownames(selected_metrics)=c(select1,select2,select3,select4)
  }
  if(community_number==5){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,], selectedbio[select3,], selectedbio[select4,], selectedbio[select5,]); rownames(selected_metrics)=c(select1,select2,select3,select4,select5)
  }
  if(community_number==6){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,], selectedbio[select3,], selectedbio[select4,], selectedbio[select5,], selectedbio[select6,]); rownames(selected_metrics)=c(select1,select2,select3,select4,select5,select6)
  }
  if(community_number==7){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,], selectedbio[select3,], selectedbio[select4,], selectedbio[select5,], selectedbio[select6,], selectedbio[select7,]); rownames(selected_metrics)=c(select1,select2,select3,select4,select5,select6,select7)
  }
  if(community_number==8){
    selected_metrics=rbind(selectedbio[select1,],selectedbio[select2,], selectedbio[select3,], selectedbio[select4,], selectedbio[select5,], selectedbio[select6,], selectedbio[select7,], selectedbio[select8,]); rownames(selected_metrics)=c(select1,select2,select3,select4,select5,select6,select7,select8)
  }
  
  transpozx=exprs(partition[[4]])
  
  
  
  if(community_number==2){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,])); colnames(gene_selected)=c(select1,select2)
  }
  if(community_number==3){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,], transpozx[select3,])); colnames(gene_selected)=c(select1,select2,select3)
  }
  if(community_number==4){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,], transpozx[select3,], transpozx[select4,])); colnames(gene_selected)=c(select1,select2,select3, select4)
  }
  if(community_number==5){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,], transpozx[select3,], transpozx[select4,], transpozx[select5,])); colnames(gene_selected)=c(select1,select2,select3,select4,select5)
  }
  if(community_number==6){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,], transpozx[select3,], transpozx[select4,], transpozx[select5,], transpozx[select6,])); colnames(gene_selected)=c(select1,select2,select3,select4,select5,select6)
  }
  if(community_number==7){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,], transpozx[select3,], transpozx[select4,], transpozx[select5,], transpozx[select6,], transpozx[select7,])); colnames(gene_selected)=c(select1,select2,select3,select4,select5,select6,select7)
  }
  if(community_number==8){
    gene_selected=t(rbind(transpozx[select1,], transpozx[select2,], transpozx[select3,], transpozx[select4,], transpozx[select5,], transpozx[select6,], transpozx[select7,], transpozx[select8,])); colnames(gene_selected)=c(select1,select2,select3,select4,select5,select6,select7,select8)
  }
  
  
  general_dataset=data.frame(partition[[4]]$type,gene_selected)
  colnames(general_dataset)[1]=c("y")
  
  #splitting train and test
  train<- general_dataset[ partition[[1]],]
  test<-  general_dataset[-partition[[1]],]
  
  set.seed(123)
  xtrain<- train[,-1]
  ytrain<- train$y
  
  require("caret")
  tuning<-tune.svm(xtrain,ytrain, kernel="radial", cost=10^(-10:10), gamma=2^(-10:10), validation.x=test[,-1], validation.y=test[,1], predict.func = predict)
  modeltrain<-svm(xtrain,ytrain, probability=TRUE, type="C-classification", kernel="radial", gamma=tuning$best.model$gamma, cost=tuning$best.model$cost)
  
  #modeltrain<-svm(xtrain,ytrain, probability=TRUE, type="C-classification")
  predict_train<-predict(modeltrain, xtrain, decision.values = TRUE, probability = TRUE)
  confusion_train<-confusionMatrix(predict(modeltrain, xtrain, decision.values = TRUE, probability = TRUE), train$y)
  accurracy_train=confusion_train$overall[[1]]
  sensitivity_train=confusion_train$byClass[[1]]
  specificity_train=confusion_train$byClass[[2]]
  
  
  predict_test<-predict(modeltrain, newdata=test[,-1], probability = TRUE, decision.values = TRUE)
  confusion_test<-confusionMatrix(predict(modeltrain, newdata=test[,-1], probability = TRUE, decision.values = TRUE), test$y)
  accurracy_test=confusion_test$overall[[1]]
  sensitivity_test=confusion_test$byClass[[1]]
  specificity_test=confusion_test$byClass[[2]]
  
  
  
  #require(ROCR)
  #svm.prob.rocr <- prediction(attr(predict_test, "probabilities")[,2], test$y)
  #svm.perf <- performance(svm.prob.rocr, "tpr","fpr")
  #plot(svm.perf, col=2)
  #slot(performance(svm.prob.rocr, "auc"),"y.values")[[1]]
  
  require(pROC)
  gc_pROC_train <- roc(response = train$y, predictor = attr(predict_train, "probabilities")[,2])
  plot(gc_pROC_train)
  gc_pROC_train$auc
  
  gc_pROC_test <- roc(response = test$y, predictor = attr(predict_test, "probabilities")[,2])
  plot(gc_pROC_test)
  gc_pROC_test$auc 
  
  output=list(partition, gc_pROC_test$auc, accurracy_test, sensitivity_test, specificity_test, confusion_test, gc_pROC_train$auc, accurracy_train, sensitivity_train, specificity_train, confusion_train, community_number, selected_metrics)
  names(output)=c("Partition", "AUC Test","Acccuracy Test","Sensitivity Test","Specificity Test", "Confusion Matrix of Test Data", "AUC Train","Acccuracy Train","Sensitivity Train","Specificity Train", "Confusion Matrix of Train Data", "Community number", "Metrics of Selected Biomarkers")
  
  
  return(output)
  detach(train_set)
  
}