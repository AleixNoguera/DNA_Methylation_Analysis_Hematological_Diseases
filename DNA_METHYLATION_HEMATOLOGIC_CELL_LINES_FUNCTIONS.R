# Colours Palette
Colours_Pal<-c("Leukaemia" = "#F3B61F",
               "Lymphoma" = "#15A390",
               "Multiple_Myeloma" = "#CE5755",
               "B_ALL" = "#c22f87",
               "Burkitt" = "#5E3203",
               "DLBCL" = "#F37612",
               "Hodgkin_Lymphoma" = "#F5ADC5",
               "Lymphoblastic_Leukemia_Others" = "#719CB9",
               "MBCL" = "#ad8f57",              
               "CML" = "#A9D39F",
               "AML" = "#4A904A",          
               "NH_B_lymphoma_Others" = "#E3E500",                         
               "T_ALL" = "#C7F7FF",
               "T_cell_Non_Hodgkin_Lymphoma" = "#3957c4",         
               "Transformed_cell_line" = "#C994F2",
               "Human" = "black",
               "Mouse" = "gray80",
               "NA" = "white",
               "B_cell" = "#F9A620",
               "Myeloblast" = "#8A6C09",
               "Myeloid_lineage"= "#8A6C09",
               "Myeloid_Lineage"= "#8A6C09",
               "T_cell" = "#A8D5E2",
               "Cancer_Cells" = "#548C2F",
               "Normal_Cells" = "#FF715B",
               "Primary_Sample" = "#5398BE",
               "Cell_Line" = "#F2CD5D")

# Preprocessing function
Preprocessing_function<-function(IDATs_dir,Organism,targets_df,Ann){
  
  Results_List<-list()
  
  if (Organism == "Homo_Sapiens"){
    bVals_temp <- openSesame(IDATs_dir, prep="QCDPB")
  } 
  
  else if (Organism == "Mus_Musculus"){
    bVals_temp <- openSesame(IDATs_dir, prep="QCDPB")
  }
  
  bVals<-bVals_temp[,match(targets_df$EPIC_Barcode,colnames(bVals_temp))]
  colnames(bVals)<-targets_df$Cell_Name
  # Remove probes that all the values are NA
  bVals<-bVals[rowSums(is.na(bVals)) != ncol(bVals), ]
  # Filter the probes from the annotation file
  bVals<-bVals[rownames(bVals) %in% Ann$CpG_ID,]
  # Remove sex chromosome probes   
  Sex_CpG<-Ann[Ann$chr %in% c("chrX","chrY"),]
  bVals<-bVals[!rownames(bVals) %in% (Sex_CpG$CpG_ID),]
  
  # Imputation of NA values
  Results_List[["Total_NA"]]<-sum(is.na(bVals))
  Results_List[["bVals_NI"]]<-bVals
  bVals<-t(apply(bVals,1,function(x) replace(x, is.na(x), median(x,na.rm=TRUE))))
  Results_List[["bVals"]]<-bVals
  
  return(Results_List)
  
}  

# Unsupervised HeatMap function
Unsupervised_HeatMap<-function(bdf,nProbes,targetsdf,BoottomAnn,Col_Split){
  
  #Select a random list of CpGs to perform the clustering and heatmap
  CpG_List_Random<-rownames(bdf)
  
  ## Set seed
  set.seed(684)
  ## Select 1% of CpGs at random to perform the clustering
  CpG_List_Random<-sample(CpG_List_Random, round(nrow(bdf)*nProbes,0))
  
  # Select the random CpGs into the Bvals dataframe
  bVals_Random<-as.data.frame(bdf)
  bVals_Random<-bVals_Random[rownames(bVals_Random) %in% CpG_List_Random,]
  
  distance <- dist(t(bdf), method="euclidean") # calcular las distancias
  clust.model <- hclust(distance, method="ward.D2") # clustering
  
  distance_r <- dist(bVals_Random, method="euclidean") # calcular las distancias
  clust.model_r <- hclust(distance_r, method="ward.D2") # clustering
  
  
  
  Heatmap<-draw(Heatmap(as.matrix(bVals_Random), name = "Methylation", 
                        col = gplots::greenred(75), 
                        bottom_annotation = BoottomAnn,
                        cluster_columns = clust.model,
                        cluster_rows = clust.model_r,
                        column_split= Col_Split, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE, 
                        show_row_dend = FALSE))
  
  return(Heatmap)
  
}

# t-SNE function
tSNE_Process_function<-function(TSNEObj,targets_df,Column_Col,Organism){
  
  colnames(targets_df)[colnames(targets_df) == Column_Col] <- "Column_Col"
  
  
  
  colnames(TSNEObj$Y)<-c("tSNE_1","tSNE_2")
  
  TSNEObj <- TSNEObj$Y %>% 
    as.data.frame() %>%
    mutate(ID=row_number())
  
  targets_df$ID<-1:nrow(targets_df)
  
  TSNEObj <- TSNEObj %>% inner_join(targets_df , by="ID")
  
  if(Organism == "Mouse_Human"){
  
    Plot<-ggplot(TSNEObj,aes(x = tSNE_2, 
                              y = tSNE_1))+
      geom_point(aes(fill = Column_Col,
                     shape = Organism),
                 size= 4)+
      scale_fill_manual(values = Colours_Pal)+
      scale_shape_manual(values = c("Human" = 21,
                                    "Mouse" = 22))+
      theme(legend.position="bottom")+
      theme_few()
  }
  
  else if(Organism == "Human"){
  
    Plot<-ggplot(TSNEObj,aes(x = tSNE_2, 
                             y = tSNE_1,
                             fill = Column_Col))+
    geom_point(size= 4,
               shape = 21)+
    scale_fill_manual(values = Colours_Pal)+
    theme(legend.position="bottom")+
    theme_few()
  }
   
  else if(Organism == "Mouse"){
    
    Plot<-ggplot(TSNEObj,aes(x = tSNE_2, 
                             y = tSNE_1,
                             fill = Column_Col))+
      geom_point(size= 4,
                 shape = 22)+
      scale_fill_manual(values = Colours_Pal)+
      theme(legend.position="bottom")+
      theme_few()
  }
   
  return(Plot)

}

# Phylogenetic Tree
Phylo_Tree_function<-function(bVals, Tree_Ann, Colours_Pal){
  
  tree <- as.phylo(hclust(dist(t(bVals), method="euclidean"), method="ward.D2")) # clustering))
  rownames(Tree_Ann)<-tree$tip.label
  Tree_ALL<-gheatmap(ggtree(tree, layout="fan"), Tree_Ann, offset=.8, width=.2,colnames = FALSE) +
    scale_fill_manual(values = Colours_Pal)
}


Dendogram_Function<-function(targets_df,bVals){
  
  targets_df$Class_Org<-ifelse(targets_df$Final_Classification == "",
                               paste0(targets_df$Disease,"_",targets_df$Organism),
                               paste0(targets_df$Final_Classification,"_",targets_df$Organism))
  
  
  res_df<-data.frame(CpG_ID = rownames(bVals))
  
  for (i in 1:length(unique(targets_df$Class_Org))) {
    Group<-unique(targets_df$Class_Org)[i]
    Samples<-targets_df[targets_df$Class_Org == Group,]
    bdf<-data.frame(CpG_ID = rownames(bVals), Group = rowMeans(bVals[,colnames(bVals) %in% Samples$Cell_Name]))
    colnames(bdf)[2]<-Group
    res_df<-merge(res_df, bdf, "CpG_ID")
  }
  
  res_df<-res_df[,!colnames(res_df) %like% "Others"]
  
  rownames(res_df)<-res_df$CpG_ID
  res_df<-res_df[-1]
  
  
  Metdata<-data.frame(Group_ID = colnames(res_df))
  Metdata$General_Disease<-ifelse(Metdata$Group_ID %like% "AML", "Leukaemia",
                           ifelse(Metdata$Group_ID %like% "B_ALL", "Leukaemia",      
                           ifelse(Metdata$Group_ID %like% "CML", "Leukaemia",
                           ifelse(Metdata$Group_ID %like% "Leukaemia", "Leukaemia",              
                           ifelse(Metdata$Group_ID %like% "Lymphoblastic_Leukemia_Others", "Leukaemia",
                           ifelse(Metdata$Group_ID %like% "T_ALL", "Leukaemia",
                           ifelse(Metdata$Group_ID %like% "Burkitt", "Lymphoma",
                           ifelse(Metdata$Group_ID %like% "Hodgkin_Lymphoma", "Lymphoma",        
                           ifelse(Metdata$Group_ID %like% "MBCL", "Lymphoma",       
                           ifelse(Metdata$Group_ID %like% "NH_B_lymphoma_Others", "Lymphoma",         
                           ifelse(Metdata$Group_ID %like% "T_cell_Non_Hodgkin_Lymphoma", "Lymphoma",
                           ifelse(Metdata$Group_ID %like% "DLBCL", "Lymphoma",
                           ifelse(Metdata$Group_ID %like% "Lymphoma", "Lymphoma",        
                           ifelse(Metdata$Group_ID %like% "Multiple_Myeloma", "Multiple_Myeloma","Transformed_cell_line"))))))))))))))
  
  Metdata$Group<-gsub("_Human","",Metdata$Group_ID)
  Metdata$Group<-gsub("_Mouse","",Metdata$Group)
  
  Metdata$Organism<-ifelse(Metdata$Group_ID %like% "_Human", "Human", "Mouse")
  
  Metdata$Group<-ifelse(Metdata$Organism == "Mouse","NA",Metdata$Group)
  
  Vardf<-data.frame(CpG_ID = rownames(res_df) ,Variance = rowVars(as.matrix(res_df)))
  Vardf<-Vardf[order(-Vardf$Variance),]
  Vardf<-Vardf[c(1:(round((nrow(Vardf)*0.01),0))),]
  res_df_var<-res_df[rownames(res_df) %in% Vardf$CpG_ID,]
  
  Tree<-Phylo_Tree_function(res_df_var,as.data.frame(Metdata[c(2:3)]),Colours_Pal)
  
  return(Tree)
}


# Supervised HeatMaps
Draw_HeatMap<-function(bVals_Cluster,bVals_HeatMap, bottom_annotation, ColSplit){
  Result_list<-list()
  distance <- dist(t(bVals_Cluster), method="euclidean") # calcular las distancias
  clust.model <- hclust(distance, method="ward.D2") # clustering
  
  distance_r <- dist(bVals_Cluster, method="euclidean") # calcular las distancias
  clust.model_r <- hclust(distance_r, method="ward.D2") # clustering
  
  Heatmap<-draw(Heatmap(as.matrix(bVals_HeatMap), name = "Methylation", 
                        col = gplots::greenred(75), 
                        bottom_annotation = bottom_annotation,
                        cluster_columns = clust.model,
                        cluster_rows = clust.model_r,
                        column_split= ColSplit, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE, 
                        show_row_dend = FALSE))
  
  
  Result_list[["HeatMap"]] <- Heatmap
  
  
  Cluster_Rows<-Cluster_Order_Rows(clust.model_r,Heatmap)
  Result_list[["Cluster_Rows"]] <- Cluster_Rows
  
  Cluster_Columns<-Cluster_Order_Columns(clust.model,Heatmap)
  Result_list[["Cluster_Columns"]] <- Cluster_Columns
  
  
  return(Result_list)
}

# Function to determine Rows order in a HeatMap
Cluster_Order_Rows<-function(clust.model,HeatMap_split){
  Labels<-as.data.frame(clust.model$labels)
  Labels$Row_order_HeatMap<-1:nrow(Labels)
  HM_Column_list<-row_order(HeatMap_split)
  
  if (length(HM_Column_list) != nrow(Labels)){
    Clust_Num<-length(HM_Column_list)
    Cluster_Data<-data.frame()
    for (i in 1:Clust_Num){
      temp<-as.data.frame(HM_Column_list[[i]])
      colnames(temp)<-"Row_order_HeatMap"
      temp$Cluster<-paste0("Cluster_",i)
      Cluster_Data<-rbind(Cluster_Data,temp)
    }
    
    Labels<-Labels[match(Cluster_Data$Row_order_HeatMap,Labels$Row_order_HeatMap),]
    Cluster_Data<-merge(Labels,Cluster_Data,"Row_order_HeatMap")
    Cluster_Data<-Cluster_Data[match(Labels$Row_order_HeatMap,Cluster_Data$Row_order_HeatMap),]
    Cluster_Data<-Cluster_Data[,c(2,1,3)]
    colnames(Cluster_Data)[1]<-"CpG_Name"
    
    return(Cluster_Data)
    
  } else {
    Clust_Num<-length(HM_Column_list)
    Cluster_Data<-data.frame()
    for (i in 1:Clust_Num){
      temp<-as.data.frame(HM_Column_list[[i]])
      colnames(temp)<-"Row_order_HeatMap"
      Cluster_Data<-rbind(Cluster_Data,temp)
    }
    Labels<-Labels[match(Cluster_Data$Row_order_HeatMap,Labels$Row_order_HeatMap),]
    Cluster_Data<-Labels
    colnames(Cluster_Data)[1]<-"CpG_Name"
    return(Cluster_Data)
  } 
}

# Function to determine Column order in a HeatMap
Cluster_Order_Columns<-function(clust.model,HeatMap_split){
  Labels<-as.data.frame(clust.model$labels)
  Labels$Column_order_HeatMap<-1:nrow(Labels)
  HM_Column_list<-column_order(HeatMap_split)
  
  if (length(HM_Column_list) != nrow(Labels)){
    Clust_Num<-length(HM_Column_list)
    Cluster_Data<-data.frame()
    for (i in 1:Clust_Num){
      temp<-as.data.frame(HM_Column_list[[i]])
      colnames(temp)<-"Column_order_HeatMap"
      temp$Cluster<-paste0("Cluster_",i)
      Cluster_Data<-rbind(Cluster_Data,temp)
    }
    
    Labels<-Labels[match(Cluster_Data$Column_order_HeatMap,Labels$Column_order_HeatMap),]
    Cluster_Data<-merge(Labels,Cluster_Data,"Column_order_HeatMap")
    Cluster_Data<-Cluster_Data[match(Labels$Column_order_HeatMap,Cluster_Data$Column_order_HeatMap),]
    Cluster_Data<-Cluster_Data[,c(2,1,3)]
    colnames(Cluster_Data)[1]<-"Sample_ID"
    
    return(Cluster_Data)
    
  } else {
    Clust_Num<-length(HM_Column_list)
    Cluster_Data<-data.frame()
    for (i in 1:Clust_Num){
      temp<-as.data.frame(HM_Column_list[[i]])
      colnames(temp)<-"Column_order_HeatMap"
      Cluster_Data<-rbind(Cluster_Data,temp)
    }
    
    Labels<-Labels[match(Cluster_Data$Column_order_HeatMap,Labels$Column_order_HeatMap),]
    Cluster_Data<-Labels
    colnames(Cluster_Data)[1]<-"Sample_ID"
    
    return(Cluster_Data)
  } 
}

# Biological DMPs calling function (Human)
Bio_Diff_Meth_Analysis_Human<-function(targets_df,Col_Comparison,bVals_mat,mvals_mat,Diff_Limma,Annotation,Ref_Variable,Tresh,fit_2){
  Levels_Study<-levels(as.factor(targets_df[[Col_Comparison]]))
  Bio<-data.frame(CpG_ID = rownames(bVals_mat))
  
  for (i in 1:length(Levels_Study)){
    x<-Levels_Study[i]
    y<-targets_df[targets_df[[Col_Comparison]]==x,]
    z<-bVals_mat[,colnames(bVals_mat) %in% y[["Cell_Name"]]]
    z<-as.data.frame(rowMeans(z))
    colnames(z)<-x
    z$CpG_ID<-rownames(z)
    Bio<-merge(Bio,z,by="CpG_ID")
  }
  
  Bio$AB<-Bio[[Levels_Study[1]]] - Bio[[Levels_Study[2]]]
  Bio$Meth_Status<-ifelse(Bio$AB<0, "Hypomethylated","Hypermethylated")
  targets_df_2<-targets_df[targets_df[[Col_Comparison]] == Levels_Study[1],]
  bVals_mat_2<-bVals_mat[,colnames(bVals_mat) %in% targets_df_2[["Cell_Name"]]]
  dfVar<-data.frame(CpG_ID = rownames(bVals_mat_2),Meth_var =rowVars(as.matrix(bVals_mat_2), useNames = TRUE))
  df<-merge(Bio,dfVar,"CpG_ID")
  temp<-as.data.frame(Diff_Limma)
  colnames(temp)<-"Meth_diff"
  temp$CpG_ID<-rownames(temp)
  df<-merge(df,temp,"CpG_ID")
  
  DMPs_df<-df[df$Meth_diff!= 0 &
                abs(df$AB) >= Tresh &
                df$Meth_var <= 0.1,]
  
  DMPs_df<-DMPs_df[DMPs_df[[Ref_Variable]] >= 0.8 | 
                     DMPs_df[[Ref_Variable]] <=0.2,]
  
  AnnotationSub <- Annotation[match(rownames(mvals_mat),Annotation$CpG_ID),]
  DMPs <- topTable(fit_2, num=Inf, coef=1, genelist=AnnotationSub)
  DMPs <- merge(DMPs,DMPs_df[1:5],"CpG_ID")
  DMPs <- DMPs[c(1:6,8,13:15,7,10,11,16)]
  DMPs <- DMPs[DMPs$adj.P.Val <= 0.01,]
  DMPs$Meth_Status<-paste(DMPs$Meth_Status,Ref_Variable,sep="_")
  DMPs<-DMPs[order(-abs(DMPs$AB), DMPs$adj.P.Val),]
  return(DMPs)
}

# Biological DMPs function (Mouse)
Bio_Diff_Meth_Analysis_Mouse<-function(targets_df,Col_Comparison,bVals_mat,mvals_mat,Diff_Limma,Annotation,Ref_Variable,Tresh,fit_2){
  
  Levels_Study<-levels(as.factor(targets_df[[Col_Comparison]]))
  Bio<-as.data.frame(rownames(bVals_mat))
  colnames(Bio)<-"CpG_ID"
  for (i in 1:length(Levels_Study)){
    x<-Levels_Study[i]
    y<-targets_df[targets_df[[Col_Comparison]]==x,]
    z<-bVals_mat[,colnames(bVals_mat) %in% y[["Cell_Name"]]]
    z<-as.data.frame(rowMeans(z))
    colnames(z)<-x
    z$CpG_ID<-rownames(z)
    Bio<-merge(Bio,z,by="CpG_ID")
  }
  Bio$AB<-Bio[[Levels_Study[1]]] - Bio[[Levels_Study[2]]]
  Bio$Meth_Status<-ifelse(Bio$AB<0, "Hypomethylated","Hypermethylated")
  targets_df_2<-targets_df[targets_df[[Col_Comparison]] == Levels_Study[1],]
  bVals_mat_2<-bVals_mat[,colnames(bVals_mat) %in% targets_df_2[["Cell_Name"]]]
  dfVar<-as.data.frame(rowVars(as.matrix(bVals_mat_2), useNames = TRUE))
  colnames(dfVar)<-"Meth_var"
  dfVar$CpG_ID<-rownames(dfVar)
  df<-merge(Bio,dfVar,"CpG_ID")
  temp<-as.data.frame(Diff_Limma)
  colnames(temp)<-"Meth_diff"
  temp$CpG_ID<-rownames(temp)
  df<-merge(df,temp,"CpG_ID")
  DMPs_df<-df[df$Meth_diff!= 0 &
                abs(df$AB) >= Tresh &
                df$Meth_var <= 0.1,]
  AnnotationSub <- Annotation[match(rownames(mvals_mat),Annotation$CpG_ID),]
  DMPs <- topTable(fit_2, num=Inf, coef=1, genelist=AnnotationSub)
  colnames(DMPs)[1]<-"CpG_ID"
  DMPs <- merge(DMPs,DMPs_df[1:5],"CpG_ID")
  DMPs$Coordinate<-paste0(DMPs$CHR,":",DMPs$MAPINFO,"-",DMPs$MAPINFO)
  DMPs <- DMPs[c(1,53,42,41,39,40,43,44,50,49,51,46,47,52)]
  DMPs <- DMPs[DMPs$adj.P.Val <= 0.05,]
  DMPs$Meth_Status<-paste(DMPs$Meth_Status,Ref_Variable,sep="_")
  DMPs<-DMPs[order(-abs(DMPs$AB), DMPs$adj.P.Val),]
  return(DMPs)
}

# Statistical DMPs analysis

Supervised_Disease_DMPs_Analysis<-function(Disease,bVals_df,mVals_df, targets_df, Tresh, Organism){
  
  targets_df$Disease_Study<-ifelse(targets_df$Disease == Disease, "Disease", "Other")
  
  Status <- factor(targets_df$Disease_Study)
  Cell_Type<-factor(targets_df$Cell_Type)
  Sentrix_Pos<-factor(targets_df$Sentrix_Position)
  
  # Design the matrix
  design <- model.matrix(~0+Status +
                           Cell_Type +
                           Sentrix_Pos,
                         data=targets_df)
  
  colnames(design) <- c(levels(Status),
                        levels(Cell_Type)[-1],
                        levels(Sentrix_Pos)[-1])
  
  
  # Fit the linear model 
  fit <- lmFit(mVals_df, design)
  
  # Create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(Disease - Other,
                              levels=design)
  
  # Fit the contrasts
  fit_2 <- contrasts.fit(fit, contMatrix)
  fit_2 <- eBayes(fit_2)
  
  Diff <- decideTests(fit_2)
  summary(Diff)
  
  if (Organism == "Human") {
    
    Tresh<-ifelse(Disease == "Multiple_Myeloma", 0.6, 0.4)
    
    DMPs<-Bio_Diff_Meth_Analysis_Human(targets_df=targets_df,
                                       Col_Comparison="Disease_Study",
                                       bVals_mat=bVals_df,
                                       mvals_mat=mVals_df,
                                       Diff_Limma=Diff,
                                       Annotation=ann850k,
                                       Ref_Variable="Disease",
                                       Tresh = Tresh,
                                       fit_2 = fit_2)
  }
  
  else if (Organism == "Mouse"){
    
    Tresh<-ifelse(Disease == "Multiple_Myeloma", 0.8, 0.4)
    
    DMPs<-Bio_Diff_Meth_Analysis_Mouse(targets_df=targets_df,
                                       Col_Comparison="Disease_Study",
                                       bVals_mat=bVals_df,
                                       mvals_mat=mVals_df,
                                       Diff_Limma=Diff,
                                       Annotation=annMous,
                                       Ref_Variable="Disease",
                                       Tresh = Tresh,
                                       fit_2 = fit_2)
    
    
  }
  
  return(DMPs)
}

# Enrichment Analysis Function
Gene_Enrichment_Analysis<-function(DMPs,Organism){
  
  Results_list<-list()
  
  if(Organism == "Human"){
    
    DMPs<-DMPs[!DMPs$UCSC_RefGene_Name=="",]
    Global_Genes<-unique(unlist(str_split(DMPs$UCSC_RefGene_Name,";")))
    
  }
  
  else if (Organism == "Mouse"){
    
    DMPs<-DMPs[!DMPs$UCSC_Gene=="",]
    Global_Genes<-unique(unlist(str_split(DMPs$UCSC_Gene,";")))
    
  }
  
  EnrichR_Table<-EnricheR_Table(Global_Genes)
  
  Results_list[["Table"]]<-EnrichR_Table
  
  EnrichR_Table$Count<-as.numeric(gsub("\\/.*","",EnrichR_Table$Overlap))
  
  Plot<-ggplot(data = EnrichR_Table, aes(x = reorder(Term,-Gene_Ratio), y = Gene_Ratio, 
                                         fill = Adjusted.P.value, size = Count)) + 
    geom_point(shape = 21) +
    scale_fill_continuous(low="#882fb5", high="#ffcf33")+
    scale_size_continuous(range = c(6,10))+
    theme_bw() + 
    ylab("GO Terms") + 
    xlab("Gene Ratio")+
    theme_few()+
    theme(axis.text.x=element_text(angle=-270, hjust=1),axis.title.x = element_blank())
  
  Results_list[["Plot"]]<-Plot
  
  
  return(Results_list)
  
}

EnricheR_Table<-function(Genes_Vector){
  Enriched_Results<-enrichr(as.vector(Genes_Vector),dbs_selection)
  
  Enriched_Table_Results<-data.frame()
  
  for (i in 1:length(dbs_selection)) {
    Database<-dbs_selection[i]
    Enriched_Table<-Enriched_Results[[Database]]
    Enriched_Table$Database<-Database
    Enriched_Table<-Enriched_Table[,-c(5:8)]
    Enriched_Table$Obs_Genes<-as.numeric(gsub("/..*","",Enriched_Table$Overlap))
    Enriched_Table$Set_Genes<-as.numeric(gsub("..*/","",Enriched_Table$Overlap))
    Enriched_Table<-Enriched_Table[Enriched_Table$Adjusted.P.value<=0.05,]
    if(!nrow(Enriched_Table)==0)
      Enriched_Table_Results<-rbind(Enriched_Table_Results,Enriched_Table)
    else{
      next 
    }
  }
  
  if(nrow(Enriched_Table_Results)>0){
    Enriched_Table_Results$Gene_Ratio<-as.numeric(Enriched_Table_Results$Obs_Genes/Enriched_Table_Results$Set_Genes)
    Enriched_Table_Results<-Enriched_Table_Results[,c(6,1:5,9)]
    Enriched_Table_Results<-Enriched_Table_Results[order(-Enriched_Table_Results$Gene_Ratio),]
    return(Enriched_Table_Results)
  } else {
    return("No differential pathways found")
  }
}


# Pie plot function
Pie_plot_function<-function(column_name,Annotation_file,Status){
  Var1 <- c("CpG_Island","CpG_Shore","CpG_Shelf","Open_Sea",
            "5_Regulatory_Region","Body","3'UTR","Intergenic",
            "B_ALL","Burkitt","DLBCL","Hodgkin_Lymphoma","Lymphoblastic_Leukemia_Others",
            "MBCL","CML","AML","NH_B_lymphoma_Others","T_ALL","T_cell_Non_Hodgkin_Lymphoma","Multiple_Myeloma")
  
  Colours<-c("gold","orchid2","bisque3","lightskyblue",
             "lightsalmon","darkolivegreen3","blueviolet","cyan3",
             "#c22f87","#5E3203","#F37612","#F5ADC5", "#719CB9", 
             "#ad8f57","#A9D39F","#4A904A","#E3E500","#C7F7FF",
             "#3957c4","#CE5755") 
                   
  To_Merge_Table<-data.frame(Var1,Colours)  
  
  if(Status=="CGI"){
    
    Annotation_file$CGI_CAT<-ifelse(Annotation_file[[column_name]] %like% "Shore","CpG_Shore",
                             ifelse(Annotation_file[[column_name]] %like% "Shelf","CpG_Shelf",
                             ifelse(Annotation_file[[column_name]] %like% "Island","CpG_Island","Open_Sea")))
    
    plot_table<-as.data.frame(table(Annotation_file[["CGI_CAT"]]))
    plot_table$Percentage <-round(((plot_table$Freq)/colSums(plot_table[2])*100),2)
    
  } else if(Status =="Gene"){
    
    Annotation_file$Gene_CAT<-ifelse(Annotation_file[[column_name]] %like% "TSS1500","5_Regulatory_Region",
                              ifelse(Annotation_file[[column_name]] %like% "TSS200","5_Regulatory_Region",       
                              ifelse(Annotation_file[[column_name]] %like% "5'UTR","5_Regulatory_Region",        
                              ifelse(Annotation_file[[column_name]] %like% "1stExon","5_Regulatory_Region",
                              ifelse(Annotation_file[[column_name]] %like% "Body","Body",
                              ifelse(Annotation_file[[column_name]] %like% "ExnBnd","Body",
                              ifelse(Annotation_file[[column_name]] %like% "3'UTR","3'UTR","Intergenic")))))))
    
    plot_table<-as.data.frame(table(Annotation_file[["Gene_CAT"]]))
    plot_table$Percentage <-round(((plot_table$Freq)/colSums(plot_table[2])*100),2)
    
  } else if (Status == "Disease"){
    
    plot_table<-as.data.frame(table(Annotation_file[["Disease"]]))
    plot_table$Percentage <-round(((plot_table$Freq)/colSums(plot_table[2])*100),2)
    
  }
  
  plot_table<-merge(plot_table,To_Merge_Table,"Var1")
  Plot<-ggplot(plot_table, aes(x = "", y = Percentage, fill = Var1)) +
    scale_fill_manual(values = c("Coding_Gene" = "orange",
                                 "Non_Coding_Gene" = "#b5b2ed",
                                 "Intergenic" = "#afdced",
                                 "CpG_Island" = "#faee84",
                                 "CpG_Shore" = "#ebaaa4",
                                 "CpG_Shelf" = "#e0a9eb",
                                 "Open_Sea" = "#ebcba4",
                                 "5_Regulatory_Region" = "#8ff2c6",
                                 "Body" = "#b5b2ed",
                                 "3'UTR" = "#c9eb9d",
                                 "Leukaemia" = "#F3B61F",
                                 "Lymphoma" = "#15A390",
                                 "Multiple_Myeloma" = "#CE5755",
                                 "B_ALL" = "#c22f87",
                                 "Burkitt" = "#5E3203",
                                 "DLBCL" = "#F37612",
                                 "Hodgkin_Lymphoma" = "#F5ADC5",
                                 "Lymphoblastic_Leukemia_Others" = "#719CB9",
                                 "MBCL" = "#ad8f57",              
                                 "CML" = "#A9D39F",
                                 "AML" = "#4A904A",          
                                 "NH_B_lymphoma_Others" = "#E3E500",                         
                                 "T_ALL" = "#C7F7FF",
                                 "T_cell_Non_Hodgkin_Lymphoma" = "#3957c4",         
                                 "Transformed_cell_line" = "#C994F2",
                                 "Human" = "black",
                                 "Mouse" = "gray80",
                                 "NA" = "white",
                                 "Multiple_Myeloma" = "#CE5755"))+
    geom_col(color = "black",size = 0.6) +
    coord_polar(theta = "y")+
    theme_few()+
    theme(legend.position = "none")
  
  return(list(Plot=Plot,table=plot_table))
}

# Function to perfrom differential analysis of cancer cell lines 
# vs normal isollated cells

Supervised_Cancer_Normal_DMPs_Analysis<-function(bVals_df,mVals_df,targets_df,Organism){
  
  Status <- factor(targets_df$Type)
  Sentrix_Pos<-factor(targets_df$Sentrix_Position)
  
  # Design the matrix
  design <- model.matrix(~0+Status +
                           Sentrix_Pos,
                         data=targets_df)
  
  colnames(design) <- c(levels(Status),
                        levels(Sentrix_Pos)[-1])
  
  # Fit the linear model 
  fit <- lmFit(mVals_df, design)
  
  # Create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(Cancer_Cells - Normal_Cells,
                              levels=design)
  
  # Fit the contrasts
  fit_2 <- contrasts.fit(fit, contMatrix)
  fit_2 <- eBayes(fit_2)
  
  Diff <- decideTests(fit_2)
  summary(Diff)
  
  if (Organism == "Human") {
    
    Tresh<-0.4
    
    DMPs<-Bio_Diff_Meth_Analysis_Human(targets_df=targets_df,
                                       Col_Comparison="Type",
                                       bVals_mat=bVals_df,
                                       mvals_mat=mVals_df,
                                       Diff_Limma=Diff,
                                       Annotation=ann850k,
                                       Ref_Variable="Cancer_Cells",
                                       Tresh = Tresh,
                                       fit_2 = fit_2)
  }
  
  else if (Organism == "Mouse"){
    
    Tresh<-0.4
    
    DMPs<-Bio_Diff_Meth_Analysis_Mouse(targets_df=targets_df,
                                       Col_Comparison="Type",
                                       bVals_mat=bVals_df,
                                       mvals_mat=mVals_df,
                                       Diff_Limma=Diff,
                                       Annotation=annMous,
                                       Ref_Variable="Cancer_Cells",
                                       Tresh = Tresh,
                                       fit_2 = fit_2)
  }
  
  return(DMPs)
  
}

# Classifier functions
limma_function<-function(targets_df){
  
  Comparison <- factor(targets_df$Comparison, levels = c("RefDisease","Others"))
  Sentrix_Pos<-factor(targets_df$Sentrix_Position)
  
  # Design the matrix
  design <- model.matrix(~ 0 + Comparison +
                           Sentrix_Pos,
                         data=targets_df)
  
  colnames(design) <- c(levels(Comparison),
                        levels(Sentrix_Pos)[-1])
  # Fit the linear model 
  fit <- lmFit(mVals_Class, design)
  # Create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(RefDisease - Others,
                              levels=design)
  
  fit_2 <- contrasts.fit(fit, contMatrix)
  fit_2 <- eBayes(fit_2)
  
  Diff <- decideTests(fit_2)
  
  Difftemp<-as.data.frame(Diff)
  Difftemp<-Difftemp[Difftemp$`RefDisease - Others` != 0,]
  
  if(length(Difftemp) < 1){
    
    return("There is no statistical DMP")
    
  }
  
  return(list(Diff_df = as.data.frame(Diff),
              fit_2 = fit_2))
}

Diff_Meth_Analysis_Classifier<-function(targets_df,Diff_Limma,Annotation,Tresh,fit_2){
  
  Levels_Study<-factor(unique(targets_df$Comparison), levels = c("RefDisease","Others"))
  
  Bio<-data.frame(CpG_ID = rownames(bVals_Class))
  
  for (i in 1:length(Levels_Study)){
    x<-Levels_Study[i]
    y<-targets_df[targets_df$Comparison == x,]
    z<-bVals_Class[,colnames(bVals_Class) %in% y$Cell_Name]
    z<-as.data.frame(rowMeans(z))
    colnames(z)<-x
    z$CpG_ID<-rownames(z)
    Bio<-merge(Bio,z,by="CpG_ID")
  }
  
  Bio$AB<-Bio$RefDisease - Bio$Others
  Bio$Meth_Status<-ifelse(Bio$AB<0, "Hypomethylated","Hypermethylated")
  targets_df_2<-targets_df[targets_df$Comparison == "RefDisease",]
  bVals_mat_2<-bVals_Class[,colnames(bVals_Class) %in% targets_df_2$Cell_Name]
  dfVar<-data.frame(CpG_ID = rownames(bVals_mat_2),Meth_var =rowVars(bVals_mat_2, useNames = TRUE))
  df<-merge(Bio,dfVar,"CpG_ID")
  temp<-as.data.frame(Diff_Limma)
  colnames(temp)<-"Meth_diff"
  temp$CpG_ID<-rownames(temp)
  df<-merge(df,temp,"CpG_ID")
  
  DMPs_df<-df[df$Meth_diff!= 0 &
                abs(df$AB) >= Tresh &
                df$Meth_var <= 0.1,]
  
  if (nrow(DMPs_df)<1){
    
    return("There is no biological DMP")
  }
  
  DMPs_df<-DMPs_df[DMPs_df$RefDisease >= 0.8 | 
                     DMPs_df$RefDisease <=0.2,]
  
  if (nrow(DMPs_df)<1){
    
    return("There is no biological DMP")
  }
  
  AnnotationSub <- Annotation[match(rownames(mVals_Class),Annotation$CpG_ID),]
  DMPs <- topTable(fit_2, num=Inf, coef=1, genelist=AnnotationSub)
  DMPs <- merge(DMPs,DMPs_df[1:5],"CpG_ID")
  DMPs <- DMPs[c(1:6,13:16,7,10,11)]
  DMPs <- DMPs[DMPs$adj.P.Val <= 0.01,]
  DMPs<-DMPs[order(-abs(DMPs$AB), DMPs$adj.P.Val),]
  return(DMPs)
}

# Functions for Primary vs Cell lines analyses
limma_function_Cell_Primaries<-function(targets_df, mVals_df){
  
  Comparison <- factor(targets_df$Type, levels = c("Primary_Sample","Cell_Line"))
  Sentrix_Pos<-factor(targets_df$Sentrix_Position)
  
  # Design the matrix
  design <- model.matrix(~ 0 + Comparison +
                           Sentrix_Pos,
                         data=targets_df)
  
  colnames(design) <- c(levels(Comparison),
                        levels(Sentrix_Pos)[-1])
  # Fit the linear model 
  fit <- lmFit(mVals_df, design)
  # Create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(Primary_Sample - Cell_Line,
                              levels=design)
  
  fit_2 <- contrasts.fit(fit, contMatrix)
  fit_2 <- eBayes(fit_2)
  
  Diff <- decideTests(fit_2)
  
  Difftemp<-as.data.frame(Diff)
  Difftemp$CpG_Name<-rownames(Difftemp)
  Difftemp<-Difftemp[Difftemp$`Primary_Sample - Cell_Line` != 0,]
  
  if(length(Difftemp) < 1){
    
    return("There is no statistical DMP")
    
  }
  
  return(list(Diff_df = as.data.frame(Diff),
              fit_2 = fit_2))
  
}



Diff_Meth_Analysis_Cells_Primaries<-function(targets_df,bVals_df,Diff_Limma,Annotation,Tresh,fit_2,mVals_df){
  
  Levels_Study<-factor(unique(targets_df$Type), levels = c("Primary_Sample","Cell_Line"))
  
  Bio<-data.frame(CpG_ID = rownames(bVals_df))
  
  for (i in 1:length(Levels_Study)){
    x<-Levels_Study[i]
    y<-targets_df[targets_df$Type == x,]
    z<-bVals_df[,colnames(bVals_df) %in% y$Cell_Name]
    z<-as.data.frame(rowMeans(z))
    colnames(z)<-x
    z$CpG_ID<-rownames(z)
    Bio<-merge(Bio,z,by="CpG_ID")
  }
  
  Bio$AB<-Bio$Primary_Sample - Bio$Cell_Line
  Bio$Meth_Status<-ifelse(Bio$AB<0, "Hypomethylated","Hypermethylated")
  targets_df_2<-targets_df[targets_df$Type == "Primary_Sample",]
  bVals_mat_2<-bVals_df[,colnames(bVals_df) %in% targets_df_2$Cell_Name]
  dfVar<-data.frame(CpG_ID = rownames(bVals_mat_2),Meth_var =rowVars(as.matrix(bVals_mat_2), useNames = TRUE))
  df<-merge(Bio,dfVar,"CpG_ID")
  temp<-as.data.frame(Diff_Limma)
  colnames(temp)<-"Meth_diff"
  temp$CpG_ID<-rownames(temp)
  df<-merge(df,temp,"CpG_ID")
  
  DMPs_df<-df[df$Meth_diff!= 0 &
                abs(df$AB) >= Tresh &
                df$Meth_var <= 0.1,]
  
  if (nrow(DMPs_df)<1){
    
    return("There is no biological DMP")
  }
  
  DMPs_df<-DMPs_df[DMPs_df$Primary_Sample >= 0.8 | 
                     DMPs_df$Primary_Sample <=0.2,]
  
  if (nrow(DMPs_df)<1){
    
    return("There is no biological DMP")
  }
  
  AnnotationSub <- Annotation[match(rownames(mVals_df),Annotation$CpG_ID),]
  DMPs <- topTable(fit_2, num=Inf, coef=1, genelist=AnnotationSub)
  DMPs <- merge(DMPs,DMPs_df[1:5],"CpG_ID")
  DMPs <- DMPs[c(1:6,13:16,7,10,11)]
  DMPs <- DMPs[DMPs$adj.P.Val <= 0.01,]
  DMPs<-DMPs[order(-abs(DMPs$AB), DMPs$adj.P.Val),]
  bVals_df<-bVals_df[rownames(bVals_df) %in% DMPs$CpG_ID,]
  
  return(list(DMPs = DMPs,
              bVals_Filt = bVals_df))
}


# Drug Analysis functions
Correlation_Function<-function(Genes,DMPs_Meth,DEG_Expr,Meth_mx,RNA_mx){
  DMPs_Meth<-DMPs_Meth[c(1,4,6)]
  DMPs_Meth<-DMPs_Meth[DMPs_Meth$UCSC_RefGene_Name != "",]
  DMPs_Meth<-unique(separate_rows(DMPs_Meth,UCSC_RefGene_Name,UCSC_RefGene_Group,sep = ";"))
  DMPs_Meth<-DMPs_Meth %>% 
    group_by(CpG_ID,UCSC_RefGene_Name) %>% 
    summarise(UCSC_RefGene_Group = paste(UCSC_RefGene_Group, collapse=";"))
  DMPs_Meth$UCSC_RefGene_Group<-ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "TSS1500","TSS1500",
                                ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "TSS200","TSS200",        
                                ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "5'UTR","5'UTR",
                                ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "1stExon","1stExon",
                                ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "Body","Body",
                                ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "ExonBnd","Body",
                                ifelse(DMPs_Meth$UCSC_RefGene_Group %like% "3'UTR","3'UTR","Intergenic")))))))
  
  DMPs_Meth<-DMPs_Meth[DMPs_Meth$UCSC_RefGene_Name %in% Genes,]
  
  if (nrow(DMPs_Meth)<2) {
    Meth_mx<-as.data.frame((Meth_mx[rownames(Meth_mx) %in% DMPs_Meth$CpG_ID,]))
    colnames(Meth_mx)<-DMPs_Meth$CpG_ID
  } else {
    Meth_mx<-as.data.frame(t(Meth_mx[rownames(Meth_mx) %in% DMPs_Meth$CpG_ID,]))  
  }
  
  Meth_mx$CELL_LINE_MOD<-toupper(gsub('[^[:alnum:] ]','', rownames(Meth_mx)))
  
  RNA_mx<-as.data.frame(t(RNA_mx[rownames(RNA_mx) %in% Genes,]))
  RNA_mx$CELL_LINE_MOD<-toupper(gsub('[^[:alnum:] ]','', rownames(RNA_mx)))
  
  Corr_df<-data.frame()
  
  for (i in 1:length(Genes)) {
    
    Gene<-Genes[i]
    DMPs_Meth_Gene<-DMPs_Meth[DMPs_Meth$UCSC_RefGene_Name == Gene,]
    Meth_mx_Gene<-Meth_mx[,colnames(Meth_mx) %in% c("CELL_LINE_MOD",DMPs_Meth_Gene$CpG_ID)]
    RNA_mx_Gene<-RNA_mx[,colnames(RNA_mx) %in% c("CELL_LINE_MOD",Gene)]
    Corr_mx<-merge(Meth_mx_Gene,RNA_mx_Gene,"CELL_LINE_MOD")
    
    for (j in 1:(ncol(Corr_mx)-2)) {
      
      CpG<-colnames(Corr_mx)[j+1]
      Group<-DMPs_Meth_Gene[DMPs_Meth_Gene$CpG_ID == CpG,]
      Group<-Group$UCSC_RefGene_Group
      Corr<-cor.test(Corr_mx[[CpG]],Corr_mx[[Gene]],method="spearman", exact = NULL)  
      df<-data.frame(CpG = CpG, Gene = Gene, CpG_Relation_Gene = Group, 
                     Corr.R = Corr$estimate, 
                     p.val = Corr$p.value, 
                     Sig = ifelse(Corr$p.value <=0.05 & abs(Corr$estimate) >=0.3, TRUE,FALSE))
      Corr_df<-rbind(Corr_df,df)
    }
    
  }
  
  return(Corr_df)
  
}


Diff_Meth_Expression_function<-function(Drug,mVals,targetsdf,Annotation,Analysis){
  
  design <- model.matrix(~+targetsdf[[Drug]],
                         data=targetsdf)
  colnames(design)[2]<-Drug
  fit <- lmFit(mVals, design)
  fit <- eBayes(fit)
  
  Diff <- decideTests(fit)
  
  if(Annotation != FALSE){
    Annotation<-Annotation[Annotation$CpG_ID %in% rownames(Diff),]
    Annotation<-Annotation[match(rownames(Diff),Annotation$CpG_ID),]
  } 
  
  Res_Meth<-topTable(fit, num=Inf, coef=2, genelist=Annotation)
  
  
  Res_Meth<-Res_Meth[Res_Meth$adj.P.Val<=0.01,]
  
  if (Analysis == "Methylation"){
    
    Res_Meth<-Res_Meth[,-c(9,12)]
    
  } else if (Analysis == "Expression"){
    
    Res_Meth<-Res_Meth[,c(2,3,5,6)]
    
  }
  
  Res_Meth<-Res_Meth[order(Res_Meth$adj.P.Val),]
  
  return(Res_Meth)
  
}

Venn_Diagram_function<-function(DMPs_Meth,DEG_Expr,Drug){
  Meth_Genes<-unique(unlist(str_split(DMPs_Meth$UCSC_RefGene_Name,pattern = ";")))
  RNA_Genes<-rownames(DEG_Expr)
  
  return(intersect(Meth_Genes,RNA_Genes))
}


Draw_HeatMap_No_Rows_Clust<-function(bVals_Cluster,bVals_HeatMap, bottom_annotation,ColSplit){
  Result_list<-list()
  distance <- dist(t(bVals_Cluster), method="euclidean") # calcular las distancias
  clust.model <- hclust(distance, method="ward.D2") # clustering
  
  Heatmap<-draw(Heatmap(as.matrix(bVals_HeatMap), name = "Methylation", 
                        col = gplots::greenred(75), 
                        bottom_annotation = bottom_annotation,
                        cluster_columns = clust.model,
                        cluster_rows = FALSE,
                        column_split= ColSplit, 
                        #row_split = 2,
                        show_row_names = TRUE, 
                        show_column_names = TRUE, 
                        show_row_dend = FALSE,
                        column_names_gp = grid::gpar(fontsize = 8)))
  
  Result_list[["HeatMap"]] <- Heatmap

  Cluster_Columns<-Cluster_Order_Columns(clust.model,Heatmap)
  Result_list[["Cluster_Columns"]] <- Cluster_Columns

  return(Result_list)
}
