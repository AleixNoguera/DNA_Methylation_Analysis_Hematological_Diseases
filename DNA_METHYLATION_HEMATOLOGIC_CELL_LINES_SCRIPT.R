# Set working directory
setwd("./")
# read functions script
source("./DNA_METHYLATION_HEMATOLOGIC_CELL_LINES_FUNCTIONS")

# Load necessary libraries
library (sesame)
library (IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library (ComplexHeatmap)
library (Rtsne)
library (dplyr)
library (ggplot2)
library (ggthemes)
library (data.table)
library (treeio)
library (ggtree)
library (ENmix)
library (limma)
library (enrichR)
library (stringr)
library (tidyr)

# Create annotation objects
ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
colnames(ann850k)[4]<-"CpG_ID"

annMous<-as.data.frame(data.table::fread("./mm10_cpg_annotation.csv"))
## Preprocess the annotation
annMous$CpG_Location_CGI<-ifelse(!annMous$CpG_Island == "", "CpG_Island",
                          ifelse(!annMous$S_Shore == "","CpG_Shore",
                          ifelse(!annMous$N_Shore == "","CpG_Shore",
                          ifelse(!annMous$N_Shelf == "","CpG_Shelf",
                          ifelse(!annMous$S_Shelf == "","CpG_Shelf","Open_Sea")))))
colnames(annMous)[c(1,16)]<-c("CpG_ID","chr")

# Open sample sheet
## Human
targets_Human<-readxl::read_xlsx("./Human_Hematologica_Cell_Lines.xlsx")

#### Create a column with cell names in capital letters and without punctuation 
#### marks for merging with other datasets
targets_Human$CELL_LINE_NAME_MOD<-toupper(gsub('[^[:alnum:] ]', '', targets_Human$Cell_Name))

## Mouse
targets_Mouse<-readxl::read_xls("./Mouse_Hematological_Cell_Lines.xlsx")
colnames(targets_Mouse)[1]<-"Cell_Name"


# SAMPLE PREPROCESSING
## Human Cell Lines
Human_Cell_Lines<-Preprocessing_function(IDATs_dir<-"./IDATs",
                                         Organism<-"Homo_Sapiens",
                                         targets_df<-targets_Human,
                                         Ann<-ann850k)

bVals_Human<-Human_Cell_Lines$bVals
bVals_Human_NI<-Human_Cell_Lines$bVals_NI

## Mouse Cell Lines
Mouse_Cell_Lines<-Preprocessing_function(IDATs_dir<-"./IDATs",
                                         Organism<-"Mus_Musculus",
                                         targets_df<-targets_Mouse,
                                         Ann<-annMous)

bVals_Mouse<-Mouse_Cell_Lines$bVals
bVals_Mouse_NI<-Mouse_Cell_Lines$bVals_NI

################################################################################
########################## Figure 1B ###########################################
################################################################################

# Unsupervised HeatMap Human

## Design top annotation
bVals_Human<-bVals_Human[,match(targets_Human$Cell_Name,colnames(bVals_Human))]

ColAnn<-HeatmapAnnotation(Disease = targets_Human[["Disease"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

Human_Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_Human,
                                         nProbes<-0.01,
                                         targetsdf<-targets_Human,
                                         BoottomAnn<-ColAnn,
                                         Col_Split<-3)
# Unsupervised HeatMap Mouse

## Design top annotation
bVals_Mouse<-bVals_Mouse[,match(targets_Mouse$Cell_Name,colnames(bVals_Mouse))]

ColAnn<-HeatmapAnnotation(Disease = targets_Mouse[["Disease"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

Mouse_Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_Mouse,
                                         nProbes<-0.01,
                                         targetsdf<-targets_Mouse,
                                         BoottomAnn<-ColAnn,
                                         Col_Split<-3)

################################################################################
##################### Supplementary Figure 1A ##################################
################################################################################

# Unsupervised HeatMap Human (with no imputation)
## Remove probes with at least 1 NA value
bVals_Human_NI<-bVals_Human_NI[complete.cases(bVals_Human_NI),]

## Design top annotation
bVals_Human_NI<-bVals_Human_NI[,match(targets_Human$Cell_Name,
                                      colnames(bVals_Human_NI))]

ColAnn<-HeatmapAnnotation(Disease = targets_Human[["Disease"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

Human_Unsupervised_NI<-Unsupervised_HeatMap(bdf<-bVals_Human_NI,
                                            nProbes<-0.01,
                                            targetsdf<-targets_Human,
                                            BoottomAnn<-ColAnn,
                                            Col_Split<-3)
# Unsupervised HeatMap Mouse (with no imputation)
## Remove probes with at least 1 NA value
bVals_Mouse_NI<-bVals_Mouse_NI[complete.cases(bVals_Mouse_NI),]

## Design top annotation
bVals_Mouse_NI<-bVals_Mouse_NI[,match(targets_Mouse$Cell_Name,
                                      colnames(bVals_Mouse_NI))]

ColAnn<-HeatmapAnnotation(Disease = targets_Mouse[["Disease"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

Mouse_Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_Mouse_NI,
                                         nProbes<-0.01,
                                         targetsdf<-targets_Mouse,
                                         BoottomAnn<-ColAnn,
                                         Col_Split<-3)


################################################################################
##################### Supplementary Figure 2A ##################################
################################################################################

# Human t-SNE
bVals_Human<-bVals_Human[,match(targets_Human$Cell_Name,colnames(bVals_Human))]

set.seed(785)
tSNE_Human <- Rtsne(t(bVals_Human), perplexity=25, check_duplicates = FALSE)

tSNE_Human<-tSNE_Process_function(TSNEObj<-tSNE_Human, 
                                  targets_df<-targets_Human, 
                                  Column_Col<-"Disease",
                                  Organism = "Human")

# Mouse t-SNE
bVals_Mouse<-bVals_Mouse[,match(targets_Mouse$Cell_Name,colnames(bVals_Mouse))]

set.seed(123)
tSNE_Mouse <- Rtsne(t(bVals_Mouse), perplexity=5, check_duplicates = FALSE)

tSNE_Mouse<-tSNE_Process_function(TSNEObj<-tSNE_Mouse, 
                                  targets_df<-targets_Mouse, 
                                  Column_Col<-"Disease",
                                  Organism = "Mouse")

# Mouse and Human t-SNE
## Load a list of already liftovered CpGs
CpG_Lift<-readxl::read_xlsx("./CpGs_Lift_to_Genome.xlsx")

## Prepare the bvals matrix
df<-as.data.frame(bVals_Human)
df$Human_CpG_ID<-rownames(df)

df2<-as.data.frame(bVals_Mouse)
df2$Mouse_CpG_ID<-rownames(df2)

bVals_HM<-as.data.frame(merge(CpG_Lift,df, "Human_CpG_ID"))
bVals_HM<-as.data.frame(merge(bVals_HM,df2,"Mouse_CpG_ID"))
rownames(bVals_HM)<-bVals_HM$Mouse_CpG_ID

bVals_HM<-bVals_HM[-c(1:2)]

# Prepare the targets data frame
df<-targets_Human[,colnames(targets_Human) %in% c("Cell_Name", "Disease", 
                                                  "Final_Classification")]
df$Organism<-"Human"

df1<-targets_Mouse[,colnames(targets_Mouse) %in% c("Cell_Name", "Disease")]
df1$Final_Classification<-""
df1$Organism<-"Mouse"

targets_HM<-rbind(df,df1)

targets_HM<-targets_HM[match(colnames(bVals_HM),targets_HM$Cell_Name),]

set.seed(123)
tSNE_HM <- Rtsne(t(bVals_HM), perplexity=5, check_duplicates = FALSE)

tSNE_HM<-tSNE_Process_function(TSNEObj<-tSNE_HM, 
                               targets_df<-targets_HM, 
                               Column_Col<-"Disease",
                               Organism = "Mouse_Human")

################################################################################
########################### Figure 2C ##########################################
################################################################################

Human_Mouse_Dend<-Dendogram_Function(bVals=bVals_HM,
                                     targets_df=targets_HM)

################################################################################
########################### Figure 2D ##########################################
################################################################################

# Differential methylation between Leukemia, Lymphoma and Multiple Myeloma
## Human
### Prepare Annotation
ann850k$Coordinate<-paste0(ann850k$chr,":",ann850k$pos,"-",ann850k$pos)
ann850k<-ann850k[,c(4,47,19,22:24)]

### Remove transformed cell lines
targets_LLM<-targets[targets$Disease != "Transformed_cell_line",]
bVals_Human_LLM<-bVals_Human[,colnames(bVals_Human) %in% targets_LLM$Cell_Name]
mVals_Human_LLM<-B2M(bVals_Human_LLM)

Human_DMPs_List<-list()

for (i in 1:length(unique(targets_LLM$Disease))) {
  
  Disease<-unique(targets_LLM$Disease)[i]
  
  print(Disease)
  
  DMPs<-Supervised_Disease_DMPs_Analysis(targets_df = targets_LLM,
                                         bVals_df = bVals_Human_LLM,
                                         mVals_df = mVals_Human_LLM,
                                         Disease = Disease,
                                         Organism = "Human")
  
  Human_DMPs_List[[Disease]]<-DMPs
  
}

# Supplementary table 2
H_Leu_DMPs<-Human_DMPs_List$Leukaemia
H_Lym_DMPs<-Human_DMPs_List$Lymphoma
H_MM_DMPs<-Human_DMPs_List$Multiple_Myeloma

# Supervised Heatmap
bVals_Disease_H<-bVals_Human_LLM[rownames(bVals_Human_LLM) %in% 
                                   c(H_Leu_DMPs$CpG_ID,
                                     H_Lym_DMPs$CpG_ID,
                                     H_MM_DMPs$CpG_ID),]

bVals_Disease_H<-bVals_Disease_H[,match(targets_LLM$Cell_Name, 
                                        colnames(bVals_Disease_H))]

ColAnn<-HeatmapAnnotation(Disease = targets_LLM[["Disease"]],
                         col = list("Disease" = Colours_Pal),
                         border = TRUE,
                         simple_anno_size = unit(0.75, "cm"))


HeatMap_Human_Disease<-Draw_HeatMap(bVals_Cluster = bVals_Disease_H,
                                    bVals_HeatMap = bVals_Disease_H, 
                                    bottom_annotation = ColAnn,
                                    ColSplit = 3)



## Mouse
### Remove transformed cell lines
targets_M_LLM<-targets_Mouse[targets_Mouse$Disease != "Transformed_cell_line",]
colnames(targets_M_LLM)[13]<-c("Sentrix_Position")
bVals_Mouse_LLM<-bVals_Mouse[,colnames(bVals_Mouse) %in% 
                               targets_M_LLM$Cell_Name]
mVals_Mouse_LLM<-B2M(bVals_Mouse_LLM)


Mouse_DMPs_List<-list()

for (i in 1:length(unique(targets_M_LLM$Disease))) {
  
  Disease<-unique(targets_M_LLM$Disease)[i]
  
  print(Disease)
  
  DMPs<-Supervised_Disease_DMPs_Analysis(targets_df = targets_M_LLM,
                                         bVals_df = bVals_Mouse_LLM,
                                         mVals_df = mVals_Mouse_LLM,
                                         Disease = Disease,
                                         Organism = "Mouse")
  
  Mouse_DMPs_List[[Disease]]<-DMPs
  
}

# Supplementary table 2
M_Leu_DMPs<-Mouse_DMPs_List$Leukaemia
M_Lym_DMPs<-Mouse_DMPs_List$Lymphoma
M_MM_DMPs<-Mouse_DMPs_List$Multiple_Myeloma

# Supervised Heatmap
bVals_Disease_M<-bVals_Mouse_LLM[rownames(bVals_Mouse_LLM) %in% 
                                   c(M_Leu_DMPs$CpG_ID,M_Lym_DMPs$CpG_ID,
                                     M_MM_DMPs$CpG_ID),]
bVals_Disease_M<-bVals_Disease_M[,match(targets_M_LLM$Cell_Name, 
                                        colnames(bVals_Disease_M))]

ColAnn<-HeatmapAnnotation(Disease = targets_M_LLM[["Disease"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))


HeatMap_Mouse_Disease<-Draw_HeatMap(bVals_Cluster = bVals_Disease_M,
                                    bVals_HeatMap = bVals_Disease_M, 
                                    bottom_annotation = ColAnn,
                                    ColSplit = 3)


################################################################################
################## Supplementary Figure 2B #####################################
################################################################################

# Bind all the DMPs in a single file
H_DMPs_Disease<-rbind(H_Leu_DMPs,H_Lym_DMPs,H_MM_DMPs)

# Summary Pie Plot
## Pie plot in relation to CpG island
Human_CGI<-Pie_plot_function("Relation_to_Island",H_DMPs_Disease,"CGI")
## Pie plot in relation to Gene location
Human_GENE<-Pie_plot_function("UCSC_RefGene_Group",H_DMPs_Disease,"Gene")

# EnrichR Analysis
setEnrichrSite("Enrichr")
websiteLive<-TRUE

# Select the databases we are interested in
dbs_selection<-c("GO_Biological_Process_2023","GO_Molecular_Function_2023",
                 "GO_Cellular_Component_2023")

GSEA_Human_Disease<-Gene_Enrichment_Analysis(H_DMPs_Disease,"Human")

################################################################################
################## Supplementary Figure 2C #####################################
################################################################################

# Bind all the DMPs in a single file
M_DMPs_Disease<-rbind(M_Leu_DMPs,M_Lym_DMPs,M_MM_DMPs)

# Summary Pie Plot

## Pie plot in relation to CpG island
Mouse_CGI<-Pie_plot_function("CpG_Location_CGI",M_DMPs_Disease,"CGI")
## Pie plot in relation to Gene location
Mouse_GENE<-Pie_plot_function("UCSC_Group",M_DMPs_Disease,"Gene")

# EnrichR Analysis
GSEA_Mouse_Disease<-Gene_Enrichment_Analysis(M_DMPs_Disease,"Mouse")

################################################################################
################## Supplementary Figure 3 ######################################
################################################################################
# Comparison between normal isolated cells and cancer cell lines

## Human
### Human isolated cell lines were obtained from GSE110554 and GSE167998.

# Load normal cell lines sample sheet
targets_H_Normal<-readxl::read_xlsx("./HUMAN_Normal_Isolated_Cell_Lines.xlsx")
colnames(targets_H_Normal)[2]<-"Cell_Name"
targets_H_Normal$Sentrix_Position<-gsub(".*\\_","",targets_H_Normal$Basename)

IDATs_dir<-"./"

Human_Normal_Cells<-Preprocessing_function(IDATs_dir,"Homo_Sapiens",
                                           targets_H_Normal,ann850k)

bVals_Human_Normal<-Human_Normal_Cells$bVals

bVals_Human_Normal<-bVals_Human_Normal[,match(targets_H_Normal$Cell_Name, 
                                              colnames(bVals_Human_Normal))]

# Merge both datasets

df<-targets_Human[targets_Human$Cell_Type %in% c("B_cell","Myeloid_lineage",
                                                 "T_cell"),]
df<-df[df$Final_Classification != "Transformed_cell_line",]

bVals_df<-bVals_Human[,colnames(bVals_Human) %in% df$Cell_Name]
bVals_df<-bVals_df[,match(df$Cell_Name,colnames(bVals_df))]

df<-df[c(1,3,4,18)]
colnames(df)[c(3)]<-c("Cell_Subtype")
df$Type<-"Cancer_Cells"

df1<-targets_H_Normal[c(2:4,7)]
colnames(df1)[1:3]<-c("Cell_Name", "Cell_Type", "Cell_Subtype")
df1$Type<-"Normal_Cells"


targets_H_CN<-rbind(df,df1)

bVals_df<-as.data.frame(bVals_df)
bVals_df$CpG_ID<-rownames(bVals_df)

bVals_df2<-as.data.frame(bVals_Human_Normal)
bVals_df2$CpG_ID<-rownames(bVals_df2)

bVals_Human_CN<-merge(bVals_df,bVals_df2,"CpG_ID")
rownames(bVals_Human_CN)<-bVals_Human_CN$CpG_ID
bVals_Human_CN<-bVals_Human_CN[-1]

bVals_Human_CN<-bVals_Human_CN[,match(targets_H_CN$Cell_Name,
                                      colnames(bVals_Human_CN))]

## Unsupervised HeatMap of all isolated normal cells and cancer cells

ColAnn<-HeatmapAnnotation(Cell_Type = targets_H_CN[["Cell_Type"]],
                          Type = targets_H_CN[["Type"]],
                          col = list("Cell_Type" = Colours_Pal,
                                     "Type" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

#Human_Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_Human_CN,
#                                         nProbes<-0.01,
#                                         targetsdf<-targets_Human_CN,
#                                         BoottomAnn<-ColAnn,
#                                         Col_Split<-2)

## Unsupervised & Supervised analysis between cancer cell lines and 
## isolated normal cells

Human_DMPs_List<-list()

for (i in 1:length(unique(targets_H_CN$Cell_Type))) {
  
  List<-list()
  
  Cell_Type<-unique(targets_H_CN$Cell_Type)[i]
  print(Cell_Type)
  
  targets_df<-targets_H_CN[targets_H_CN$Cell_Type == Cell_Type,]
  bVals_df<-bVals_Human_CN[,colnames(bVals_Human_CN) %in% targets_df$Cell_Name]
  mVals_df<-B2M(bVals_df)
  
  bVals_df<-bVals_df[,match(targets_df$Cell_Name,colnames(bVals_df))]
  
  ColAnn<-HeatmapAnnotation(Type = targets_df[["Type"]],
                            col = list("Type" = Colours_Pal),
                            border = TRUE,
                            simple_anno_size = unit(0.5, "cm"))
  
  Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_df,
                                     nProbes<-0.01,
                                     targetsdf<-targets_df,
                                     BoottomAnn<-ColAnn,
                                     Col_Split<-2)
  
  List[["Unuspervised_HeatMap"]]<-Unsupervised
  
  DMPs<-Supervised_Cancer_Normal_DMPs_Analysis(targets_df = targets_df,
                                               bVals_df = bVals_df,
                                               mVals_df = mVals_df,
                                               Organism = "Human")
  
  List[["DMPs"]]<-DMPs
  
  bVals_df<-bVals_df[rownames(bVals_df) %in% DMPs$CpG_ID,]
  
  
  Supervised<-Unsupervised_HeatMap(bdf<-bVals_df,
                                   nProbes<-0.1,
                                   targetsdf<-targets_df,
                                   BoottomAnn<-ColAnn,
                                   Col_Split<-2)
  
  List[["Supervised_HeatMap"]]<-Supervised
  
  Human_DMPs_List[[Cell_Type]]<-List
  
}

# Mouse 
### Human isolated cell lines were obtained from GSE184410.
# Load normal cell lines sample sheet
targets_M_Normal<-readxl::read_xlsx("./MOUSE_Isolated_Cell_Lines.xlsx")
colnames(targets_M_Normal)[1]<-"Cell_Name"

IDATs_dir<-"./"

Mouse_Normal_Cells<-Preprocessing_function(IDATs_dir,"Mus_Musculus",
                                           targets_M_Normal,annMous)

bVals_Mouse_Normal<-Mouse_Normal_Cells$bVals

bVals_Mouse_Normal<-bVals_Mouse_Normal[,match(targets_M_Normal$Cell_Name, 
                                              colnames(bVals_Mouse_Normal))]

# Merge both datasets

bVals_df<-as.data.frame(bVals_Mouse)
bVals_df$CpG_Name<-rownames(bVals_df)

bVals_df2<-as.data.frame(bVals_Mouse_Normal)
bVals_df2$CpG_Name<-rownames(bVals_df2)

bVals_Mouse_CN<-merge(bVals_df,bVals_df2,"CpG_Name")

rownames(bVals_Mouse_CN)<-bVals_Mouse_CN$CpG_Name
bVals_Mouse_CN<-bVals_Mouse_CN[-1]

df<-targets_Mouse[c(1,4,13,6)]
colnames(df)[3]<-"Sentrix_Position"
df<-df[df$Disease != "Transformed_cell_line",][-ncol(df)]

df$Cell_Type<-ifelse(df$Cell_Type == "Myeloid", "Myeloid_Lineage",
              ifelse(df$Cell_Type == "Hematopoietic_Cell_Progenitor", "Myeloid_Lineage",
              ifelse(df$Cell_Type == "Myeloid_lineage", "Myeloid_Lineage",       
              ifelse(df$Cell_Type == "Myeloblast", "Myeloid_Lineage",df$Cell_Type))))

df$Type<-"Cancer_Cells"

targets_Mouse_CN<-df

df<-targets_M_Normal[c(1:3)]
colnames(df)[2]<-"Cell_Subtype"
df$Cell_Type<-ifelse(df$Cell_Subtype %like% "B_Cell","B_cell",
              ifelse(df$Cell_Subtype %like% "T_Cell","T_cell","Myeloid_Lineage"))

df$Type<-"Normal_Cells"
df$Sentrix_Position<-gsub(".*\\_","",df$EPIC_Sentrix_ID)

targets_Mouse_CN<-rbind(targets_Mouse_CN,df[c(1,4,6,5)])

bVals_Mouse_CN<-bVals_Mouse_CN[,colnames(bVals_Mouse_CN) %in% targets_Mouse_CN$Cell_Name,]

bVals_Mouse_CN<-bVals_Mouse_CN[,match(targets_Mouse_CN$Cell_Name,colnames(bVals_Mouse_CN))]


## Unsupervised HeatMap of all isolated normal cells and cancer cells

ColAnn<-HeatmapAnnotation(Cell_Type = targets_Mouse_CN[["Cell_Type"]],
                          Type = targets_Mouse_CN[["Type"]],
                          col = list("Cell_Type" = Colours_Pal,
                                     "Type" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

Mouse_Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_Mouse_CN,
                                         nProbes<-0.01,
                                         targetsdf<-targets_Mouse_CN,
                                         BoottomAnn<-ColAnn,
                                         Col_Split<-2)

## Unsupervised & Supervised analysis between cancer cell lines and 
## isolated normal cells
Mouse_DMPs_List<-list()

for (i in 1:length(unique(targets_Mouse_CN$Cell_Type))) {
  
  List<-list()
  
  Cell_Type<-unique(targets_Mouse_CN$Cell_Type)[i]
  print(Cell_Type)
  
  targets_df<-targets_Mouse_CN[targets_Mouse_CN$Cell_Type == Cell_Type,]
  bVals_df<-bVals_Mouse_CN[,colnames(bVals_Mouse_CN) %in% targets_df$Cell_Name]
  mVals_df<-B2M(bVals_df)
  
  bVals_df<-bVals_df[,match(targets_df$Cell_Name,colnames(bVals_df))]
  
  ColAnn<-HeatmapAnnotation(Type = targets_df[["Type"]],
                            col = list("Type" = Colours_Pal),
                            border = TRUE,
                            simple_anno_size = unit(0.5, "cm"))
  
  Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_df,
                                     nProbes<-0.01,
                                     targetsdf<-targets_df,
                                     BoottomAnn<-ColAnn,
                                     Col_Split<-2)
  
  List[["Unuspervised_HeatMap"]]<-Unsupervised
  
  DMPs<-Supervised_Cancer_Normal_DMPs_Analysis(targets_df = targets_df,
                                               bVals_df = bVals_df,
                                               mVals_df = mVals_df,
                                               Organism = "Mouse")
  
  List[["DMPs"]]<-DMPs
  
  bVals_df<-bVals_df[rownames(bVals_df) %in% DMPs$CpG_ID,]
  
  
  Supervised<-Unsupervised_HeatMap(bdf<-bVals_df,
                                   nProbes<-0.1,
                                   targetsdf<-targets_df,
                                   BoottomAnn<-ColAnn,
                                   Col_Split<-2)
  
  List[["Supervised_HeatMap"]]<-Supervised
  
  Mouse_DMPs_List[[Cell_Type]]<-List
  
}


################################################################################
########################### Figure 2A ##########################################
################################################################################

# Prepare targets file
# Remove Transformed cell lines, NH_B_Lymphoma_Others and Lymphoblastic_Leukemia_Others
targets_Class<-targets_Human[!targets_Human$Final_Classification %in% 
                              c("Transformed_cell_line",
                                "NH_B_lymphoma_Others",
                                "Lymphoblastic_Leukemia_Others"),]

bVals_Class<-bVals_Human[,colnames(bVals_Human) %in% targets_Class$Cell_Name]
bVals_Class<-bVals_Class[,match(targets_Class$Cell_Name,colnames(bVals_Class))]
mVals_Class<-B2M(bVals_Class)

# Differential methylation analysis

Classifier_DMPs_df<-data.frame()
Classifier_DMPs_list<-list()

for (i in 1:length(unique(targets_Class$Final_Classification))) {
  
  Ref_Disease<-unique(targets_Class$Final_Classification)[i]
  print(paste0("Calculating data from ",Ref_Disease))
  
  sink(file="./Summary_DMPs_Classifier.txt",append=TRUE)
  cat("\n",Ref_Disease,"\n")
  sink()
  
  targets_df<-targets_Class
  targets_df$Comparison<-ifelse(targets_df$Final_Classification == Ref_Disease,"RefDisease","Others")
  
  print(paste0("Calculating stats DMPs from ",Ref_Disease))
  
  limma_res<-limma_function(targets_df)
  
  if (class(limma_res) != "list"){
    
    sink(file="./Summary_DMPs_Classifier.txt",append=TRUE)
    cat("\n There is no statistical DMP\n")
    sink()
    
    next
    
  }
  
  print(paste0("Calculating bio DMPs from ",Ref_Disease))
  
  DMPs<-Diff_Meth_Analysis_Classifier(targets_df,limma_res$Diff_df,ann850k,0.60,limma_res$fit_2)
  
  if (class(DMPs) != "data.frame"){
    
    sink(file="./Summary_DMPs_Classifier.txt",append=TRUE)
    cat("\n There is no biological DMP\n")
    sink()
    
    next
    
  }
  
  sink(file="./Summary_DMPs_Classifier.txt",append=TRUE)
  cat("\n Number of BIO DMPs\n")
  print(table(DMPs$Meth_Status))
  sink()
  
  DMPs$Disease<-Ref_Disease
  
  writexl::write_xlsx(DMPs,paste0("./DMPs_",Ref_Disease,".xlsx"))
  
  Classifier_DMPs_df<-rbind(Classifier_DMPs_df,DMPs)
  Classifier_DMPs_list[[Ref_Disease]]<-DMPs
  
  print(paste0(Ref_Disease," finished"))
  
}

DMPs_Classifier<-Classifier_DMPs_df

bVals_Class<-bVals_Class[rownames(bVals_Class) %in% DMPs_Classifier$CpG_ID,]

bVals_Class<-bVals_Class[,match(targets_Class$Cell_Name,colnames(bVals_Class))]

ColAnn<-HeatmapAnnotation(Disease = targets_Class[["Final_Classification"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))


HeatMap_Classificator<-Draw_HeatMap(bVals_Cluster = bVals_Class,
                                    bVals_HeatMap = bVals_Class, 
                                    bottom_annotation = ColAnn,
                                    ColSplit = 8)

################################################################################
##################### Supplementary Figure 4 ##################################
################################################################################

# Summary Pie Plot
## Pie plot disease representation

Class_PP_Disease<-Pie_plot_function("Disease",DMPs_Classifier,"Disease")
## Pie plot in relation to CpG island
Class_PP_CGI<-Pie_plot_function("Relation_to_Island",DMPs_Classifier,"CGI")
## Pie plot in relation to Gene location
Class_PP_Gene<-Pie_plot_function("UCSC_RefGene_Group",DMPs_Classifier,"Gene")

#t-SNE analysis
targets_Class<-targets_Class[match(colnames(bVals_Class),targets_Class$Cell_Name),]

set.seed(214)
tSNE_Class <- Rtsne(t(bVals_Class), perplexity=30, check_duplicates = FALSE)

tSNE_Class<-tSNE_Process_function(TSNEObj<-tSNE_Class, 
                                  targets_df<-targets_Class, 
                                  Column_Col<-"Final_Classification",
                                  Organism = "Human")
# EnrichR Analysis
GSEA_Human_Disease_Class<-Gene_Enrichment_Analysis(DMPs_Classifier,"Human")

################################################################################
###########################  Figure 2B #########################################
################################################################################

# Validation of the classifier in Primary samples
## AML samples are from GSE153349  and DLBCL samples are from GSE255869, 
##T-ALL and B-ALL samples are from in-house methylation data

## Load sample sheet
targets_Primary<-readxl::read_xlsx("./Primary_Samples_850K.xlsx")
colnames(targets_Primary)[1]<-"Cell_Name"

# Preprocess IDATs
IDATs_dir<-"./"

Primary_Samples<-Preprocessing_function(IDATs_dir,"Homo_Sapiens",
                                        targets_Primary,ann850k)

bVals_Primary<-Primary_Samples$bVals

bVals_Primary<-bVals_Primary[,match(targets_Primary$Cell_Name,
                                    colnames(bVals_Primary))]



# Validate the classifier in cell lines in tumor samples
bVals_Primary_Class<-bVals_Primary[rownames(bVals_Primary) %in% 
                                     DMPs_Classifier$CpG_ID,]

bVals_Primary_Class<-bVals_Primary_Class[,match(targets_Primary$Cell_Name,
                                                colnames(bVals_Primary_Class))]

ColAnn<-HeatmapAnnotation(Disease = targets_Primary[["Final_Classification"]],
                          col = list("Disease" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))


HeatMap_Classificator_Primary<-Draw_HeatMap(bVals_Cluster = bVals_Primary_Class,
                                            bVals_HeatMap = bVals_Primary_Class, 
                                            bottom_annotation = ColAnn,
                                            ColSplit = 3)

################################################################################
##################### Supplementary  Figure 5 ##################################
################################################################################

# Join cancer cell lines and primary samples datasets
## bVals
temp<-as.data.frame(bVals_Human)
temp$CpG_ID<-rownames(temp)

temp2<-as.data.frame(bVals_Primary)
temp2$CpG_ID<-rownames(temp2)

bVals_CP<-merge(temp, temp2, "CpG_ID")
rownames(bVals_CP) <- bVals_CP$CpG_ID
bVals_CP<-bVals_CP[-1]

## Sample sheets
targets_CP<-targets_Human[targets_Human$Final_Classification %in% 
                          targets_Primary$Final_Classification,][c(1,24,18)]

targets_CP$Type<-"Cell_Line"

temp<-targets_Primary[c(1,4,5)]
colnames(temp)[c(3)]<-c("Sentrix_Position")
temp$Sentrix_Position<-gsub(".*\\_","",temp$Sentrix_Position)
temp$Type<-"Primary_Sample"

targets_CP<-rbind(targets_CP,temp)

bVals_CP<-bVals_CP[,colnames(bVals_CP) %in% targets_CP$Cell_Name,]
bVals_CP<-bVals_CP[,match(targets_CP$Cell_Name,colnames(bVals_CP))]

nrow(targets_CP)==ncol(bVals_CP)

## Unsupervised HeatMap
## Design top annotation
bVals_CP<-bVals_CP[,match(targets_CP$Cell_Name,colnames(bVals_CP))]

ColAnn<-HeatmapAnnotation(Disease = targets_CP[["Disease"]],
                          Type = targets_CP[["Type"]],
                          col = list("Disease" = Colours_Pal,
                                     "Type" = Colours_Pal),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

Human_Unsupervised<-Unsupervised_HeatMap(bdf<-bVals_CP,
                                         nProbes<-0.01,
                                         targetsdf<-targets_CP,
                                         BoottomAnn<-ColAnn,
                                         Col_Split<-2)

# Disease speciffic-analysis

Result_List<-list()

for (i in 1:length(unique(targets_CP$Final_Classification))) {
  
  Disease<-unique(targets_CP$Final_Classification)[i]
  
  print(Disease)
  
  targets_df<-targets_CP[targets_CP$Final_Classification == Disease,]
  bVals_df<-bVals_CP[,colnames(bVals_CP) %in% targets_df$Cell_Name]
  bVals_df<-bVals_df[,match(targets_df$Cell_Name,colnames(bVals_df))]
  mVals_df<-B2M(bVals_df)
  
  #HeatMap_Un<-Unsupervised_HeatMap(bVals_df,0.01,targets_df)
  
  Result<-limma_function_Cell_Primaries(targets_df, mVals_df)
  Result<-Diff_Meth_Analysis_Cells_Primaries(targets_df,bVals_df,Result$Diff_df,as.data.frame(ann850k),0.4,Result$fit_2,mVals_df)
  
  ColAnn<-HeatmapAnnotation(Type = targets_df[["Type"]],
                            col = list("Type" = Colours_Pal),
                            border = TRUE,
                            simple_anno_size = unit(0.50, "cm"))
  
  HeatMap<-Draw_HeatMap(bVals_Cluster=Result$bVals_Filt,
                        bVals_HeatMap=Result$bVals_Filt, 
                        bottom_annotation = ColAnn,
                        ColSplit=2)
  
  Result[["HeatMap_Un"]]<-HeatMap_Un
  Result[["HeatMap"]]<-HeatMap
  
  Result_List[[Disease]]<-Result
  
}

# Supervised HeatMap of all diseases

DMPs_DLBCL<-Result_List$DLBCL$DMPs
DMPs_TALL<-Result_List$T_ALL$DMPs
DMPs_BALL<-Result_List$B_ALL$DMPs
DMPs_AML<-Result_List$AML$DMPs

bVals_CP_Sup<-bVals_CP[rownames(bVals_CP) %in% c(DMPs_DLBCL$CpG_Name,
                                                 DMPs_TALL$CpG_Name,
                                                 DMPs_BALL$CpG_Name,
                                                 DMPs_AML$CpG_Name),]

bVals_CP_Sup<-bVals_CP_Sup[,match(targets_CP$Cell_Name,colnames(bVals_CP_Sup))]

HeatMap_Supervised_Cell_Primary<-Draw_HeatMap(bVals_Cluster=bVals_CP_Sup,
                                              bVals_HeatMap=bVals_CP_Sup, 
                                              bottom_annotation = ColAnn,
                                              ColSplit=3)

################################################################################
############################# Figure 2G ########################################
################################################################################

# Study of methylation - expression - Drug IC50 associations
## Data preparation
### Drugs approved to treat hematological malignancies
Drugs_Approved<-toupper(c("Vinblastine","Cytarabine","Methotrexate",
                          "Vorinostat","Nilotinib","Bosutinib","Lenalidomide",
                          "Dasatinib","Crizotinib","Bortezomib",
                          "Cyclophosphamide","Ibrutinib","Mitoxantrone",
                          "Bleomycin","Fludarabine","Nelarabine","Dacarbazine",
                          "Romidepsin","5azacytidine"))


# GDSC2
GDSC2<-readxl::read_xlsx("./GDSC2_fitted_dose_response_24Jul22.xlsx")
GDSC2$DRUG_NAME_MOD<-toupper(gsub('[^[:alnum:] ]', '', GDSC2$DRUG_NAME))
GDSC2$CELL_LINE_NAME_MOD<-toupper(gsub('[^[:alnum:] ]', '', GDSC2$CELL_LINE_NAME))

# Select the drugs of the study
GDSC2_Drug<-GDSC2[GDSC2$DRUG_NAME_MOD %in% Drugs_Approved,]
# Remove a cell line that in its modified form is the same
GDSC2<-GDSC2[GDSC2$CELL_LINE_NAME != "KMH-2",]
GDSC2<-GDSC2[GDSC2$CELL_LINE_NAME_MOD %in% targets_Human$CELL_LINE_NAME_MOD,]

targets_Drug<-targets_Human[targets_Human$CELL_LINE_NAME_MOD %in% GDSC2_Drug$CELL_LINE_NAME_MOD,]

# RNA DATA
RNA<-as.data.frame(data.table::fread("./Cell_line_RMA_proc_basalExp.txt"))
RNA<-RNA[-2]
colnames(RNA)<-gsub("DATA.","",colnames(RNA))
RNA<-RNA[RNA$GENE_SYMBOLS != "",]
rownames(RNA)<-RNA$GENE_SYMBOLS
RNA<-RNA[-1]

targets_RNA<-readxl::read_xlsx("./Clinical_Data_CELL_LINES_EXPRESSION.xlsx")
targets_RNA$CELL_LINE_NAME_MOD<-toupper(gsub('[^[:alnum:] ]', '', targets_RNA$Sample_Name))

targets_RNA_METH<-targets_Human[targets_Human$CELL_LINE_NAME_MOD %in% targets_RNA$CELL_LINE_NAME_MOD,]

RNA<-RNA[,colnames(RNA) %in% targets_RNA$COSMIC_identifier]
RNA<-RNA[,match(targets_RNA$COSMIC_identifier, colnames(RNA))]
colnames(RNA)<-targets_RNA$Sample_Name

# Analysis

for (i in 1:length(unique(GDSC2_Drug$DRUG_NAME_MOD))) {
  Drug<-unique(GDSC2_Drug$DRUG_NAME_MOD)[i]
  df<-GDSC2[GDSC2$DRUG_NAME_MOD == Drug,]
  targetsdf<-merge(targets_Human[c(1,23)],df[c(21,16)],"CELL_LINE_NAME_MOD")
  colnames(targetsdf)[3]<-Drug
  mVals<-mVals_Human[,colnames(mVals_Human) %in% targetsdf$Cell_Name]
  mVals<-mVals[,match(targetsdf$Cell_Name,colnames(mVals))]
  bVals<-bVals_Human[,colnames(bVals_Human) %in% targetsdf$Cell_Name]
  bVals<-bVals[,match(targetsdf$Cell_Name,colnames(bVals))]
  DMPs<-Diff_Meth_Expression_function(Drug,mVals,targetsdf,as.data.frame(ann850k),"Methylation")
  
  sink(file="./Summary_Linear_Regression.txt",append=TRUE)
  cat(paste0("\n\n",Drug,"\n"))
  cat("\nNumber of DMPs (CpGs)\n")
  cat(length(unique(DMPs$CpG_ID)))
  sink()
  
  if (nrow(DMPs) == 0) {
    next
  }
  
  writexl::write_xlsx(DMPs, paste0(".DMPs_",Drug,".xlsx"))
  
  bVals<-bVals[rownames(bVals) %in% DMPs$CpG_ID,]
  
  targets_RNA_temp<-targets_RNA[targets_RNA$CELL_LINE_NAME_MOD %in% df$CELL_LINE_NAME_MOD,]
  targets_RNA_temp<-merge(targets_RNA_temp,df[c(21,16)],"CELL_LINE_NAME_MOD")
  colnames(targets_RNA_temp)[8]<-Drug
  RNAdf<-RNA[,colnames(RNA) %in% targets_RNA_temp$Sample_Name]
  RNAdf<-RNAdf[,match(targets_RNA_temp$Sample_Name,colnames(RNAdf))]
  DEG<-Diff_Meth_Expression_function(Drug,RNAdf,targets_RNA_temp,FALSE,"Expression")
  
  sink(file="./Summary_Linear_Regression.txt",append=TRUE)
  cat("\nNumber of DEG (Genes)\n")
  cat(length(unique(rownames(DEG))))
  sink()
  
  if (nrow(DEG) == 0) {
    next
  }
  
  DEG$Gene_ID<-rownames(DEG)
  DEG<-DEG[,c(ncol(DEG),1:(ncol(DEG)-1))]
  writexl::write_xlsx(DEG, paste0("./DEGs_",Drug,".xlsx"))
  
  Common_Genes<-Venn_Diagram_function(DMPs,DEG,Drug)
  
  sink(file="./Summary_Linear_Regression.txt",append=TRUE)
  cat("\nNumber of Common Genes\n")
  cat(length(unique(Common_Genes)))
  sink()
  
  
  if (length(Common_Genes) == 0) {
    next
  }
  
  Corr_table<-Correlation_Function(Common_Genes,DMPs,DEG,bVals,RNAdf)
  
  writexl::write_xlsx(Corr_table, paste0("./Correlation_DMPs_DEGs_",Drug,".xlsx"))
  
  sink(file="./Summary_Linear_Regression.txt",append=TRUE)
  cat("\n\nNumber of CpGs that correlate with gene expression (CpGs)\n")
  print(as.matrix(table(Corr_table$Sig)))
  sink()
  
  Corr_table_Sig<-Corr_table[Corr_table$Sig == TRUE,]
  
  if (nrow(Corr_table_Sig) == 0) {
    next
  }
  
  writexl::write_xlsx(Corr_table_Sig, paste0("./Correlation_DMPs_DEGs_",Drug,"_Sig.xlsx"))
  
  sink(file="./Summary_Linear_Regression.txt",append=TRUE)
  cat("\nNumber of unique CpGs that correlate with gene expression (CpGs):",length(unique(Corr_table_Sig$CpG)),"\n") 
  cat("\nNumber of unique Genes that correlate with gene expression (Genes):",length(unique(Corr_table_Sig$Gene)),"\n\n\n") 
  sink()
}


## Intersect the common CpG within nucleotide analogues

CYT<-readxl::read_xlsx("./Correlation_DMPs_DEGs_CYTARABINE_Sig.xlsx")
FLU<-readxl::read_xlsx("./Correlation_DMPs_DEGs_FLUDARABINE_Sig.xlsx")
NEL<-readxl::read_xlsx("./Correlation_DMPs_DEGs_NELARABINE_Sig.xlsx")

Common_CpGs<-intersect(CYT$CpG,intersect(FLU$CpG,NEL$CpG))
length(Common_CpGs)

bVals_Drug<-bVals_Human[,colnames(bVals_Human) %in% targets_Drug$Cell_Name]
bVals_Drug_AN<-bVals_Drug[rownames(bVals_Drug) %in% Common_CpGs,]

Drugs_Study<-GDSC2[GDSC2$DRUG_NAME %in% c("Cytarabine", "Fludarabine", "Nelarabine"),][c(9,16,21)]

df_res<-data.frame(CELL_LINE_NAME_MOD = unique(Drugs_Study$CELL_LINE_NAME_MOD))

for (i in 1:length(unique(Drugs_Study$DRUG_NAME))) {
  Drug<-unique(Drugs_Study$DRUG_NAME)[i]
  df<-Drugs_Study[Drugs_Study$DRUG_NAME == Drug,][2:3]
  colnames(df)[1]<-Drug
  df_res<-merge(df_res,df,"CELL_LINE_NAME_MOD", all.x=TRUE)
  
}

targets_Drug<-merge(targets_Drug,df_res,"CELL_LINE_NAME_MOD")

# Perform the heatmap to see if they cluster by drug sensitivity
targets_Drug<-targets_Drug[!targets_Drug$Final_Classification %in% c("Transformed_cell_line",
                                                                     "NH_B_lymphoma_Others",
                                                                     "Lymphoblastic_Leukemia_Others"),]

bVals_Drug_AN<-bVals_Drug_AN[,colnames(bVals_Drug_AN) %in% targets_Drug$Cell_Name]

bVals_Drug_AN<-bVals_Drug_AN[,match(targets_Drug$Cell_Name,colnames(bVals_Drug_AN))]

palCYT <-colorRamp2::colorRamp2(seq(from=min(na.omit(targets_Drug[["Cytarabine"]])),to=max(na.omit(targets_Drug[["Cytarabine"]])),length.out=10),
                                rev(viridis::inferno(10)))

palFLU <-colorRamp2::colorRamp2(seq(from=min(na.omit(targets_Drug[["Fludarabine"]])),to=max(na.omit(targets_Drug[["Fludarabine"]])),length.out=10),
                                rev(viridis::viridis(10)))

palNEL <-colorRamp2::colorRamp2(seq(from=min(na.omit(targets_Drug[["Nelarabine"]])),to=max(na.omit(targets_Drug[["Nelarabine"]])),length.out=10),
                                rev(viridis::mako(10)))

ColAnn<-HeatmapAnnotation(Final_Classification =as.matrix(targets_Drug[["Final_Classification"]]),
                          Cyterabine = (targets_Drug[["Cytarabine"]]),
                          Fludarabine = (targets_Drug[["Fludarabine"]]),
                          Nelarabine = (targets_Drug[["Nelarabine"]]),                           
                          col = list("Final_Classification" = Colours_Pal,
                                     Cyterabine = palCYT,
                                     Fludarabine = palFLU,
                                     Nelarabine = palNEL),
                                     border = TRUE,
                                     simple_anno_size = unit(0.75, "cm"))

HeatMap_Drug<-Draw_HeatMap(bVals_Drug_AN,bVals_Drug_AN, ColAnn,2)


temp<-HeatMap_Drug$Cluster_Columns[c(1,3)]
colnames(temp)<-c("Cell_Name","Cluster_General")
temp$Cluster_General<-ifelse(temp$Cluster_General == "Cluster_1", "Sensitive_Cluster","Resistant_Cluster")
targets_Drug<-merge(targets_Drug,temp,"Cell_Name")

targets_Drug<-targets_Drug[match(colnames(bVals_Drug_AN),targets_Drug$Cell_Name),]

ColAnn_TOP<-HeatmapAnnotation(Cluster =as.matrix(targets_Drug[["Cluster_General"]]),                 
                              col = list("Cluster" = c("Sensitive_Cluster" = "#CAF2C2FF",
                                                       "Resistant_Cluster" = "#1E1E32FF")),
                               border = TRUE,
                               simple_anno_size = unit(0.75, "cm"))

distance <- dist(t(bVals_Drug_AN), method="euclidean")
clust.model <- hclust(distance, method="ward.D2")

distance_r <- dist(bVals_Drug_AN, method="euclidean")
clust.model_r <- hclust(distance_r, method="ward.D2")

HeatMap_Drug<-draw(Heatmap(as.matrix(bVals_Drug_AN), name = "Methylation", 
                           col = gplots::greenred(75), 
                           bottom_annotation = ColAnn,
                           top_annotation = ColAnn_TOP,  
                           cluster_columns = clust.model,
                           cluster_rows = clust.model_r,
                           column_split= 2,
                           show_row_names = FALSE, 
                           show_column_names = FALSE, 
                           show_row_dend = FALSE))

# Barplots

# IC50 between resistant and sensitive clusters
Drugs<-c("Cytarabine","Nelarabine","Fludarabine")

res_df<-data.frame()
res_stat<-data.frame()

for (i in 1:length(Drugs)){
  Drug<-Drugs[i]
  df<-targets_Drug[,colnames(targets_Drug) %in% c("Cell_Name","Cluster_General",Drug)]
  df$Drug<-Drug
  colnames(df)[2]<-"IC50"
  df<-df[complete.cases(df$IC50),]
  res_df<-rbind(res_df,df)
  test<-wilcox.test(IC50 ~ Cluster_General, data = df)
  res_stat<-rbind(res_stat, data.frame(Drug = Drug, Stat = test$p.value))
  df
}

res_df$Cluster_General<-factor(res_df$Cluster_General, levels = c("Sensitive_Cluster","Resistant_Cluster"))

ggplot(res_df, aes(x=Cluster_General, y=IC50, fill=Cluster_General)) +
  scale_fill_manual(values=c("Resistant_Cluster" = "#1E1E32FF",
                             "Sensitive_Cluster" = "#CAF2C2FF")) +
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(shape=21, position=position_jitter(0.2), size = 3, color = "black", fill = "gray80", alpha = 0.75)+
  facet_wrap(~Drug,nrow = 3,scales = "fixed")+
  ylim(-10,10)+
  theme_few()

################################################################################
#################### Supplementary Figure 2G ###################################
################################################################################

# To obtain the CpGs associated to AZA sensitivity you have to repeat the 
# previous analysis changing the threshold into the Diff_Meth_Expression_function 
# function, instead of 0.01, it has to be changed to 0.5

# And repeat the previous graphs:

DMPs_AZA<-readxl::read_xlsx("./AZA_ANALYSIS_0.05.xlsx")
DMPs_AZA


bVals_Aza<-as.data.frame(bVals_Human)
bVals_Aza<-as.data.frame(t(bVals_Aza[rownames(bVals_Aza) %in% DMPs_AZA$CpG,]))
bVals_Aza$Cell_Name<-rownames(bVals_Aza)

Drug_Aza<-GDSC2[GDSC2$DRUG_NAME_MOD == "5AZACYTIDINE",][c(16,21)]
targets_AZA<-merge(targets_Human,Drug_Aza,"CELL_LINE_NAME_MOD")

bVals_Aza<-bVals_Aza[rownames(bVals_Aza) %in% targets_AZA$Cell_Name,]
bVals_Aza<-bVals_Aza[,match(targets_AZA$Cell_Name,colnames(bVals_Aza))]

palAZA <-colorRamp2::colorRamp2(seq(from=min(na.omit(targets_AZA[["LN_IC50"]])),to=max(na.omit(targets_AZA[["LN_IC50"]])),length.out=10),
                                rev(viridis::magma(10)))

ColAnn<-HeatmapAnnotation(Final_Classification =as.matrix(targets_AZA[["Final_Classification"]]),
                          Azacytidine = (targets_AZA[["LN_IC50"]]),
                          col = list("Final_Classification" = Colours_Pal,
                          Azacytidine = palAZA),
                          border = TRUE,
                          simple_anno_size = unit(0.75, "cm"))

HeatMap_AZA<-Draw_HeatMap_No_Rows_Clust(t(bVals_Aza[1]),t(bVals_Aza[1]), ColAnn, 2)

temp<-HeatMap_AZA$Cluster_Columns[c(1,3)]
colnames(temp)<-c("Cell_Name","Cluster_General")
temp$Cluster_General<-ifelse(temp$Cluster_General == "Cluster_1", "Sensitive_Cluster","Resistant_Cluster")
targets_AZA<-merge(targets_AZA,temp,"Cell_Name")

# Barplots

# IC50 between resistant and sensitive clusters
targets_AZA$Cluster_General<-factor(targets_AZA$Cluster_General, levels = c("Sensitive_Cluster","Resistant_Cluster"))
AZA_IC50_PLOT<-ggplot(targets_AZA, aes(x=Cluster_General, y=LN_IC50, fill=Cluster_General)) +
                      scale_fill_manual(values=c("Resistant_Cluster" = "#1E1E32FF",
                                                 "Sensitive_Cluster" = "#CAF2C2FF")) +
                      geom_boxplot(outlier.alpha = 0)+
                      geom_jitter(shape=21, position=position_jitter(0.2), size = 3, color = "black", fill = "gray80", alpha = 0.75)+
                      ylim(-10,10)+
                      theme_few()
wilcox.test(LN_IC50~Cluster_General, data = targets_AZA)

# Expression Plot
RNA_Aza<-as.data.frame(t(RNA[rownames(RNA) %in% DMPs_AZA$Gene,]))
RNA_Aza$CELL_LINE_NAME_MOD<-toupper(gsub("[^[:alnum:] ]","",rownames(RNA_Aza)))
targets_AZA<-merge(targets_AZA,RNA_Aza,"CELL_LINE_NAME_MOD", all.x=TRUE)

AZA_EXP_PLOT<-ggplot(targets_AZA, aes(x=Cluster_General, y=TNFAIP3, fill=Cluster_General)) +
                     scale_fill_manual(values=c("Resistant_Cluster" = "#1E1E32FF",
                                                "Sensitive_Cluster" = "#CAF2C2FF")) +
                     geom_boxplot(outlier.alpha = 0)+
                     geom_jitter(shape=21, position=position_jitter(0.2), size = 3, color = "black", fill = "gray80", alpha = 0.75)+
                     ylim(0,15)+
                     theme_few()

wilcox.test(TNFAIP3~Cluster_General, data = targets_AZA)
