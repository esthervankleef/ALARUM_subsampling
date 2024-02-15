
rm(list=ls())


library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(scales)
library(stringi)
library(MASS)
library(stargazer)
library(reporttools)
library(epitools)
library(gdata)
library(car)
library(plyr)
library(dplyr)
library(data.table)
library(tibble)
library(psych)
library(tidyr)
library(janitor)
library(psych)
library(plotrix)
#library(slopegraph)
library(Lock5Data)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(treemap)
library (treemapify)
library(ggraph)
library(igraph)



###The code, uses the following "INPUT" files:

# CARD_read_count_specific.tsv                # Can be found in input_files
# ResPipe_CARD-3.0.3.meta.tsv                 # Can be found in input_files
# CARD_lateral_coverage_specific.tsv          # Can be found in input_files
# CARD_read_lengths_specific.tsv              # Can be found in input_files
# Metagenomics_Metadata.csv                   # Could not be found, took from code original paper
# Infection_Data_For_Bayesian_Model.csv       # Could not be found, took from code original paper
# bracken_combined_reads.tsv                  # Could not be found, took from code original paper

# to produce the corrected resistance gene counts, plus a matrix that links each resistance gene and antibiotic based on the "Confers_Resistance_to_Antibiotic" relationship ontology term.
    # in the matrix, 1 == the gene is associated with clear experimental evidence of elevated MIC for that antibiotic but the "Confers_Resistance_to_Antibiotic" relationship ontology term is missing
    # in the matrix, 2 == the gene is associated with demonstrably elevated MIC for that antibiotic and is known to confer or contribute to clinically relevant resistance to that antibiotic ("Confers_Resistance_to_Antibiotic" relationship ontology term is present)

#OUTPUT FILES:

# Corrected_Counts.csv
# AB_Matrix_1_or_2.csv

# AMR_DEF.csv
# AMR_ALL.csv

# And to produce the final dataset for the bayesian modelling (i.e. "OUTPUT" file)

# Dataset_For_Bayesian_Model.csv

# The code for the non-metric multidimensional scaling (NMDS) ordination method is also presented (see supplementary results for the validation of the pooling using 30-sample pools)


##################PRODUCE THE DATASET WITH CORRECTED GENE COUNTS#############################


CARD_read_count_specific<-read.csv("./Input Files/CARD_read_count_specific.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
attach(CARD_read_count_specific)
fix(CARD_read_count_specific)


Gene_Lenghts<-read.csv("./Input Files/ResPipe_CARD-3.0.3.meta.tsv", header=T,check.names = F, row.names = 1, sep = "\t")
MT<-Gene_Lenghts
attach(MT)
GL<-subset(MT, select=("SeqLength"))
fix(GL)


Spec_Lat_Cov<-read.csv("./Input Files/CARD_lateral_coverage_specific.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
SLC<-Spec_Lat_Cov
attach(SLC)
fix(SLC)

Read_Lenghts<-read.csv("./Input Files/CARD_read_lengths_specific.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
ARL<-Read_Lenghts
attach(ARL)
fix(ARL)



memory.limit(size = 1000000000)
names(SLC) = gsub(pattern = "CARD_", replacement = "", x = names(SLC))
names(ARL) = gsub(pattern = "CARD_", replacement = "", x = names(ARL))
names(CARD_read_count_specific) = gsub(pattern = "CARD_", replacement = "", x = names(CARD_read_count_specific))

identical(rownames(GL), rownames(ARL))
identical(rownames(ARL), rownames(SLC))
identical(rownames(SLC), rownames(CARD_read_count_specific))
identical(colnames(SLC), colnames(ARL))
identical(colnames(SLC), colnames(CARD_read_count_specific))

GSL<-mapply("*", SLC,GL)
GLARL=GSL/ARL 
GLARL[is.na(GLARL)] <- 0
Corr_Count=CARD_read_count_specific/GLARL 
Corr_Count[is.na(Corr_Count)] <- 0


##########################MAP CORRECTED COUNTS TO AROs################################################################

newnamesd1<-names(Corr_Count)

metadata = read.csv(file = "./Input Files/Original_paper/Metagenomics_Metadata.csv", header = T, row.names = 1, check.names = F)
attach(metadata)

fix(metadata)
z = metadata[match(newnamesd1, rownames(metadata)),]


#REMOVE 24 FOLLOW-UP SAMPLES FROM PARALLEL STUDY


#data(Corr_Count) # Seems to not work (Evk 15 Feb 2024)
dropList <- c("NCS-ST-0047",
              "NCS-ST-0049",
              "NCS-ST-0051",
              "NCS-ST-0065",
              "NCS-ST-0073",
              "NCS-ST-0079",
              "NCS-ST-0110",
              "NCS-ST-0175",
              "NCS-ST-0176",
              "NCS-ST-0191",
              "NCS-ST-0247",
              "NCS-ST-0251",
              "NCS-ST-0257",
              "NCS-ST-0380",
              "NCS-ST-0389",
              "NCS-ST-0413",
              "NCS-ST-0416",
              "NCS-ST-0417",
              "NCS-ST-0423",
              "NCS-ST-0426",
              "NCS-ST-0429",
              "NCS-ST-0438",
              "NCS-ST-0454",
              "NCS-ST-0457")

Int_1 <- Corr_Count[, !colnames(Corr_Count) %in% dropList]
Corr_Count<-Int_1

#####merge gene data with ARO labels#######################


df<-MT
df[,2:12]<- list(NULL)
labels_RG<-df

All_Data <- merge(Corr_Count, labels_RG, by=0, all=FALSE)


#remove Row names

All_Data$Row.names<-NULL

## move ARO accession number column to first position

Int_1 <- All_Data %>%
  select(ARO_accession, everything()) 

####a few aro accession numbers map up to more than one R_XXX UID. Hence we need a dataset that
#### sums counts across all rows that map to the same aro number


Int_2<-aggregate(.~ARO_accession, data=Int_1, FUN=sum) 

##convert aroaccession to row names

row.names(Int_2) <- Int_2[,1]
Int_2["ARO_accession"]<-NULL 


##remove target genes where the sum across ALL SAMPLES AND POOLS is zero. 


total_col<-apply(Int_2[,], 1, sum)
Int_3<-as.data.frame(cbind(Int_2,total_col))
Int_4<-Int_3[Int_3$total_col!=0, ]

nrow(Int_4) ##Total number of ARO accession IDs detected in at least one of the samples

ARO_Counts<-Int_4
ARO_Counts$total_col<-NULL

#ARCHIVE ARO DATA

write.csv(ARO_Counts, "./Output Files/Corrected Gene Counts.csv") 




##############################MAP AROS TO ANTIBIOTICS AND CALCULATE TOTAL CORRECTED COUNTS FOR EACH ANTIBIOTIC##################################### 

##FIRST CREATE A FILE MAPPING WHAT ANTIBIOTICS CONTAIN THE "CONFERS RESISTANCE TO ANTIBIOTIC" ONTOLOGY TERM AND WHICH DON'T

da<-MT
da$Antibiotic <- gsub('_', '', da$Antibiotic)
da$Antibiotic <- gsub('-', '', da$Antibiotic)
Ab_b<-subset(da, select=c("ARO_accession", "Antibiotic"))
Ab_c <- unique(Ab_b[ , ])


Ab_Drug <- Ab_c %>%
  separate_rows(Antibiotic) %>%
  mutate(Value = case_when(grepl("confersresistancetodrug",Antibiotic)~ 2, !grepl("confersresistancetodrug",Antibiotic)~1)) %>%
  spread(Antibiotic, Value, fill = 0) %>%
  rename_at(vars(-ARO_accession), funs(paste0("AB_", .)))
Ab_Drug  

Ab_Drug$AB_Unassigned<-NULL

Ab_Drug<-as.data.frame(Ab_Drug)
Ab_Drug = setNames(data.frame(t(Ab_Drug[,-1])), Ab_Drug[,1])


Ab_Drug<-rownames_to_column(Ab_Drug, var = "ABS")
Ab_Drug$ABS <- gsub('confersresistancetodrug', '', Ab_Drug$ABS)


Gr <- Ab_Drug %>% group_by(ABS) %>% summarise_all(sum)

Gr<-as.data.frame(Gr)
Gr2 = setNames(data.frame(t(Gr[,-1])), Gr[,1])

####map aro counts to antibiotics

Antibiotics<-Gr2

Gr2<-rownames_to_column(Gr2, var = "ARO")


#ARCHIVE ABS CLASSIFICATION (1 OR 2)

write.csv(Gr2, "./Output Files/AB_Matrix_1_or_2.csv") 



##Group antibiotics that are the same/similar:

names(Antibiotics) = gsub(pattern = "AB_", replacement = "", x = names(Antibiotics))


Antibiotics$Bacitracin_All<-apply(Antibiotics[,c("bacitracinA","bacitracinB","bacitracinF")], 1, max)
Antibiotics$Bleomycin_All<-apply(Antibiotics[,c("bleomycinA2","bleomycinB2","bleomycinicacid")], 1, max)
Antibiotics$Colistin_All<-apply(Antibiotics[,c("colistinA","colistinB")], 1, max)
Antibiotics$Edeine_All<-apply(Antibiotics[,c("edeineA","edeineB","edeineD", "edeineF")], 1, max)
Antibiotics$Gentamicin_All<-apply(Antibiotics[,c("gentamicinA","gentamicinB","gentamicinC")], 1, max)
Antibiotics$Lividomycin_All<-apply(Antibiotics[,c("lividomycinA","lividomycinB")], 1, max)
Antibiotics$Penicillin_All<-apply(Antibiotics[,c("penicillinN","benzylpenicillin","isopenicillinN","phenoxymethylpenicillin","penicillin")], 1, max)
Antibiotics$Pristinamycin_All<-apply(Antibiotics[,c("pristinamycinIA","pristinamycinIB","pristinamycinIIA")], 1, max)
Antibiotics$Polymyxin_All<-apply(Antibiotics[,c("polymyxinB1", "polymyxinB2", "polymyxinB3", "polymyxinB4")], 1, max)
Antibiotics$Vernamycin_All<-apply(Antibiotics[,c("vernamycinBgamma", "vernamycinC")], 1, max)
Antibiotics$Sulfonamides_All<-apply(Antibiotics[,c("sulfacetamide",
                                                  "sulfadiazine",
                                                  "sulfadimidine",
                                                  "sulfadoxine",
                                                  "sulfamethizole",
                                                  "sulfamethoxazole",
                                                  "sulfasalazine",
                                                  "sulfisoxazole",
                                                  "mafenide")], 1, max)

rem_pat <- c("bacitracin","bleomycin","colistin","edeine","gentamicin","lividomycin","penicillin","pristinamycin","polymyxin","vernamycin")
Antibiotics2<-Antibiotics[, -grep(paste(rem_pat,collapse = "|"),colnames(Antibiotics),ignore.case = FALSE)]
Antibiotics2[c("sulfacetamide","sulfadiazine","sulfadimidine","sulfadoxine","sulfamethizole","sulfamethoxazole","sulfasalazine","sulfisoxazole", "mafenide")]<- list(NULL)  


#####MERGE ARO DATA WITH ANTIBIOTICS DATA


All_Antibiotic_Data <- merge(ARO_Counts, Antibiotics2, by=0, all=FALSE)


#2/remove antibiotics where sum of 0/1/2 = zero across all AROs. This removes all antibiotics for which resistance was not detected in any of the samples 


All_Antibiotic_Data2 = All_Antibiotic_Data[,!sapply(All_Antibiotic_Data, function(col) all(col == 0))]


###GIVE THE SAME PATTERN TO ALL ANTIBIOTIC COLUMNS SO THAT THE SAME FUNCTION CAN BE DONE FOR ALL THESE COLUMNS

####the following will work if the last sample before the first antibiotic is "UK_POP_POOL"

All_Antibiotic_Data2$Row.names<-NULL

#First_Index <- grep("UK_POP_POOL", colnames(All_Antibiotic_Data2))+1 ####AMEND THIS DEPENDING ON NEED TO START ON FIRST AB
First_Index <- grep("UK_POP_POOL_50M", colnames(All_Antibiotic_Data2))+1 ####AMEND THIS DEPENDING ON NEED TO START ON FIRST AB : FOR SUBANALYSES CHANGE THE GREP AS LAST COLUMN WITH SAMPLING COUNTS IS NOW WITH _50 ADDED TO IT
Second_Index<-ncol(All_Antibiotic_Data2)                             ####AMEND THIS DEPENDING ON NEED TO FINISH ON LAST AB

colnames(All_Antibiotic_Data2)[First_Index:Second_Index] <- paste("AB", colnames(All_Antibiotic_Data2[,c(First_Index:Second_Index)]), sep = "_")


###########AGGREGATE COUNTS BY ANTIBIOTIC


Long_AB <- gather(All_Antibiotic_Data2, key = "Antibiotic", value = value, First_Index:Second_Index) #####ADD HERE THE VALUES OF FIRST_INDEX AND SECOND_INDEX, WHICH
#####WILL BE DIFFERENT DEPENDING ON THE NUMBER OF SAMPLES IN THE
#####ANALYSIS

####CREATE THE RISK ADJUSTED DATASET, WHERE WE ONLY CONSIDER COUNTS OF GENES WHERE ANTIBIOTIC SCORE ==2 (I.E. CONFERS RESISTANCE TO ACCORDING TO CARD DATABASE)
####THIS DATASET WILL BE CALLED "RES_SCORE_TWO"AMR_DEF"

Long_AB1<-Long_AB                                                                             
Long_AB1$value<-as.factor(Long_AB1$value)
Long_AB1$value <- revalue(Long_AB1$value, c("0"="Zero", "1"="One", "2"="Two"))
names(Long_AB1)[names(Long_AB1)=="value"] <- "Resistance_Status"
Aggregated_Antibiotics1 <- Long_AB1 %>% group_by(Antibiotic,Resistance_Status) %>% summarise_all(sum)
Res_Score_Two<-filter(Aggregated_Antibiotics1, Resistance_Status == "Two")
Res_Score_Two$Resistance_Status<-NULL



####CREATE THE OVERALL DATASET, WHERE WE CONSIDER COUNTS OF GENES WHERE ANTIBIOTIC SCORE ==1 |==2 
####THIS DATASET WILL BE CALLED "AMR_ALL"

Long_AB2<-Long_AB                                                                             
Long_AB2$value<-as.factor(Long_AB2$value)
Long_AB2$value <- revalue(Long_AB2$value, c("0"="Zero", "1"="One", "2"="Two"))
names(Long_AB2)[names(Long_AB2)=="value"] <- "Resistance_Status"
levels(Long_AB2$Resistance_Status)[levels(Long_AB2$Resistance_Status)=="One"] <- "OneorTwo"
levels(Long_AB2$Resistance_Status)[levels(Long_AB2$Resistance_Status)=="Two"] <- "OneorTwo"
Aggregated_Antibiotics2 <- Long_AB2 %>% group_by(Antibiotic,Resistance_Status) %>% summarise_all(sum)
Res_Score_One_Two<-filter(Aggregated_Antibiotics2, Resistance_Status == "OneorTwo")
Res_Score_One_Two$Resistance_Status<-NULL




#ARCHIVE Antibiotic Datase == One and Antibiotic Dataset == One or Two 

write.csv(Res_Score_Two, "./Output Files/AMR_def.csv",row.names=F)
write.csv(Res_Score_One_Two, "./Output Files/AMR_all.csv",row.names=F)


##################CREATE THE DATASET FOR BAYESIAN MODELLING################################################################################################################################

Rcgc12<-Res_Score_One_Two

#Rcgc12<-Rcgc12[,c("Antibiotic","CAMBODIA_POP_POOL", "KENYA_POP_POOL", "UK_POP_POOL")]
Rcgc12<-Rcgc12[,c("Antibiotic","KENYA_30_SAMPLE_POOL_20M", "KENYA_30_SAMPLE_POOL_50M", "KENYA_30_SAMPLE_POOL",
                  "KENYA_POP_POOL_20M","KENYA_POP_POOL_50M", "UK_30_SAMPLE_POOL_20M", "UK_30_SAMPLE_POOL_50M",
                  "UK_30_SAMPLE_POOL","UK_POP_POOL_20M","UK_POP_POOL_50M")]

Rcgc12$Antibiotic <- gsub('AB_', '', Rcgc12$Antibiotic)

Tot_AMR12<-Rcgc12 %>%
  adorn_totals("row") %>%
  mutate_at(vars(-Antibiotic), list(~replace(., row_number() < n(),
                                       (.[-n()]/.[n()])))) 

Tot_AMR12<- Tot_AMR12 %>%
  filter(Antibiotic=="amikacin" | Antibiotic=="ampicillin" | Antibiotic=="cefotaxime" | Antibiotic=="cefoxitin" | Antibiotic=="cefpodoxime" | 
           Antibiotic=="ceftazidime" | Antibiotic=="ceftriaxone" | Antibiotic=="cefuroxime" | Antibiotic=="chloramphenicol" | Antibiotic=="ciprofloxacin" | 
           Antibiotic=="Gentamicin_All" | Antibiotic=="imipenem" | Antibiotic=="meropenem" | Antibiotic=="nalidixicacid" | Antibiotic=="nitrofurantoin" | Antibiotic=="trimethoprimsulfamethoxazole") %>% 
  droplevels()

#long_Tot_AMR12 <- Tot_AMR12 %>% gather(Setting, Rcgc_1_2,2:4)
long_Tot_AMR12 <- Tot_AMR12 %>% gather(Setting, Rcgc_1_2,2:11) # Needs to be changed to 11 as more datasets, i.e. columns (EvK 15 Feb 2024)

Rcgc2<-Res_Score_Two

#Rcgc2<-Rcgc2[,c("Antibiotic","CAMBODIA_POP_POOL", "KENYA_POP_POOL", "UK_POP_POOL")]
Rcgc2<-Rcgc2[,c("Antibiotic","KENYA_30_SAMPLE_POOL_20M", "KENYA_30_SAMPLE_POOL_50M", "KENYA_30_SAMPLE_POOL",
                  "KENYA_POP_POOL_20M","KENYA_POP_POOL_50M", "UK_30_SAMPLE_POOL_20M", "UK_30_SAMPLE_POOL_50M",
                  "UK_30_SAMPLE_POOL","UK_POP_POOL_20M","UK_POP_POOL_50M")]



Rcgc2$Antibiotic <- gsub('AB_', '', Rcgc2$Antibiotic)

Tot_AMR2<-Rcgc2 %>%
  adorn_totals("row") %>%
  mutate_at(vars(-Antibiotic), list(~replace(., row_number() < n(),
                                             (.[-n()]/.[n()])))) 

Tot_AMR2<- Tot_AMR2 %>%
  filter(Antibiotic=="amikacin" | Antibiotic=="ampicillin" | Antibiotic=="cefotaxime" | Antibiotic=="cefoxitin" | Antibiotic=="cefpodoxime" | 
           Antibiotic=="ceftazidime" | Antibiotic=="ceftriaxone" | Antibiotic=="cefuroxime" | Antibiotic=="chloramphenicol" | Antibiotic=="ciprofloxacin" | 
           Antibiotic=="Gentamicin_All" | Antibiotic=="imipenem" | Antibiotic=="meropenem" | Antibiotic=="nalidixicacid" | Antibiotic=="nitrofurantoin" | Antibiotic=="trimethoprimsulfamethoxazole") %>% 
  droplevels()

#long_Tot_AMR2 <- Tot_AMR2 %>% gather(Setting, Rcgc_2,2:4)
long_Tot_AMR2 <- Tot_AMR2 %>% gather(Setting, Rcgc_2,2:11)

#ByMo1<-merge(long_Tot_AMR12,long_Tot_AMR2,by=c("Antibiotic","Setting"),all=TRUE)
ByMo1<-merge(long_Tot_AMR12,long_Tot_AMR2,all=TRUE)


ByMo2<-ByMo1 %>% 
  mutate(Rcgc2_zero_replaced = case_when(Rcgc_2==0 ~ Rcgc_1_2, 
                                         is.na(Rcgc_2)~Rcgc_1_2, 
                                         TRUE ~ Rcgc_2))

ByMo2$Rcgc_2[is.na(ByMo2$Rcgc_2)] <- 0


Inf1<-read.csv("./Input Files/Original_paper/Infection_Data_For_Bayesian_Model.csv", header=T)
attach(Inf1)
fix(Inf1)

ByMo3<-merge(ByMo2,Inf1,by=c("Antibiotic","Setting"),all=TRUE)


######################################################################
# WAITING FOR THE RIGHT FILE TO BE SHARED BY KEVIN (EVK 15 February 2024)
Tx1<-read.csv("./Input Files/Original_paper/bracken_combined_reads.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
attach(Tx1)
fix(Tx1)

Tx2<-Tx1[,grepl("POOL",colnames(Tx1))]
Tx2 <- add_rownames(Tx2, "Taxonomy")
Tx3<- Tx2 %>%
  filter(grepl("k__Bacteria",Taxonomy)) %>% 
  droplevels()

d3<-filter(Tx3 [,2:7])
d3<-d3 %>% replace(is.na(.), 0)
AllBacteria<-as.data.frame(t(colSums(d3 [,])))


Tx3<-Tx3 %>% replace(is.na(.), 0)
ecol<-Tx3[grep("g__Escherichia;s__Escherichia coli", Tx3$Taxonomy), ]
ecol<-filter(ecol [,2:7])
kle<-Tx3[grep("g__Klebsiella;s__Klebsiella pneumoniae", Tx3$Taxonomy), ]
kle<-filter(kle [,2:7])

eb<-Tx3[grep("g__Enterobacter;", Tx3$Taxonomy), ]
eb<-filter(eb [,2:7])
eb<-as.data.frame(t(colSums(eb [,])))
sal<-Tx3[grep("g__Salmonella;", Tx3$Taxonomy), ]
sal<-filter(sal [,2:7])
sal<-as.data.frame(t(colSums(sal [,])))
f_E<-Tx3[grep("f__Enterobacteriaceae;", Tx3$Taxonomy), ]
f_E<-filter(f_E [,2:7])
f_E<-as.data.frame(t(colSums(f_E [,])))
o_E<-Tx3[grep("o__Enterobacterales;", Tx3$Taxonomy), ]
o_E<-filter(o_E [,2:7])
o_E<-as.data.frame(t(colSums(o_E [,])))

ClinicalGroups<-rbind(ecol,kle,eb,sal)
SumClinicalGroups<-as.data.frame(t(colSums(ClinicalGroups)))
Rtaxcolkleebsal<-SumClinicalGroups/AllBacteria
Rtaxf_E<-f_E/AllBacteria
Rtaxo_E<-o_E/AllBacteria

Rtaxcolkleebsal<-t(Rtaxcolkleebsal)
Rtaxcolkleebsal<-as.data.frame(Rtaxcolkleebsal)
Rtaxcolkleebsal<- Rtaxcolkleebsal %>% rownames_to_column("SampleID1")
colnames(Rtaxcolkleebsal)[2]<-"Rtax_colkleebsal_Bracken"

Rtaxf_E<-t(Rtaxf_E)
Rtaxf_E<-as.data.frame(Rtaxf_E)
Rtaxf_E<- Rtaxf_E %>% rownames_to_column("SampleID2")
colnames(Rtaxf_E)[2]<-"Rtax_f_E_Bracken"

Rtaxo_E<-t(Rtaxo_E)
Rtaxo_E<-as.data.frame(Rtaxo_E)
Rtaxo_E<- Rtaxo_E %>% rownames_to_column("SampleID3")
colnames(Rtaxo_E)[2]<-"Rtax_o_E_Bracken"


dataBracken<-cbind(Rtaxcolkleebsal,Rtaxf_E,Rtaxo_E)

names(dataBracken)[names(dataBracken) == "SampleID1"] <- "Setting"
dataBracken$SampleID2<-NULL
dataBracken$SampleID3<-NULL

dataBracken<- dataBracken %>%
  filter(grepl("POP_POOL",Setting)) %>% 
  droplevels()


ByMo4<-merge(ByMo3,dataBracken,by="Setting", all=TRUE)

ByMo4<-ByMo4[,c(1,2,3,4,5,8,9,10,6,7)]

#ARCHIVE DATA FOR BAYESIAN MODEL

write.csv(ByMo4, "./Output Files/Dataset_For_Bayesian_Model.csv",row.names=F)



#########################################################non-metric multidimensional scaling (NMDS) ############################

Intrg_3<-ARO_Counts


Intrg_4<-Intrg_3[, -grep("POOL", colnames(Intrg_3))]

patterns <- unique(substr(names(Intrg_4), 1, 2))
new <- sapply(patterns, function(xx) rowSums(Intrg_4[,grep(xx, names(Intrg_4)), drop=FALSE]))

new<-as.data.frame(new)
new$FC<-as.numeric(new$FC)
new$PK<-as.numeric(new$PK)
new$FCPK <- new$FC + new$PK
new$FC <- NULL
new$PK<-NULL
colnames(new)[colnames(new) == "FCPK"] <- "KENYA_Sum_Ind"
colnames(new)[colnames(new) == "NC"] <- "CAMBODIA_sum_Ind"
colnames(new)[colnames(new) == "UK"] <- "UK_Sum_Ind"

Intrg_5 <- Intrg_3[ , grepl( "POOL" , names( Intrg_3 ) ) ]
Intrg_6 <- merge(Intrg_5, new, by="row.names", all=TRUE)


#/transpose the matrix

Intrg_7 <- Intrg_6[,-1]##Intrg_7 <- Intrg_6[,-1]
rownames(Intrg_7) <- Intrg_6[,1]

counts_lat_cov<-t(Intrg_7)
fix(counts_lat_cov)


counts_lat_cov <- counts_lat_cov[ order(row.names(counts_lat_cov)), ]

classification=c(rep("Cambodia",3),rep("Kenya",3),rep("UK",3))

counts_lat_cov<-as.data.frame(counts_lat_cov)

total_col<-apply(counts_lat_cov[,], 1, sum)

relative_abundance = lapply(counts_lat_cov[,], function(x) {
  (x / total_col)*100
})

relative_abundance<-as.data.frame(relative_abundance)


set.seed(10)
nmds1 <- metaMDS(counts_lat_cov, binary = F,k = 2,try=1000)


set.seed(10)
nmds2 <- metaMDS(relative_abundance, binary = F,k = 2,try=1000)





