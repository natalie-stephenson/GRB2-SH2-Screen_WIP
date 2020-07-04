#########################################################################################################
# 													                                                                            #
# Script to identify proteins that contain SH3 domains and Proline Rich Motifs (subsetting TM proteins)	#
#													                                                                              #
#########################################################################################################
#													                                                                              # 
# AUTHOR LIST:												                                                                  #
#	Natalie Stephenson,   										                                                            #	
#													                                                                              #
# CURRENT STATUS OF SCRIPT:												                                                      #
#	                                				                                                 							#
# This script creates lists of proteins containing SH3 domains based on their GO:ID and PxxP motifs     #
# based on the sequences described in J. Teyra et al. (2017), Structure, 25(10).  A list of all proteins#
# containing the PxxP motif that also harbours a transmembrane domain, again based on GO:IDs, was       #
# produced to narrow this down to portential recepter interactions. There is also a section, currently  #
# blocked out, that looks at receptor tyrosine kinases (TKs) that contain these PxxP motifs.            #
#                                                                                                       #
#########################################################################################################


# Intro text.
cat("[INSERT PAPER DETAILS HERE] \n\n",
    "This script identifies proteins within the ensembl database that contain SH3 domains and proline rich motifs \n",
    sep="")


# Install required packages.
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
#biocLite(c("UniProt.ws", "RTCGA", "RTCGA.CNV"))
install.packages("dplyr")
#install.packages("stringr")
#install.packages("cgdsr")


# Load required libraries.
#library(UniProt.ws)
library(biomaRt)
library(dplyr)
#library(stringr)
#library(cgdsr)
#library(RTCGA)


# Load the ensembl data, to be used for the analysis later. NB: ACCESSED FEB2017 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# Create a gene list of proteins that have domain (pfam) information, genetic locations and sequences. NB: ACCESSED FEB2017
gene_list <- getBM(attributes = c('uniprot_gn', 'uniprotswissprot', 'go_id', 'pfam', 'chromosome_name', 'start_position', 'end_position'), filters = 'with_pfam', values = list(TRUE), mart = ensembl)
full_gene_sequences <- getSequence(id=gene_list$uniprotswissprot, type = "uniprotswissprot", seqType = "peptide", mart = ensembl)
complete_gene_list <- merge(gene_list, full_gene_sequences,by="uniprotswissprot")
write.csv(complete_gene_list, file = "complete_gene_list.csv")


####### Identify all proteins containing SH3 domains within the created gene list.
####### This list is filtered to produce a list of proteins containing SH3 domains regardless of the number of SH3 domains present.
####### PFAM accession codes used are PF00018 (SH3_1, SH3 domain), PF07653 (SH3_2, variant SH3 domain), PF12913 (SH3_6, SH3 domain - SH3b1 type), PF12914 (SH3_7, SH3 domain - SH3b2 type), PF13457 (SH3_8, SH3-like domain), PF14604 (SH3_9, variant SH3 domain). Not included PF08239, PF06347 and PF08460 (SH3_3-6), which are described as Bacterial SH3 domains).
###complete_gene_list <- read.csv(file = "complete_gene_list.csv")
###all_SH3_domains <- filter(complete_gene_list, complete_gene_list$pfam=="PF00018" | complete_gene_list$pfam=="PF07653" | complete_gene_list$pfam=="PF12913" | complete_gene_list$pfam=="PF12914" | complete_gene_list$pfam=="PF13457"| complete_gene_list$pfam=="PF14604")
###SH3_Domains_filtered <- all_SH3_domains[!duplicated(all_SH3_domains$uniprot_swissprot),]


#### SH3 Domains ####
#Identify all proteins containing SH3 domains based on their GO:ID (GO:0017124)
complete_gene_list <- read.csv(file = "complete_gene_list.csv")
GO_SH3_domains <- filter(complete_gene_list, complete_gene_list$go_id=="GO:0017124")
GO_SH3_Domains_filtered <- GO_SH3_domains[!duplicated(GO_SH3_domains$uniprot_genename),]
GO_SH3_Domains_filtered$X <- NULL

###SH3_domain_proteins <- subset (SH3_Domains_filtered, select=c("uniprot_gn"))
###SH3_domain_proteins <- filter(SH3_domain_proteins, !SH3_domain_proteins$uniprot_genename=="")
###write.csv(SH3_domain_proteins, file = "SH3_domain_proteins_GO.csv")

#Identify all proteins containing SH3 domains based on their PFAM ID (PF00018)
PFAM_SH3_domains <- filter(complete_gene_list, complete_gene_list$pfam=="PF00018")
PFAM_SH3_Domains_filtered <- PFAM_SH3_domains[!duplicated(PFAM_SH3_domains$uniprot_genename),]
PFAM_SH3_Domains_filtered$X <- NULL

#Merge SH3 domain lists based on uniprot genename and save files.
SH3_domains_merged <- rbind(PFAM_SH3_Domains_filtered, GO_SH3_Domains_filtered)
SH3_domains_merged_dedup <- SH3_domains_merged[!duplicated(SH3_domains_merged$uniprot_genename),]
write.csv(SH3_domains_merged_dedup, file="SH3_domain_containing_proteins.csv")
SH3_domain_merged_dedup_list <- subset (SH3_domains_merged_dedup, select=c("uniprot_genename"))
write.csv(SH3_domain_merged_dedup_list, file = "SH3_domain_containing_proteins_list.csv")



#### PxxP Motifs ####
# Identify all proteins containing Proline Rich (PR) motifs through sequence searching for PxxP motifs. 
PR_containing <- complete_gene_list[grep("P[A-Z]{2}P", complete_gene_list$peptide), ]                       
PR_containing_filtered <- PR_containing[!duplicated(PR_containing$uniprot_genename),]  
PR_containing_filtered_named <- filter(PR_containing_filtered, !PR_containing_filtered$uniprot_genename=="")
write.csv(PR_containing_filtered_named, file = "PR_containing_proteins.csv")
PR_containing_filtered_list <- subset (PR_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_containing_filtered_list, file="PR_containing_proteins_list.csv")


# Paper identifying different classes of PxxP motif recognised by Sh3 domains: Teyra et al., 2017, Structure, 25, 1598-1610
# Identify all proteins containing Class 1 Proline Rich (PR) motifs through sequence searching for RxxPxxP motifs. 
PR_Class1_containing <- complete_gene_list[grep("R[A-Z]{2}P[A-Z]{2}P", complete_gene_list$peptide), ]                       
PR_Class1_containing_filtered <- PR_Class1_containing[!duplicated(PR_Class1_containing$uniprot_genename),]  
PR_Class1_containing_filtered_named <- filter(PR_Class1_containing_filtered, !PR_Class1_containing_filtered$uniprot_genename=="")
write.csv(PR_Class1_containing_filtered_named, file = "PR_Class1_containing_proteins.csv")
PR_Class1_containing_filtered_list <- subset (PR_Class1_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class1_containing_filtered_list, file="PR_Class1_containing_proteins_list.csv")


# Identify all proteins containing Class 2 Proline Rich (PR) motifs through sequence searching for PxxPxR motifs. 
PR_Class2_containing <- complete_gene_list[grep("P[A-Z]{2}P[A-Z]{1}R", complete_gene_list$peptide), ]                       
PR_Class2_containing_filtered <- PR_Class2_containing[!duplicated(PR_Class2_containing$uniprot_genename),]  
PR_Class2_containing_filtered_named <- filter(PR_Class2_containing_filtered, !PR_Class2_containing_filtered$uniprot_genename=="")
write.csv(PR_Class2_containing_filtered_named, file = "PR_Class2_containing_proteins.csv")
PR_Class2_containing_filtered_list <- subset (PR_Class2_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class2_containing_filtered_list, file="PR_Class2_containing_proteins_list.csv")


# Merging Class 1 and 2 Proline Rich motif sequence lists to produce a canonical list.
PR_Class1and2_containing_filtered_named <- rbind(PR_Class1_containing_filtered_named, PR_Class2_containing_filtered_named)
PR_Class1and2_containing_filtered_named_dedup <- PR_Class1and2_containing_filtered_named[!duplicated(PR_Class1and2_containing_filtered_named),]
write.csv(PR_Class1and2_containing_filtered_named_dedup, file="PR_Class1and2_containing_filtered_named_dedup.csv")
PR_Class1and2_containing_filtered_list <- subset (PR_Class1and2_containing_filtered_named_dedup, select=c("uniprot_genename")) 
write.csv(PR_Class1and2_containing_filtered_list, file="PR_Class1and2_containing_proteins_list.csv")

# Identify all proteins containing Class 3 Proline Rich (PR) motifs through sequence searching for RxxxxxP motifs. 
PR_Class3_containing <- complete_gene_list[grep("R[A-Z]{5}P", complete_gene_list$peptide), ]                       
PR_Class3_containing_filtered <- PR_Class3_containing[!duplicated(PR_Class3_containing$uniprot_genename),]  
PR_Class3_containing_filtered_named <- filter(PR_Class3_containing_filtered, !PR_Class3_containing_filtered$uniprot_genename=="")
write.csv(PR_Class3_containing_filtered_named, file = "PR_Class3_containing_proteins.csv")
PR_Class3_containing_filtered_list <- subset (PR_Class3_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class3_containing_filtered_list, file="PR_Class3_containing_proteins_list.csv")


# Identify all proteins containing Class 4 Proline Rich (PR) motifs through sequence searching for PxxxxR motifs. 
PR_Class4_containing <- complete_gene_list[grep("P[A-Z]{4}P", complete_gene_list$peptide), ]                       
PR_Class4_containing_filtered <- PR_Class4_containing[!duplicated(PR_Class4_containing$uniprot_genename),]  
PR_Class4_containing_filtered_named <- filter(PR_Class4_containing_filtered, !PR_Class4_containing_filtered$uniprot_genename=="")
write.csv(PR_Class4_containing_filtered_named, file = "PR_Class4_containing_proteins.csv")
PR_Class4_containing_filtered_list <- subset (PR_Class4_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class4_containing_filtered_list, file="PR_Class4_containing_proteins_list.csv")


# Identify all proteins containing Class 5 Proline Rich (PR) motifs through sequence searching for RxxPxxx motifs. 
PR_Class5_containing <- complete_gene_list[grep("R[A-Z]{2}P[A-Z]{3}", complete_gene_list$peptide), ]                       
PR_Class5_containing_filtered <- PR_Class5_containing[!duplicated(PR_Class5_containing$uniprot_genename),]  
PR_Class5_containing_filtered_named <- filter(PR_Class5_containing_filtered, !PR_Class5_containing_filtered$uniprot_genename=="")
write.csv(PR_Class5_containing_filtered_named, file = "PR_Class5_containing_proteins.csv")
PR_Class5_containing_filtered_list <- subset (PR_Class5_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class5_containing_filtered_list, file="PR_Class5_containing_proteins_list.csv")


# Identify all proteins containing Class 6 Proline Rich (PR) motifs through sequence searching for xxxPxR motifs. 
PR_Class6_containing <- complete_gene_list[grep("[A-Z]{3}P[A-Z]{1}R", complete_gene_list$peptide), ]                       
PR_Class6_containing_filtered <- PR_Class6_containing[!duplicated(PR_Class2_containing$uniprot_genename),]  
PR_Class6_containing_filtered_named <- filter(PR_Class6_containing_filtered, !PR_Class6_containing_filtered$uniprot_genename=="")
write.csv(PR_Class6_containing_filtered_named, file = "PR_Class6_containing_proteins.csv")
PR_Class6_containing_filtered_list <- subset (PR_Class6_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class6_containing_filtered_list, file="PR_Class6_containing_proteins_list.csv")


# Identify all proteins containing Class 7 Proline Rich (PR) motifs through sequence searching for PxKP motifs. 
PR_Class7_containing <- complete_gene_list[grep("P[A-Z]{1}KP", complete_gene_list$peptide), ]                       
PR_Class7_containing_filtered <- PR_Class7_containing[!duplicated(PR_Class7_containing$uniprot_genename),]  
PR_Class7_containing_filtered_named <- filter(PR_Class7_containing_filtered, !PR_Class7_containing_filtered$uniprot_genename=="")
write.csv(PR_Class7_containing_filtered_named, file = "PR_Class7_containing_proteins.csv")
PR_Class7_containing_filtered_list <- subset (PR_Class7_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class7_containing_filtered_list, file="PR_Class7_containing_proteins_list.csv")


# Identify all proteins containing Class 8 Proline Rich (PR) motifs through sequence searching for PxxPxxP motifs. 
PR_Class8_containing <- complete_gene_list[grep("P[A-Z]{2}P[A-Z]{2}P", complete_gene_list$peptide), ]                       
PR_Class8_containing_filtered <- PR_Class8_containing[!duplicated(PR_Class8_containing$uniprot_genename),]  
PR_Class8_containing_filtered_named <- filter(PR_Class8_containing_filtered, !PR_Class8_containing_filtered$uniprot_genename=="")
write.csv(PR_Class8_containing_filtered_named, file = "PR_Class8_containing_proteins.csv")
PR_Class8_containing_filtered_list <- subset (PR_Class8_containing_filtered_named, select=c("uniprot_genename")) 
write.csv(PR_Class8_containing_filtered_list, file="PR_Class8_containing_proteins_list.csv")

#Identify all transmembrane proteins containing PR motifs
PR_TM_containing <- filter(PR_containing, PR_containing$go_id=="GO:0016021")
PR_TM_containing_nodup <- PR_TM_containing[!duplicated(PR_TM_containing$uniprot_genename),]
PR_TM_containing_nodup_named <- filter(PR_TM_containing_nodup, !PR_TM_containing_nodup$uniprot_genename=="")
write.csv(PR_TM_containing_nodup_named, file = "PR_TM_containing_proteins.csv")
PR_TM_proteins <- subset(PR_TM_containing_nodup_named, select=c("uniprot_genename"))
write.csv(PR_TM_proteins, file = "PR_TM_containing_proteins_list.csv")


####Identifying just tyrosine kinases that contain PR motifs
##PR_YK_containing <- filter(PR_containing, PR_containing$pfam=="PF07714")
##write.csv(PR_TM_containing, file="PR_YK_containing.csv")
##PR_YK_proteins <- subset(PR_YK_containing, select=c("uniprot_genename"))
##PR_YK_proteins <- filter(PR_YK_proteins, !PR_YK_proteins$uniprot_genename=="")
##write.csv(PR_YK_proteins, file="PR_YK_proteins.csv")
### Join SH3 and PR_YK files
##SH3_PRTM_proteins <- rbind(SH3_domain_proteins, PR_TM_proteins)
##SH3_PRTM_proteins <- filter(SH3_PRTM_proteins, !SH3_PRTM_proteins$uniprot_genename=="")
##write.csv(SH3_PRTM_proteins, file = "SH3_PRTM_proteins.csv")


# Prints Session information for this code.
sessionInfo()
