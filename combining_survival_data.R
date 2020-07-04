## Practice script - getting patient survival data into CNV/genenames scripts.

#Install dependencies
library(data.table)

#read in required files
LUAD_clinical_patientsurvival <- read.csv("LUAD_clinical_patient_survival.csv")
LUAD_CNV_10000 <- read.csv("LUAD_CNV_genenames_bind.csv")

#Splitting Sample ID into information about the patient - Patient ID and sample type.
LUAD_CNV_10000_samplesplit <- within(LUAD_CNV_10000, {
  Patient_ID <- substr(LUAD_CNV_10000$Sample, 1, 12)
  Sample_Type_ID <- substr(LUAD_CNV_10000$Sample, 14, 15)
})

#Reducing the data down into components required for next step (#### FOR GENENAMES ADD "uniprot_genename" AT END ####)
LUAD_CNV_10000_samplesplit <- LUAD_CNV_10000_samplesplit[, c("Patient_ID", "Sample_Type_ID", "Sample", "Segment_Mean")]

#Characterising the Sample Type ID into Sample Type (tumour vs normal vs control)
LUAD_CNV_10000_samplesplit <- LUAD_CNV_10000_samplesplit %>%
  mutate(Sample_Type=Sample_Type_ID)

for (i in c("01", "02", "03", "04", "05", "06", "07", "08", "09")){
  LUAD_CNV_10000_samplesplit$Sample_Type <- ifelse (LUAD_CNV_10000_samplesplit$Sample_Type == i, "Tumour", LUAD_CNV_10000_samplesplit$Sample_Type)
}

for (j in c("10", "11", "12", "13", "14", "15", "16", "17", "18", "19")){
  LUAD_CNV_10000_samplesplit$Sample_Type <- ifelse (LUAD_CNV_10000_samplesplit$Sample_Type == j, "Normal", LUAD_CNV_10000_samplesplit$Sample_Type)
}

for (k in c("20", "21", "22", "23", "24", "25", "26", "27", "28", "29")){
  LUAD_CNV_10000_samplesplit$Sample_Type <- ifelse (LUAD_CNV_10000_samplesplit$Sample_Type == k, "Control", LUAD_CNV_10000_samplesplit$Sample_Type)
}

LUAD_clinical_patientsurvival <- data.table(LUAD_clinical_patientsurvival)
LUAD_CNV_10000_samplesplit <- data.table(LUAD_CNV_10000_samplesplit)
colnames(LUAD_clinical_patientsurvival)[1] <- "Patient_ID" 

merge <- merge(LUAD_CNV_10000_samplesplit, LUAD_clinical_patientsurvival, by="Patient_ID")
write.csv(merge, file="LUAD_CNV_genenames_suvival.csv")
