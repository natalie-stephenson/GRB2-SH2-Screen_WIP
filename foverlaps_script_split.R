library(data.table)
LUAD_CNV_dt <- read.csv(file = "LUAD_CNV_dt.csv")
complete_gene_list_dt <- read.csv(file = "complete_gene_list_dt.csv")
LUAD_CNV_dt$X <- NULL
LUAD_CNV_dt <- data.table(LUAD_CNV_dt)
complete_gene_list_dt <- data.table(complete_gene_list_dt)

LUAD_CNV_5000 <- LUAD_CNV_dt[1:5000,]
LUAD_CNV_10000 <- LUAD_CNV_dt[5001:10000,]
LUAD_CNV_15000 <- LUAD_CNV_dt[10001:15000,]
LUAD_CNV_20000 <- LUAD_CNV_dt[15001:20000,]
LUAD_CNV_25000 <- LUAD_CNV_dt[20001:25000,]
LUAD_CNV_30000 <- LUAD_CNV_dt[25001:30000,]
LUAD_CNV_35000 <- LUAD_CNV_dt[30001:35000,]
LUAD_CNV_40000 <- LUAD_CNV_dt[35001:40000,]
LUAD_CNV_45000 <- LUAD_CNV_dt[40001:45000,]
LUAD_CNV_50000 <- LUAD_CNV_dt[45001:50000,]
LUAD_CNV_55000 <- LUAD_CNV_dt[50001:55000,]
LUAD_CNV_60000 <- LUAD_CNV_dt[55001:60000,]
LUAD_CNV_65000 <- LUAD_CNV_dt[60001:65000,]
LUAD_CNV_70000 <- LUAD_CNV_dt[65001:70000,]
LUAD_CNV_75000 <- LUAD_CNV_dt[70001:75000,]
LUAD_CNV_80000 <- LUAD_CNV_dt[75001:80000,]
LUAD_CNV_85000 <- LUAD_CNV_dt[80001:85000,]
LUAD_CNV_90000 <- LUAD_CNV_dt[85001:90000,]
LUAD_CNV_95000 <- LUAD_CNV_dt[90001:95000,]
LUAD_CNV_100000 <- LUAD_CNV_dt[95001:100000,]
LUAD_CNV_105000 <- LUAD_CNV_dt[100001:105000,]
LUAD_CNV_110000 <- LUAD_CNV_dt[105001:110000,]
LUAD_CNV_115000 <- LUAD_CNV_dt[110001:115000,]
LUAD_CNV_120000 <- LUAD_CNV_dt[115001:120000,]
LUAD_CNV_125000 <- LUAD_CNV_dt[120001:125000,]
LUAD_CNV_130000 <- LUAD_CNV_dt[125001:130000,]
LUAD_CNV_135000 <- LUAD_CNV_dt[130001:135000,]
LUAD_CNV_140000 <- LUAD_CNV_dt[135001:140000,]
LUAD_CNV_145000 <- LUAD_CNV_dt[140001:145000,]
LUAD_CNV_150000 <- LUAD_CNV_dt[145001:150000,]
LUAD_CNV_155000 <- LUAD_CNV_dt[150001:155000,]
LUAD_CNV_160000 <- LUAD_CNV_dt[155001:160000,]
LUAD_CNV_165000 <- LUAD_CNV_dt[160001:165000,]
LUAD_CNV_170000 <- LUAD_CNV_dt[165001:170000,]
LUAD_CNV_175000 <- LUAD_CNV_dt[170001:175000,]
LUAD_CNV_180000 <- LUAD_CNV_dt[175001:180000,]
LUAD_CNV_185000 <- LUAD_CNV_dt[180001:185000,]
LUAD_CNV_190000 <- LUAD_CNV_dt[185001:190000,]
LUAD_CNV_195000 <- LUAD_CNV_dt[190001:191418,]


setkey(LUAD_CNV_5000, Chromosome, Start, End)
LUAD_CNV_genenames1 <- foverlaps(complete_gene_list_dt, LUAD_CNV_5000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_10000, Chromosome, Start, End)
LUAD_CNV_genenames2 <- foverlaps(complete_gene_list_dt, LUAD_CNV_10000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_15000, Chromosome, Start, End)
LUAD_CNV_genenames3 <- foverlaps(complete_gene_list_dt, LUAD_CNV_15000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_20000, Chromosome, Start, End)
LUAD_CNV_genenames4 <- foverlaps(complete_gene_list_dt, LUAD_CNV_20000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_25000, Chromosome, Start, End)
LUAD_CNV_genenames5 <- foverlaps(complete_gene_list_dt, LUAD_CNV_25000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_30000, Chromosome, Start, End)
LUAD_CNV_genenames6 <- foverlaps(complete_gene_list_dt, LUAD_CNV_30000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_35000, Chromosome, Start, End)
LUAD_CNV_genenames7 <- foverlaps(complete_gene_list_dt, LUAD_CNV_35000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_40000, Chromosome, Start, End)
LUAD_CNV_genenames8 <- foverlaps(complete_gene_list_dt, LUAD_CNV_40000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_45000, Chromosome, Start, End)
LUAD_CNV_genenames9 <- foverlaps(complete_gene_list_dt, LUAD_CNV_45000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_50000, Chromosome, Start, End)
LUAD_CNV_genenames10 <- foverlaps(complete_gene_list_dt, LUAD_CNV_50000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_55000, Chromosome, Start, End)
LUAD_CNV_genenames11 <- foverlaps(complete_gene_list_dt, LUAD_CNV_55000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_60000, Chromosome, Start, End)
LUAD_CNV_genenames12 <- foverlaps(complete_gene_list_dt, LUAD_CNV_60000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_65000, Chromosome, Start, End)
LUAD_CNV_genenames13 <- foverlaps(complete_gene_list_dt, LUAD_CNV_65000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_70000, Chromosome, Start, End)
LUAD_CNV_genenames14 <- foverlaps(complete_gene_list_dt, LUAD_CNV_70000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_75000, Chromosome, Start, End)
LUAD_CNV_genenames15 <- foverlaps(complete_gene_list_dt, LUAD_CNV_75000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_80000, Chromosome, Start, End)
LUAD_CNV_genenames16 <- foverlaps(complete_gene_list_dt, LUAD_CNV_80000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_85000, Chromosome, Start, End)
LUAD_CNV_genenames17 <- foverlaps(complete_gene_list_dt, LUAD_CNV_85000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_90000, Chromosome, Start, End)
LUAD_CNV_genenames18 <- foverlaps(complete_gene_list_dt, LUAD_CNV_90000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_95000, Chromosome, Start, End)
LUAD_CNV_genenames19 <- foverlaps(complete_gene_list_dt, LUAD_CNV_95000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_100000, Chromosome, Start, End)
LUAD_CNV_genenames20 <- foverlaps(complete_gene_list_dt, LUAD_CNV_100000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_105000, Chromosome, Start, End)
LUAD_CNV_genenames21 <- foverlaps(complete_gene_list_dt, LUAD_CNV_105000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_110000, Chromosome, Start, End)
LUAD_CNV_genenames22 <- foverlaps(complete_gene_list_dt, LUAD_CNV_110000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_115000, Chromosome, Start, End)
LUAD_CNV_genenames23 <- foverlaps(complete_gene_list_dt, LUAD_CNV_115000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_120000, Chromosome, Start, End)
LUAD_CNV_genenames24 <- foverlaps(complete_gene_list_dt, LUAD_CNV_120000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_125000, Chromosome, Start, End)
LUAD_CNV_genenames25 <- foverlaps(complete_gene_list_dt, LUAD_CNV_125000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_130000, Chromosome, Start, End)
LUAD_CNV_genenames26 <- foverlaps(complete_gene_list_dt, LUAD_CNV_130000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_135000, Chromosome, Start, End)
LUAD_CNV_genenames27 <- foverlaps(complete_gene_list_dt, LUAD_CNV_135000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_140000, Chromosome, Start, End)
LUAD_CNV_genenames28 <- foverlaps(complete_gene_list_dt, LUAD_CNV_140000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_145000, Chromosome, Start, End)
LUAD_CNV_genenames29 <- foverlaps(complete_gene_list_dt, LUAD_CNV_145000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_150000, Chromosome, Start, End)
LUAD_CNV_genenames30 <- foverlaps(complete_gene_list_dt, LUAD_CNV_150000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_155000, Chromosome, Start, End)
LUAD_CNV_genenames31 <- foverlaps(complete_gene_list_dt, LUAD_CNV_155000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_160000, Chromosome, Start, End)
LUAD_CNV_genenames32 <- foverlaps(complete_gene_list_dt, LUAD_CNV_160000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_165000, Chromosome, Start, End)
LUAD_CNV_genenames33 <- foverlaps(complete_gene_list_dt, LUAD_CNV_165000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_170000, Chromosome, Start, End)
LUAD_CNV_genenames34 <- foverlaps(complete_gene_list_dt, LUAD_CNV_170000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_175000, Chromosome, Start, End)
LUAD_CNV_genenames35 <- foverlaps(complete_gene_list_dt, LUAD_CNV_175000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_180000, Chromosome, Start, End)
LUAD_CNV_genenames36 <- foverlaps(complete_gene_list_dt, LUAD_CNV_180000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_185000, Chromosome, Start, End)
LUAD_CNV_genenames37 <- foverlaps(complete_gene_list_dt, LUAD_CNV_185000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_190000, Chromosome, Start, End)
LUAD_CNV_genenames38 <- foverlaps(complete_gene_list_dt, LUAD_CNV_190000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

setkey(LUAD_CNV_195000, Chromosome, Start, End)
LUAD_CNV_genenames39 <- foverlaps(complete_gene_list_dt, LUAD_CNV_195000, by.x = c("chromosome_name", "start_position", "end_position"), by.y = c("Chromosome", "Start", "End"),  maxgap=0, minoverlap=1, type = "any", mult = "all", nomatch=0, which=FALSE)

write.csv(LUAD_CNV_genenames1, file="LUAD_CNV_genenames1.csv")
write.csv(LUAD_CNV_genenames2, file="LUAD_CNV_genenames2.csv")
write.csv(LUAD_CNV_genenames3, file="LUAD_CNV_genenames3.csv")
write.csv(LUAD_CNV_genenames4, file="LUAD_CNV_genenames4.csv")
write.csv(LUAD_CNV_genenames5, file="LUAD_CNV_genenames5.csv")
write.csv(LUAD_CNV_genenames6, file="LUAD_CNV_genenames6.csv")
write.csv(LUAD_CNV_genenames7, file="LUAD_CNV_genenames7.csv")
write.csv(LUAD_CNV_genenames8, file="LUAD_CNV_genenames8.csv")
write.csv(LUAD_CNV_genenames9, file="LUAD_CNV_genenames9.csv")
write.csv(LUAD_CNV_genenames10, file="LUAD_CNV_genenames10.csv")
write.csv(LUAD_CNV_genenames11, file="LUAD_CNV_genenames11.csv")
write.csv(LUAD_CNV_genenames12, file="LUAD_CNV_genenames12.csv")
write.csv(LUAD_CNV_genenames13, file="LUAD_CNV_genenames13.csv")
write.csv(LUAD_CNV_genenames14, file="LUAD_CNV_genenames14.csv")
write.csv(LUAD_CNV_genenames15, file="LUAD_CNV_genenames15.csv")
write.csv(LUAD_CNV_genenames16, file="LUAD_CNV_genenames16.csv")
write.csv(LUAD_CNV_genenames17, file="LUAD_CNV_genenames17.csv")
write.csv(LUAD_CNV_genenames18, file="LUAD_CNV_genenames18.csv")
write.csv(LUAD_CNV_genenames19, file="LUAD_CNV_genenames19.csv")
write.csv(LUAD_CNV_genenames20, file="LUAD_CNV_genenames20.csv")
write.csv(LUAD_CNV_genenames21, file="LUAD_CNV_genenames21.csv")
write.csv(LUAD_CNV_genenames22, file="LUAD_CNV_genenames22.csv")
write.csv(LUAD_CNV_genenames23, file="LUAD_CNV_genenames23.csv")
write.csv(LUAD_CNV_genenames24, file="LUAD_CNV_genenames24.csv")
write.csv(LUAD_CNV_genenames25, file="LUAD_CNV_genenames25.csv")
write.csv(LUAD_CNV_genenames26, file="LUAD_CNV_genenames26.csv")
write.csv(LUAD_CNV_genenames27, file="LUAD_CNV_genenames27.csv")
write.csv(LUAD_CNV_genenames28, file="LUAD_CNV_genenames28.csv")
write.csv(LUAD_CNV_genenames29, file="LUAD_CNV_genenames29.csv")
write.csv(LUAD_CNV_genenames30, file="LUAD_CNV_genenames30.csv")
write.csv(LUAD_CNV_genenames31, file="LUAD_CNV_genenames31.csv")
write.csv(LUAD_CNV_genenames32, file="LUAD_CNV_genenames32.csv")
write.csv(LUAD_CNV_genenames33, file="LUAD_CNV_genenames33.csv")
write.csv(LUAD_CNV_genenames34, file="LUAD_CNV_genenames34.csv")
write.csv(LUAD_CNV_genenames35, file="LUAD_CNV_genenames35.csv")
write.csv(LUAD_CNV_genenames36, file="LUAD_CNV_genenames36.csv")
write.csv(LUAD_CNV_genenames37, file="LUAD_CNV_genenames37.csv")
write.csv(LUAD_CNV_genenames38, file="LUAD_CNV_genenames38.csv")
write.csv(LUAD_CNV_genenames39, file="LUAD_CNV_genenames39.csv")
