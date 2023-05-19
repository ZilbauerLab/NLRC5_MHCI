# SAMPLE FILTERING OF QCed, BATCH CORRECTED DATASET PRODUCED BY RACHEL EDGAR
# QCed dataset: "threebatch_combined_organoids_combatted.Rdata"
# Current metadata: "AllEPICOrganoid_FullMetadata_27Sep22.RData"
# 07Oct22


# SET ENVIRONMENT, LOAD DATA & RENAME INPUT:

library(plyr)

load("~/Documents/NLRC5_PAPER/METHYLATION_ARRAY_DATA/threebatch_combined_organoids_combatted.Rdata")
load("~/Documents/NLRC5_PAPER/METHYLATION_ARRAY_DATA/AllEPICOrganoid_FullMetadata_27Sep22.RData")

all.beta <- combat_organoid_Beta
re.meta <- epic.organoid_combined
all.ok.samples <- subset(epic.samples, array.id %in% colnames(all.beta))

dim(all.beta)
nrow(all.ok.samples)


# SAMPLE WITH MISSING UP TO DATE METADATA:
# Missing as it was removed from the older dataset before we realised that the T202 & T203 samples had been switched
# For more detail on this, see: 07_ALL-ORGANOID_PaediatricSCTI_DatasetUpdates.v1.0.R
# T202 (real): CD, age 12, F
# T203 (real): CD, age 14, M

missing <- subset(colnames(all.beta), !(colnames(all.beta) %in% all.ok.samples$array.id))
missing <- subset(re.meta, array.id==missing)

missing

colnames(epic.samples)
head(subset(epic.samples, sentrix.id=="205605870061"),n=1)[,c("sentrix.id","cgs.submission.id")]
unique(subset(epic.samples, case.no=="T203"))[,c("case.no","DUO.Inflammation","TI.Inflammation","SC.Inflammation","Biologics","Surgery","AZA","Perianal.Disease","Disease.Severity")]

# Update "epic.samples" with the missing sample
# NB, although this sample is labelled as T202, it's actually T203 (M, age=14, CD), see notes at top of methylation section

update <- data.frame(case.no="T203",array.id=missing$array.id,diagnosis="CD",sex="M",age=14,tissue="TI",sample.type="organoid",sampling.time.point="original",passage.or.rescope.no="P6",fraction=NA,age.group="paediatric",alternative.id=NA,pdata.id="T203 TI P6",treatment="NT",array.type="EPIC",sentrix.id=missing$sentrix_ID,sentrix.pos=missing$sentrix.pos,notes=NA,cgs.submission.id="METH-K.NAYAK-20035_RE",sample.id="T203_TI",mean.detectionP=NA,sample.id.3="T203_TI_NT",DUO.Inflammation="N",TI.Inflammation="Macroscopic",SC.Inflammation="Macroscopic",Biologics="Y",Surgery="N",AZA="Y",Perianal.Disease="N",Disease.Severity="severe",aza.trial.id="T203_NT",RE.batch="RE")

ncol(update)
dim(epic.samples)

epic.samples <- rbind(epic.samples,update)

nrow(epic.samples)

# Quick sanity check of T202 & T203:

subset(epic.samples, case.no=="T202" | case.no=="T203")[,c("case.no","sample.id","sample.id.3","sex","age","passage.or.rescope.no","notes","pdata.id")]


# METADATA SIMPLIFICATION & UPDATES:

epic.samples$rescope <- NA
epic.samples$rescope[epic.samples$sampling.time.point=="rescope"] <- epic.samples$passage.or.rescope.no[epic.samples$sampling.time.point=="rescope"]
epic.samples$rescope <- gsub(".P1","",epic.samples$rescope)
epic.samples$rescope <- gsub(".P2","",epic.samples$rescope)
epic.samples$rescope <- gsub(".P4","",epic.samples$rescope)
epic.samples <- rename(epic.samples, c("passage.or.rescope.no"="passage"))
epic.samples$passage <- gsub("RE1.","",epic.samples$passage)
epic.samples$passage <- gsub("RE3.","",epic.samples$passage)

re.stuff <- re.meta[,c("array.id","sample_ID","det_pval","passage.or.rescope.no_numeric","Sample_Name","Wnt.type","condition","Biobank.Rachel.Replicates","batch")]
colnames(re.stuff) <- c("array.id","RE.sample.id","RE.mean.detectionP","numeric.passage","RE.sample.name","Wnt.type","RE.condition","RE.Biobank.Replicates","RE.batch.2")

epic.samples <- join(epic.samples, re.stuff, type="left", match="all")

table(epic.samples$RE.batch,epic.samples$RE.batch.2)
table(epic.samples$treatment,epic.samples$RE.condition)
table(epic.samples$rescope)
table(epic.samples$sample.type,epic.samples$sampling.time.point)
table(epic.samples$passage,epic.samples$numeric.passage)
subset(epic.samples, sampling.time.point=="rescope")[,c("sample.id","sampling.time.point","passage")]

# There are a couple of samples with mislabelled passage in the new metadata:

subset(epic.samples, passage=="P1" & numeric.passage==2 | passage=="P2" & numeric.passage==1)[,c("case.no","sample.id","passage","numeric.passage","pdata.id","array.id")]

# Correct this:

epic.samples$passage[epic.samples$array.id=="205605880094_R06C01"] <- "P2"
epic.samples$passage[epic.samples$array.id=="205605880096_R08C01"] <- "P1"

table(epic.samples$passage,epic.samples$numeric.passage)

# Remove the redundant batch column & update "all.ok.samples":

epic.samples$RE.batch.2 <- NULL
all.ok.samples <- subset(epic.samples, array.id %in% colnames(all.beta))

dim(all.beta)
nrow(re.meta)
nrow(epic.samples)
nrow(all.ok.samples)


# SAMPLE FILTERING
# Based on Rachel's filtering:
#	- exclude any samples with diagnosis "Other.GI"
#	- remove treated/differentiated samples, but leave in frozen biopsies
#	- where there are at diagnosis & rescope versions of the same sample, retain the at diagnosis sample
#	- where there are duplicate samples, retain the one with the lowest passage
#	- remove duplicates (preferentially removing frozen biopsies), retaining the newest sample (NEW > RE > OLD)
# 	- exclude foetal & neonatal samples
# 	- check rescopes

# EXCLUDE "Other GI":

samples <- subset(all.ok.samples, diagnosis!="Other.GI")

table(samples$diagnosis)
nrow(all.ok.samples)
nrow(samples)


# REMOVE TREATED/DIFFERENTIATED SAMPLES:
# Retaining frozen biopsies
# Check & compare the up to date metadata with Rachel's metadata to make sure there are no discrepancies

nrow(subset(samples, treatment=="NT"))
table(samples$treatment,samples$RE.condition)

# "T074 RE3 DUO P2 FB" isn't labelled as a frozen biopsy even though it should be
# Add "differentiated" & "frozen.biopsy" columns & drop the "RE.condition" column:

epic.samples$differentiated <- "N"
epic.samples$differentiated[epic.samples$RE.condition=="D"] <- "Y"
epic.samples$frozen.biopsy <- "N"
epic.samples$frozen.biopsy[epic.samples$RE.condition=="FB"] <- "Y"
epic.samples$frozen.biopsy[epic.samples$notes=="frozen biopsy"] <- "Y"
epic.samples$RE.condition <- NULL

table(epic.samples$differentiated)
table(epic.samples$frozen.biopsy)

all.ok.samples <- subset(epic.samples, array.id %in% all.ok.samples$array.id)
samples <- subset(epic.samples, array.id %in% samples$array.id)

nrow(epic.samples)
nrow(all.ok.samples)
nrow(samples)

# Drop the treated & differentiated samples:

samples <- subset(samples, treatment=="NT" & differentiated=="N")

nrow(samples)


# CHECK RESCOPES:
# Where there are at diagnosis & rescope versions of the same sample, drop the rescope(s)

rescopes <- subset(samples, sampling.time.point=="rescope")

subset(samples, sample.id.3 %in% rescopes$sample.id.3)[,c("case.no","sample.id.3","sampling.time.point","passage","rescope","RE.batch","notes")]

# T091 & T116 both have at diagnosis & rescope samples

dups <- subset(samples, sample.id %in% rescopes$sample.id)
dups <- subset(dups, duplicated(sample.id))

subset(samples, case.no=="T091" & (sample.id %in% dups$sample.id))[,c("case.no","sample.id","sampling.time.point","passage","rescope","RE.batch","notes","array.id")]
subset(samples, case.no=="T116" & (sample.id %in% dups$sample.id))[,c("case.no","sample.id","sampling.time.point","passage","rescope","RE.batch","notes","array.id")]

# Drop the T091 & T116 rescopes:

remove <- c("203548970061_R01C01","203548970061_R06C01")

samples <- subset(samples, !(array.id %in% remove))

nrow(samples)


# FILTER BY PASSAGE:
# Where there are duplicate samples from multiple passages, retain the sample with the lowest passage & remove the rest

ids <- unique(samples$sample.id)

samples <- lapply(1:length(ids), function(x){
    samp <- samples[which(samples$sample.id==ids[x]),]
    samp[which(samp$numeric.passage==min(samp$numeric.passage)),]
})

samples <- do.call(rbind, samples)

nrow(samples)


# FINAL DUPLICATE FILTERING:
# Preferentially remove the frozen biopsy duplicates
# Preferentially retain the newest duplicate (RE.batch: NEW > RE > OLD)

dups <- subset(samples, duplicated(sample.id))
dups <- subset(samples, sample.id %in% dups$sample.id)

nrow(dups)

dups[,c("sample.id","passage","frozen.biopsy","RE.batch","RE.sample.name","notes","array.id")]

samples <- subset(samples, !(RE.sample.name %in% c("T421 DUO P2 FB","T074 RE3 DUO P2 FB")))

nrow(samples)

# Remove "OLD" duplicates:

remove <- subset(dups, RE.batch=="OLD")$array.id
samples <- subset(samples, !(array.id %in% remove))

nrow(samples)

dups <- subset(samples, duplicated(sample.id))
dups <- subset(samples, sample.id %in% dups$sample.id)

dups[,c("sample.id","passage","frozen.biopsy","RE.batch","RE.sample.name","notes","array.id")]

# Remove "RE" duplicates:

remove <- subset(dups, RE.batch=="RE")$array.id
samples <- subset(samples, !(array.id %in% remove))

nrow(samples)

# Verify no duplicates remaining:

nrow(subset(samples, duplicated(sample.id)))


# REMOVE NEONATAL & FOETAL SAMPLES:

table(samples$age.group)

samples <- subset(samples, age.group=="paediatric")

nrow(samples)


# CHECK REMAINING RESCOPES:

subset(samples, sampling.time.point=="rescope")[,c("sample.id","sampling.time.point","passage","pdata.id","array.id")]

table(samples$tissue)


# SAVE UPDATED DATASET:

readme.epic <- data.frame(Dataset=c("epic.samples","all.beta","re.meta","all.ok.samples","samples"),Description=c("All EPIC array methylation samples (10Oct22)","All RE QCed EPIC methylation array beta values (04Oct22, RE name: combat_organoid_Beta)","RE's (out of date) meta data matching all.beta (RE name: epic.organoid_combined)","epic.samples subset to match all.beta","Final filtered sample info table (07Oct22)"))
save(epic.samples,all.beta,re.meta,all.ok.samples,samples,readme.epic,file="METHYLATION_ARRAY_DATA/AllEPICOrganoid_AnalysisReadyDataset_10Oct22.RData")


# SESSION INFO:

sessionInfo()

