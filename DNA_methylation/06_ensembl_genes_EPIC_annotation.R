#'---
#'title: Add promoter intragenic and 3'UTR features to Ensembl transcript file generated from biomart
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---
#'

library(here)
options(stringsAsFactors = FALSE)


#' The emsembl human genes were downloaded using biomart
ensembl_transcripts<-read.csv(here("DNAm/data","ensembl_biomart_GRCh37.p13_human_genes.csv")) # 215,170 transcripts

#' remove transcripts on chromosome scaffolds
#' no EPIC CpGs on MT
ensembl_transcripts<-ensembl_transcripts[-grep("GL000|PATCH|HSCHR|HG|MT",ensembl_transcripts$Chromosome.scaffold.name),] 
table(ensembl_transcripts$Chromosome.scaffold.name) # 196,317 transcripts


## add coordinates for promoter, gene body, 3'UTR
ensembl_transcripts_features<-lapply(1:nrow(ensembl_transcripts), function(x){
  
  if(ensembl_transcripts$Strand[x]==1){
    
    data.frame(
      total_start = ensembl_transcripts$Transcription.start.site..TSS.[x]-1500,
      total_end = ensembl_transcripts$Transcript.end..bp.[x]+300,
      
      promoter_start = ensembl_transcripts$Transcription.start.site..TSS.[x]-1500,
      promoter_end = ensembl_transcripts$Transcription.start.site..TSS.[x]+300,
      
      intragenic_start = ensembl_transcripts$Transcription.start.site..TSS.[x]+300,
      intragenic_end = ensembl_transcripts$Transcript.end..bp.[x]-300,
      
      UTR3_start = ensembl_transcripts$Transcript.end..bp.[x]-300,
      UTR3_end = ensembl_transcripts$Transcript.end..bp.[x]+300)
    
  }else{ # here I am flipping start and end so that everything is 5'to 3' and the start of a feature is alway upstream of the end regradless of gene direction. 
    #This means the CpG to gene association can be starnd agnostic 
    
    data.frame(
      total_start = ensembl_transcripts$Transcript.start..bp.[x]-300,
      total_end = ensembl_transcripts$Transcription.start.site..TSS.[x]+1500,
      
      promoter_start = ensembl_transcripts$Transcription.start.site..TSS.[x]-300,
      promoter_end = ensembl_transcripts$Transcription.start.site..TSS.[x]+1500,
      
      intragenic_start = ensembl_transcripts$Transcript.start..bp.[x]+300,
      intragenic_end = ensembl_transcripts$Transcription.start.site..TSS.[x]-300,
      
      UTR3_start = ensembl_transcripts$Transcript.start..bp.[x]-300,
      UTR3_end = ensembl_transcripts$Transcript.start..bp.[x]+300)
    
  }
  
})

ensembl_transcripts_features<-do.call(rbind, ensembl_transcripts_features)
ensembl_transcripts<-cbind(ensembl_transcripts, ensembl_transcripts_features)
head(ensembl_transcripts)

save(ensembl_transcripts, file=paste(here("DNAm/data"),"/ensembl_transcripts.RData",sep=""))
