## Please uncomment the first time you run this and re-install packages

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(MultiAmplicon)
library(ggplot2)
library(data.table)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE

doMultiAmp <- TRUE

doTax <- TRUE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run raccoon#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/2018_22_raccoon_complete/" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 


samples <- gsub("_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[215]])

#Creation of a folder for filtrated reads 

filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/filtered/"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(250,250), minLen=c(250,250), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
}


names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file 

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/primer_file_raccoon.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MAR <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/Stratified_files"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MAR <- sortAmplicons(MAR, n=1e+05, filedir=filedir)
  
  errF <-  learnErrors(unlist(getStratifiedFilesF(MAR)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MAR)), nbase=1e8,
                      verbose=0, multithread = 12)
  
  MAR <- derepMulti(MAR, mc.cores=12) 
  #MAR.1 <- dadaMulti(MAR[, which(grepl("^P1", colnames(MAR)))], Ferr=errF, Rerr=errR,  pool=FALSE,
    #               verbose=0, mc.cores=12)
  #MAR.2 <- dadaMulti(MAR[, which(grepl("^P2", colnames(MAR)))], Ferr=errF, Rerr=errR,  pool=FALSE,
    #                 verbose=0, mc.cores=12)
  #MAR.3 <- dadaMulti(MAR[, which(grepl("^P3", colnames(MAR)))], Ferr=errF, Rerr=errR,  pool=FALSE,
   #                  verbose=0, mc.cores=12)
  #MA <- concatenateMultiAmplicon(MAR.1, MAR.2, what = "samples") 
  #MA <- concatenateMultiAmplicon(MA, MAR.3, what = "samples")
  
  MAR <- dadaMulti(MAR, Ferr=errF, Rerr=errR,  pool=FALSE,
                     verbose=0, mc.cores=12)
  
  MAR <- mergeMulti(MAR, mc.cores=12) 
  
  propMerged <- MultiAmplicon::calcPropMerged(MAR)
  
  MAR <- mergeMulti(MAR, mc.cores=12) 
  
  MAR <- makeSequenceTableMulti(MAR, mc.cores=12) ## FIXME in package!!!
  
  MAR <- removeChimeraMulti(MAR, mc.cores=12)
  
  saveRDS(MAR, "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR_complete.RDS")
} else{
  MAR <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR_complete.RDS") ###START from here now! 
}

trackingF <- getPipelineSummary(MAR) 
## doesn't work for now

plotAmpliconNumbers(MAR)

## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()

###New taxonomic assignment

## When loading an old MA object that lacks sample data, simply:
MAR <- addSampleData(MAR)

if (newTax){
  MAR <- blastTaxAnnot(MAR,
                       infasta="/SAN/Zebra/all_in.fasta",
                       outblast="/SAN/Zebra/blast_out.fasta", 
                       taxonSQL="/SAN/db/taxonomy/taxonomizr.sql")
  saveRDS(MA7, file="/SAN/Zebra/MA7.Rds")
} else {
  MA7 <- readRDS(file="/SAN/Zebra/MA7.Rds")
}


#Sys.setenv("BLASTDB" = "/SAN/db/blastdb/") To make the annotation work, boss will fix this in the package
#library("vctrs", lib.loc="/usr/local/lib/R/site-library")
MAR <- blastTaxAnnot(MAR,  dataBaseDir = Sys.getenv("BLASTDB"), negative_gilist = "/SAN/db/blastdb/uncultured.gi", num_threads = 15)

###Extract sequences to do taxonomic assignment 

STNCR <- getSequenceTableNoChime(MAR)
#save(STNCR, file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/STNCR_table_combi.Rdata")

sequences <- unlist(lapply(STNCR, colnames))
names(sequences) <- paste0("asv_", 1:length(sequences))

if(doTax){
  library(taxonomizr)
  library(taxize)
  
  Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
                              "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/RaccRun_seq_combined.fasta")
  
  clusters <- plotAmpliconNumbers(MAR) 
  
  ###BLAST
  ## blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query /SAN/Victors_playground/Metabarcoding/AA_Raccoonn/RaccRun_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -out /SAN/Victors_playground/Metabarcoding/AA_Raccoon/asv_vs_nt_raccfinal.asn
  
  ## blast_formatter -archive /SAN/Victors_playground/Metabarcoding/AA_Racc/asv_vs_nt_raccfinal.asn -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" > /SAN/Victors_playground/Metabarcoding/AA_Raccoon/asv_nt_raccfinal.blttax
  
  ###Read blast result 
  ## we read that ouput into R blast <-
  blast <- read.csv("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/asv_vs_nt_RaccComb.blttax", header=FALSE)
  
  names(blast) <- c("query", "subject", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                    "bitscore", "staxid")
  blast <- as.data.table(blast)
  blast$staxid <- as.character(blast$staxid)
  
  read.nodes.sql("/SAN/db/taxonomy/nodes.dmp",
                 "/SAN/db/taxonomy/taxonomizr.sql")
  read.names.sql("/SAN/db/taxonomy/names.dmp",
                 "/SAN/db/taxonomy/taxonomizr.sql")
  
  blast.tax <- getTaxonomy(unique(blast$staxid),
                           "/SAN/db/taxonomy/taxonomizr.sql")
  
  blast.tax <- as.data.table(blast.tax, keep.rownames="staxid")
  blast.tax$staxid <- gsub("\\s*", "", blast.tax$staxid)
  
  blt <- merge(blast, blast.tax, by="staxid", all=TRUE)
  
  ## ## ## We need to be more clever if we want to use multiple
  ## ## ## hsps, this does not work for whole genome subjects eg.
  ### blt <- blt[,.(bitsum=sum(bitscore),
  ###              superkingdom, phylum, class, order, family, genus, species),
  ###           by=c("query", "subject")]
  
  ###    blt <- unique(blt)
  
  blt <- blt[,.(bitdiff= bitscore - max(bitscore),
                superkingdom, phylum, class, order, family, genus, species),
             by=c("query")]
  
  get.unique.or.na <- function (x){
    ## unique taxa at that level excluding potential NA's 
    ux <- unique(as.character(x[!is.na(x)]))
    ## but return NA if they are not unique
    if(length(ux)==1){return(ux)} else {as.character(NA)}
  }
  
  species <- blt[bitdiff>-2, .(species=get.unique.or.na(species)),
                 by=query]
  
  genus <- blt[bitdiff>-2, .(genus=get.unique.or.na(genus)),
               by=query]
  
  family <- blt[bitdiff>-7, .(family=get.unique.or.na(family)),
                by=query]
  
  order <- blt[bitdiff>-12, .(order=get.unique.or.na(order)),
               by=query]
  
  class <- blt[bitdiff>-20, .(class=get.unique.or.na(class)),
               by=query]
  
  phylum <- blt[bitdiff>-30, .(phylum=get.unique.or.na(phylum)),
                by=query]
  
  superkingdom <- blt[bitdiff>-50, .(superkingdom=get.unique.or.na(superkingdom)),
                      by=query]
  
  annot <- cbind(superkingdom[,c("query", "superkingdom")],
                 phylum[,"phylum"],
                 class[,"class"],
                 order[,"order"],
                 family[,"family"],
                 genus[,"genus"],
                 species[,"species"])
  
  seqnametab <- as.data.table(cbind(query=names(sequences), sequences))
  seqnametab <- merge(seqnametab, annot)
  
  dupseq <- seqnametab$sequences[duplicated(seqnametab$sequences)]
  
  seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]
  
  annot.list <- lapply(STNCR, function (x) {
    setkey(seqnametab, sequences)
    seqnametab[colnames(x),
               c("superkingdom", "phylum", "class", "order", "family", "genus", "species")]
  })
  
  saveRDS(annot.list, file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/Racc_blast_tax_all.Rds")
} else{
  annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/Racc_blast_tax_all.Rds")
}


## ## Not needed anymore
## keepr <- unlist(lapply(annot.list, nrow))>0
## annot.list <- annot.list[keep]
## STNC <- STNC[keep]


## name the annotation lists to have the names of the taxa 
annot.list <- lapply(seq_along(annot.list), function (i){
  an <- as.matrix(annot.list[[i]])
  rownames(an) <- colnames(STNCR[[i]])
  an
})

names(annot.list) <- names(STNCR)

phylalist <- lapply(annot.list, function (x) {
  if(nrow(x)>0){
    table(x[, "phylum"])
  }
})


tabulate.taxa <- function(taxtab, taxon, phylumsubset){
  if(nrow(taxtab)>0){
    t <- taxtab[taxtab[, "phylum"]%in%phylumsubset, ]
    if(!is.null(ncol(t))){
      table(t[, taxon])
    } else {NULL} 
  }else {NULL} 
}


## Tabulate by specific phylum
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Cestoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "species", "Nematoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "species", "Apicomplexa"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Ascomycota"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x, "species", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x, "phylum", "Ascomycota"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Amphibia"))


### Add sample information
library(phyloseq)
sample.dataR <- read.csv("~/AA_Raccoon/Raccoon_sample_data.csv")
rownames(sample.dataR) <- sample.dataR$plate_ID


PS.l <- lapply(seq_along(STNCR)[keepr], function(i){
  phyloseq(otu_table(STNCR[[i]], taxa_are_rows=FALSE),
           sample_data(sample.dataR[rownames(STNCR[[i]]), ]),
           tax_table(annot.list[[i]]))
})

#PS.l <- lapply(seq_along(STNCR)[keep], function(i){
#        phyloseq(otu_table(STNCR[[i]], taxa_are_rows=FALSE),
#           tax_table(annot.list[[i]]))
#}) ###Just used without sample information


#samples <- as.data.frame(list.files("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/2018_22_Raccoon/"))
#colnames(samples) <- "sample_id" 
#samples <- gsub(pattern = "_R\\d+_001.fastq.gz", "", samples$sample_id)
#samples <- as.data.frame(unique(samples))
#colnames(samples) <- "sample_id" 
#write.csv(samples, "/localstorage/victor/AA_Primer_evaluation/Sample_raccoon.csv")


#samples <- gsub(pattern = "S\\d+_\\d+_", "", samples$`unique(samples)`)

sumSeqByTax <- function (Phy, tax) {
  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
  counts$asvCount <- as.numeric(as.character(counts$asvCount))
  tapply(counts$asvCount, counts[, tax], sum)
}

readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
names(readNumByPhylum) <- names(STNCR)[keepr]


readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByGenus) <- names(STNCR)[keepr]


readNumByfamily <- lapply(PS.l, sumSeqByTax, "family") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByfamily) <- names(STNCR)[keepr]

readNumByspecies <- lapply(PS.l, sumSeqByTax, "species") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByspecies) <- names(STNCR)[keepr]

####
makeSequenceTableMulti(MAR)

##funtion fillSampleTables is not working

fill <- fillSampleTables(MAR)
MAR@sequenceTableFilled <- fill@sequenceTableFilled

MAR@sampleData

## Analyse all at once for now
ALL <- Reduce(cbind, fill@sequenceTableFilled[keepr])

## Problem: over all amplicons some ASVs are identical...
table(duplicated(colnames(ALL)))

## sum up same reads over amplicons
ALL.u <- do.call(rbind, by(t(ALL), rownames(t(ALL)), colSums))

## same for tax
all.tax <- Reduce(rbind, annot.list[rownames(MAR)[keep]])
all.tax <- all.tax[rownames(ALL.u), ]


PS <- phyloseq(otu_table(ALL.u, taxa_are_rows=TRUE),
               sample_data(sample.data[rownames(ALL), ]),
               tax_table(all.tax))

prune_both_zero <- function (ps) {
  p <- prune_samples(sample_sums(ps) > 0 , ps)
  prune_taxa(taxa_sums(p) > 0 , p)
}

PS <- prune_both_zero(PS)
PS.l <- lapply(PS.l, prune_both_zero)
names(PS.l) <- names(STNCR)[keepr]


################# ## HOW TO GO ON FROM HERE ## ######################
#### PS is now a single Phyloseq object over all amplicons. 

## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/PhyloSeqList_Raccall.Rds") ###For primer analysis (Victor)
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqCombi.Rds") ###For Fox analysis (Caro and Sophia)

