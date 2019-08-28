## Please uncomment the first time you run this and re-install packages

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(taxize)
library(parallel)

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
  filt.track <- lapply(seq_along(fastqF),  function (i) {
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
## When loading an old MA object that lacks sample data, simply:
MAR <- addSampleData(MAR)

###New taxonomic assignment

if (doTax){ ## simply save the blast files, that's even faster than
            ## setting doTax to FALSE and re-loading the object
    MAR2 <- blastTaxAnnot(MAR,  
                          negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                          db = "/SAN/db/blastdb/nt/nt",
                          infasta = "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/in.fasta",
                          outblast = "/SAN/Victors_playground/Metabarcoding/AA_Raccoon/out.blt",
                          num_threads = 22)
    saveRDS(MAR2, file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR2.Rds")
} else {
  MAR2 <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR2.Rds")
}

### Bugging probably an edge case with an empty amplicon
## trackingF <- getPipelineSummary(MAR2) 

plotAmpliconNumbers(MAR2, cluster_cols= F)

## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()

lapply(getTaxonTable(MAR2), function (x) table(as.vector(x[, "phylum"])))
lapply(getTaxonTable(MAR2), function (x) table(as.vector(x[, "genus"])))
lapply(getTaxonTable(MAR2), function (x) table(as.vector(x[, "species"])))

PS.l <- toPhyloseq(MAR2, samples=colnames(MAR2), multi2Single=FALSE)
PS <- toPhyloseq(MAR2, samples=colnames(MAR2), multi2Single=TRUE)


## ## just commented out everything that will work much more nicely
## ## now, but needs to be rewritten...

### ### 

## getNumReadsByTax <- function (Phy, tax){
##     sumSeqByTax <- function (Phy, tax) {
##         counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
##         counts$asvCount <- as.numeric(as.character(counts$asvCount))
##         tapply(counts$asvCount, counts[, tax], sum)
##     }
##     readNumByPhylum <- lapply(Phy, getNumReadsByTax, tax)
##     names(readNumByPhylum) <- rownames()

##     mergeDf <- function(x, y) {
##         m <- merge(x, y,  by=0, all=TRUE)
##         rownames(m) <- m$Row.names
##         m$Row.names <- NULL
##         m
##     }

##     foo <- Reduce(mergeDf, readNumByPhylum)
##     colnames(foo) <- rownames(MA8)
##     foo[is.na(foo)] <- 0
##     foo
## }


## ###Add sample data 
## sample.dataR <- read.csv("~/AA_Raccoon/Raccoon_sample_data.csv")
## rownames(sample.dataR) <- sample.dataR$plate_ID

## MARsample <- addSampleData(MAR2, sample.dataR)

## clust <- plotAmpliconNumbers(MARsample[, which(colnames(MARsample)%in%
##                                        sample.dataR$plate_ID)])

## clusters.row <- cutree(clust$tree_row, k=12) ##Eliminate the 5 primer pairs that didn't work
## clusters.col <- cutree(clust$tree_col, k=3) ##Eliminate negative controls and samples didn't work

## keep.prime <- names(clusters.row)[clusters.row!=5]
## keep.sample <- names(clusters.col)[clusters.col!=3]


## MAR3 <- MARsample[which(rownames(MARsample)%in%keep.prime),
##                 which(colnames(MARsample)%in%
##                         keep.sample &
##                         colnames(MARsample)%in% sample.dataR$plate_ID)]

## plotAmpliconNumbers(MAR3)


## ## name the annotation lists to have the names of the taxa 
## annot.list <- lapply(seq_along(annot.list), function (i){
##   an <- as.matrix(annot.list[[i]])
##   rownames(an) <- colnames(STNCR[[i]])
##   an
## })

## names(annot.list) <- names(STNCR)

## phylalist <- lapply(annot.list, function (x) {
##   if(nrow(x)>0){
##     table(x[, "phylum"])
##   }
## })


## tabulate.taxa <- function(taxtab, taxon, phylumsubset){
##   if(nrow(taxtab)>0){
##     t <- taxtab[taxtab[, "phylum"]%in%phylumsubset, ]
##     if(!is.null(ncol(t))){
##       table(t[, taxon])
##     } else {NULL} 
##   }else {NULL} 
## }


## ## Tabulate by specific phylum
## lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Cestoda"))
## lapply(annot.list, function (x) tabulate.taxa(x,  "species", "Nematoda"))
## lapply(annot.list, function (x) tabulate.taxa(x,  "species", "Apicomplexa"))
## lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
## lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))
## lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Ascomycota"))
## lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Chordata"))
## lapply(annot.list, function (x) tabulate.taxa(x, "species", "Chordata"))
## lapply(annot.list, function (x) tabulate.taxa(x, "phylum", "Ascomycota"))
## lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Amphibia"))


## ### Add sample information
## library(phyloseq)
## sample.dataR <- read.csv("~/AA_Raccoon/Raccoon_sample_data.csv")
## rownames(sample.dataR) <- sample.dataR$plate_ID


## PS.l <- lapply(seq_along(STNCR)[keepr], function(i){
##   phyloseq(otu_table(STNCR[[i]], taxa_are_rows=FALSE),
##            sample_data(sample.dataR[rownames(STNCR[[i]]), ]),
##            tax_table(annot.list[[i]]))
## })

## #PS.l <- lapply(seq_along(STNCR)[keep], function(i){
## #        phyloseq(otu_table(STNCR[[i]], taxa_are_rows=FALSE),
## #           tax_table(annot.list[[i]]))
## #}) ###Just used without sample information


## #samples <- as.data.frame(list.files("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/2018_22_Raccoon/"))
## #colnames(samples) <- "sample_id" 
## #samples <- gsub(pattern = "_R\\d+_001.fastq.gz", "", samples$sample_id)
## #samples <- as.data.frame(unique(samples))
## #colnames(samples) <- "sample_id" 
## #write.csv(samples, "/localstorage/victor/AA_Primer_evaluation/Sample_raccoon.csv")


## #samples <- gsub(pattern = "S\\d+_\\d+_", "", samples$`unique(samples)`)

## sumSeqByTax <- function (Phy, tax) {
##   counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
##   counts$asvCount <- as.numeric(as.character(counts$asvCount))
##   tapply(counts$asvCount, counts[, tax], sum)
## }

## readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
## names(readNumByPhylum) <- names(STNCR)[keepr]


## readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
## names(readNumByGenus) <- names(STNCR)[keepr]


## readNumByfamily <- lapply(PS.l, sumSeqByTax, "family") ## Change "text" in order to get counts per a different taxonomic level
## names(readNumByfamily) <- names(STNCR)[keepr]

## readNumByspecies <- lapply(PS.l, sumSeqByTax, "species") ## Change "text" in order to get counts per a different taxonomic level
## names(readNumByspecies) <- names(STNCR)[keepr]

## ####
## makeSequenceTableMulti(MAR)

## ##funtion fillSampleTables is not working

## fill <- fillSampleTables(MAR)
## MAR@sequenceTableFilled <- fill@sequenceTableFilled

## MAR@sampleData

## ## Analyse all at once for now
## ALL <- Reduce(cbind, fill@sequenceTableFilled[keepr])

## ## Problem: over all amplicons some ASVs are identical...
## table(duplicated(colnames(ALL)))

## ## sum up same reads over amplicons
## ALL.u <- do.call(rbind, by(t(ALL), rownames(t(ALL)), colSums))

## ## same for tax
## all.tax <- Reduce(rbind, annot.list[rownames(MAR)[keep]])
## all.tax <- all.tax[rownames(ALL.u), ]


## PS <- phyloseq(otu_table(ALL.u, taxa_are_rows=TRUE),
##                sample_data(sample.data[rownames(ALL), ]),
##                tax_table(all.tax))

## prune_both_zero <- function (ps) {
##   p <- prune_samples(sample_sums(ps) > 0 , ps)
##   prune_taxa(taxa_sums(p) > 0 , p)
## }

## PS <- prune_both_zero(PS)
## PS.l <- lapply(PS.l, prune_both_zero)
## names(PS.l) <- names(STNCR)[keepr]


## ################# ## HOW TO GO ON FROM HERE ## ######################
## #### PS is now a single Phyloseq object over all amplicons. 

## ## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


## saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/PhyloSeqList_Raccall.Rds") ###For primer analysis (Victor)
## saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqCombi.Rds") ###For Fox analysis (Caro and Sophia)

