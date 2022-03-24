#!/usr/bin/env Rscript
rm(list=ls())

## extract the path to the script itself
myarg <- commandArgs(trailingOnly = FALSE)
scriptName <- gsub("--file=", "", grep("--file=", myarg, value = TRUE))
scriptDir <- dirname(scriptName)

## arguments from commandline
argPos <- grep("--args", myarg, fixed = TRUE)
sequenceFile <- as.character(myarg[argPos+1])
geneFile <- as.character(myarg[argPos+2])
repeatmaskerfile <- as.character(myarg[argPos+3])
outfileName <- as.character(myarg[argPos+4])

## load common functions
source(file.path(scriptDir, "commonFunctions.R"))

################################################################################################
### load data, first get all possible sequence IDs from the fasta
allPossibleIDs <- f.extraxt.fasta.IDs(sequenceFile) # see commonFunctions
allPossibleIDs <- paste0("chr", allPossibleIDs)
seqAnno <- f.load.seq.annotation(geneFile, repeatmaskerfile, allPossibleIDs) # see commonFunctions
write.csv(seqAnno, outfileName)

cat("Loaded", length(allPossibleIDs), "sequences.\n")
cat("Annotated as gene:", sum(seqAnno$gene == "yes"), "\n")
cat("Annotated as repeat:", sum(seqAnno$repMasVerySimple == "repeat"), "\n")
cat("Annotated as transposon:", sum(seqAnno$repMasVerySimple == "transposon"), "\n")
cat("Unannotated:", sum((seqAnno$gene == "no")&(seqAnno$repMasVerySimple=="nothing")), "\n")





