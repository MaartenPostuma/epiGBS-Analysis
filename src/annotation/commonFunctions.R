#'@title a function for printing
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

#' a color gradient that I like
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.pinkredyellowblueblackNICE <- function (x) {
  r <- approx(c(0, 0.25, 0.5, 0.75, 1), c(154/255, 205/255, 238/255, 24/255, 0), n = x)$y
  g <- approx(c(0, 0.25, 0.5, 0.75, 1), c(50/255, 38/255, 201/255, 116/255, 0), n = x)$y
  b <- approx(c(0, 0.25, 0.5, 0.75, 1), c(205/255, 38/255, 0/255, 205/255, 0), n = x)$y
  return(rgb(r, g, b))
}

#' a color gradient that I like
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blackblueyellowredpinkNICE <- function(x) {
  return(rev(f.pinkredyellowblueblackNICE(x)))
}

#' add transparancy to an existing color
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.add.alpha <- function(col, alpha=1){
  if(missing(col)){stop("Please provide a vector of colours.")}
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

#'@title summarize a table according to groups and names
#'@param x <dataframe, matrix>: numeric, colnames as in byTab
#'@param byTab <dataframe>: at least two columns: sample and group
#'@param summaryFunction: a function like mean, median, sum
#'@return the summarized table (a matrix)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.summarize.columns <- function(x, byTab, summaryFunction) {
  byTab$sample <- as.character(byTab$sample)
  byTab$group <- as.character(byTab$group)
  groupnames <- unique(byTab$group)
  out <- matrix(0, nrow = nrow(x), ncol = length(groupnames), dimnames = list(rownames(x), groupnames))
  for (gn in groupnames) { 
    toSummarize <- byTab$sample[byTab$group == gn]
    if (length(toSummarize) > 1) {out[,gn] <- apply(x[,toSummarize], 1, summaryFunction)}
    else { out[,gn] <- x[,toSummarize] }
  }
  return(out)
}

#'@title simplify detailed repeatmasker annotations
#'@param seqAnno annotation data (see \code{\link{f.load.seq.annotation}})
#'@return the simplified annotation data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.simplify.repeat.masker <- function(seqAnno) {
  featureSimplification <- list(
    DNA_Ginger = c("DNA_Ginger"),
    DNA_hAT = c("DNA_hAT", "DNA_hAT_Ac", "DNA_hAT_Charlie", "DNA_hAT_Tag1", "DNA_hAT_Tip100"),
    DNA_TcMar = c("DNA_TcMar_Pogo", "DNA_TcMar_Mariner", "DNA_TcMar_Tc1", "DNA_TcMar_Stowaway"),
    DNA_MuLE = c("DNA_MULE_MuDR", "DNA_MULE"),
    DNA_PIF_Harbinger = c("DNA_PIF_Harbinger"),
    DNA_CMC = c("DNA_CMC_EnSpm"),
    DNA_Sola = c("DNA_Sola_1"),
    DNA = c("DNA"),
    LTR_Cassandra = c("LTR_Cassandra"),
    LTR_Caulimovirus = c("LTR_Caulimovirus"),
    LTR_Copia = c("LTR_Copia"),
    LTR_Gypsy = c("LTR_Gypsy"),
    LTR = c("LTR"),
    RC_Helitron = c("RC_Helitron"),
    LINE_Penelope = c("LINE_Penelope"),
    LINE_CRE = c("LINE_CRE"),
    LINE_RTE = c("LINE_RTE_BovB"),
    LINE_L1 = c("LINE_L1", "LINE_L1_Tx1"),
    LINE_I = c("LINE_I"),
    SINE_tRNA = c("SINE_tRNA", "SINE_tRNA_RTE"),
    SINE = c("SINE"),
    retroposon = c("Retroposon"),
    lowComplexity = c("Low_complexity"),
    uknRepeats = c("unClassRep", "Unknown"),
    simple_repeat = c("Simple_repeat"),
    satellite_centromeric = c("Satellite_centr"),
    satellite = c("Satellite"),
    tRNA = c("tRNA"),
    rRNA = c("rRNA"),
    snRNA = c("snRNA"),
    others = c("Other_Composite", "ARTEFACT", "nothing", "unAnnotated")
  )
  extendedToSimple <- list()
  for (simple in names(featureSimplification)) {
    for (entry in featureSimplification[[simple]]){
      extendedToSimple[[entry]] <- simple
    }
  }
  extendedToSimple <- unlist(extendedToSimple)
  # check what's missing (should be added manually)
  missing <- unique(seqAnno$repMas[!(seqAnno$repMas %in% names(extendedToSimple))])
  if (length(missing) > 0) {
    f.print.message("Missing repeatmasker features (add in f.simplify.repeat.masker):")
    cat(paste0(missing, collapse = '\n'), '\n')
  }
  seqAnno$repMas <- extendedToSimple[seqAnno$repMas]
  return(seqAnno)
}

#'@title read gene, transposon, and repeat annotations
#'@param geneFile path to the file with the genes
#'@param repMasFile path to the file with the repeats/transposons
#'@param allSeqs (optional) a vector with all sequence numbers
#'@param simplifyRepMasker attempt to simplify the detailed repeatmasker annotation (see \code{\link{f.simplify.repeat.masker}})
#'@return a merged annotation
#'@note 
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.seq.annotation <- function(geneFile, repMasFile, allSeqs = c(), simplifyRepMasker = FALSE) {
  if (grepl("\\.gz$", geneFile)) {
    geneAnno <- read.table(gzfile(geneFile), sep = '\t', stringsAsFactors = FALSE, header = FALSE, col.names = c("seqID", "refID", "score"))
  } else {
    geneAnno <- read.table(geneFile, sep = '\t', stringsAsFactors = FALSE, header = FALSE, col.names = c("seqID", "refID", "score"))
  }
  if (grepl("\\.gz$", repMasFile)) {
    repMasAnno <- read.table(gzfile(repMasFile), sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  } else {
    repMasAnno <- read.table(repMasFile, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  }
  geneAnno$seqID <- paste0("chr", geneAnno$seqID)
  repMasAnno$queryID <- paste0("chr", repMasAnno$queryID)
  repMasAnno$repeatClassOrFamilySimple <- sapply(repMasAnno$repeatClassOrFamily, function(x) unlist(strsplit(x, "\\?|\\/"))[1])
  repMasAnno$teOrRep <- "transposon"
  repMasAnno$teOrRep[repMasAnno$repeatClassOrFamilySimple %in% c("Low_complexity", "Other", "rRNA", "Satellite", "Simple_repeat", "snRNA", "tRNA", "Unknown")] <- "repeat"
  if (length(allSeqs) == 0) { allSeqs <- union(geneAnno$seqID, repMasAnno$queryID) }
  out <- data.frame(gene = rep("no", length(allSeqs)),
                    repMas = rep("nothing", length(allSeqs)),
                    repMasSimple = rep("nothing", length(allSeqs)),
                    repMasVerySimple = rep("nothing", length(allSeqs)),
                    geneScore = rep(0, length(allSeqs)),
                    repMasScore = rep(0, length(allSeqs)),
                    stringsAsFactors = FALSE, row.names = allSeqs)
  out[geneAnno$seqID, "gene"] <- "yes"
  out[repMasAnno$queryID, "repMas"] <- repMasAnno$repeatClassOrFamily
  out[repMasAnno$queryID, "repMasSimple"] <- repMasAnno$repeatClassOrFamilySimple
  out[repMasAnno$queryID, "repMasVerySimple"] <- repMasAnno$teOrRep
  out[geneAnno$seqID, "geneScore"] <- geneAnno$score
  out[repMasAnno$queryID, "repMasScore"] <- repMasAnno$swScore
  # replace some odd characters
  out$repMas <- gsub("-", "_", out$repMas, fixed = TRUE)
  out$repMas <- gsub("/", "_", out$repMas, fixed = TRUE)
  out$repMas <- gsub("?", "", out$repMas, fixed = TRUE)
  if (simplifyRepMasker) {
    out <- f.simplify.repeat.masker(out)
  }
  return(out)
}

#'@title extract all IDs from a fasta file
#'@param geneFile path to the fasta file (can be .gz)
#'@return vector with IDs
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.extraxt.fasta.IDs <- function(infileName) {
  if (grepl("\\.gz$", infileName)) {
    infile <- gzcon(file(infileName, "rb"))
  }
  else {
    infile <- file(infileName, open = "rt")
  }
  out <- rep(NA, 1e6)
  counter <- 0
  while (length(oneLine <- readLines(infile, n = 1, warn = FALSE)) > 0) {
    if (substr(oneLine, 1, 1) != ">") { next }
    counter <- counter + 1
    seqID <- gsub(">|[[:space:]]", "", oneLine)
    if (counter > length(out)) {
      out <- c(out, rep(NA, 1e6))
    }
    out[counter] <- seqID
  } 
  close(infile)
  out <- out[!is.na(out)]
  return(out)
}

#'@title read the design table and check if the sample names are ok
#'@param designTablePath path to the design table file
#'@param colNamesForGrouping a vector with the column names used for grouping, if given, will check and add a column "group" with all of them pasted together
#'@return sampleTab
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.sampleTable <- function(designTablePath, colNamesForGrouping = c()) {
  sampleTab <- read.table(designTablePath, sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  # avoid using anything else than underscores, letters, and numbers in sample names
  rownames(sampleTab) <- gsub("-|\\.", "_", rownames(sampleTab))
  # and don't start sample names with a number
  if (sum(grepl("^[[:digit:]]", rownames(sampleTab))) > 0) {
    rownames(sampleTab) <- paste0("sample_", rownames(sampleTab))
  }
  # check if some samples need to be removed
  if ("sampleRemovalInfo" %in% colnames(sampleTab)) {
    f.print.message("Removing", sum(sampleTab$sampleRemovalInfo == "ERASE"), "samples due to the sampleRemovalInfo column")
    sampleTab <- subset(sampleTab, sampleRemovalInfo != "ERASE")
  }
  # check for the grouping variables
  if (length(colNamesForGrouping) > 0) {
    if (sum(colNamesForGrouping %in% colnames(sampleTab)) != length(colNamesForGrouping)) {
      f.print.message("The column names for grouping don't exist!")
      cat(paste0(setdiff(colNamesForGrouping, colnames(sampleTab)), collapse = '\n'), '\n')
      f.print.message("Columns in the design table:")
      cat(paste0(colnames(designTable), collapse = '\n'), '\n')
      stop("STOPPING!")
    }
    if (length(colNamesForGrouping) > 1) {
      sampleTab$group <- apply(sampleTab[,colNamesForGrouping], 1, function(x) paste(x, collapse = '_'))
    } else {
      sampleTab$group <- sampleTab[[colNamesForGrouping]]
    }
  }
  return(sampleTab)
}

#'@title read a vcf file and correct names
#'@param infileName path to the vcf file
#'@return a list with two entries, genoTypes and snpCoverage
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.vcf.file <- function(infileName) {
  myData <- read.vcfR(infileName)
  genoTypes <- extract.gt(myData, element = "GT", as.numeric = FALSE)
  snpCoverage <- extract.gt(myData, element = "DP", as.numeric = TRUE)
  # correct Fallopia
  colnames(genoTypes) <- gsub("^FALLOPIA_", "", colnames(genoTypes))
  colnames(snpCoverage) <- gsub("^FALLOPIA_", "", colnames(snpCoverage))
  # avoid using anything else than underscores, letters, and numbers in sample names
  colnames(genoTypes) <- gsub("-|\\.", "_", colnames(genoTypes)) 
  colnames(snpCoverage) <- gsub("-|\\.", "_", colnames(snpCoverage))
  # and don't start sample names with a number
  if (sum(grepl("^[[:digit:]]", colnames(genoTypes))) > 0) {
    colnames(genoTypes) <- paste0("sample_", colnames(genoTypes))
    colnames(snpCoverage) <- paste0("sample_", colnames(snpCoverage))
  }
  out <- list(
    genoTypes = genoTypes,
    snpCoverage = snpCoverage
  )
  return(out)
}
 
#'@title read a methylation.bed(.gz)
#'@param infileName path to the bed file
#'@param removeMulticontextPositions if TRUE, positions with more than one context are removed
#'@param percentages if TRUE, the methylated and total columns are converted into percentages
#'@return table with the methylation data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.methylation.bed <- function(infileName, removeMulticontextPositions = TRUE, percentages = FALSE) {
  require("data.table")
  if (grepl("\\.gz$", infileName)) {
    myData <- fread(cmd=paste("zcat", infileName), sep = '\t', data.table = FALSE, fill = TRUE)
  } else {
    myData <- fread(infileName, sep = '\t', data.table = FALSE, fill = TRUE)
  }
  # correct Fallopia
  colnames(myData) <- gsub("^FALLOPIA_", "", colnames(myData))
  # columns with sampleInfo
  sampleInfoColNumbers <- grep("_total$|_methylated$", colnames(myData), value = FALSE)
  # avoid using anything else than underscores, letters, and numbers in sample names
  colnames(myData)[sampleInfoColNumbers] <- gsub("-|\\.", "_", colnames(myData)[sampleInfoColNumbers]) 
  # and don't start sample names with a number
  if (sum(grepl("^[[:digit:]]", colnames(myData)[sampleInfoColNumbers])) > 0) {
    colnames(myData)[sampleInfoColNumbers] <- paste0("sample_", colnames(myData)[sampleInfoColNumbers])
  }
  # remove positions with multiple contexts if requested
  if (removeMulticontextPositions) {
    temp <- unique(myData[,c("chr", "pos", "context")])
    ctxtCount <- table(paste0(temp$chr, '_', temp$pos))
    multipleContexts <- names(ctxtCount)[ctxtCount>1]
    if (length(multipleContexts) > 0) {
      f.print.message("removing", length(multipleContexts), "positions with multiple contexts.")
      myData <- subset(myData, !(paste0(chr, '_', pos) %in% multipleContexts))
    }
  }
  # convert to percentages if requested
  if (!percentages) {
    out <- myData
  } else {
    out <- myData[,c("chr", "pos", "context")]
    methCols <- grep("_methylated$", colnames(myData), value = TRUE)
    totCols <- grep("_total$", colnames(myData), value = TRUE)
    mePerc <- round(myData[,methCols]/myData[,totCols]*100, 3)
    colnames(mePerc) <- gsub("_methylated$", "", colnames(mePerc))
    for (toAdd in colnames(mePerc)) {
      out[[toAdd]] <- mePerc[[toAdd]]
    }
  }
  return(out)
}

#'@title read a filtered SNP file
#'@param infileName path to the filtered snp file
#'@return table with the snp data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.filtered.SNPs <- function(infileName) {
  if (grepl("\\.gz$", infileName)) {
    myData <- read.table(gzfile(infileName), sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  } else {
    myData <- read.table(infileName, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  }
  # correct Fallopia
  colnames(myData) <- gsub("^FALLOPIA_", "", colnames(myData))
  # avoid using anything else than underscores, letters, and numbers in sample names
  colnames(myData) <- gsub("-|\\.", "_", colnames(myData)) 
  # and don't start sample names with a number
  if (sum(grepl("^[[:digit:]]", colnames(myData))) > 0) {
    colnames(myData) <- paste0("sample_", colnames(myData))
  }
  return(myData)
}

#'@title read a merged annotation file and subset for a feature if requested
#'@param infileName path to the merged annotation (see the script mergeAnnotation.R)
#'@param selectedFeature all, gene, transposon, repeat, nothing 
#'@param simplifyRepMasker attempt to simplify the detailed repeatmasker annotation (see \code{\link{f.simplify.repeat.masker}})
#'@return table with the snp data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.merged.annotation <- function(infileName, selectedFeature = "all", simplifyRepMasker = FALSE) {
  mergedAnno <- read.csv(infileName, row.names = 1, stringsAsFactors = FALSE, header = TRUE)
  # subset feature
  if (selectedFeature != "all") {
    if (selectedFeature == "gene") {
      out <- subset(mergedAnno, gene == "yes")
    } else if (selectedFeature == "repeat") {
      out <- subset(mergedAnno, repMasVerySimple == "repeat")
    } else if (selectedFeature == "transposon") {
      out <- subset(mergedAnno, repMasVerySimple == "transposon")
    } else if (selectedFeature == "nothing") {
      out <- subset(mergedAnno, (repMasVerySimple == "nothing") & (gene == "no"))
    } else {
      cat("ERROR, feature", feature, "does not exist.\n")
      quit("no", 1)
    }
  } else {
    out <- mergedAnno
  }
  if (simplifyRepMasker) {
    out <- f.simplify.repeat.masker(out)
  }
  return(out)
}

#'@title read and combine snp and methylation data
#'@param filteredSnpFile path to the filtered snp file (.txt(.gz))
#'@param filteredMethFile path to the filtered methylation file (.bed(.gz))
#'@param replaceNASNPsWithRefVar if TRUE, NA SNPs (..) are replaced with the reference variant (00). Otherwise, the NA is carried on (contig will be NA for this individual)
#'@return table with the snp data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.and.combine.snp.and.methylation.data <- function(filteredSnpFile, filteredMethFile, replaceNASNPsWithRefVar = FALSE) {
  #filteredSnpFile <- "/media/mwschmid/myData/MWSchmid/ChristinaRichards_FallopiaEAGER/GitIgnore_filteredData/filteredSNPs.txt.gz"
  #filteredMethFile <- "/media/mwschmid/myData/MWSchmid/ChristinaRichards_FallopiaEAGER/GitIgnore_filteredData/filteredMETH.bed.gz"
  snpDat <- f.load.filtered.SNPs(filteredSnpFile)
  meDat <- f.load.methylation.bed(filteredMethFile, removeMulticontextPositions = TRUE, percentages = TRUE)
  #################################################################################################### SNP data
  # paste the SNPs into a single-string genotype per contig
  # not sure how to deal with the NA, I replaced them with the reference: reason, if there is another SNP, then it will still get a unique genotype if false (the factor still works). 
  # If the NA position is the only thing that's different, I think that it's bad treating is as it's own variant. In case the NAs have a technical reason, that might also influence the methylation.
  # With Fallopia I saw that this increased the number of "zero variation" SNPs from 734 to 35119
  # one could also just set the entire contig as NA for individuals with NA SNPs... (that happens in case of the SNPs)
  # there are both versions available here (replaceNASNPsWithRefVar)
  sampleColumnsSnps <- colnames(snpDat)
  if (replaceNASNPsWithRefVar) {
    f.print.message("Replacing", sum((snpDat == "..")|(snpDat == "./.")), "NA SNPs with the reference variant (see comment in commonFunction.R, f.read.and.combine.snp.and.methylation.data()).")
    snpDat[snpDat == ".."] <- "00"
    snpDat[snpDat == "./."] <- "00"
    # remove all lines that are identical in all samples
    varCount <- apply(snpDat, 1, function(x) length(unique(x)))
    f.print.message("Removing", sum(varCount==1), "SNPs without variation.")
    snpDat <- snpDat[varCount > 1,]
    # paste the SNPs into a single-string genotype per contig
    snpRefs <- sapply(rownames(snpDat), function(x) unlist(strsplit(x, '_'))[1])
    snpPerRef <- aggregate(snpDat, by = list(chr = snpRefs), function(x) paste(x, collapse = '|'))
  } else {
    # use a special paste function
    f.paste.genotype.if.not.NA <- function(x) {
      if (sum(is.na(x)) > 0) {
        out <- NA
      } else {
        out <- paste0(x, collapse = "|")
      }
      return(out)
    }
    # paste the SNPs into a single-string genotype per contig
    snpDat[snpDat == ".."] <- NA
    snpDat[snpDat == "./."] <- NA
    snpRefs <- sapply(rownames(snpDat), function(x) unlist(strsplit(x, '_'))[1])
    snpPerRef <- aggregate(snpDat, by = list(chr = snpRefs), f.paste.genotype.if.not.NA)
    # now replace the ones with all equal (except for SNPs)
    varCount <- apply(snpPerRef[,sampleColumnsSnps], 1, function(x) length(unique(x[!is.na(x)])))
    f.print.message("Removing", sum(varCount==1), "contigs without variation.")
    snpPerRef <- snpPerRef[varCount > 1,]
  }
  rownames(snpPerRef) <- paste0("chr", snpPerRef$chr)
  f.print.message("Created", nrow(snpPerRef), "SNP contigs.")
  #################################################################################################### Methylation data
  # average methylation data per contig and context
  sampleColumnsMethylation <- setdiff(colnames(meDat), c("chr", "pos", "context"))
  mePerRef <- aggregate(meDat[,sampleColumnsMethylation], by = list(chr = meDat$chr, context = meDat$context), function(x) mean(x, na.rm = TRUE))
  # there are a few "." contigs. Remove them
  f.print.message("Removing", sum(mePerRef$context == "."), "contig/context combinations with a context == '.' (a dot).")
  mePerRef <- subset(mePerRef, context != ".")
  mePerRef$chr <- paste0("chr", mePerRef$chr)
  f.print.message("Created", length(unique(mePerRef$chr)), "methylation contigs.")
  #################################################################################################### 
  # intersect them
  commonContigs <- intersect(rownames(snpPerRef), mePerRef$chr)
  commonSamples <- intersect(sampleColumnsMethylation, sampleColumnsSnps)
  f.print.message("Common contigs:", length(commonContigs), "and common samples:", length(commonSamples))
  mePerRef <- subset(mePerRef, chr %in% commonContigs)
  availableContexts <- split(mePerRef$context, mePerRef$chr)
  temp <- split(mePerRef[,c("chr", commonSamples)], mePerRef$context)
  for (cc in names(temp)) {
    rownames(temp[[cc]]) <- temp[[cc]]$chr
    temp[[cc]] <- temp[[cc]][,commonSamples]
  }
  out <- list(
    snpPerRef = snpPerRef[commonContigs,commonSamples],
    mePerRef = temp,
    availableContexts = availableContexts
  )
  #out <- list(
  #  snpPerRef = snpPerRef[commonContigs,commonSamples],
  #  mePerRef = mePerRef[,commonSamples],
  #  meRef =  mePerRef$chr,
  #  meContext = mePerRef$context
  #)
  return(out)
}

#'@title remove col/rows with NA from a matrix
#'@return a distance matrix without NA
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.remove.NA.from.distance.matrix <- function(x) {
  toRemove <- which(is.na(x), arr.ind = TRUE)
  toKeep <- setdiff(rownames(x), rownames(toRemove))
  x <- x[toKeep,toKeep]
  return(x)
}

#'@title calculate distance between two genotypes (at a single SNP)
#'@param stringA genotype of sample A (e.g. 0/1)
#'@param stringB genotype of sample B (e.g. 0/0)
#'@param stringSeparator the separator (e.g. / is default for VCF, nothing is default for the filtered tables)
#'@return 0 if both all alleles are the same, 0.5 if one or more but not all the same, 1 if none the same
#'@notes This function is slow, if you want to process many SNPs, create first a lookup table (see filterSNPs.R)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.genotype.distance <- function(stringA, stringB, stringSeparator = '/') {
  if (is.numeric(stringA) | is.numeric(stringB)) { 
    f.print.message("ERROR in f.get.genotype.distance: numeric genotype.")
    return(NA)
  }
  partA <- unlist(strsplit(stringA, stringSeparator, TRUE))
  partB <- unlist(strsplit(stringB, stringSeparator, TRUE))
  if ((length(partA) != 2) | (length(partB) != 2)) {
    f.print.message("ERROR in f.get.genotype.distance: really not diploid data?.")
    print(stringA)
    print(partA)
    print(stringB)
    print(partB)
    return(NA)
  }
  numSame <- length(intersect(partA, partB))
  numTotal <- length(union(partA, partB))
  numEqualToTot <- numSame/numTotal
  if (numEqualToTot == 1) {
    out <- 0 # change from similarity to 
  } else if (numEqualToTot == 0) {
    out <- 1
  } else {
    out <- 0.5
  }
  return(out)
}

#'@title read a bunch of distance files and check if the names are alright
#'@param allDistanceFiles a vector with path to distance files
#'@param imputeNA TRUE/FALSE if NAs should be imputed
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.many.distance.files <- function(allDistanceFiles, imputeNA = TRUE) {
  require("ape")
  allDistances <- list()
  for (curFile in allDistanceFiles) {
    curFileName <- basename(curFile)
    curName <- gsub("\\.csv$", "", curFileName)
    temp <- read.csv(curFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    rownames(temp) <- gsub("-|\\.", "_", rownames(temp)) # avoid using anything else than underscores, letters, and numbers in sample names
    colnames(temp) <- gsub("-|\\.", "_", colnames(temp)) # avoid using anything else than underscores, letters, and numbers in sample names
    if (sum(grepl("^[[:digit:]]", colnames(temp))) > 0) {# and don't start with a number
      rownames(temp) <- paste0("sample_", rownames(temp))# and don't start with a number
      colnames(temp) <- paste0("sample_", colnames(temp))# and don't start with a number
    }# and don't start with a number
    temp <- as.matrix(temp)
    diag(temp) <- 0
    temp[lower.tri(temp)] <- t(temp)[lower.tri(t(temp))]
    if (imputeNA) {
      if (mean(is.na(temp)) > 0.2) {
        cat("Skipping data set because more than 20% are NA.\n")
        next
      }
      if (sum(is.na(temp)) > 0) {
        prevNames <- dimnames(temp)
        temp <- additive(temp)
        dimnames(temp) <- prevNames
      }
      if (sum(temp < 0) > 0) {
        cat("Detected negative distances. Removing some samples, trying again\n")
        percNeg <- rowMeans(temp<0)
        toKeep <- setdiff(rownames(temp), names(percNeg)[percNeg>0.05]) # remove the samples with more than 5 % negative distances, then set neg distances to NA and remove them
        temp <- temp[toKeep, toKeep]
        temp[temp<0] <- NA
        temp <- f.remove.NA.from.distance.matrix(temp)
        if (nrow(temp) < 4) {
          cat("Skipping table with less than 4 samples after NA removal.\n")
          next
        }
      }
    } else {
      temp <- f.remove.NA.from.distance.matrix(temp)
      if (nrow(temp) < 4) {
        cat("Skipping table with less than 4 samples after NA removal.\n")
        next
      }
    }
    allDistances[[curName]] <- temp
  }
  return(allDistances)
}

#'@title plot an image without text
#'@param x x-coords
#'@param y y-coords
#'@param z the image value
#'@param xLabel label for x-axis
#'@param yLabel label for y-axis
#'@param mainLabel a title
#'@param useLog TRUE/FALSE if data should be log2(x+1) transformed
#'@param ... other parameters forwarded to \code{\link{image}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.image.without.text <- function(x, y, z, xLabel, yLabel, mainLabel, useLog = FALSE, ...) {
  xChars <- as.character(x)
  yChars <- as.character(y)
  xPos <- 1:length(xChars)
  yPos <- 1:length(yChars)
  if (useLog) {
    toPlot <- log2(z+1)
  } else {
    toPlot <- z
  }
  image(xPos, yPos, toPlot, xlab = xLabel, ylab = yLabel, main = mainLabel, yaxt = "n", xaxt = "n", ...)
  axis(1, at = xPos, labels = xChars, outer = FALSE, las = 2)
  axis(2, at = yPos, labels = yChars, outer = FALSE, las = 1)
  return(NULL)
}

#'@title get the allele frequency of the second most abundant allele
#'@param snpDataRow a row of the snpData, see format for adegenet (rows: SNPs, cols: samples, entries: 00 or 11 or 01, etc)
#'@return a vector with the allele frequencies as percentage
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.allele.frequency.second.most.abundant <- function(snpDataRow) {
  snpDataRow[snpDataRow == ".."] <- NA
  snpDataRow[snpDataRow == "./."] <- NA
  temp <- snpDataRow[!is.na(snpDataRow)]
  if (sum(grepl("/", snpDataRow)) == length(snpDataRow)) {
    separator <- "/"
  } else {
    separator <- ""
  }
  temp <- unlist(strsplit(unlist(temp), separator))
  alleleCounts <- sort(table(temp), decreasing = TRUE)
  if (length(alleleCounts) < 2) {
    out <- 0
  } else {
    out <- round(alleleCounts[2]/sum(alleleCounts)*100, 2)
  }
  return(out)
}

#'@title Formulate a simple contrast using two sets of regular expressions.
#'@param plusTermsRE a vector of regular expressions or strings for the terms which will be summed up.
#'@param minusTermsRE a vector of regular expressions or strings for the terms which will be substracted.
#'@param design a design matrix obtained with \code{model.matrix}
#'@param invertPlus see note below
#'@param invertMinus see note below
#'@param fixed set to TRUE to disable regular expression mode of \code{\link{grep}}
#'@return a string containing the contrast definition (to use in \code{\link{makeContrast}})
#'@note The function will use \code{\link{grep}} to retrieve all \code{colnames(design)} 
#'which fit to either the plusTermsRE or the minusTermsRE. The contrast is then formulated as
#'plusTerms/length(plusTerms) - minusTerms/length(minusTerms). If invertPlus/invertMinus are set
#'to TRUE, the terms which do NOT match the search criteria are extracted (e.g. if one would like
#'to test a certain set A versus all non-A, one can do that with:
#'\code{f.formulate.simple.contrast(setA, setA, design, invertMinus = TRUE)})
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.formulate.simple.contrast <- function(plusTermsRE, minusTermsRE, design, invertPlus = FALSE, invertMinus = FALSE, fixed = FALSE) {
  plusTerms <- sapply(plusTermsRE, function(x) grep(x, colnames(design), fixed = fixed, value = TRUE, invert = invertPlus))
  minusTerms <- sapply(minusTermsRE, function(x) grep(x, colnames(design), fixed = fixed, value = TRUE, invert = invertMinus))
  if (length(plusTermsRE) > 1) {
    if (invertPlus) {
      plusTerms <- Reduce(intersect, plusTerms)
    } else { 
      plusTerms <- Reduce(union, plusTerms)
    }
  } else {
    plusTerms <- unique(unlist(plusTerms))
  }
  if (length(minusTermsRE) > 1) {
    if (invertMinus) {
      minusTerms <- Reduce(intersect, minusTerms)
    } else { 
      minusTerms <- Reduce(union, minusTerms)
    }
  } else { 
    minusTerms <- unique(unlist(minusTerms))
  }
  out <- paste0(
    "(", paste(plusTerms, collapse = "+"), ")/", length(plusTerms), "-",
    "(", paste(minusTerms, collapse = "+"), ")/", length(minusTerms)
  )
  return(out)
}
################################################################################################
################################################################################################
################################################################################################
### functions to check terms and remove the ones with too few levels (lm throws an error otherwise)
#'@title check how many levels a term has
#'@param dataToCheck a vector with data
#'@return number of levels (unique entries) in case of factor and character, 1 in case of numeric
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.number.of.term.levels <- function(dataToCheck) {
  if (is.numeric(dataToCheck)) {
    out <- 1
  } else {
    numLevels <- length(unique(dataToCheck))
    out <- ifelse(numLevels > 1, numLevels, NA)
  }
  return(out)
}

#'@title check whether a single term in the model has enough levels
#'@param toCheck the term (either a single column name or an interaction separated by :)
#'@param dataForTest the data frame that will be used for the LM
#'@return NA if not enough levels (numeric: at least two unique values, others: at least two levels, interactions: TODO
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.check.one.term <- function(toCheck, dataForTest) {
  out <- NA
  if (grepl(":", toCheck)) {
    individualTerms <- unlist(strsplit(toCheck, ":", fixed = TRUE))
    individualTerms <- sapply(individualTerms, function(x) f.check.one.term(x, dataForTest))
    individualTerms <- individualTerms[!is.na(individualTerms)]
    if (length(individualTerms) > 1) {
      # there should be a minimal number of values occuring for each combination
      levelsPerTerm <- sapply(individualTerms, function(x) f.get.number.of.term.levels(dataForTest[[x]]))
      if (sum(levelsPerTerm == 1) > 0) { # means that at least one term is numeric
        f.print.message("WARNING: interactions with numeric variables assume that most of the numeric values should be non-identical.")
        averageOccurences <- mean(table(dataForTest[,individualTerms])) # this should not be too high because that means that many values are identical
        if (averageOccurences < 2) {
          out <- paste0(individualTerms, collapse = ":")
        }
      } else {
        averageOccurences <- mean(table(dataForTest[,individualTerms])) # combinations should occur at least twice on average, otherwise it's too much 1-against-1
        if (averageOccurences >= 2) {
          out <- paste0(individualTerms, collapse = ":")
        }
      }
      #numLevelsCombined <- length(unique(apply(dataForTest[, individualTerms], 1, function(x) paste0(x, collapse = "|"))))
      #if (numLevelsCombined > max(levelsPerTerm)) { # there must be at least one combination more - note that this does not mean that it's balanced
      #  out <- paste0(individualTerms, collapse = ":")
      #}
    }
  } else {
    numLevels <- f.get.number.of.term.levels(dataForTest[[toCheck]])
    out <- ifelse(!is.na(numLevels), toCheck, NA)
  }
  return(out)
}

#'@title check whether all terms in a model can be fit (works only with simple formulas, only + and : operators)
#'@param formulaString the original formula as a string
#'@param dataForTest the data frame that will be used for the LM
#'@return the formula with terms removed that would make lm fail
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.check.and.clean.formula <- function(formulaString, dataForTest) {
  formulaString <- gsub(" ", "", formulaString, fixed = TRUE) # remove whitespaces
  formulaParts <- unlist(strsplit(formulaString, "~", fixed = TRUE))
  allTerms <- unlist(strsplit(formulaParts[2], "+", fixed = TRUE))
  termsToKeep <- sapply(allTerms, function(x) f.check.one.term(x, dataForTest))
  termsToKeep <- termsToKeep[!is.na(termsToKeep)]
  out <- paste0(formulaParts[1], "~", paste0(termsToKeep, collapse = "+"))
  return(out)
}
################################################################################################
################################################################################################
################################################################################################

####################
# A function from the BayeScan software
# 	 This file is used to plot figures for the software Bayescan in R.

#    This program, BayeScan, aims at detecting genetics markers under selection,
#	 based on allele frequency differences between population. 
#    Copyright (C) 2010  Matthieu Foll
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Arguments:
# - file is the name of your file ex: "output_fst.txt"
# - the q-value threshold corresponding to the target False Discovery Rate (FDR)
# - size is the size of the points and text labels for outliers
# - pos is the distance between the points and the labels 
# - highlight is a optional list of marker indices to display in red.
# - name_highlighted alows to write the indices of highlighted markers instead of using a point like the other markers
# - add_text adds the indices of the outlier markers

# Output:
# This function returns different paremeters in a list
# - outliers: the list of outliers
# - nb_outliers: the number of outliers

# Typical usage: 
# - load this file into R (file/source R code)
# - in R, go to the directory where "output_fst.txt" is (file/change current dir)
# - at the R prompt, type 
# > plot_bayescan("output_fst.txt",0,FDR=0.05)
# if you save the output in a variable, you can recall the different results:
# results<-plot_bayescan("output_fst.txt",0,FDR=0.05)
# results$outliers
# results$nb_outliers

#
# plotting posterior distribution is very easy in R with the output of BayeScan:
# first load the output file *.sel produced by BayeScan
# > mydata=read.table("bi.sel",colClasses="numeric")
# choose the parameter you want to plot by setting for example:
# parameter="Fst1"
# then this line will make the plot for:
# > plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
# you can plot population specific Fst coefficient by setting
# parameter="Fst1"
# if you have non-codominant data you can plot posterior for Fis coefficients in each population:
# parameter="Fis1"
# if you test for selection, you can plot the posterior for alpha coefficient for selection:
# parameter="alpha1"
# you also have access to the likelihood with:
# parameter="logL"
# if you have the package "boa" installed, you can very easily obtain Highest Probability 
# Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
# > boa.hpd(mydata[[parameter]],0.05)
f.plot.bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers), "n"=nrow(res)))
}
####################


