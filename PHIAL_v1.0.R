## PHIAL version 1.0
## Eliezer Van Allen
## 4/2013

## The Broad Institute of MIT and Harvard / Cancer program.
## phial-help@broadinstitute.org

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## It is only to be distributed to editors and reviewers confidentially.

###############################################################################

#Input data
suppressPackageStartupMessages(require(optparse))

option_list <- list(
                    make_option(c("-i", "--individual"), action="store", type="character", help= "Individual patient to be analyzed"),
                    make_option(c("-t", "--tumor_type"), action="store", type="character", help= "Tumor type of patient (i.e. Lung, Colon)"),
                    make_option(c("-o", "--output_dir"), action="store", type="character", help="Output directory for results of analysis", default=getwd()),
                    make_option("--mut.path", action="store", type="character", help="maf_file_oxoG3_capture"),
                    make_option("--mut.dumpster", action="store", type="character", help="maf_file_capture_no_filter"),
                    make_option("--indel.path", action="store", type="character", help="indel_maf_file_capture"),
                    make_option("--dranger.path", action="store", type="character", help="dRanger capture results"),
                    make_option("--segfile.path", action="store", type="character", help="CapSeg segementation file"),
                    make_option("--actdb.mini", action="store", type="character", help="Actionable genes rationales file"),
                    make_option("--actdb.large", action="store", type="character", help="Large actionable genes file"),
                    make_option("--current_panel", action="store", type="character", help="Current clinical panel file"),
                    make_option("--cosmic", action="store", type="character", help="COSMIC file"),
                    make_option("--gsea.pathways", action="store", type="character", help="GSEA Pathways file"),
                    make_option("--gsea.overlap", action="store", type="character", help="GSEA overlap file"),
                    make_option("--gsea.modules", action="store", type="character", help="GSEA modules file"),
                    make_option("--pertinent.negs", action="store", type="character", help="pertinent negs file"),
                    make_option("--refseq", action="store", type="character", help="RefSeq file")
                    )

opt <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)
save(opt, file="debug.RData")
print(opt)

#Load patient data
individual <- opt$individual
tumor_type <- opt$tumor_type
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings=FALSE)
mut.path <- opt$mut.path
indel.path <- opt$indel.path
segfile.path <- opt$segfile.path
dranger.path <- opt$dranger.path

#Sanity check on input
if(is.na(mut.path)&is.na(indel.path)&is.na(segfile.path)&is.na(dranger.path)) {
	print("All input from patient is NA...Goodbye")
	quit()
}
	
print("Loading patient data")
if(is.na(mut.path)) patient.mut <- NA else patient.mut <- read.delim(mut.path, header=TRUE, as.is=TRUE, comment.char="#")
if(is.na(indel.path)) patient.indel <- NA else patient.indel <- read.delim(indel.path, header=TRUE, as.is=TRUE, comment.char="#")
if(is.na(segfile.path)) seg <- NA else seg <- read.delim(segfile.path, header=TRUE, as.is=TRUE, comment.char="#")
if(is.na(dranger.path)) dranger <- NA else dranger <- read.delim(dranger.path, header=TRUE, as.is=TRUE, comment.char="#")

##Load input databases
print("Loading input databases")
actdb_mini <- read.delim(as.character(opt$actdb.mini), header=TRUE, as.is=TRUE)
actdb_large <- read.delim(as.character(opt$actdb.large), header=TRUE, as.is=TRUE)
current_panel <- read.delim(as.character(opt$current_panel), header=TRUE, as.is=TRUE)
COSMIC <- read.delim(as.character(opt$cosmic), header=TRUE, as.is=TRUE)
CGC <- subset(COSMIC, COSMIC$Cancer_census. == "y")
gsea_pathways <- read.delim(as.character(opt$gsea.pathways), header=FALSE, as.is=TRUE)
gsea_overlap_actdb <- read.delim(as.character(opt$gsea.overlap), header=FALSE, as.is=TRUE)
gsea_modules <- read.delim(as.character(opt$gsea.modules), header=FALSE, as.is=TRUE)
refseq <- read.delim(as.character(opt$refseq), header=FALSE, as.is=TRUE)

#---------------------Score Mutations and Indels
print("Scoring somatic alterations")
#Define variant classifications that get immediately demoted to the bottom and are appended at the end
reject_subset <- function(patient, is.reject) {
	rejects <- c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region", "De_novo_Start_InFrame")
	if(is.reject==TRUE) {patient <- subset(patient, patient$Variant_Classification %in% rejects)}
	else patient <- {patient <- subset(patient, !(patient$Variant_Classification %in% rejects))}
	return(patient)
}

#Prepare MAFs for subsequent analyses
build_patient_bin <- function(patient, type) {
  hugo_length <- length(patient$Hugo_Symbol)
  big.bin <- rep("Filtered Calls", hugo_length)
  bin <- rep(11, hugo_length)
  rationale <- rep("", hugo_length)
  notes <- rep("", hugo_length)
  type <- rep(type, hugo_length)

  patient_bin <- cbind(patient, bin, big.bin, rationale, notes, type)
  patient_bin$big.bin <- as.character(patient_bin$big.bin)
  patient_bin$rationale <- as.character(patient_bin$rationale)
  patient_bin$notes <- as.character(patient_bin$notes)

  return(patient_bin)
}

#Complements scoring funciton, finds relevant pathways within MSigDB data to annotate as part of the scoring function
build_pathway_string <- function(db, temp_gene) {
  relevant_col <- grep(temp_gene, db)
  pathway.string <- ""
  
  for (j in seq_along(relevant_col)) {
    temp.col <- relevant_col[j]
    temp.row <- grep(temp_gene, t(db[temp.col]))
    for (k in seq_along(temp.row)) { ##FIX ME??
      if(length(grep(db[temp.row[k],1], pathway.string)) < 1)  {
        pathway.string <- paste(pathway.string, db[temp.row[k], 1], " | ", collapse = " ")
      }
    }
  }
  return (pathway.string)  
}

## Mutations and Indels Scoring function - input is annotated and modified MAF, output is a scored version of that data frame
mut_indel.score.list <- function(An_patient) {
  An_patient[is.na(An_patient)] <- ""
  
  for (i in seq_along(An_patient$Hugo_Symbol)) {
    temp.gene <- An_patient$Hugo_Symbol[i]
    temp.alt <- An_patient$Protein_Change[i]

    if (temp.gene %in% actdb_mini$Gene) {
      temp.loc <- grep(temp.gene, actdb_mini$Gene)
      An_patient$rationale[i] <- actdb_mini$Rationale[temp.loc]
      relevant.subdb <- subset(actdb_large, actdb_large$Gene == temp.gene)
      if ((length(relevant.subdb$Gene) > 0) && (temp.alt %in% relevant.subdb$Alteration)) {
        An_patient$bin[i] <- 0
        An_patient$big.bin[i] <- "Actionable"
      } else if (grepl(temp.alt, An_patient$COSMIC_overlapping_mutations[i], fixed=TRUE)) {			
        ## recurrent alteration in COSMIC
        An_patient$bin[i] <- 1 
        An_patient$big.bin[i] <- "Investigate Actionability"
      } else {
        An_patient$bin[i] <- 2
        An_patient$big.bin[i] <- "High Priority"
      }
    } else if (temp.gene == "") {
      ## Special case for empty gene      
      An_patient$bin[i] <- 9
      An_patient$big.bin[i] <- "VUS"
    } else if (temp.gene %in% CGC$COSMIC_GENE_NAME & grepl(temp.alt, An_patient$COSMIC_overlapping_mutations[i], fixed=TRUE)) {	
      ## CGC gene and COSMIC alteration 
      An_patient$bin[i] <- 3
      An_patient$big.bin[i] <- "Investigate Actionability"
    } else if (temp.gene %in% CGC$COSMIC_GENE_NAME) {
      ## CGC gene, alteration unknown
      An_patient$bin[i] <- 4
      An_patient$big.bin[i] <- "Cancer Gene Census"
       pathway.string <- build_pathway_string(gsea_overlap_actdb, temp.gene)
       An_patient$notes[i] <- pathway.string
    } else if (length(grep(temp.gene, gsea_overlap_actdb)) > 0) {
      ## Gene in cancer pathway that overlaps with MSigDB cancer pathway
      pathway.string <- build_pathway_string(gsea_overlap_actdb, temp.gene)
      An_patient$bin[i] <- 5
      An_patient$big.bin[i] <- "Cancer Pathway"
      An_patient$notes[i] <- pathway.string
    } else if (length(grep(temp.gene, gsea_pathways)) > 0) {
      ## Gene is MSigDB cancer pathway list
      pathway.string <- build_pathway_string(gsea_pathways, temp.gene)
      An_patient$bin[i] <- 6
      An_patient$big.bin[i] <- "Cancer Pathway"
      An_patient$notes[i] <- pathway.string
    } else if (length(grep(temp.gene, gsea_modules)) > 0) {
      ## Gene in MSigDB cancer module list
      An_patient$bin[i] <- 7
      An_patient$big.bin[i] <- "Cancer Module"
    } else if (is.element(temp.gene, COSMIC$COSMIC_GENE_NAME)) {
      ## Gene in COSMIC, alteration = VUS
      An_patient$bin[i] <- 8	
      An_patient$big.bin[i] <- "Cancer Gene"
    } else {
      ## SNPS
      An_patient$bin[i] <- 10
      An_patient$big.bin[i] <- "SNP"
    }
  }

  return(An_patient)
}

#----UniProt Adjustment Function: If mutation/indel is in kinase domain, make it slightly more important within it's own subcategory
sub_sort_uniprot <- function(An_patient) {
	kinase.domain <- which(grepl("kinase", An_patient$UniProt_Region))
	An_patient$bin[kinase.domain] <- An_patient$bin[kinase.domain]-0.25
	return(An_patient)
}
#-----

if(!is.na(mut.path)) {
	print("Scoring somatic mutations")
	patient.mut.bin <- build_patient_bin(patient.mut, "Mutation")
	
	lowest_prob.mut <- reject_subset(patient.mut.bin, TRUE)
	An_patient.mut <- reject_subset(patient.mut.bin, FALSE)
	
	if(nrow(An_patient.mut) == 0) {patient.mut.scored <- patient.mut.bin}
		else {
			patient.mut.scored <- mut_indel.score.list(An_patient.mut)
			patient.mut.scored <- sub_sort_uniprot(patient.mut.scored)
			patient.mut.scored <- rbind(patient.mut.scored, lowest_prob.mut[1:nrow(lowest_prob.mut),])
		}
}

if(!is.na(indel.path)) {
	print("Scoring somatic indels")
	patient.indel.bin <- build_patient_bin(patient.indel, "Indel")
	
	lowest_prob.indel <- reject_subset(patient.indel.bin, TRUE)
	An_patient.indel <- reject_subset(patient.indel.bin, FALSE)
	             
	if (nrow(An_patient.indel) == 0) { patient.indel.scored <- patient.indel.bin } 
		else {
  			patient.indel.scored <- mut_indel.score.list(An_patient.indel)
  			patient.indel.scored <- sub_sort_uniprot(patient.indel.scored)
  			patient.indel.scored <- rbind(patient.indel.scored, lowest_prob.indel[1:nrow(lowest_prob.indel), ])
	}
}

#--------Score Copy Number----------
print("Scoring copy number alterations")
#Function for getting a list of genes in the region defined by the segment
scna.list <- function(seg.sub, refseq, boundary) {
	
	Gene <- c()
	Chromosome <- c()
	Gene_Start <- c()
	Gene_End <- c()
	Segment_Start <- c()
	Segment_End <- c()
	Num_Probes <- c()
	Segment_Mean <- c()
	Tumor_Sample_Barcode <- c()
	
	for (i in 1:nrow(seg.sub)) {
		gene.list <- c()
		seg.line <- seg.sub[i,]
		temp.low <- seg.line$Start[1] - boundary
		temp.high <- seg.line$End[1] + boundary
	
		refseq.sub <- subset(refseq, refseq$V3 == seg.line$Chromosome[1])
		refseq.sub.start <- subset(refseq.sub, (refseq.sub$V5 >= temp.low) & (refseq.sub$V5 <= temp.high))
		refseq.sub.end <- subset(refseq.sub, (refseq.sub$V6 >= temp.low) & (refseq.sub$V6 <= temp.high))
		refseq.sub.merge <- rbind(refseq.sub.start, refseq.sub.end)
		
		refseq.sub <- refseq.sub.merge[!duplicated(refseq.sub.merge$V13),]
		
		if(nrow(refseq.sub)>0) {
		Gene <- c(Gene, refseq.sub$V13)
		Chromosome <- c(Chromosome, refseq.sub$V3)
		
		Gene_Start <- c(Gene_Start, refseq.sub$V5)
		Gene_End <- c(Gene_End, refseq.sub$V6)
		
		Segment_Start <- c(Segment_Start, rep(seg.line$Start[1], nrow(refseq.sub)))
		Segment_End <- c(Segment_End, rep(seg.line$End[1], nrow(refseq.sub)))
			
		Num_Probes <- c(Num_Probes, rep(seg.line$Num_Probes[1], nrow(refseq.sub)))
		Segment_Mean <- c(Segment_Mean, rep(seg.line$Segment_Mean[1], nrow(refseq.sub)))
		Tumor_Sample_Barcode <- c(Tumor_Sample_Barcode, rep(seg.line$Tumor_Sample_Barcode[1], nrow(refseq.sub)))
		}
	}

	if(length(Gene) >0) {
	
		master.df <- data.frame(Gene=Gene, 
			Chromosome=Chromosome, 
			Gene_Start=Gene_Start,
			Gene_End=Gene_End,
			Segment_Start = Segment_Start, 
			Segment_End=Segment_End, 
			Num_Probes=Num_Probes, 
			Segment_Mean=Segment_Mean,
			Tumor_Sample_Barcode=Tumor_Sample_Barcode
		)
	
		In_Segment <- c()
		for (i in 1:nrow(master.df)) {
			temp.gene <- master.df[i,]
			if ((temp.gene$Gene_Start >= temp.gene$Segment_Start) & (temp.gene$Gene_End <= temp.gene$Segment_End))
			In_Segment <- c(In_Segment, TRUE)
			else In_Segment <- c(In_Segment, FALSE)
				
		}
	master.df <- cbind(master.df, In_Segment)
	}
	else master.df <- c()
	master.df
}

build_scna_data <- function(seg, refseq) {
  build_seg_sort_genes <- function(boundary, seg.sort, class) {
    genes <- NULL
    if (nrow(seg.sort) > 0) {
      genes <- scna.list(seg.sort, refseq, boundary)
    } 
    if ((nrow(seg.sort) == 0) || length(genes) == 0) {
      genes <- seg.sort[0, ]
    }
    genes <- cbind(genes, Class=rep(class, nrow(genes)))
    return(genes)
  }
  seg.sort <- seg[order(-seg$Segment_Mean), ]
  ## Start with gains - threshold to Segment_Mean > 2
  seg.sort.gain <- subset(seg.sort, seg.sort$Segment_Mean > 2)
  boundary <- 0  
  amplified.genes <- build_seg_sort_genes(boundary, seg.sort.gain, "Amplified")
  
  ##Then losses - threshold to Segment_Mean < 0.5 ?
  seg.sort.loss <- subset(seg.sort, seg.sort$Segment_Mean < -1)  
  boundary <- 15000 #approximating for one gene away, better solution would be to get to the next discrete gene over
  deletion.genes <- build_seg_sort_genes(boundary, seg.sort.loss, "Deleted")

  return(rbind(amplified.genes, deletion.genes))
}

scna.directionality <- function(seg.mean, temp.type) {
	score <- FALSE
	if ((seg.mean > 0) && grepl("Amplification", temp.type)) {
    score <- TRUE
	} 
  if ((seg.mean < 0) && (grepl("Deletion", temp.type) || grepl("Biallelic Inactivation", temp.type))) {
		score <- TRUE
	}
  return(score)
}

#Copy Number modified heuristics - gene name only
scna.score.list <- function(An_patient) {	
  An_patient[is.na(An_patient)] <- ""
  for (i in seq_len(nrow(An_patient))) {
	temp.gene <- An_patient$Gene[i]
	
    if (temp.gene %in% actdb_mini$Gene) {
  		temp.loc <- grep(temp.gene, actdb_mini$Gene)
	  	An_patient$rationale[i] <- actdb_mini$Rationale[temp.loc]
		  temp.type <- actdb_mini$Types_of_recurrent_alterations[temp.loc]
		
  		if(scna.directionality(An_patient$Segment_Mean[i], temp.type)) {
	  		An_patient$bin[i] <- 0
		  	An_patient$big.bin[i] <- "Investigate Actionability"	
  		} else {
			  An_patient$bin[i] <- 4
			  An_patient$big.bin[i] <- "Investigate Biological Significance"
		  }	
  	} else if (temp.gene == "") {
      ## Special case for empty gene
	  	An_patient$bin[i] <- 9
	  	An_patient$big.bin[i] <- "VUS"
  	} else if (temp.gene %in% CGC$COSMIC_GENE_NAME) {
  	  ## CGC gene, alteration unknown
	  	if (An_patient$Segment_Mean[i] < 0) {
        ## less since low level deletions = garbage
        An_patient$bin[i] <- 7 
	  	} else {
	  		pathway.string <- build_pathway_string(gsea_overlap_actdb, temp.gene)
       		 An_patient$bin[i] <- 4
       		 An_patient$notes[i] <- pathway.string
	  	}
		  An_patient$big.bin[i] <- "Cancer Gene Census"
	  } else if (length(grep(temp.gene, gsea_overlap_actdb)) > 0) {
        ## Gene in cancer pathway that overlaps with MSigDB cancer pathway
        pathway.string <- build_pathway_string(gsea_overlap_actdb, temp.gene)
     		An_patient$bin[i] <- 5
  	  	An_patient$big.bin[i] <- "Cancer Pathway"
  		  An_patient$notes[i] <- pathway.string
	  } else if (length(grep(temp.gene, gsea_pathways)) > 0) {
  	    ## Gene is MSigDB cancer pathway list
        pathway.string <- build_pathway_string(gsea_pathways, temp.gene)
    		An_patient$bin[i] <- 6
  	  	An_patient$big.bin[i] <- "Cancer Pathway"
  	  	An_patient$notes[i] <- pathway.string
    } else if (length(grep(temp.gene, gsea_modules)) > 0) {
  	    ## Gene in MSigDB cancer module list
    		An_patient$bin[i] <- 7
  	  	An_patient$big.bin[i] <- "Cancer Module"
  	} else if (is.element(temp.gene, COSMIC$COSMIC_GENE_NAME)) {
	      ## Gene in COSMIC, alteration = VUS
  	  	An_patient$bin[i] <- 8	
  		  An_patient$big.bin[i] <- "Cancer Gene"
	  } else {
	      ## VUS that is not a SNP
	  	  An_patient$bin[i] <- 9
	    	An_patient$big.bin[i] <- "VUS"
	  }	
  }
  return(An_patient)	
}

if(!is.na(segfile.path)) {
	seg$Tumor_Sample_Barcode <- rep(patient.mut$Tumor_Sample_Barcode[1], nrow(seg)) #hack to match up Tumor_Sample_Barcode naming conventions
	scna.data <- build_scna_data(seg, refseq)
	
	scna.data <- cbind(scna.data, big.bin=rep(c("Filtered Calls"), nrow(scna.data)))
	scna.data$big.bin <- as.character(scna.data$big.bin)
	scna.data <- cbind(scna.data, bin=rep(11, nrow(scna.data)))
	scna.data <- cbind(scna.data, rationale=rep(c(""), nrow(scna.data)))
	scna.data$rationale <- as.character(scna.data$rationale)
	scna.data <- cbind(scna.data, notes=rep(c(""), nrow(scna.data)))
	scna.data$notes <- as.character(scna.data$notes)

	if(nrow(scna.data) > 0) scna.data.scored <- scna.score.list(scna.data)
	if(nrow(scna.data) == 0) scna.data.scored <- scna.data

	scna.data.scored$Gene <- as.character(scna.data.scored$Gene)
	#temp.false.pos <- grep(FALSE, scna.data.scored$In_Segment)
	#scna.data.scored$Gene[temp.false.pos] <- paste(scna.data.scored$Gene[temp.false.pos], "-", sep="")

#scna.data.sort <- scna.data.scored[order(scna.data.scored$bin),]
}

#-------------Score Rearrangements------------------
print("Scoring rearrangements")
#Modified heuristics for rearrangements
dranger.score.list <- function(An_patient) {
  An_patient[is.na(An_patient)] <- ""
  
  for (i in seq_len(nrow(An_patient))) {
    temp.gene1 <- An_patient$gene1[i]
    temp.gene2 <- An_patient$gene2[i]
    temp_genes <- c(temp.gene1, temp.gene2)
    
    if (any(temp_genes %in% actdb_mini$Gene)) {
      temp.gene = ifelse(temp_genes[1] %in% actdb_mini$Gene, temp_genes[1], temp_genes[2])
      temp.loc <- grep(temp.gene, actdb_mini$Gene)
      An_patient$rationale[i] <- actdb_mini$Rationale[temp.loc]
      temp.type <- actdb_mini$Types_of_recurrent_alterations[temp.loc]
      
      if(is.element("Rearrangement", temp.type)) {
        An_patient$bin[i] <- 0
        An_patient$big.bin[i] <- c("Investigate Actionability")
      }
      else {
        An_patient$bin[i] <- 4
        An_patient$big.bin[i] <- c("Investigate Actionability")
      }      
    } else if (any(temp_genes == "")) {
      ## Special case for empty gene
      An_patient$bin[i] <- 9
      An_patient$big.bin[i] <- "VUS"
    } else if (any(temp_genes %in% CGC$COSMIC_GENE_NAME)) {
      #CGC gene, alteration unknown      
      An_patient$bin[i] <- 4
      An_patient$big.bin[i] <- "High Priority"
    } else if (length(grep(temp.gene1, gsea_overlap_actdb)) > 0  ||  length(grep(temp.gene2, gsea_overlap_actdb)) > 0) {
      ## FIXME: can we simplify the || on the grep() calls?
      #Gene in cancer pathway that overlaps with MSigDB cancer pathway
      An_patient$bin[i] <- 5
      An_patient$big.bin[i] <- "Cancer Pathway"
    } else if (length(grep(temp.gene1, gsea_pathways)) > 0  ||  length(grep(temp.gene2, gsea_pathways)) > 0) {
      ## FIXME: can we simplify the || on the grep() calls?
      #Gene is MSigDB cancer pathway list      
      An_patient$bin[i] <- 6
      An_patient$big.bin[i] <- "Cancer Pathway"
      An_patient$notes[i] <- pathway.string
    } else if (length(grep(temp.gene1, gsea_modules)) > 0  ||  length(grep(temp.gene2, gsea_modules)) > 0) {
      ## FIXME: can we simplify the || on the grep() calls?      
      #Gene in MSigDB cancer module list      
      An_patient$bin[i] <- 7
      An_patient$big.bin[i] <- "Cancer Module"
    } else if (any(temp_genes %in% COSMIC$COSMIC_GENE_NAME)) {
      #Gene in COSMIC, alteration = VUS      
      An_patient$bin[i] <- 8  
      An_patient$big.bin[i] <- "Cancer Gene"
    } else {
      #VUS that is not a SNP      
      An_patient$bin[i] <- 9
      An_patient$big.bin[i] <- "VUS"
    } 
  }
  
  return(An_patient)  
}

if(nrow(dranger) > 0) {
	dranger <- cbind(dranger, Gene_fusion = paste(dranger$gene1, ":", dranger$gene2, sep=""))
	dranger$bin <- rep(11, nrow(dranger))
	dranger$big.bin <- rep(c("Filtered Calls"), nrow(dranger))
	dranger$rationale <- rep(c(""), nrow(dranger))
	dranger.scored <- dranger.score.list(dranger)
}
if(nrow(dranger) == 0) dranger.scored <- dranger


#-----Create data frames from each of the four outputs so that they can be merged into one master df, knowing that some data will be in only a subset

print("Merging data sets")
if(nrow(dranger.scored) == 0) {
	dranger.sorted = NULL	
} else {
	dranger.sorted <- data.frame(Gene=dranger.scored$Gene_fusion,
		Variant_Classification = rep("Rearrangement", nrow(dranger.scored)),
		Alteration = paste(dranger.scored$site1, " : ", dranger.scored$site2, sep=""),
		Tumor_allele=rep(NA, nrow(dranger.scored)),
		Reference_allele = rep(NA, nrow(dranger.scored)),
		dbSNP_RS = rep(NA, nrow(dranger.scored)),
		bin=dranger.scored$bin,
		Score_bin = dranger.scored$big.bin,
		Chromosome = paste(dranger.scored$chr1, ":", dranger.scored$chr2, sep=""),
		Start_position = dranger.scored$pos1,
		End_position = dranger.scored$pos2,
		Coverage = rep(NA, nrow(dranger.scored)),
		Allelic_fraction = rep(NA, nrow(dranger.scored)),
		Number_of_Probes = rep(NA, nrow(dranger.scored)),
		Segment_Mean = rep(NA, nrow(dranger.scored)),
		In_Segment=rep(NA, nrow(dranger.scored)),
		Pathways=rep(NA, nrow(dranger.scored)),
		COSMIC_overlapping_mutations=rep(NA, nrow(dranger.scored)),
		COSMIC_total_alterations_in_gene=rep(NA, nrow(dranger.scored)),
		UniProt_Region = rep(NA, nrow(dranger.scored)),
		Rationale= rep(NA, nrow(dranger.scored)),
		Tumor_Sample_Barcode = rep(patient.mut.scored$Tumor_Sample_Barcode[1], nrow(dranger.scored)), stringsAsFactors=F
	)
}

if (nrow(patient.mut.scored) == 0) {
  patient.mut.sorted <- NULL
} else {
	patient.mut.sorted <- data.frame(Gene=patient.mut.scored$Hugo_Symbol, 
		Variant_Classification=patient.mut.scored$Variant_Classification, 
		Alteration=patient.mut.scored$Protein_Change,
		Tumor_allele=patient.mut.scored$Tumor_Seq_Allele1,
		Reference_allele = patient.mut.scored$Reference_Allele,
		dbSNP_RS=patient.mut.scored$dbSNP_RS, 
		bin=patient.mut.scored$bin,
		Score_bin=as.character(patient.mut.scored$big.bin),
		Chromosome=patient.mut.scored$Chromosome, 
		Start_position=patient.mut.scored$Start_position, 
		End_position=patient.mut.scored$End_position, 
		Coverage = (patient.mut.scored$t_ref_count+patient.mut.scored$t_alt_count), 
		Allelic_fraction = (patient.mut.scored$t_alt_count/(patient.mut.scored$t_ref_count+patient.mut.scored$t_alt_count)),
		Number_of_Probes=rep(NA, nrow(patient.mut.scored)),
		Segment_Mean=rep(NA, nrow(patient.mut.scored)),
		In_Segment=rep(NA, nrow(patient.mut.scored)),
		Pathways=patient.mut.scored$notes,
		COSMIC_overlapping_mutations=patient.mut.scored$COSMIC_overlapping_mutations,
		COSMIC_total_alterations_in_gene=patient.mut.scored$COSMIC_total_alterations_in_gene, 
		UniProt_Region = patient.mut.scored$UniProt_Region,
		Rationale = patient.mut.scored$rationale,
		Tumor_Sample_Barcode=patient.mut.scored$Tumor_Sample_Barcode, stringsAsFactors=F
		)
}


if (nrow(patient.indel.scored) == 0) {
  patient.indel.sorted <- NULL
} else {
	patient.indel.sorted <- data.frame(Gene=patient.indel.scored$Hugo_Symbol, 
		Variant_Classification=patient.indel.scored$Variant_Classification, 
		Alteration=patient.indel.scored$Protein_Change, 
		Tumor_allele=patient.indel.scored$Tumor_Seq_Allele1,
		Reference_allele = patient.indel.scored$Reference_Allele,
		dbSNP_RS = patient.indel.scored$dbSNP_RS,
		bin=patient.indel.scored$bin,
		Score_bin=as.character(patient.indel.scored$big.bin),
		Chromosome=patient.indel.scored$Chromosome, 
		Start_position=patient.indel.scored$Start_position, 
		End_position=patient.indel.scored$End_position, 
		Coverage = (patient.indel.scored$t_ref_count+patient.indel.scored$t_alt_count), 
		Allelic_fraction = (patient.indel.scored$t_alt_count/(patient.indel.scored$t_ref_count+patient.indel.scored$t_alt_count)),
		Number_of_Probes=rep(NA, nrow(patient.indel.scored)),
		Segment_Mean=rep(NA, nrow(patient.indel.scored)),
		In_Segment=rep(NA, nrow(patient.indel.scored)),
		Pathways=patient.indel.scored$notes,
		COSMIC_overlapping_mutations=patient.indel.scored$COSMIC_overlapping_mutations,
		COSMIC_total_alterations_in_gene=patient.indel.scored$COSMIC_total_alterations_in_gene,
		UniProt_Region = patient.indel.scored$UniProt_Region,
		Rationale=patient.indel.scored$rationale,
		Tumor_Sample_Barcode=patient.indel.scored$Tumor_Sample_Barcode, stringsAsFactors=F
		) 
}

if (nrow(scna.data.scored) == 0) {
  patient.scna.sorted <- NULL 
} else {
	patient.scna.sorted <- data.frame(Gene=scna.data.scored$Gene, 
		Variant_Classification=rep(c("Copy Number"), nrow(scna.data.scored)), 
		Alteration=scna.data.scored$Class,
		Tumor_allele=rep(NA, nrow(scna.data.scored)),
		Reference_allele = rep(NA, nrow(scna.data.scored)),
		dbSNP_RS= rep(NA, nrow(scna.data.scored)),
		bin=scna.data.scored$bin,
		Score_bin=as.character(scna.data.scored$big.bin),
		Chromosome=scna.data.scored$Chromosome, 
		Start_position=scna.data.scored$Segment_Start, 
		End_position=scna.data.scored$Segment_End, 
		Coverage = rep(NA, nrow(scna.data.scored)), 
		Allelic_fraction = rep(NA, nrow(scna.data.scored)),
		Number_of_Probes=scna.data.scored$Num_Probes,
		Segment_Mean=scna.data.scored$Segment_Mean,
		In_Segment=scna.data.scored$In_Segment,
		Pathways=scna.data.scored$notes,
		COSMIC_overlapping_mutations= rep(NA, nrow(scna.data.scored)),
		COSMIC_total_alterations_in_gene=rep(NA, nrow(scna.data.scored)), 
		UniProt_Region = rep(NA, nrow(scna.data.scored)),
		Rationale = scna.data.scored$rationale,
		Tumor_Sample_Barcode = scna.data.scored$Tumor_Sample_Barcode, stringsAsFactors=F
)
}

patient.merged <- rbind(patient.mut.sorted, patient.indel.sorted, patient.scna.sorted, dranger.sorted)

if (is.null(patient.merged)) {
  ## Handling an extremely pathological case, basically the case where there's
  ## absolutely no data to push downstream
  stop("No data for individual")
}

# Then re-loop around pathway findings to see if any are in a gene with somatic actionable event and elevate those
patient.sorted <- patient.merged
pathway.links <- unique(as.character(subset(patient.sorted$Gene, patient.sorted$bin < 4)))
patient.sorted$Score_bin <- as.character(patient.sorted$Score_bin)
patient.sorted$Pathways <- as.character(patient.sorted$Pathways)
if(length(pathway.links) > 0) {
for (j in 1:nrow(patient.sorted)) {
	scored <- FALSE
	linked.genes <- c()
	if (patient.sorted$bin[j] == 4 | patient.sorted$bin[j] == 5 | patient.sorted$bin[j] == 6) {
		#temp.gene <- patient.sorted$Gene[j]
		temp.pathway <- patient.sorted$Pathways[j]
		for (k in 1:nrow(gsea_pathways)) {
			if(grepl(gsea_pathways[k,1], temp.pathway) & !scored) {
				temp.row <- gsea_pathways[k,]
				for (l in 1:length(pathway.links)) {
					temp.gene <- pathway.links[l]					
					if (is.element(temp.gene, temp.row[3:ncol(temp.row)])) {
						if (!scored) patient.sorted$bin[j] <- (patient.sorted$bin[j] - 0.5)
						linked.genes <- c(linked.genes, temp.gene)
						#patient.sorted$Score_bin[j] <- as.character(paste("Cancer Pathway Linked to ", temp.gene, sep=""))
						scored <- TRUE
					}
				}
			}		
		}
	if(length(linked.genes) > 1) {
		temp.link <- linked.genes[1]
		patient.sorted$bin[j] <- patient.sorted$bin[j] - (0.5-(0.5/length(linked.genes)))
		for(w in 2:length(linked.genes)) temp.link <- paste(temp.link, linked.genes[w], sep = "; ")
		linked.genes <- temp.link		
	}
	#linked.genes <- paste(linked.genes, sep = "; ")
	if(scored) patient.sorted$Score_bin[j] <- paste("Cancer Pathway Linked to: ", linked.genes, sep="")		
	}
}
}

##---Find links function-------
## Appends relevant HTML links for subsets of alterations that will be used for reporting downstream
find.links <- function(patient.df) {
	mutation.class <- c("Missense_Mutation", "Nonsense_Mutation")
	
	patient.df$Clinical_trials <- rep("", nrow(patient.df))
	patient.df$mutation.assessor <- rep("", nrow(patient.df))
	patient.df$Variant_Classification <- as.character(patient.df$Variant_Classification)
	
	for (i in 1:nrow(patient.df)) {
		#tmp.disease	<- paste(tumor_type, sep="")
		trial.link <- paste("<a href=\"http://clinicaltrials.gov/ct2/results?term=", tumor_type, "+AND+cancer+AND+", patient.df$Gene[i], "&recr=Open\">Click here</a>", sep="")
		patient.df$Clinical_trials[i] <- trial.link		
		
		if(is.element(patient.df$Variant_Classification[i], mutation.class)) {
			patient.df$mutation.assessor[i] <- paste("<a href=\"http://mutationassessor.org/?cm=var&var=hg19,", patient.df$Chromosome[i], ",", patient.df$Start_position[i], ",", patient.df$Reference_allele[i], ",", patient.df$Tumor_allele[i], "&fts=all\">", patient.df$Alteration[i], "</a>", sep="")
		}		
	}
	return(patient.df)
}

#--------ORGANIZATION FOR OUTPUT----------
patient.sorted <- find.links(patient.sorted)
patient.sorted <- patient.sorted[order(patient.sorted$bin, patient.sorted$Gene),]

temp.filename <- paste(output_dir, "/", individual, "_complete_muts_indels_scna_detailed.txt", sep="")
write.table(patient.sorted, file=temp.filename, sep="\t", row.names=FALSE, quote=FALSE)

#Complete table 
complete.sorted <- patient.sorted
complete.sorted[is.na(complete.sorted)] <- ""
temp.filename <- paste(output_dir, "/", individual, "_complete_muts_indels_scna.txt", sep="")
write.table(complete.sorted, file=temp.filename, sep="\t", row.names=FALSE, quote=FALSE)


# First take highly actionable genes
top.sorted <- subset(patient.sorted, patient.sorted$bin <= 1)
top.sorted.t <- data.frame(Gene=top.sorted$Gene,
	Alteration=top.sorted$Alteration,
	Variant=top.sorted$Variant_Classification,
	Coverage=top.sorted$Coverage,
	Allelic_fraction=top.sorted$Allelic_fraction,
	Number_of_Probes=top.sorted$Number_of_Probes,
	Segment_Mean=top.sorted$Segment_Mean,
	UniProt_Region=top.sorted$UniProt_Region,
	Trials = top.sorted$Clinical_trials,
	Rationale = top.sorted$Rationale,
	Mutation_assessor = top.sorted$mutation.assessor, stringsAsFactors=F
)

top.sorted.t[is.na(top.sorted.t)] <- ""
temp.filename <- paste(output_dir, "/", individual, "_investigate_clinical_relevance_high.txt", sep="")
write.table(top.sorted.t, file=temp.filename, sep="\t", row.names=FALSE, quote=FALSE)

# Now take the actionable/possibly actionable and tier them seperately for output in Nozzle doc
act.sorted <- subset(patient.sorted, (patient.sorted$bin >1)&(patient.sorted$bin <3.5))
act.sorted.t <- data.frame(Gene=act.sorted$Gene,
	Alteration=act.sorted$Alteration,
	Variant=act.sorted$Variant_Classification,
	Coverage=act.sorted$Coverage,
	Allelic_fraction=act.sorted$Allelic_fraction,
	Number_of_Probes=act.sorted$Number_of_Probes,
	Segment_Mean=act.sorted$Segment_Mean,
	UniProt_Region=act.sorted$UniProt_Region,
	Trials = act.sorted$Clinical_trials,
	Rationale = act.sorted$Rationale,
	Mutation_assessor = act.sorted$mutation.assessor, stringsAsFactors=F
)

act.sorted.t[is.na(act.sorted.t)] <- ""
temp.filename <- paste(output_dir, "/", individual, "_investigate_clinical_relevance_low.txt", sep="")
write.table(act.sorted.t, file=temp.filename, sep="\t", row.names=FALSE, quote=FALSE)

# Now grab everything above 'junk in COSMIC' OR cancer module - may use score or score_bin to highlight some of these...
additional.sorted <- subset(patient.sorted, (patient.sorted$bin >=3.5)&(patient.sorted$bin <7))
additional.sorted.t <- data.frame(Gene=additional.sorted$Gene,
	Alteration=additional.sorted$Alteration,
	Variant=additional.sorted$Variant_Classification,
	Score = additional.sorted$bin, 
	Score_Bin=additional.sorted$Score_bin,
	Coverage=additional.sorted$Coverage,
	Allelic_fraction=additional.sorted$Allelic_fraction,
	Pathways_involved=additional.sorted$Pathways,
	UniProt_Region=additional.sorted$UniProt_Region,
	Mutation_assessor = additional.sorted$mutation.assessor, stringsAsFactors=F
)

additional.sorted.t[is.na(additional.sorted.t)] <- ""
temp.filename <- paste(output_dir, "/", individual,"_investigate_biological_relevance.txt", sep="")
write.table(additional.sorted.t, file=temp.filename, sep="\t", row.names=FALSE, quote=FALSE)

#--------Make PHIAL Gel-----
print("Making PHIAL gel")

make_patient_gel <- function(maf, png.path) {

suppressPackageStartupMessages(require(gplots))

par(mar=c(6, 2, 2, 2))
png(png.path, width=5, height=9, res=300, units="in")
	maf$bin2 <- -maf$bin
	maf$pos <- rep(1, nrow(maf))
	maf$labels <- rep("", nrow(maf))
	for(j in 1:nrow(maf)) {
			maf$labels[j] <- ifelse(maf$bin[j] < 5, paste("   ", maf$Gene[j], "   ", sep=""), "")
	}
	x.stuff <- ifelse(maf$bin2>=-1, jitter(maf$pos, factor=10), 
					ifelse(maf$bin2>=-4, jitter(maf$pos, factor=30), jitter(maf$pos, factor=40)))
						#ifelse(maf$bin2 == 4.5, jitter))
	y.stuff <- ifelse(maf$bin2>=-1, jitter(maf$bin2, factor=5), 
					ifelse(maf$bin2>=-3, jitter(maf$bin2, factor=7), 
						ifelse(maf$bin2>=-4, jitter(maf$bin2, factor=5), jitter(maf$bin2, factor=5))))
						
	myplot <- plot(x.stuff, y.stuff, col=(ifelse(maf$bin2>-3, "red", ifelse(maf$bin2>-6, "orange",
														ifelse(maf$bin2>-8, "gold", 
														ifelse(maf$bin2 >-11, "yellow", "lightgray"))))), 
	pch=20, axes=F, xlab="", ylab="", xlim=c(0,2), ylim=c(-12,1), cex=3)
		
	for(k in 1:nrow(maf)) {
		if (maf$bin[k] >=3 & maf$bin[k] <5) {
			if (grepl("Cancer Pathway Linked to: ", maf$Score_bin[k])) {
				tmp.line <- unlist(strsplit(maf$Score_bin[k], "Cancer Pathway Linked to: "))
				tmp.line <- unlist(strsplit(tmp.line[2], "; "))
				for(l in 1:length(tmp.line)) {
					higher.gene.loc <- grep(tmp.line[l], maf$Gene, fixed=T)
					for(m in 1:length(higher.gene.loc)) {
						tmp.m <- higher.gene.loc[m]
						if(maf$bin[tmp.m] < 4) {
						segments(x.stuff[k], y.stuff[k], x.stuff[tmp.m], y.stuff[tmp.m], 
						col=ifelse(maf$bin[tmp.m]<2, "lightcoral", "lightgoldenrod"), lty=2)	
						}
					}					
				}
			}
		}
	}
	
	par(font=2)
	for(j in 1:length(x.stuff)) {
		text(x.stuff[j], y.stuff[j], labels=maf$labels[j], cex=1, offset = 1, adj=ifelse(j%%2==0, 0, 1))
	}
	axis(1, at=c(1), labels=" ", tick=F, cex.axis=1, font.axis=2)
	par(las=0)
	smartlegend(x="right", y="top", c("Clinical", "Biological", "Pathway", "COSMIC", "Syn."), c("red", "orange", "gold", "yellow", "lightgray"), cex=1, inset=0.02, bty="n")
	box(which="plot", lty=1, col="black", lwd=2.5)
dev.off()
return(TRUE)
}
tmp.path=paste(output_dir, "/", individual, "_phial_gel.png", sep="")
make_patient_gel(patient.sorted, tmp.path)

#---Make Nozzle HTML report----
print("Creating Nozzle HTML report")
source("./Nozzle_template.R")
make_nozzle_report(individual, tumor_type, output_dir)

print("Somatic PHIAL completed.")
quit()
