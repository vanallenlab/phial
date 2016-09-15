## Eliezer Van Allen
## Nozzle exome report generator
## 4/2013

## The Broad Institute of MIT and Harvard / Cancer program.
## phial-help@broadinstitute.org

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## It is only to be distributed to editors and reviewers confidentially.

###############################################################################


make_nozzle_report <- function(individual, tumor_type, output_dir) {
	
	suppressPackageStartupMessages(require(gplots))
	suppressPackageStartupMessages(require(Nozzle.R1))

#Standard text that is used periodically
empty.section <- newParagraph("There were no alterations that qualified for this section.")
spacer <- newParagraph(" ")

#Set up definition mini-reports for downstream usage
actionable.definition <- newResult("INVESTIGATE CLINICAL RELEVANCE", isSignificant=TRUE)

#-----Patient Information Section
patient.ID <- newParagraph("ID: ", individual)
tumor.type <- newParagraph("Tumor type: ", tumor_type, " adenocarcinoma")
patient.info.list <- newList(patient.ID, tumor.type, isNumbered=FALSE)

patient.info <- newSection("Clinical Information")
patient.info <- addTo(patient.info, patient.info.list)

##-----Sequencing Metrics 
#AF histogram as purity QC
complete.scored.path <- paste(output_dir, "/", individual, "_complete_muts_indels_scna_detailed.txt", sep="")
complete.scored <- read.delim(complete.scored.path, header=T, as.is=T, comment.char="#")
png.file.path <- paste(output_dir, "/", individual, "_allelicfx_hist.png", sep="")
png(png.file.path, width=10, height=10, res=300, units="in")
af.hist <- hist(complete.scored$Allelic_fraction, xlab="Allelic Fraction", ylab="Frequency", main="Histogram of allelic fraction for mutations or insertion/deletions", xlim=c(0,1), col="cornflowerblue")
text(0.5, max(af.hist$counts)/2, paste("Median: ", signif(median(complete.scored$Allelic_fraction, na.rm=TRUE), digits=2)))
text(0.5, max(af.hist$counts)/2+5, paste("Mean: ", signif(mean(complete.scored$Allelic_fraction, na.rm=TRUE), digits=2)))
dev.off()
png.file.path <- paste("./", individual, "_allelicfx_hist.png", sep="")
af.fig <- newFigure(png.file.path, "Allelic fraction histogram")

summary.metrics <- newSection("Sequencing Metrics")
summary.metrics <- addTo(summary.metrics,  af.fig)

##-----Clinically relevant gene/alteration highlights
actionable.p <- newParagraph("The following alterations were classified as ", actionable.definition, ":")
actionable.path.high <- paste(output_dir, "/", as.character(individual), "_investigate_clinical_relevance_high.txt", sep="")
actionable.path.low <- paste(output_dir, "/", as.character(individual), "_investigate_clinical_relevance_low.txt", sep="")
actionable.table.file.high <- read.delim(actionable.path.high, header=TRUE, as.is=TRUE)
actionable.table.file.low <- read.delim(actionable.path.low, header=TRUE, as.is=TRUE)
actionable.table.file <- rbind(actionable.table.file.high, actionable.table.file.low)

actionable.mini <- actionable.table.file[,1:3] #keep it simple for this display
if(nrow(actionable.mini)>0) actionable.nozzle.table <- newTable(actionable.mini)
if(nrow(actionable.mini)==0) actionable.nozzle.table <- empty.section

#Patient gel addition
tmp.path=paste("./", individual, "_phial_gel.png", sep="")
somatic.gel <- newFigure(tmp.path, fileHighRes=NA, "PHIAL Somatic Patient Gel")

actionable.highlights <- newSection("Clinical Highlights")
actionable.highlights <- addTo(actionable.highlights, actionable.p, actionable.nozzle.table, somatic.gel)


##-----Somatic Alterations
#--Assign somatic links function
assign.som.links <- function(som.df) {
	mutation.class <- c("Missense_Mutation", "Nonsense_Mutation")
	if(nrow(som.df) > 0) {
	for(i in 1:nrow(som.df)) {
		if(som.df$Variant[i] %in% mutation.class) {
			som.df$Alteration[i] <- som.df$Mutation_assessor[i]
			}
		}
		}
	return(som.df)
}

actionable.table.file$Allelic_fraction <- as.character(round(as.numeric(actionable.table.file$Allelic_fraction), digits=2))
actionable.table.file[is.na(actionable.table.file)] <- ""
if(actionable.table.file$Number_of_Probes[1] == "" & length(unique(actionable.table.file$Number_of_Probes == ""))==1) actionable.table.file <- subset(actionable.table.file, select = -c(Number_of_Probes, Segment_Mean))
if(actionable.table.file$Coverage[1] == "" & length(unique(actionable.table.file$Coverage == ""))==1) actionable.table.file <- subset(actionable.table.file, select = -c(Coverage, Allelic_fraction))

actionable.table.file <- assign.som.links(actionable.table.file)
actionable.table.file <- subset(actionable.table.file, select = -c(Mutation_assessor))
actionable.table <- newTable(actionable.table.file, "Clinically relevant findings with details, sorted by PHIAL score", significantDigits=2)

actionable.section <- newSubSection("Investigate Clinical Relevance")
actionable.section <- addTo(actionable.section, actionable.table)

##Investigate Biological Relevance
possible.path <- paste(output_dir, "/", as.character(individual), "_investigate_biological_relevance.txt", sep="")
possible.table.file <- read.delim(possible.path, header=TRUE, as.is=TRUE)
possible.section <- newSubSection("Investigate Biological Relevance")
if (nrow(possible.table.file)>0) {
	possible.table.file[is.na(possible.table.file)]<-""
	tmp.score <- possible.table.file$Score
	possible.table.file <- assign.som.links(possible.table.file)
	possible.table.file <- subset(possible.table.file, select= -c(Score, Pathways_involved, Mutation_assessor))
	possible.table <- newTable(possible.table.file, file=possible.path, "Biologically relevant findings, sorted by PHIAL score. Highest scoring alterations or those with mutations in the kinase domain are labeled red.")
		
	for (i in 1:nrow(possible.table.file)) {
		if(grepl("kinase", possible.table.file$UniProt_Region[i], ignore.case=TRUE)) {
			result2 <- newResult("", isSignificant=TRUE)		
			possible.table <- addTo(possible.table, result2, row=i, column=1)
		}
	}
	possible.section <- addTo(possible.section, possible.table)
	}
if(nrow(possible.table.file)==0) possible.section <- addTo(possible.section, empty.section)

#Complete table - to get empty table with link to full table just pass column names
complete.table.path <- paste(output_dir, "/", as.character(individual), "_complete_muts_indels_scna_detailed.txt", sep="")
complete.table <- read.delim(complete.table.path, header=TRUE, as.is=TRUE)
complete.table.path <- paste("./", as.character(individual), "_complete_muts_indels_scna_detailed.txt", sep="")
complete.table.names <- colnames(complete.table)
complete.table.header <- append(complete.table.names[1:3], "...")
complete.table.header <- newTable(complete.table.header, file=complete.table.path, "A link to complete table of somatic alterations variants")
complete.section <- newSubSection("Complete listing of somatic alterations")
complete.section <- addTo(complete.section, complete.table.header)

##Organize this section
somatic.alts <- newSection("Somatic Analysis")
somatic.p <- newParagraph("This section investigates somatic mutations, insertion/deletions, copy number alterations, and rearrangements across the exome.")
somatic.alts <- addTo(somatic.alts, somatic.p, spacer, actionable.section, possible.section, complete.section)

#-----References
analysis.references.header <- newParagraph("References pertaining to these analyses:")
Nozzle.citation <- newCitation(authors="Nils Gehlenborg, et al.", title="The Nozzle Report Package", year="2011", url="http://gdac.broadinstitute.org/nozzle")
cga.citation <- newCitation(authors="Broad Cancer Genome Analaysis Team", title="Analysis Tools", year="2012", url="http://www.broadinstitute.org/cancer/cga")
heuristic.citation <- newCitation(authors="Van Allen, Wagle, et al.", title="A heuristic platform for clinical interpretation of cancer genome sequencing data.", year="2012", url="http://abstract.asco.org/AbstView_114_97268.html")

references.section <- newSection("References")
references.section <- addTo(references.section, analysis.references.header, Nozzle.citation, cga.citation, heuristic.citation)

#-----Assemble Report
patient.report <- newCustomReport("Cancer Genome Report")
patient.report <- addTo(patient.report, patient.info, summary.metrics, actionable.highlights,somatic.alts, references.section)
temp.filename <- paste(output_dir, "/", as.character(individual), "_cancer_genome_report", sep="")
writeReport(patient.report, filename=temp.filename)
}