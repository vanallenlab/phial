## Wrapper script to run somatic PHIAL
## Eliezer Van Allen
## 4/2013

## The Broad Institute of MIT and Harvard / Cancer program.
## phial-help@broadinstitute.org

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


suppressPackageStartupMessages(require(optparse))

option_list <- list(
	make_option(c("-i", "--individual_input_data"), action="store", type="character", help= "Individual patient to be analyzed"),
	make_option(c("-d", "--databases_paths"), action="store", type="character", help="Databases to be used")
)

# Load files

opt <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

print ("Loading relevant files...")

patient_data <- read.delim(opt$individual_input_data, header=T, as.is=T, comment.char="#")
databases <- read.delim(opt$databases_paths, header=T, as.is=T, comment.char="#")

#Somatic PHIAL
for(i in 1:nrow(patient_data)) {
	bsub.command <- paste("Rscript PHIAL_v1.0.R", 
		" -i ", patient_data$individual[i],
		" -t ", patient_data$tumor_type[i],
		" -o ", patient_data$output_dir[i],
		" --mut.path ", patient_data$mutation_data[i],
		" --indel.path ", patient_data$indel_data[i],
		" --dranger.path ", patient_data$rearrangement_data[i],
		" --segfile.path ", patient_data$copynumber_data[i],
		" --actdb.mini ", databases$actdb.mini,
		" --actdb.large ", databases$actdb.large,
		" --current_panel ", databases$current_panel,
		" --cosmic ", databases$cosmic,
		" --gsea.pathways ", databases$gsea.pathways,
		" --gsea.overlap ", databases$gsea.overlap,
		" --gsea.modules ", databases$gsea.modules,
		" --refseq ", databases$refseq,
		sep="")

	system(bsub.command, ignore.stdout=FALSE)
}
