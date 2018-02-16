# PHIAL (Research use only)
Precision Heuristics for Interpreting the Alteration Landscape (PHIAL) is a heuristic-based clinical interpretation algorithm for genomic tumor data that sorts somatic variants by clinical and biological relevance. 

## Getting PHIAL
The codebase is available for download through this Github repository, the [v1.0.0 release](https://github.com/vanallenlab/phial/releases), and [Dockerhub](https://hub.docker.com/r/vanallenlab/phial/). The code is also available to run as a method on [FireCloud](https://portal.firecloud.org/#methods/vanallenlab/phial/1).

To download via Github
```
git clone https://github.com/vanallenlab/phial.git
```

To download via Dockerhub
```
docker pull vanallenlab/phial:1.0.0
```

## Using PHIAL
PHIAL can be run using the wrapper script, `PHIAL_wrapper.R` by preparing a patient input file, such as [paitent_demo_1_input](https://github.com/vanallenlab/phial/blob/master/patient_demos/patient_demo_1/phial_patient_demo_1_input.txt), and running the following
```
Rscript PHIAL_wrapper.R -i patient_demos/patient_demo_1/phial_patient_demo_1_input.txt -d databases/databases.txt
```
To run PHIAL without the wrapper and without having to prepare an input file,
```
Rscript PHIAL_v1.0.R -i ${individual} -t ${tumorType} -o ${outDir} --mut.path ${snvHandle} --indel.path ${indelHandle} --segfile.path ${segHandle} --dranger.path ${drangerHandle} --actdb.mini databases/Actionable_genes_rationales_4.29.13.txt --actdb.large databases/CAdb_large_1_2012.txt --current_panel databases/current_panel_2011.txt --cosmic databases/CosmicHGNC_v56_151111.tsv --gsea.pathways databases/GSEA_cancer_gene_sets.txt --gsea.overlap databases/GSEA_cancer_actDB_overlap_sets.txt --gsea.modules databases/GSEA_cancer_modules.txt --refseq databases/refGene.hg19.20100825.sorted.txt
```

An example run command can be found in the [WDL](https://github.com/vanallenlab/phial/blob/master/firecloud.wdl) used to run [PHIAL on FireCloud](https://portal.firecloud.org/#methods/vanallenlab/phial/1).

## References
1. [Van Allen EM, Wagle N, Stojanov P, et al. Whole-exome sequencing and clinical interpretation of formalin-fixed, paraffin-embedded tumor samples to guide precision cancer medicine. Nat Med. 2014;20(6):682-8.](https://www.nature.com/articles/nm.3559)

# Disclaimer - For Research use only 
Please consult PHIAL_demo_README.pdf for requirements, input details, example usage, and description of outputs.

"This corresponds to the code and tables freeze linked to the Van Allen, Wagle et al Nature Medicine 2014 publication. Please stay tuned for future updates to this algorithm. DIAGNOSTIC AND CLINICAL USE PROHIBITED.  THE BROAD INSTITUTE and DFCI MAKE NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT or VALIDITY OF ANY INTELLECTUAL PROPERTY RIGHTS OR CLAIMS, WHETHER ISSUED OR PENDING, AND THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  In no event shall Broad or DFCI or their Trustees, Directors, Officers, Employees, Students, Affiliates, Core Faculty, Associate Faculty and Contractors, be liable for incidental, punitive, consequential or special damages, including economic damages or injury to persons or property or lost profits, regardless of whether the party was advised, had other reason to know or in fact knew of the possibility of the foregoing, regardless of fault, and regardless of legal theory or basis. You may not download or use any portion of this program for any use not expressly authorized by the Broad. You further agree that the program shall not be used as the basis of a commercial product and that the program shall not be rewritten or otherwise adapted to circumvent the need for obtaining permission for use of the program other than as specified herein."
