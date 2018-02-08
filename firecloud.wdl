workflow PHIAL {
    String individual
    String tumorType
    File? snvHandle
    File? indelHandle
    File? segHandle
    File? drangerHandle

    call RunPHIAL {
        input:
	        individual=individual,
	        tumorType=tumorType,
	        snvHandle=snvHandle,
	        indelHandle=indelHandle,
	        segHandle=segHandle,
	        drangerHandle=drangerHandle
    }
}

task RunPHIAL {
    String individual
    String tumorType
    File? snvHandle
    File? indelHandle
    File? segHandle
    File? drangerHandle
    
    command <<<
        actdbMini="/databases/Actionable_genes_rationales_4.29.13.txt"
        actdbLarge="/databases/CAdb_large_1_2012.txt"
        currentPanel="/databases/current_panel_2011.txt"
        cosmic="/databases/CosmicHGNC_v56_151111.tsv"
        gseaPathways="/databases/GSEA_cancer_gene_sets.txt"
        gseaOverlap="/databases/GSEA_cancer_actDB_overlap_sets.txt"
        gseaModules="/databases/GSEA_cancer_modules.txt"
        refSeq="/databases/refGene.hg19.20100825.sorted.txt"

        Rscript /PHIAL_v1.0.R -i ${individual} -t ${tumorType} -o . \
        ${"--mut.path " + snvHandle} ${"--indel.path " + indelHandle} \
        ${"--segfile.path " + segHandle} ${"--dranger.path " + drangerHandle} \
        --actdb.mini $actdbMini --actdb.large $actdbLarge --current_panel $currentPanel --cosmic $cosmic \
        --gsea.pathways $gseaPathways --gsea.overlap $gseaOverlap --gsea.modules $gseaModules \
        --refseq $refSeq
    >>>

    output {
        File phialGel = "${individual}_phial_gel.png"
	    File phialClinicalRelevanceLow = "${individual}_investigate_clinical_relevance_low.txt"
	    File phialClinicalRelevanceHigh = "${individual}_investigate_clinical_relevance_high.txt"
	    File phialBiologicalRelevance = "${individual}_investigate_biological_relevance.txt"
	    File phialScoredDetailed = "${individual}_complete_muts_indels_scna_detailed.txt"
	    File phialScored = "${individual}_complete_muts_indels_scna.txt"
	    File phialReport = "${individual}_cancer_genome_report.html"
	    File phialReportRData = "${individual}_cancer_genome_report.RData"
	    File phialAFHistogram = "${individual}_allelicfx_hist.png"
    }

    runtime {
        docker: "vanallenlab/phial:1.0.0"
        memory: "4 GB"
    }
}
