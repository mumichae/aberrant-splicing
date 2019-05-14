blood_prokisch <- readRDS("/s/project/crg_seq_data/processed_results/v29_overlap/datasets/blood_prokisch/counts_blood_prokisch.Rds")
bt <- data.table(sampleID = colnames(blood_prokisch), condition = colnames(blood_prokisch))
bt[, bamFile := paste0("/s/project/mitoMultiOmics/raw_data/helmholtz/", sampleID, "/RNAout/paired-endout/stdFilenames/", sampleID, ".bam")]

