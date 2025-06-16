# Bioinformatics-Lung-Adenocarcinoma

**Background**

Wanted to build don the results from the 2012 study, "The transcriptional landscape and mutational profile of lung adenocarcinoma' my seo et al. This study looked at differences in mutations across the patients in the study based on factors such as smoking status. The goal of my work was to use the genomic data from this study to explore the differences between SNPs on known cancer driving genes based on age of diagnoses (above or below 60 years old at diagnosis). For more information, refer to the 'Methods' section below or take a look at the .pdf file in this repository which is the report that was completed for this project.

**Data**

The patient data from the 2012 study, “The transcriptional landscape and mutational profile of lung adenocarcinoma” by Seo et al. including patient demographic information, the unhealthy tissue transcriptomes, and adjacent healthy tissue transcriptome data for each patient can be found from the SRA database through NCBI.

**pull_trim_map.sh**

The healthy and unhealthy transcriptome data for 2 patients above the age of 60 at diagnosis and 2 patients below the age of 60 at diagnosis are downloaded. The patients utilized were LC-C1, LC-C7, LC-S14, and LC-S20 who’s SRA values for healthy and unhealthy transcriptome were found from NCBI Gene Expression Omnibus (GEO) (http://www.ncbi.nlm.nih.gov/geo/) under accession number GSE40419. Following this, the reads from each transcriptome were trimmed to remove bases with quality scores <10 and length =>80. Following trimming, the reads for the healthy and unhealthy tissue trimmed transcriptomes were mapped with and aligned to the GRCh38.p14 transcriptome assembly as a reference using BBMAP .sh. We then sorted and indexed the resulting .bam file and ran analysis using bcftools to identify SNPs and short indel variations on the known cancer driving genes of interest (EGFR, KRAS, NRAS, PIK3CA, BRAF , CTNNB1, and MET) for the healthy and unhealthy transcriptomes. To help distinguish the somatic point mutations vs. germline mutations for each patient, we filtered out the mutations that appeared on the gene of interest in the healthy and unhealthy tissue transcriptomes. 

**variant_analysis.py**

Following this an organization of the mutational profiles on the genes of interest for each patient was done using the pyVCF library.
This data was organized to create the visual representations of the data shown in the figures within the results.
