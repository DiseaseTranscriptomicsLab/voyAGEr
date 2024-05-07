This repository contains the corrected logCPM valoes for each tissue analysed in voyAGEr. Data preprocessing followed the methods outlined in voyAGEr's original paper:

> Arthur L. Schneider, Rita Martins-Silva, Alexandre Kaizeler, Nuno Saraiva-Agostinho, Nuno L. Barbosa-Morais. (2023). voyAGEr: free web interface for the analysis of age-related gene expression alterations in human tissues. eLife, 12:RP88623.

In summary, the RNA-seq read count matrix for each gene in GTEx v8 samples was obtained from the project’s data portal (https://www.gtexportal.org/). 
Out of the 54 available tissues, five with fewer than 50 samples (kidney medulla, fallopian tube, bladder, ectocervix, and endocervix) were excluded.

Read count data for each tissue were pre-processed individually. Genes with very low expression across samples (fewer than 1 CPM in fewer than 40% of samples) were filtered out, 
resulting in the final number of genes used for the analysis.

Normalisation factors were calculated using ```edgeR```'s ```calcNormFactors``` function, implementing the trimmed mean of M-values, and read counts were subsequently normalised and log-transformed 
with the ```voom``` function from the ```limma``` package.

Principal component analysis was conducted for each tissue, and potential batch effects were identified based on empirically defined thresholds. The **COHORT** variable was identified as a primary batch 
effect. Additionally, **SMRIN**, **DTHHRDY**, **MHSMKYRS**, and the number of detected genes per sample contributed significantly to data variance. 
These conditions were corrected for, on a tissue-by-tissue basis using linear models and an adaptation of the ```removeBatchEffect``` function from the limma package. Imputation was performed for missing 
values using the ```mice``` package.

The metadata associated with this dataset is available for download from the GTEx portal 
(https://gtexportal.org/home/downloads/adult-gtex/metadata), or can be accessed through dbGaP with restricted access. Requests to access this data can be made following the guidelines provided on the 
GTEx portal (https://gtexportal.org/home/protectedDataAccess).

