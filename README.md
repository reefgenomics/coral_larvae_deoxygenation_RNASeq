# coral_larvae_deoxygenation_RNASeq

Project Summary
We experimentally exposed planula larave of Acropora selago, to 12 h continuous deoxygenation (<2 mg O2 L-1) followed by 12 h re-oxygenation (~6 mg O2 L-1) via their natural night-day light cycle and followed gene expression changes.

Published manuscript (open access): https://onlinelibrary.wiley.com/doi/10.1111/mec.16259

RNA-Seq data @ NCBI: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA723188 


## Workflow
1. The script `coral_larvae_deoxygenation_RNASeq.sh` was used to produce transcript count tables
2. Differential expression analysis was done using the script `coral_larvae_deoxygenation_DESeq.R`
3. PCA plots were produced using the script `coral_larvae_deoxygenation_PCA_plots.R`

The scripts are taken from: https://github.com/ajcardenasb/coral_deoxygenation_RNASeq
