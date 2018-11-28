##Meth_zebrafish_germline

This repository is intended to provide the code used for the bioinformatic analysis for the project Zebrafish preserve germline epigenetic memory globally but demethylate and amplify sex-linked rDNA during feminisation.

The code is implemented in R (v3.4.4) and the following libraries:
- ggbio v1.26.1
- rtracklayer v1.38.3
- biovizBase v 1.26.0
- BSgenome.Drerio.UCSC.danRer10 v1.4.2
- dplyr 0.7.6
- reshape 0.8.7
- ggplot2 3.0.0

#Components of the analysis

- Grand plot
Probes of 1 Mb through the zebrafish genome were created in SeqMonk (v1.43.0)to quantify overrepresentated regions thought the genome. This code plots the results obtained.

- Ideogram
This code plots an ideogram highlighting the peak obtained previously.

- Coverage plot
BAM files obtained by mapping reads wit Bowtie2 v(2.3.2) as mention in the manuscript were used for this analysis. This code plots the reads mapped to the genomic region of interest.

- Wilson data fem-rDNA
This code plots the location of the amplified region and its relation with the sex linked SNPs obtained in a previous study (Wilson et al, 2014)

- Relation amplification - methylation
This code 

- Summary amplification