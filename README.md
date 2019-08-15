# Cross-Validated Rank-based Marker Selection (CVRMS)

CVRMS is an R package designed to extract marker subset from ranked marker dataset induced by genome-wide association study (GWAS) or marker effect analysis for genome-wide prediction.

## Requirement
 * R (>= 3.5.1)
 * pre-installed library in R
   - argparse
   - ggplot2
   - rrBLUP
   - ggpmisc
   - caret
 * How to install libraries (in R console)
   - install.packages(c("argparse", "ggplot2", "rrBLUP", "ggpmisc", "caret"))
   
 ## Quick example (in linux, ios, windows command console)
  * Rscript CVRMS_v1.3.R -g example_geno.txt -p example_pheno.txt -pn 1 -gw gwas_output.txt -gw_snp_cn 1 -gw_chr_cn 2 -gw_pos_cn 3 -gw_pv_cn 4 -min 10 -max 300 -cv 5 -a 0.9 -d 0.001 -ss 1 -m rrblup -t 1
  
## Detail options
 -g : genotypes input (row : markers, column : samples) - first row is the header, the first column includes marker names and the sample names should be included from the second column. !!Caution - genotypes are encoding to 0, 1, 2!!
 
 -p : phenotypes input (row : samples, column : phenotypes) - the first column includes the samples names and the traits of interests are from the second to end column)
 
 - 