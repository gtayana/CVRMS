#!/data/lovemun/src_packages/R-3.4.0/bin/Rscript
#
#
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("SNP_info", nargs = 1, help = "SNP information table")

args <- parser$parse_args()
snpfile <- args$SNP_info
s_markers <- scan("Final_selected_markers.txt", what = "character")
snpinfo <- read.table(snpfile, header=T, sep='\t')

cat("Reading input argument", "\n")

load("Marker_selection_rrBLUP.RData")
cat("Loading Working Space File", "\n")
library(rrBLUP)
library(ggplot2)
library(ggpmisc)
library(progress)
val_geno <- readRDS("Validation_genotype.rds")
val_pheno <- readRDS("Validation_phenotype.rds")
geno2 <- geno2[names(phenotype1),s_markers]
cmd <- paste0("mkdir -p ./Validation")
system(cmd)
setwd(paste0("./Validation/"))
all_train_acc <- NULL
e_mat <- NULL
beta <- NULL

for (j in 1:50){
  cat("Iteration Number ", j, "\n")
  train_samples <- sample(1:nrow(geno2), size = nrow(geno2)*0.5, replace = FALSE)
  train_geno <- geno2[train_samples,]
  test_geno <- geno2[-train_samples,]
  train_pheno <- phenotype1[train_samples]
  test_pheno <- phenotype1[-train_samples]
  train_BLUP <- mixed.solve(y = train_pheno, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
  e_mat <- cbind(e_mat, train_BLUP$u)
  beta <- c(beta, train_BLUP$beta)
  train_e = rowMeans(e_mat)
  train_beta <- mean(beta)
  pheno_valid = data.matrix(test_geno) %*% as.matrix(train_e)
  pred_pheno <- pheno_valid + c(train_beta)
  train_acc <- cor(pred_pheno, test_pheno, use = "complete")
  all_train_acc <- c(all_train_acc, as.vector(train_acc))
  cat("Correlation of j-th cross validation step is ", train_acc, "\n")
}
cat("Final model is generated", "\n")
pheno_valid = data.matrix(val_geno[,s_markers]) %*% train_e
pred_pheno <- pheno_valid + c(train_beta)
df1 <- data.frame(Phenotype = val_pheno, GEBV = pred_pheno)
my.formula <- y~x
g0 <- ggplot(df1, aes(x=Phenotype, y = GEBV)) + stat_smooth(method = "lm", se = FALSE, color = "black", formula = my.formula) + 
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + geom_point() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
ggsave(filename = "Test_prediction_GEBV.png", plot = g0, dpi = 600)
cat("Correlation for test dataset is ", cor(pred_pheno, val_pheno, use = "complete"), "\n")
cat("Average correlation for 50 replications is ", mean(all_train_acc), "\n")
cat("Standard deviation for 50 replications is ", sd(all_train_acc), "\n")
rownames(e_mat) <- s_markers
colnames(e_mat) <- 1:50
sel_ix <- match(s_markers, snpinfo$Name)
df2 <- data.frame(SNP = s_markers, Chromosome = snpinfo$Chromosome[sel_ix], Position = snpinfo$Position[sel_ix])
df3 <- as.data.frame(table(factor(df2$Chromosome, levels = c(1:max(snpinfo$Chromosome)))))
colnames(df3) <- c("Chromosome", "freq")
g2 <- ggplot(df3, aes(Chromosome, freq)) +
  geom_bar(stat = "identity", width = .6, fill = "tomato2") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
  geom_text(aes(label = freq), vjust = -0.3, size = 3.5) + theme_classic()
ggsave(filename = "Final_selected_markers_chr_dist.png", plot = g2, dpi = 600)
write.csv(df2, file = "Final_selected_markers_info.csv", row.names=F, quote=F)

write.table(all_train_acc, file = "Training_accuracy.txt", row.names=F)
write.csv(rowMeans(e_mat), file = "e_matrix.csv", row.names=F, quote=F)
write.table(beta, file = "beta.txt", row.names=F, sep='\t')
