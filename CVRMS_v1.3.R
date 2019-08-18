#!/data/lovemun/src_packages/R-3.4.0/bin/Rscript
# created by Seongmun Jeong
# Contact : likemun@gmail.com

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-g", "--geno", default = NULL, help = "genotype file")
parser$add_argument("-p", "--pheno", default = NULL, help = "phenotype file")
parser$add_argument("-pn", "--pheno_num", type = "integer", default = 1, help="Phenotype column number")
parser$add_argument("-gw", "--gwas", default = NULL, help = "GWAS output")
parser$add_argument("-gw_snp_cn", "--gw_snp_column_number", type = "integer", default = 2, help = "Column number of SNP in GWAS results")
parser$add_argument("-gw_chr_cn", "--gw_chr_column_number", type = "integer", default = 1, help = "Column number of chromosome of SNP in GWAS results")
parser$add_argument("-gw_pos_cn", "--gw_position_column_number", type = "integer", default = 3, help = "Column number of position of SNP in GWAS results")
parser$add_argument("-gw_pv_cn", "--gw_pvalue_column_number", type = "integer", default = 9, help = "Column number of p values in GWAS results")
parser$add_argument("-min", "--min_snps", type = "integer", default = 10, help = "Minimum number of SNPs")
parser$add_argument("-max", "--max_snps", type = "integer", default = 20, help = "Maximum number of SNPs")
parser$add_argument("-cv", "--cv", type = "integer", default = 10, help = "Cross validation")
parser$add_argument("-a", "--acc", type = "double", default = 0.9, help = "Goal of correlation")
parser$add_argument("-d", "--delta", type = "double", default = 0.001, help = "Increasing rate for cut-off")
parser$add_argument("-ss", "--sel_snps", type = "integer", default = 1, help = "# of selected snps for each epoch")
parser$add_argument("-m", "--method", default = "rrblup", help = "prediction method (rrblup, rf)")
parser$add_argument("-t", "--time", type = "double", default = 1, help = "Runtime")

## Load input argument
args <- parser$parse_args()
if (is.null(args$geno)){
    stop("Genotype input file is null")
}
if (is.null(args$pheno)){
    stop("phenotype input file is null")
}
genofile <- args$geno
phenofile <- args$pheno
gwasfile <- args$gwas
snp_cn <- args$gw_snp_column_number
chr_cn <- args$gw_chr_column_number
pos_cn <- args$gw_position_column_number
pv_cn <- args$gw_pvalue_column_number
pheno_num_input <- as.numeric(args$pheno_num)
min_num <- as.numeric(args$min_snps)
max_num <- as.numeric(args$max_snps)
cv <- as.numeric(args$cv)
delta <- as.numeric(args$delta)
sel_snps <- as.numeric(args$sel_snps)
acc1 <- as.numeric(args$acc)
mm <- args$method
time_cutoff <- args$time

cat("Reading input argument", "\n")
cat("Genotype file is ", genofile, "\n")
cat("Phenotype file is ", phenofile, "\n")
cat("Column number of phenotype is ", pheno_num_input, "\n")
cat("GWAS result file is ", gwasfile, "\n")
cat("Mininum number of selected SNPs is ", min_num, "\n")
cat("Cross validation k is ", cv, "\n")
cat("Accuracy goal is ", acc1, "\n")
cat("Delta is ", delta, "\n")
cat("The number of selected SNPs for each step is ", sel_snps, "\n")
cat("Limit of runtime is ", time_cutoff, "\n")

## Load input datasets
geno <- read.table(genofile, header=T, check.names=F)
rownames(geno) <- geno[,1]
geno <- data.matrix(geno[,-1])
geno <- t(geno)
pheno <- read.table(phenofile, header=T, check.names=F)
rownames(pheno) <- pheno[,1]
pheno <- data.matrix(pheno[,-1])
common_samples <- intersect(rownames(geno), rownames(pheno))
geno1 <- geno[common_samples,]
pheno1 <- pheno[common_samples,]
rm(geno); rm(pheno)
pheno_num <- pheno_num_input
cat("Column number of phenotype is ", pheno_num, '\n')

## Load associated packages
library(rrBLUP)
library(caret)
library(ggplot2)
library(ggpmisc)
library(progress)

## Separate validation samples from whole datasets
phenotype <- as.vector(pheno1[,pheno_num])
names(phenotype) <- rownames(pheno1)
val_samples <- sample(1:nrow(geno1), size = max(nrow(geno1)*0.1, 100), replace = FALSE)
val_geno <- geno1[val_samples,]
val_pheno <- phenotype[val_samples]
geno2 <- geno1[-val_samples,]
phenotype1 <- phenotype[-val_samples]
all_train_acc <- NULL
cat("Trait name is ", colnames(pheno1)[pheno_num], '\n')
if (is.null(gwasfile)){
    cmd <- paste0("mkdir -p ", colnames(pheno1)[pheno_num], "_wo_gwas")
    w_dir <- paste0(colnames(pheno1)[pheno_num], "_wo_gwas/")
} else {
    cmd <- paste0("mkdir -p ", colnames(pheno1)[pheno_num], "_with_gwas")
    w_dir <- paste0(colnames(pheno1)[pheno_num], "_with_gwas/")
}
system(cmd)

if (is.null(gwasfile)){
    if (mm == "rrblup"){
        BLUP = mixed.solve(y = phenotype1, Z = geno2, K = NULL, SE = FALSE, return.Hinv = FALSE)
        marker_effects = BLUP$u
    }
    if (mm == "rf"){
        ttMod <- train(x = geno2, y = phenotype1, method = "rf") 
        marker_effects = varImp(ttMod)$importance
    }
    gwas_results <- data.frame(SNP = names(marker_effects), P = marker_effects)
    pv_cn <- 2
    write.csv(gwas_results, file = "Marker_effects.csv", quote=F, row.names=F)
} else {
    if (grepl(pattern = ".csv", x = gwasfile)){
        gwas_results <- read.csv(gwasfile, header=T, check.names=F)
    } else {
        gwas_results <- read.table(gwasfile, header=T, check.names=F)
    }
}
gwas_results <- gwas_results[order(gwas_results[,pv_cn]),]
ix <- as.vector(gwas_results[,snp_cn])

setwd(w_dir)
saveRDS(val_geno, file = "Validation_genotype.rds", compress = "gzip")
saveRDS(val_pheno, file = "Validation_phenotype.rds", compress = "gzip")
#write.csv(val_geno, file = paste0("Validation_genotype.csv"), quote=F)
#write.csv(val_pheno, file = paste0("Validation_phenotype.csv"), quote=F)

## Cross validation : divide samples
cv_samples <- sample(1:cv, nrow(geno2), replace = TRUE)

## Set progress bar
#pb <- progress_bar$new(total = cv)

for (j in 1:cv){
    o_delta <- 1
    cat("CV Number ", j, "\n")
    #pb$tick()
    
    train_samples <- which(cv_samples != j)
    train_geno <- geno2[train_samples,]
    test_geno <- geno2[-train_samples,]
    train_pheno <- phenotype1[train_samples]
    test_pheno <- phenotype1[-train_samples]
    if (mm == "rrblup"){
        train_BLUP <- mixed.solve(y = train_pheno, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
        train_e = as.matrix(train_BLUP$u)
        pheno_valid = test_geno %*% train_e
        pred_pheno <- pheno_valid + c(train_BLUP$beta)
    }
    if (mm == "rf"){
        ttMod <- train(x = train_geno, y = train_pheno, method = "rf")
        pred_pheno <- predict(ttMod, newdata = test_geno)
    }
    train_acc <- cor(pred_pheno, test_pheno, use = "complete")
    all_train_acc <- c(all_train_acc, as.vector(train_acc))
    cat("Training samples correlation for all markers is ", train_acc, "\n")

    ix1 <- ix
    if (sel_snps > 1){
        init_set <- ix1[1:sel_snps]
        cor_mat <- cor(train_geno[,init_set], method = "spearman")
        cor_mat[upper.tri(cor_mat)] <- 0
        diag(cor_mat) <- 0
        train_ix <- init_set[!apply(cor_mat, 2, function(x) any(x > 0.9))]
    } else {
        train_ix <- ix1[1:2]
    }
    except_ix <- NULL
    add_ix <- NULL
    result_acc <- 0 
    ms_out <- list()
    train_itr_geno <- train_geno[,train_ix]
    test_itr_geno <- test_geno[,train_ix]
    train_itr_BLUP <- mixed.solve(y = train_pheno, Z = train_itr_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
    train_itr_e = as.matrix(train_itr_BLUP$u)
    pheno_itr_valid = test_itr_geno %*% train_itr_e
    pred_itr_pheno = pheno_itr_valid + c(train_itr_BLUP$beta)
    train_itr_acc <- cor(pred_itr_pheno, test_pheno, use = "complete")
    if (sel_snps == 1){
    	ms_out[[1]] <- c(0, train_ix[1])
    	ms_out[[2]] <- c(train_itr_acc, train_ix[2])
    	k=3
    } else {
    	ms_out[[1]] <- c(train_itr_acc, train_ix)
    	k=2
    }
    i=1

    cat("Correlation of ", j, "-th iteration training is ", train_itr_acc, "\n")
    result_acc <- c(result_acc, train_itr_acc)
    #pb <- progress_bar$new(total = ncol(geno2))
    max_noid <- 0
    s_time <- Sys.time()
    cat("Start marker selection for ", j, "-th fold at ", format(s_time, usetz = T), "\n")
    ld_count <- 0
    while (train_itr_acc < min(train_acc*5, acc1) | length(train_ix) < min_num | o_delta < delta){
        if (length(train_ix) == max_num | ld_count == 20 | difftime(Sys.time(), s_time, units = "hours") > time_cutoff){
            cat("# of selected snps is ", length(train_ix), ", and # of low delta count is ", ld_count, "\n")
            cat("Time difference is ", Sys.time() - s_time, "\n")
            break
        }
        #pb$tick()
        if (sel_snps == 1){
        	  add_snps <- ix1[i+2]
            new_train_ix <- c(train_ix, add_snps)
        } else {
            init_set <- ix1[(sel_snps*i+1):(sel_snps*(i+1))]
            cor_mat <- cor(train_geno[,init_set], method = "spearman")
            cor_mat[upper.tri(cor_mat)] <- 0
            diag(cor_mat) <- 0
            add_snps <- init_set[!apply(cor_mat, 2, function(x) any(x > 0.9))]
            new_train_ix <- c(train_ix, add_snps)
        }
        train_itr_geno <- train_geno[,new_train_ix]
        test_itr_geno <- test_geno[,new_train_ix]
        if (mm == "rrblup"){
            train_itr_BLUP <- mixed.solve(y = train_pheno, Z = train_itr_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
            train_itr_e = as.matrix(train_itr_BLUP$u)
            pheno_itr_valid = test_itr_geno %*% train_itr_e
            pred_itr_pheno = pheno_itr_valid + c(train_itr_BLUP$beta)
        }
        if (mm == "rf"){
            ttMod <- train(x = train_itr_geno, y = train_pheno, method = "rf")
            pred_itr_pheno <- predict(ttMod, newdata = test_itr_geno)
        }
        acc_test <- cor(pred_itr_pheno, test_pheno, use = "complete")
        if (round(acc_test - train_itr_acc, digits = 4) <= 0){
            except_ix <- c(except_ix, add_snps)
            max_noid <- max_noid + 1
        } else {
            ms_out[[k]] <- c(acc_test, add_snps)
            train_ix <- c(train_ix, add_snps)
            train_itr_acc <- acc_test
            add_ix <- c(add_ix, i)
            result_acc <- c(result_acc, train_itr_acc)
            o_delta <- result_acc[length(result_acc)] - result_acc[length(result_acc)-1]
            cat("Add ", train_itr_acc,", o_delta: ", o_delta, ", time: ", format(difftime(Sys.time(), s_time), usetz = TRUE), "\n")
            max_noid <- 0
        }
        if (i == length(ix1)/sel_snps){
            ix1 <- ix1[-c(add_ix, except_ix)]
            add_ix <- NULL
            except_ix <- NULL
            i <- 1
        } else {
            i <- i+1
        }
        k = k+1
    }
    cat("# of selected markers is ", length(train_ix), " and correlation is ", train_itr_acc, "for ", j, "-th iteration", "\n")
    save.image("temp.RData")

    a <- plyr::ldply(ms_out, rbind)
    count <- length(which(!is.na(as.vector(a[1,-1]))))
    index <- count
    for (i in 2:nrow(a)){
        count <- count + length(which(!is.na(as.vector(a[i,-1]))))
        index <- c(index, count)
    }
    a <- data.frame(Index = index, a)
    colnames(a)[c(2,3)] <- c("Corr", "Marker")
    a$Corr <- as.numeric(as.vector(a$Corr))
    write.table(a, file = paste0("Iteration_", j, ".txt"), row.names = F, sep='\t', quote=F)

    df1 <- data.frame(Phenotype = test_pheno, GEBV = pred_itr_pheno)
    if (!is.null(gwasfile)){
    	sel_ix <- match(train_ix, gwas_results[,snp_cn])
    	df2 <- data.frame(SNP = train_ix, Chromosome = gwas_results[sel_ix, chr_cn], Position = gwas_results[sel_ix, pos_cn],
    										logP = log10(gwas_results[sel_ix, pv_cn]), effects = train_itr_BLUP$u[train_ix])
    	df3 <- as.data.frame(table(factor(df2$Chromosome, levels = c(1:max(gwas_results[,chr_cn])))))
    	colnames(df3) <- c("Chromosome", "freq")
    	g2 <- ggplot(df3, aes(Chromosome, freq)) +
    		geom_bar(stat = "identity", width = .6, fill = "tomato2") +
    		theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
    		geom_text(aes(label = freq), vjust = -0.3, size = 3.5) + theme_classic()
    	ggsave(filename = paste0("Iteration_", j, "_Chrom_dist.png"), plot = g2, dpi = 600)
    }
    my.formula <- y~x
    g0 <- ggplot(df1, aes(x=Phenotype, y = GEBV)) + stat_smooth(method = "lm", se = FALSE, color = "black", formula = my.formula) + 
        stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + geom_point() +
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())
    ggsave(filename = paste0("Iteration_", j, "_GEBV.png"), plot = g0, dpi = 600)

    cutoff <- data.frame(x=c(-Inf, Inf), y = train_acc, cutoff = factor(train_acc))
    g1 <- ggplot(a, aes(x=Index, y = Corr)) + geom_line() + geom_smooth() + xlab("Index") + ylab("Correlation") + ylim(c(0, 0.9)) +
        theme_bw() + geom_line(aes(x, y), cutoff) +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) 
    ggsave(filename = paste0("Iteration_", j, "_correlation.png"), plot = g1, dpi = 600)
}
save.image(file = paste0("Marker_selection_rrBLUP.RData"))
files <- list.files(path = "./", pattern = "Iteration_")
files <- files[grepl(pattern = ".txt", x = files)]
marker_list <- list()
for (i in 1:length(files)){
    tt <- read.table(files[i], header=T, sep='\t')
    marker_list[[i]] <- tt$Marker
}
at <- unlist(marker_list)
at_table <- table(at)[order(table(at), decreasing = TRUE)]
write(names(at_table)[1:ceiling(length(at)/cv)], file = "Final_selected_markers.txt")
write(all_train_acc, file = "Final_markers_accuracy.txt")
