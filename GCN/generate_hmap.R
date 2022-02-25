library(optparse)
library(pheatmap)


option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              metavar="character"),
  make_option(c("-c","--cohort"), type="character", default=NULL, 
              metavar="character"),
  make_option(c("-n","--clust_num"), type="character", default=NULL, 
              metavar="character"),
  make_option(c("-o","--out_dir"), type="character", default=NULL, 
              metavar="character")
);

# option_list = list(
#   make_option(c("-n", "--name"), type="character", default=NULL, 
#               metavar="character")
# );

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# print(opt$file)
# print(opt$directory)
# print(opt$clust_num)
# print(opt$cohort)
# setwd('square_sum_results/results_skcm/interaction_figures')
print(file.path(opt$directory, opt$file))
pair_df <- read.csv(file.path(opt$directory, opt$file))
pair_mat <- as.matrix(pair_df[,2:dim(pair_df)[2]])
rownames(pair_mat) <- pair_df$X
setwd(opt$out_dir)
pair_hmap <- pheatmap(pair_mat,cluster_cols = F, cluster_rows = F, main = paste(opt$cohort, 'Cluster', opt$clust_num),
                      filename = paste0(opt$cohort,'_Cluster',opt$clust_num,'_Interactions.pdf'))



