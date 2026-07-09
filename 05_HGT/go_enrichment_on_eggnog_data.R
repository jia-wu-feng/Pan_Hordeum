if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
install.packages("here")

library(here)
library(clusterProfiler)
library(enrichplot)

# Load data (one orthogroup - GO term per row format)
ortho_dir <- here("..", "..", "..", "annot", "orthogroups")
eggnog_groups <- read.csv(here(ortho_dir, "all_groups_eggnog_go_exploded"), sep='\t', col.names = c("GO_TERM", "GENE"))
eggnog_hgt <- read.csv(here(ortho_dir, "hgt_genes_only_eggnog_exploded"), sep='\t', col.names = c("GO_TERM", "GENE"))

# Build maps
eggnog_map <- buildGOmap(eggnog_groups)
eggnog_hgt_map <- buildGOmap(eggnog_hgt)

# Different types of plots
eggnog_res <- enricher(gene=eggnog_hgt_map$GENE, TERM2GENE = eggnog_map)
# barplot(eggnog_res, showCategory = 10)
dotplot(eggnog_res, showCategory = 10)
emapplot(eggnog_res, showCategory = 10)
cnetplot(eggnog_res, showCategory=10, circular=TRUE)

# Try translating the GO numbers to GO terms
egg_res_df <- as.data.frame(eggnog_res)
# res_filtered <- egg_res_df[egg_res_df$p.adjust < 0.05, ] # not necessary, already filtered
egg_res_df$Description <- go2term(egg_res_df$ID)

# Getting an error when missing a term, try a custom function
safe_go2term <- function(go_ids) {
  sapply(go_ids, function(go_id) {
    tryCatch(
      go2term(go_id),  # Attempt to get the term
      error = function(e) "Unmatched"  # Assign "Unmatched" for errors
    )
  })
}
egg_res_df$Description <- safe_go2term(egg_res_df$ID)

go_ids <- egg_res_df$ID
go_terms <- go2term(go_ids)
go_ids[!go_ids %in% go_terms$go_id]
res_go_translated <- merge(egg_res_df, go_terms, by.x = "ID", by.y = "go_id", all.x = T)
res_go_translated$Description <- res_go_translated$Term
res_go_translated$Term <- NULL
