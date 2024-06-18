#Install packages
install.packages('BiocManager')
BiocManager::install('clusterProfiler')
BiocManager::install("org.Dr.eg.db") #genome wide annotation for zebrafish

#GO Analysis using clusterprofiler package.
library(clusterProfiler)
library(org.Dr.eg.db)

#Get gene symbol ID from previous DEG analysis.
significant_genes <- subset(DEG2, p_val_adj < 0.05) #Get significant DEGs.
genes <- rownames(significant_genes)

#GO enrichment analysis
enrich.go <- enrichGO(gene = genes,
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',
                      ont = 'ALL',
                      pAdjustMethod = 'fdr',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = FALSE)
#Output results
write.csv(enrich.go, "wt_40hpf_earcells_GO.csv", row.names = FALSE)


#clusterProfiler default plotting
barplot(enrich.go)
dotplot(enrich.go)
cnetplot(enrich.go) #Network diagram showing the inclusion relationship between enriched functions and genes
emapplot(enrich.go) #Network diagram showing shared gene relationships between enriched functions
heatplot(enrich.go) #Heatmap showing the inclusion relationship between enriched functions and genes

