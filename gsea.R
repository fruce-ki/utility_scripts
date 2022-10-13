library(fgsea)
library(grid)
library(ggplot2)
library(data.table)


# Broad Institute Database
msigdb <- "/Volumes/groups/busslinger/Bioinf/src/cbe/broad/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs"
xorg <- "/Volumes/groups/busslinger/Bioinf/src/cbe/broad/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.5.1.chip"

dir(msigdb)
SETDB <- fgsea::gmtPathways(file.path(msigdb, 'c5.go.mf.v7.5.1.symbols.gmt'))


# Human orthologs
XREF <- fread(xorg)
setnames(XREF, c('mm', 'hs', 'desc'))
setkey(XREF, mm)


# Data
DE <- fread('/Volumes/groups/busslinger/Kimon/tanja/R13870_RNAseq_timecourse/results/DE/exon_genecounts.CellType_actB.Condition_Ikaros_00_5phIAA_vs_Aiolos_00_5phIAA.deseq2.nolab.spotfire.txt')
setorder(DE, -log2FoldChange.shrink, na.last=TRUE)
genes.rnk <- DE$log2FoldChange.shrink
names(genes.rnk) <- XREF[DE$name, hs]
names(genes.rnk)[is.na(names(genes.rnk))] <- DE$name[is.na(names(genes.rnk))]


# Wrapper for {fgsea}, https://bioinformaticsbreakdown.com/how-to-gsea/
GSEA = function(gene_list, myGO, pval=0.05) {
	set.seed(54321)
	library(dplyr)
	library(fgsea)
	
	if ( any( duplicated(names(gene_list)) )  ) {
		warning("Duplicates in gene names")
		gene_list = gene_list[!duplicated(names(gene_list))]
	}
	if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
		warning("Gene list not sorted")
		gene_list = sort(gene_list, decreasing = TRUE)
	}
	
	fgRes <- fgsea::fgsea(pathways = myGO,
												stats = gene_list,
												minSize=15, ## minimum gene set size
												maxSize=400 ) %>% # , ## maximum gene set size
												# nperm=10000) %>%
		as.data.frame() %>% 
		dplyr::filter(padj < !!pval) %>% 
		arrange(desc(NES))
	message(paste("Number of signficant gene sets =", nrow(fgRes)))
	
	message("Collapsing Pathways -----")
	concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
																			pathways = myGO,
																			stats = gene_list)
	fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
	message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
	
	fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
	filtRes = rbind(head(fgRes, n = 10),
									tail(fgRes, n = 10 ))
	
	total_up = sum(fgRes$Enrichment == "Up-regulated")
	total_down = sum(fgRes$Enrichment == "Down-regulated")
	header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
	
	colos = setNames(c("firebrick2", "dodgerblue2"),
									 c("Up-regulated", "Down-regulated"))
	
	g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
		geom_point( aes(fill = Enrichment, size = size), shape=21) +
		scale_fill_manual(values = colos ) +
		scale_size_continuous(range = c(2,10)) +
		geom_hline(yintercept = 0) +
		coord_flip() +
		labs(x="Pathway", y="Normalized Enrichment Score",
				 title=header)
	
	output = list("Results" = fgRes, "Plot" = g1)
	return(output)
}


# GSEA
res <- GSEA(genes.rnk, SETDB, 0.05)
res$Results <- as.data.table(res$Results)
print(dim(res$Results))
print(res$Plot)


# Cluster gene sets, https://bioinformaticsbreakdown.com/how-to-gsea/
plot_geneset_clusters = function(gs_results, myGO, min.sz = 4, main="GSEA clusters"){
	library(ggplot2)
	library(ggrepel)
	library(stringr)
	
	df = matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
	rownames(df) = colnames(df) = gs_results$pathway
	
	for ( i in 1:nrow(gs_results)) {
		genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
		for (j in 1:nrow(gs_results)) {
			genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
			## Jaccards distance  1 - (intersection / union )
			overlap = sum(!is.na(match(genesI, genesJ )))
			jaccards = overlap / length(unique(c(genesI, genesJ) ))
			df[i,j] = 1-jaccards
		}
	}
	
	## Cluster nodes using dynamic tree cut, for colors
	distMat = as.dist(df)
	dendro = hclust(distMat, method = "average" )
	clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
	## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it
	
	gs_results$Cluster = clust
	gs_results = gs_results[gs_results$Cluster != 0, ]
	
	## select gene sets to label for each clusters
	bests = gs_results %>%  
		group_by( Cluster ) %>% 
		top_n(wt = abs(size), n = 1) %>% 
		.$pathway
	## determine cluster order for plotting
	clust_ords = gs_results %>% 
		group_by( Cluster ) %>% 
		summarise("Average" = NES ) %>% 
		arrange(desc(Average)) %>% 
		.$Cluster %>% 
		unique
	
	gs_results$Cluster = factor(gs_results$Cluster, levels = clust_ords)
	
	gs_results$Label = ""
	gs_results$Label[gs_results$pathway %in% bests ] = gs_results$pathway[gs_results$pathway %in% bests ]
	gs_results$Label = str_remove(gs_results$Label, "GO_")
	gs_results$Label = tolower(gs_results$Label)
	
	g1 = ggplot(gs_results, aes(x = Cluster, y = NES, label = Label )) +
		geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
		scale_size_continuous(range = c(0.5,5)) +
		geom_text_repel( force = 2, max.overlaps = Inf) +
		ggtitle(main)
	
	return(g1)
}


#
print( plot_geneset_clusters( gs_results = res$Results[res$Results$NES > 0, ], 
											 main = "Up-regulated GSEA clusters",
											 GO_file = GO,
											 min.sz = 4 ) )

plot_geneset_clusters( gs_results = res$Results[res$Results$NES < 0, ], 
											 main = "Down-regulated GSEA clusters",
											 GO_file = GO,
											 min.sz = 2 )


# Enrichment plot of one cluster
plotEnrichment(SETDB[["GOMF_SINGLE_STRANDED_DNA_BINDING"]],
							 genes.rnk) + labs(title="GOMF_SINGLE_STRANDED_DNA_BINDING")


# Table of enrichment plots for a bunch of selected pathways
topPathwaysUp <- res$Results[ES > 0, ][head(order(pval), n=10), pathway]
topPathwaysDown <- res$Results[ES < 0, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
grid.newpage(); par(mfg=c(1,1)) # Prevent overlay with whatever previous plot
plotGseaTable(SETDB[topPathways], genes.rnk, res$Results, gseaParam=0.5)


# Select only independent pathways
collapsedPathways <- collapsePathways(res$Results[order(pval)][padj < 0.01], SETDB, genes.rnk)
mainPathways <- res$Results[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
grid.newpage(); par(mfg=c(1,1)) # Prevent overlay with whatever previous plot
plotGseaTable(SETDB[mainPathways], genes.rnk, res$Results, gseaParam = 0.5)



