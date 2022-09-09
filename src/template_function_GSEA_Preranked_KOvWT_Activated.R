GSEA_Preranked_KOvWT_Activated <- function(DEG_Analysis, msigdb_v6_2_with_orthologs) {

graphics.off()

# LIBRARIES ====

library(dplyr); library(fgsea); library(grid); library(gridExtra); library(gtable); library(patchwork); library(data.table)

# INPUT PARAMS ====

## dataset

geneset_db = msigdb_v6_2_with_orthologs # gene set collection is pinned (can be changed but assumes presence of columns named: "species", "collection", "gene_set_name", "gene_symbol")
genescore_df = DEG_Analysis
    
## ranking
genescore = "_tstat"
genescore_alternative = c()
geneid = "Gene"
filterInput_byContrast = 'none'
contrasts = c()    
    
## GSEA
collections_to_include = c("H: hallmark gene sets") # object Pathway/Collection lookup
FDR_correction_mode = "within each collection"
species = "Mouse"
minimum_gene_set_size = 15
maximum_gene_set_size = 500
number_of_permutations = 5000
        
## output
collapse_redundant_pathways = FALSE
sortBy = c()
sortDecreasing = FALSE
       
## visualization
image_width = 2500
image_height = 2500
image_resolution = 300       
     
    
# INPUT HANDLING AND FILTERING ====    
    
## pathway collection    
geneset_db <- geneset_db %>% filter(geneset_db[["species"]]==species) %>% filter(`%in%` (geneset_db[["collection"]], collections_to_include)) 
db_unique <- geneset_db %>% dplyr::select("gene_set_name") %>% dplyr::distinct()
db_selected <- geneset_db %>% dplyr::select("collection","gene_set_name") %>% dplyr::distinct()
db_isDuplicated <- dplyr::count(db_selected) > dplyr::count(db_unique)

if (db_isDuplicated) {

    db_selected <- dplyr::collect(db_selected)
    within_collection <-  db_selected %>% dplyr::group_by(collection) %>% dplyr::filter(duplicated(gene_set_name)) %>% dplyr::ungroup() %>% dplyr::count()
     
     if (within_collection == 0) {        
        stop("ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter")
    
    } else if (within_collection > 0) {        
        between_collection <- Reduce("intersect", split(db_selected$gene_set_name, db_selected$collection)) %>% length()
        
        if (between_collection == 0) {
            stop("ERROR: duplicated gene set names found in the 'Gene set database' due to not unique gene set names within a collection")
        
        } else if (between_collection > 0) {
            "ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter and not unique gene set names within a collection"
        } 
    } 
        
} else if (!db_isDuplicated) {    
        geneset_db <- dplyr::collect(geneset_db)     
}

geneset_db <- geneset_db %>% dplyr::group_by(collection, gene_set_name) %>% dplyr::summarize(gene_symbol = as.list(strsplit(paste0(unique(gene_symbol), collapse = " "), " "))) %>% dplyr::ungroup()
geneset_list = geneset_db$gene_symbol
names(geneset_list) = geneset_db$gene_set_name

## ranking
if (!is.null(genescore_alternative)) { genescore <- genescore_alternative }
rank_columns = colnames(genescore_df)[grepl(paste0("\\Q", genescore, "\\E$"), colnames(genescore_df))]
rank_contrasts = unlist(strsplit(rank_columns, genescore))

if ( filterInput_byContrast == "remove" ) { 

        if ( ! is.null(contrasts) ) {

            all_contrasts = rank_contrasts
            index = match(contrasts, rank_contrasts)
            rank_columns = rank_columns[-index]
            rank_contrasts = rank_contrasts[-index]
            removed = setdiff(all_contrasts, rank_contrasts)
            if (length(removed) < 1) {
                cat(sprintf('WARNING:contrast(s) to remove (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
            
            } else {
                cat(sprintf("Removed contrast(s): %s\n", paste(removed, collapse=", ")))
                cat(sprintf("Kept contrast(s): %s\n", paste(rank_contrasts, collapse=", ")))
            }

        } else if (is.null(contrasts) ) {
            cat(sprintf('WARNING:contrast(s) to remove (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
        }

    } else if ( filterInput_byContrast == "keep" ) {

        if ( ! is.null(contrasts) ) {

            all_contrasts = rank_contrasts
            index = match(contrasts, rank_contrasts)
            rank_columns = rank_columns[index]
            rank_contrasts = rank_contrasts[index]
            removed = setdiff(all_contrasts, rank_contrasts)
            if (length(rank_contrasts) < 1) {
                cat(sprintf('WARNING:contrast(s) to keep (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
            
            } else {
                cat(sprintf("Removed contrast(s): %s\n", paste(removed, collapse=", ")))
                cat(sprintf("Kept contrast(s): %s\n", paste(rank_contrasts, collapse=", ")))
            }

        } else if (is.null(contrasts) ) {
            cat(sprintf('WARNING:contrast(s) to keep (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
        }

    } else if ( filterInput_byContrast == "none" ) {

        if ( ! is.null(contrasts) ) {
            cat(sprintf('WARNING:contrast filter not specified correctly; filter not applied\nIdentified contrast(s) used: %s\n', paste(rank_contrasts, collapse=", ")))
        } else {
            cat(sprintf('Filter contrast ("none"); Identified contrast(s) used: %s\n', paste(rank_contrasts, collapse=", ")))
        }
    }

genescore_df <- genescore_df %>% dplyr::select(geneid, rank_columns) %>% tidyr::pivot_longer(!geneid, names_to="contrast", values_to="genescores", values_drop_na=TRUE) %>% dplyr::rename("gene_id"=geneid) %>% dplyr::mutate(contrast=sub(genescore, "", contrast))
duplicates <- dplyr::group_by(genescore_df, contrast, gene_id) %>% dplyr::filter(dplyr::n()>1)
if (nrow(duplicates) > 0 ) {
    genescore_grouped <- genescore_df %>% dplyr::group_by(contrast,gene_id) %>% dplyr::summarize(genescores=mean(genescores, na.rm=TRUE))
    cat(sprintf("WARNING: duplicated gene names found of %g gene(s), duplicated values of gene scores were averaged per gene", length(unique(duplicates$gene_id))))
} else {
    genescore_grouped <- dplyr::group_by(genescore_df, contrast)    
}
        

# ANALYSIS ====

## GSEA    
gsea <- dplyr::group_modify( genescore_grouped, ~run.gsea(., db=geneset_db, collections=geneset_list, minimum_size=minimum_gene_set_size, maximum_size=maximum_gene_set_size, number_perms=number_of_permutations, mode=FDR_correction_mode, organism=species) ) 

# OUTPUT ====

## visualization

#image: png
png(filename="GSEA_Preranked_KOvWT_Activated.png", width=image_width, height=image_height, units="px", pointsize=4, bg="white", res=image_resolution, type="cairo") 

tab <- dplyr::group_modify(gsea, ~table.pvalue(.x)) %>% dplyr::ungroup()
ltab <- split(tab, tab$contrast) %>% lapply(function(x) plot.table(x, score=genescore) ) %>% wrap_plots()
print(ltab + plot_annotation(title="Cumulative number of significant calls (GSEA)", subtitle = sprintf("*p value adjusted %s by the method of Benjamini and Hochberg (1995)", FDR_correction_mode), tag_levels = 'A', theme=theme(plot.title = element_text(size=20, face='bold', hjust=0.5, margin = margin(t = 0)), plot.subtitle =  element_text(size=11, face='italic', hjust=0.5,  margin = margin(t = 10, b=20)))) )

          
## logs
cat("\nThe number of tested gene sets per each collection and contrast\n")
N <- count(gsea, collection); print(N, n=nrow(N))
cat(sprintf("\nCumulative number of significant calls\np-value adjusted for the false discovery rate %s by the method of Benjamini and Hochberg (1995)\n", FDR_correction_mode))
tab %>% print(n=nrow(tab))

## collapse redundant?

if (collapse_redundant_pathways == TRUE) {

    gsea_grouped <- gsea %>% dplyr::mutate(group_contrast=contrast) %>% dplyr::group_by(group_contrast)
    gsea <- dplyr::group_modify(gsea_grouped, ~collapse.gsea(., dx=genescore_df, collections=geneset_list)) %>% dplyr::ungroup() %>% dplyr::select(-group_contrast) %>% dplyr::group_by(contrast)    
    cat("\nThe number of non-redundant gene sets per each collection and contrast\n")
    N <- count(gsea, collection); print(N, n=nrow(N))
}

## return dataset   
 
if (sortDecreasing) { sortBy = sapply( sortBy, function(x) sprintf("desc(%s)", x) ) }

gsea <- gsea %>% dplyr::arrange_(.dots = sortBy) %>% tibble::add_column(geneScore = genescore, .after="contrast") %>% tibble::add_column(fdr_correction_mode = FDR_correction_mode, .after='geneScore')

return(gsea)    

}

########################################################
# Functions defined here will be available to call in
# the code for any table.
########################################################

## run gsea function

run.gsea <- function(dx, mode, collections, db, minimum_size, maximum_size, number_perms, organism) {

    # compute gsea stats

    ranked = dx$genescores        
    names(ranked) = dx$gene_id
    db$inPathway = sapply(db$gene_symbol, function(x) paste(sort(x[x %in% names(ranked)]), collapse=","))
    
    if (mode == "over all collections") {
            
        set.seed(246642)
        gsea <- fgsea(pathways=collections, stats=ranked, minSize=minimum_size, maxSize=maximum_size, nperm=number_perms) 
        gsea$size_leadingEdge <- sapply(gsea$leadingEdge, length)
        gsea$fraction_leadingEdge <- gsea$size_leadingEdge/gsea$size
        gsea$leadingEdge <- sapply(gsea$leadingEdge, function(x) paste(x, collapse=","))
        gsea <- dplyr::inner_join(gsea, select(db, collection, gene_set_name, inPathway), by=c("pathway"="gene_set_name")) %>% dplyr::select(collection, dplyr::everything())         
        
    } else {

        included_collections <- setNames(unique(db$collection), unique(db$collection))
        gsea <- lapply(included_collections, function(x) {
        set.seed(246642)
        gsea_collection = fgsea(pathways = collections[names(collections) %in% dplyr::filter(db, collection==x)$gene_set_name], stats=ranked, minSize=minimum_size, maxSize=maximum_size, nperm=number_perms)
        gsea_collection$size_leadingEdge <- sapply(gsea_collection$leadingEdge, length)
        gsea_collection$fraction_leadingEdge <- gsea_collection$size_leadingEdge/gsea_collection$size
        gsea_collection$leadingEdge <- sapply(gsea_collection$leadingEdge, function(x) paste(x, collapse=","))
        return(dplyr::inner_join(gsea_collection, select(db, collection, gene_set_name, inPathway) %>% filter(collection == x), by=c("pathway"="gene_set_name")) %>% dplyr::select(collection, dplyr::everything()) )                
            }) %>% dplyr::bind_rows()
        }
    gsea$species = organism
    return(gsea)        
    }

    ## collapse GSEA (from Matt Angel code)

    collapse.gsea <- function(grp, dx, collections){

    # filter ranked variable
    temp = dx %>% filter(contrast %in% grp$contrast)
    ranked = temp$genescores
    names(ranked) = temp$gene_id    
    
    # collapse function    
    run.collapse <- function(cp.input, pvalue, collections, ranked) {
        collapsedPathways <- collapsePathways(as.data.table(cp.input), collections, ranked, pval.threshold=pvalue) # requires the data.table library
        return(data.frame(pathway=collapsedPathways$mainPathways))            
    }

    filter_gsea = grp %>% dplyr::filter(pval < 0.05) %>% dplyr::arrange(pval) %>% dplyr::group_by(collection)
    collapsedResults = dplyr::group_modify(filter_gsea , ~run.collapse(., pvalue = 0.05, collections = collections, ranked=ranked)) %>% dplyr::ungroup()
    out <- grp %>% dplyr::inner_join(collapsedResults, by=c('pathway'='pathway','collection'='collection'))
        
    return(out)
        
}

## pvalue cutoffs table functions
    
table.pvalue <- function(gsea) {
        
    cuts <- c(-Inf, 1e-04, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 1)
    cutsLab <- paste("<",cuts[-1], sep="")
    p= cumsum(table(cut(gsea$pval, breaks = cuts, labels=cutsLab, include.lowest = FALSE, right=TRUE)))
    q= cumsum(table(cut(gsea$padj, breaks = cuts, labels=cutsLab,include.lowest = FALSE, right=TRUE)))
    tab <- data.frame(cutsLab, p, q)
    colnames(tab) <- c("alpha","p-value","*adjusted\np-value")
    rownames(tab) <- NULL
    return(tab)
}       
       
plot.table <- function(dtab, score) {

    title <- textGrob(paste0(unique(dtab$contrast), score),gp=gpar(fontsize=10))
    tab <- as.data.frame.matrix(dtab %>% dplyr::select(-contrast))
    table <- tableGrob(tab, theme = ttheme_default(core=list(fg_params=list(cex=0.9)), colhead = list(fg_params=list(cex = 0.9, parse=FALSE)), rowhead=list(fg_params=list(cex= 0.6))))
    table <- gtable_add_rows(table, heights = grobHeight(title) + unit(2,"line"), pos=0)
    table <- gtable_add_grob(table, list(title), t=c(1), l=c(1), r=ncol(table)) 
    wrap_elements( table)        
}      

print("template_function_GSEA_Preranked_KOvWT_Activated.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Analysis<-readRDS(paste0(rds_output,"/var_DEG_Analysis.rds"))
var_DEG_Analysis<-as.data.frame(var_DEG_Analysis)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_msigdb_v6_2_with_orthologs<-readRDS(paste0(rds_output,"/var_msigdb_v6_2_with_orthologs.rds"))
var_msigdb_v6_2_with_orthologs<-as.data.frame(var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
var_GSEA_Preranked_KOvWT_Activated<-GSEA_Preranked_KOvWT_Activated(var_DEG_Analysis,var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
saveRDS(var_GSEA_Preranked_KOvWT_Activated, paste0(rds_output,"/var_GSEA_Preranked_KOvWT_Activated.rds"))
