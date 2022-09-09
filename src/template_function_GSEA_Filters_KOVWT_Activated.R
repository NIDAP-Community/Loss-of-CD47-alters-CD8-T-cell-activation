GSEA_Filters_KOVWT_Activated <- function(GSEA_Preranked_KOvWT_Activated, msigdb_v6_2_with_orthologs ) {

graphics.off()

# LIBRARIES ====

    library(dplyr); library(ggplot2); library(plotly); library(htmlwidgets); library(RColorBrewer)
    

# INPUT PARAMETERS ====

## dataset
    gsea <- GSEA_Preranked_KOvWT_Activated
    pathway_db <- msigdb_v6_2_with_orthologs
    species <- "Mouse"

# GSEA table filters
    filterInput_by_Score = 'NES (Normalized Enrichment Score)'
    filterInput_byScore = switch(filterInput_by_Score, 'NES (Normalized Enrichment Score)'='NES', 'ES (Enrichment Score)'='ES')
    filterInput_cutoffScore = 0
    filterInput_signScore = "+/-"
    filterInput_by_Pvalue = 'adjusted p-value'
    filterInput_byPvalue = switch(filterInput_by_Pvalue, 'p-value'='pval', 'adjusted p-value'='padj')
    filterInput_cutoffPvalue = 0.05
    filterInput_bySize = "Pathway size"
    filterInput_cutoffSize = 0

## collection filters
    collections =  c()
    
## pathway filters 
    pathways = c()

## gene filters
    filterInput_bySet = "Leading Edge (LE)"
    genes = c()

## contrast filters
    filterInput_byContrast = "none"
    contrasts = c()
    
## legacy parameters    
    testedContrast = c()
    geneScore = c()

## top significant filter
    filterInput_byTop = Inf

    
# OUTPUT PARAMETERS ====

## plot visuals
    bubble_color = "collection" # collection
    bubble_opacity = 0.95
    bubble_max_size = 2
    xmin = c(); xmax = c()
    ymin = c(); ymax = c()

## returned dataset
    sortOutput_decreasing = FALSE
    sortOutput_by = c()
    output_type = "Foundry dataset"
        

# RUN ====

## adjust old Preranked GSEA output (release v67 or lower)

gsea <- adjust.v67(input=gsea, gS=geneScore, tC=testedContrast, sp=species, db=pathway_db)

## apply filters

    cat("\n\nFiltering steps\n\n")

    gsea_filtered <- gsea %>% dplyr::filter(get(filterInput_byPvalue) <= filterInput_cutoffPvalue)
    filter.message(filter=filterInput_by_Pvalue, condition=filterInput_cutoffPvalue, output=gsea_filtered)

    if(filterInput_signScore == "+"){
            gsea_filtered <- gsea_filtered %>% dplyr::filter(get(filterInput_byScore) >= filterInput_cutoffScore)
            filter.message(filter="value of GSEA score", condition=sprintf("%s > %g", filterInput_byScore, filterInput_cutoffScore), output=gsea_filtered)
        
    } else if (filterInput_signScore == "-") {
            gsea_filtered <- gsea_filtered %>% dplyr::filter(get(filterInput_byScore) <= -filterInput_cutoffScore)
            filter.message(filter="value of GSEA score", condition=sprintf("%s < %s%g", filterInput_byScore, ifelse(filterInput_cutoffScore==0,"","-"), filterInput_cutoffScore), output=gsea_filtered)
        
    } else {
            gsea_filtered <- gsea_filtered %>% dplyr::filter(abs(get(filterInput_byScore)) >= filterInput_cutoffScore)
            filter.message(filter="value of GSEA score", condition=sprintf("|%s| > %g", filterInput_byScore, filterInput_cutoffScore), output=gsea_filtered)
    }

    if (filterInput_cutoffSize > 0 ) {

        if ( filterInput_bySize == "Pathway size") {
            gsea_filtered <- gsea_filtered %>% dplyr::filter(size >= filterInput_cutoffSize)
            filter.message(filter=filterInput_bySize, condition=filterInput_cutoffSize, output=gsea_filtered)
        
        } else {
            gsea_filtered <- gsea_filtered %>% dplyr::filter(size_leadingEdge >= filterInput_cutoffSize)
            filter.message(filter=filterInput_bySize, condition=filterInput_cutoffSize, output=gsea_filtered)
        }
    
    } else {
        filter.message(filter=filterInput_bySize, condition=paste0(">",filterInput_cutoffSize), output=gsea_filtered)
    }

    if ( ! is.null(collections) ) {

        gsea_filtered <- gsea_filtered %>% dplyr::filter( collection %in% collections ) 
        filter.message(filter="Collection filter ON", condition="Collections to include", output=gsea_filtered)
        cat(sprintf("    Found collections: %s\n", paste(unique(gsea_filtered$collection), collapse=", ")))
    }

    if ( ! is.null(pathways) ) {

        gsea_filtered <- gsea_filtered %>% dplyr::filter( pathway %in% pathways ) 
        filter.message(filter="Pathway filter ON", condition="Pathways to include", output=gsea_filtered)
        cat(sprintf("    Found pathways: %s\n", paste(unique(gsea_filtered$pathway), collapse=", ")))
    }

    if ( ! is.null(genes) ) {

        if (filterInput_bySet == 'Leading Edge (LE)') { 
            
            index_pathway = lapply(genes, function(x) grep(paste0("\\Q", x, "\\E"), gsea_filtered$leadingEdge))
            found = sapply(index_pathway, function(x) {length(x) > 0})
            index_pathway = Reduce(unique, index_pathway)
            gsea_filtered <- gsea_filtered %>% dplyr::slice(index_pathway) 
                
            filter.message(filter="Gene filter ON", condition="Leading Edge", output=gsea_filtered)
            cat(sprintf("    Found genes: %s\n", paste(genes[which(found)], collapse=", ")))
            if (sum(!found) > 0) { cat(sprintf("    Missing genes: %s\n", paste(genes[which(!found)], collapse=", "))) }
                
        } else if (filterInput_bySet == "Pathway") { 
                
            index_pathway = lapply(genes, function(x) grep(paste0("\\Q", x, "\\E"), gsea_filtered$inPathway))
            found = sapply(index_pathway, function(x) {length(x) > 0})
            index_pathway = Reduce(unique, index_pathway)
            gsea_filtered <- gsea_filtered %>% dplyr::slice(index_pathway) 
                
            filter.message(filter="Gene filter ON", condition="in Pathway", output=gsea_filtered)
            cat(sprintf("    Found genes: %s\n", paste(genes[which(found)], collapse=", ")))
            if (sum(!found) > 0) { cat(sprintf("    Missing genes: %s\n", paste(genes[which(!found)], collapse=", "))) }
        
        }
    }

    if ( filterInput_byContrast == "remove" ) { 

        if ( ! is.null(contrasts) ) {
            
            all_contrasts = unique(gsea_filtered$contrast)
            gsea_filtered <- gsea_filtered %>% dplyr::filter( ! contrast %in% contrasts ) 
            removed = setdiff(all_contrasts, unique(gsea_filtered$contrast))
            if (length(removed) < 1) {
                filter.message(warn=TRUE, filter="Contrast filter", condition=sprintf("remove %s missing", paste(contrasts, collapse=", ")), output=gsea_filtered)
            } else {
                filter.message(filter="Contrast filter", condition=filterInput_byContrast, output=gsea_filtered)
                cat(sprintf("    Removed contrast(s): %s\n", paste(removed, collapse=", ")))
                cat(sprintf("    Keep contrast(s): %s\n", paste(unique(gsea_filtered$contrast), collapse=", ")))
            }

        } else if (is.null(contrasts) ) {
            filter.message(warn=TRUE, filter="Contrast filter", condition=class(contrasts), output=gsea_filtered)
        }

    } else if ( filterInput_byContrast == "keep" ) {

        if ( ! is.null(contrasts) ) {

            all_contrasts = unique(gsea_filtered$contrast)
            gsea_filtered <- gsea_filtered %>% dplyr::filter( contrast %in% contrasts ) 
            kept = intersect(all_contrasts, unique(gsea_filtered$contrast))
            removed = setdiff(all_contrasts, unique(gsea_filtered$contrast))
            if (length(kept) < 1) {
                filter.message(warn=TRUE, filter="Contrast filter", condition=sprintf("keep %s missing", paste(contrasts, collapse=", ")), output=gsea_filtered)
            } else {
                filter.message(filter="Contrast filter", condition=filterInput_byContrast, output=gsea_filtered)
                cat(sprintf("    Removed contrast(s): %s\n", paste(removed, collapse=", ")))
                cat(sprintf("    Kept contrast(s): %s\n", paste(unique(gsea_filtered$contrast), collapse=", ")))
            }

        } else if (is.null(contrasts) ) {
            filter.message(warn=TRUE, filter="Contrast filter", condition=class(contrasts), output=gsea_filtered)
        }

    } else if ( filterInput_byContrast == "none" ) {

        if ( ! is.null(contrasts) ) {
            filter.message(warn=TRUE, filter="Contrast filter", condition=sprintf("none; %s", paste(contrasts, collapse=", ")), output=gsea_filtered)
        }
    }

    if (filterInput_byTop == "Inf") { filterInput_byTop = Inf } else { filterInput_byTop = as.numeric(filterInput_byTop) }
    
    gsea_filtered <- gsea_filtered %>% dplyr::group_by(contrast, collection) %>% dplyr::mutate(p_rank = rank(pval, ties.method="min")) %>% filter(p_rank <= filterInput_byTop) %>% dplyr::select(-p_rank)    
    filter.message(filter="Top significant filter", condition=sprintf("up to p-value rank of %s per contrast and collection", filterInput_byTop), output=gsea_filtered)

# OUTPUT ====

## sort output
    if (sortOutput_decreasing) { sortOutput_by = sapply( sortOutput_by, function(x) sprintf("desc(%s)", x) ) }
    gsea_filtered <- gsea_filtered %>% dplyr::arrange_(.dots = sortOutput_by)
    
## do plot and return dataset

    if ( nrow(gsea_filtered) == 0 ) { 
        stop("ERROR: filtering returned 0 pathways")

    } else {

        cat("\n\nFiltered pathways\n")
        tab = table(gsea_filtered$collection, gsea_filtered$contrast) %>% addmargins(margin=c(1,2))
        print(tab)
        
        cat("\n\nGSEA statistics (filtered pathways)\n\n")
        print(tibble(gsea_filtered))

        df <- gsea_filtered %>% dplyr::select(-leadingEdge,-inPathway) %>% dplyr::mutate(textContrast=sprintf("%s: NES = %g, P-value = %g", contrast, signif(NES,2), signif(pval,1))) %>% dplyr::group_by(pathway, collection) %>% 
            summarize(mean_pval=mean(pval), n_contrast = length(contrast), Pathway_size=mean(size), mean_NES=mean(NES), individual_values = paste(textContrast, collapse="\n")) %>% 
                dplyr::mutate(textPathway=sprintf("%s<br>%s<br>mean NES = %g, mean P-value = %g<br>Pathway size = %g<br><br>N contrasts = %g\n%s", collection, pathway, signif(mean_NES,2), signif(mean_pval,1), Pathway_size, n_contrast, individual_values)) %>%
                    dplyr::select(collection, pathway, Pathway_size, n_contrast, mean_NES, mean_pval, individual_values, textPathway)
        
        find_sort = grepl( paste(sortOutput_by, collapse="|"), colnames(df) )        
        if ( sum(find_sort) == 0 ) { sort_by = c("n_contrast", "collection", "mean_pval") } else { sort_by = colnames(df)[which(find_sort)] }
          
        if (sortOutput_decreasing) { sort_by = sapply( sort_by, function(x) sprintf("desc(%s)", x) ) }
        df <- df %>% dplyr::arrange_(.dots = sort_by)

        #cat("\n\nGSEA contrast summary (filtered pathways)\n\n")
        #print(tibble(df %>% dplyr::select(-textPathway) %>% dplyr::rename("mean_size"="Pathway_size")))
                
        if ( bubble_color == "pathway size" ) {

                ggp <- plot_ly( data= df, x = ~mean_NES, y= ~-log10(mean_pval), color= ~Pathway_size, size= ~n_contrast, text= ~textPathway, hoverinfo= "text", opacity=bubble_opacity, marker = list(sizeref = ~2*max(n_contrast)/bubble_max_size**2)) %>% 
                    
                    layout(
                        xaxis = list(title="Normalized Enrichment Score (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE),
                            yaxis = list(title="-log10 P-value (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE), 
                                legend= list(itemsizing='constant'), title=list(text="Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway", x=0, font=list(size=15)), showlegend=TRUE
                    )
                    
        } else {

                ggp <- plot_ly( data= df, x = ~mean_NES, y= ~-log10(mean_pval), color= ~collection, size= ~n_contrast, text= ~textPathway, hoverinfo= "text", opacity=bubble_opacity, marker = list(sizeref = ~2*max(n_contrast)/bubble_max_size**2)) %>% 
                    
                    layout(
                        xaxis = list(title="Normalized Enrichment Score (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE),
                            yaxis = list(title="-log10 P-value (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE), 
                                legend= list(itemsizing='constant'), title=list(text="Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway", x=0, font=list(size=15)), showlegend=TRUE
                    )
        } 
# custom axis range?

if ( !(is.null(xmin) & is.null(xmax)) ) {
    
    if (is.null(xmin)) xmin = floor(min(df$mean_NES))
    if (is.null(xmax)) xmax = ceiling(max(df$mean_NES))
    ggp <- ggp %>% layout(xaxis = list(range = list(xmin, xmax)))
}

if ( !(is.null(ymin) & is.null(ymax)) ) {
    
    if (is.null(ymin)) ymin = floor(min(-log10(df$mean_pval)))
    if (is.null(ymax)) ymax = ceiling(max(-log10(df$mean_pval)))
    ggp <- ggp %>% layout(yaxis = list(range = list(ymin, ymax)))
}

# raw file output silenced as of now   
        if (output_type == 'raw files') {

            cat(sprintf("\n\nRaw files\n\n"))

            html_by = "collection"

# auto removed:             output <- new.output()
# auto removed:             output_fs <- output$fileSystem()         
            
            fileName <- "gsea_filtered_statistics.csv"; cat(sprintf("%s\n",fileName))
            write.csv(gsea_filtered, file(fileName, 'w'))

            fileName <- "gsea_filtered_contrastSummary.csv"; cat(sprintf("%s\n",fileName))
            write.csv(df %>% dplyr::select(-textPathway) %>% dplyr::rename("mean_size"="Pathway_size") %>% dplyr::mutate(individual_values=gsub("\n", "; ", individual_values)), file(fileName, 'w'))
                    
            fileName <- "gsea_filtered.html"; cat(sprintf("%s\n\n",fileName))
            htmlwidgets::saveWidget(ggp, fileName)
            output_fs$upload(fileName, fileName) 

            df_grouped <- df %>% dplyr::group_by_at(html_by)            
            keys = dplyr::group_keys(df_grouped) %>%  tidyr::unite(html_by)    
            df_list <- df_grouped %>% dplyr::group_split() %>% setNames(keys$html_by)
            
            for ( i in keys$html_by)  {
            
                if ( bubble_color == "pathway size" ) {

                    fig <- plot_ly( data= df_list[[i]], x = ~mean_NES, y= ~-log10(mean_pval), color= ~Pathway_size, size= ~n_contrast, text= ~textPathway, hoverinfo= "text", opacity=bubble_opacity, marker = list(sizeref = ~2*max(n_contrast)/bubble_max_size**2)) %>% 
                        
                        layout(
                            xaxis = list(title="Normalized Enrichment Score (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE),
                                yaxis = list(title="-log10 P-value (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE), 
                                    legend= list(itemsizing='constant'), title=list(text="Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway", x=0, font=list(size=15)), showlegend=TRUE
                        )
                        
                } else {

                    fig <- plot_ly( data= df_list[[i]], x = ~mean_NES, y= ~-log10(mean_pval), color= ~collection, size= ~n_contrast, text= ~textPathway, hoverinfo= "text", opacity=bubble_opacity, marker = list(sizeref = ~2*max(n_contrast)/bubble_max_size**2)) %>% 
                        
                        layout(
                            xaxis = list(title="Normalized Enrichment Score (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE),
                                yaxis = list(title="-log10 P-value (mean)", tickfont=list(size=15), titlefont=list(size=15), showgrid=TRUE), 
                                    legend= list(itemsizing='constant'), title=list(text="Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway", x=0, font=list(size=15)), showlegend=TRUE
                        )
                } 

                fileName = make.names(sprintf("%s.html", i)); cat(sprintf("%s\n",fileName))
                htmlwidgets::saveWidget(fig, fileName)
                output_fs$upload(fileName, fileName)

            }
        
        } else if (output_type == "Foundry dataset" ) {

                print(ggp)                
                return( gsea_filtered )
        }
        
        print(ggp)
    }

}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

## HELPER FUNCTION

    filter.message <- function(filter, warn=FALSE, condition, output) {    
        
        n = length(unique(output$pathway))
        if ( ! warn ){
            if (n == 0) { 
                cat(sprintf("ERROR: Filter by %s (%s) returned %g unique pathway(s)\n", filter, condition, n))
                stop("Filter condition error\n")  
            } else { 
                cat(sprintf("OK: Filter by %s (%s) returned %g unique pathway(s)\n", filter, condition, n))
            }
        } else {
            cat(sprintf("WARNING: Filter by %s (%s) not specified correctly; this filter is not applied\n", filter, condition))
        }
    }

# adjust old GSEA table (v67 or lower)
    adjust.v67 <- function(input, gS, tC, sp, db, required_columns=c("contrast","geneScore","fdr_correction_mode","collection","pathway","pval","padj","ES","NES","nMoreExtreme","size","leadingEdge","size_leadingEdge","inPathway")) {

    missing_columns = setdiff(required_columns, colnames(input))

    if ( length(missing_columns) > 0 ) {

        cat("WARNING: Output from outdated 'Preranked GSEA [CCBR]' detected.\n")
        cat(sprintf("\nWARNING: Missing columns are added (%s)\n\t 'inPathway' column includes all genes annotated to a gene set - run the latest released version of 'Preranked GSEA [CCBR]'\n\t to return only genes mapped in the dataset.\n", paste(missing_columns, collapse=", ")))

        if ( is.null(gS) ) {
            gS = "NULL"
            cat("\nWARNING: 'Gene score' parameter not provided therefore 'geneScore' column is assigned the default NULL value\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it will fail due to this specification of the 'Gene score' parameter.\n")        
        } else {
            if ( length(gS) > 1 ) {
                gS = gS[1]
                cat(sprintf("\nWARNING: too many values for 'Gene score' provided; only the first one used, '%s'\n", gS))
            }
            has.underscore = substring(gS, 1, 1) == "_" 
            if ( !has.underscore ) {
                gS = paste0("_",gS[1])
                cat(sprintf("\nWARNING: 'Gene score' should start with underscore to match column naming convention in a DEG table used for GSEA run (e.g. '_tstat');\n\t column 'geneScore' is assigned the provided value with '_' added ('%s');\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it may fail if this is not an adequate specification of the 'Gene score' parameter.\n", gS))
            }
        }

        if ( is.null(tC) ) {
            tC = "NULL"
            cat("\nWARNING: 'Tested contrast' parameter not provided therefore 'contrast' column is assigned the default NULL value;\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it will fail due to this specification of the 'Tested contrast' parameter.\n")        
        } else {
            if ( length(tC) > 1 ) {
                tC = tC[1]
                cat(sprintf("\nWARNING: too many values for 'Tested contrast' provided; only the first one used, '%s'\n",tC))
            }
            is.contrast = any(grepl("-", tC))
            if ( !is.contrast ) {
                cat(sprintf("\nWARNING: 'Tested contrast' should be specified exactly the same way as when Preranked GSEA was run (e.g. treated-control);\n\t column 'contrast' is assigned the provided '%s' value;\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it may fail if this is not an adequate specification of the 'Tested contrast' parameter.\n", tC[1]))
            }   
        } 

        input <- input %>% mutate(contrast=tC, geneScore=gS, fdr_correction_mode="over all collections", size_leadingEdge=sapply(strsplit(leadingEdge, ","), length))
        input <- input[, match(required_columns[required_columns!='inPathway'], colnames(input))]    
        db <- db %>% filter(db[["species"]]==sp) %>% filter(`%in%`(db[["collection"]] ,unique(input$collection))) %>% filter(`%in%`(db[["gene_set_name"]], unique(input$pathway))) %>% collect()
        db <- db %>% dplyr::group_by(collection, gene_set_name, species) %>% dplyr::summarize(inPathway = paste(unique(gene_symbol), collapse = ",")) %>% dplyr::ungroup() %>% data.frame()
        input <- input %>% left_join(db, by=c("collection"="collection","pathway"="gene_set_name"))
        input <- input[, na.omit(match(c(required_columns,"species"),colnames(input)))]
        return(input)

    } else {

        cat("Filtering Preranked GSEA table.\n")
        return(input)
    }
}

print("template_function_GSEA_Filters_KOVWT_Activated.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_GSEA_Preranked_KOvWT_Activated<-readRDS(paste0(rds_output,"/var_GSEA_Preranked_KOvWT_Activated.rds"))
var_GSEA_Preranked_KOvWT_Activated<-as.data.frame(var_GSEA_Preranked_KOvWT_Activated)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_msigdb_v6_2_with_orthologs<-readRDS(paste0(rds_output,"/var_msigdb_v6_2_with_orthologs.rds"))
var_msigdb_v6_2_with_orthologs<-as.data.frame(var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
var_GSEA_Filters_KOVWT_Activated<-GSEA_Filters_KOVWT_Activated(var_GSEA_Preranked_KOvWT_Activated,var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
saveRDS(var_GSEA_Filters_KOVWT_Activated, paste0(rds_output,"/var_GSEA_Filters_KOVWT_Activated.rds"))
