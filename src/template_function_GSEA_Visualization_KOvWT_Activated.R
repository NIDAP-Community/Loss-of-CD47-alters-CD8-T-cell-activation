GSEA_Visualization_KOvWT_Activated <- function(DEG_Analysis, GSEA_Filters_KOVWT_Activated, Normalized_Counts, ccbr804_metadata_for_NIDAP, msigdb_v6_2_with_orthologs) {

    graphics.off()

## LIBRARIES ====
    library(tibble); library(dplyr); library(tidyr)
    library(patchwork); library(ComplexHeatmap); library(colorspace); library(RColorBrewer)
    library(fgsea)

## PARAMETERS

    which_plot = "ES+RNK+LE"
    plot_limit = "10"
    topline = 'coordinate'
    pathways = c()
    geneid = "Gene"
    geneid_gex = "Gene"
    contrasts = c()
    species = "Mouse"

    genescore_df <- DEG_Analysis
    gex_df <- Normalized_Counts
    metadata_df <- ccbr804_metadata_for_NIDAP
    gsea_df <- GSEA_Filters_KOVWT_Activated
    pathway_db <- msigdb_v6_2_with_orthologs
        
    gex_transformation = 'median centering'
    reference_phenotype = FALSE
    drop_ref = FALSE
    
    pdfType = "common PDF"
    pageWidth=11
    pageHeight=8

    geneScore = c()
    testedContrast = c()  
    

#..GSEA

## adjust old Preranked GSEA output (release v67 or lower)
gsea_df <- adjust.v67(input=gsea_df, gS=geneScore, tC=testedContrast, sp=species, db=pathway_db)   

# gene score name
genescore = unique(gsea_df$geneScore)

# pathway set
if ( !is.null(pathways) ) { gsea_df <- gsea_df %>% dplyr::filter(pathway %in% pathways) }

# gsea stats
gsea_df <- gsea_df %>% dplyr::select( colnames(gsea_df)[colnames(gsea_df) %in% c("contrast", "collection", "pathway", "ES", "NES", "pval", "padj", "leadingEdge","size_leadingEdge", "inPathway")] )
if ( !is.null(contrasts) ) { 
    gsea_df <- gsea_df %>% dplyr::filter(contrast %in% contrasts)
} 

#..GSEA ES  

if (grepl("ES", which_plot)) {
    
    # gene scores    
    rank_columns = colnames(genescore_df)[grepl(paste0("\\Q", genescore, "\\E$"), colnames(genescore_df))]
    rank_contrasts = unlist(strsplit(rank_columns, genescore))
    if ( length(rank_columns)==0 ) stop("ERROR: 'Gene score' not specified correctly")
    if ( length(rank_contrasts)==0 ) stop("ERROR: 'Tested contrasts' not specified correctly")
    
    
    if (!is.null(contrasts)) {  

        index = match(contrasts, rank_contrasts)
        rank_columns = rank_columns[index]
        rank_contrasts = rank_contrasts[index]
        groups_from_contrasts = unique(unlist(strsplit(rank_contrasts,"-"))) 
    
    } else if (is.null(contrasts)) {         
        rank_contrasts = unique(gsea_df$contrast)
    }
    
    genescore_df <- genescore_df %>% dplyr::select(geneid, rank_columns) %>% tidyr::pivot_longer(!geneid, names_to="contrast", values_to="genescores", values_drop_na=TRUE) %>% dplyr::rename("geneid"=geneid) %>% dplyr::mutate(contrast=sub(genescore, "", contrast)) %>% tidyr::drop_na() %>% dplyr::arrange(desc(genescores))

}

#..LE HEATMAP

if (grepl("LE", which_plot)) {
    
    # samples in gene expression dataset
    samples_to_include = setdiff(c("WT_UT_1","WT_CD3_CD28_1","CD47KO_UT_1","CD47KO_CD3_CD28_1","WT_UT_2","WT_CD3_CD28_2","CD47KO_UT_2","CD47KO_CD3_CD28_2","WT_UT_3","WT_CD3_CD28_3","CD47KO_UT_3","CD47KO_CD3_CD28_3","WT_UT_4","WT_CD3_CD28_4","CD47KO_UT_4","CD47KO_CD3_CD28_4"), geneid)
        
    # gene expression
    le_genes <- gsea_df %>% dplyr::select(leadingEdge) %>% tidyr::separate_rows(leadingEdge, sep=",") %>% dplyr::distinct()
    gex_df = gex_df %>% dplyr::select(geneid_gex, samples_to_include) %>% dplyr::rename("Gene"=geneid_gex) %>% inner_join(le_genes, by=c("Gene"="leadingEdge"))
}    

# DO PLOT ====

cat(sprintf("Saving files in the workbook-output:\n\nGSEA-RunningES.csv\n\n"))

# retain only top significant pathways if top rank requested

plot_limit = tolower(plot_limit)

if ( plot_limit != "all" ) {    
    
    n_input = nrow(gsea_df)
    
    plot_limit = as.numeric(plot_limit)
    
    if (is.na(plot_limit)) { 
        stop("ERROR in Top rank filter; enter ALL (case insensitive) or a numeric rank.\n")
    } else if (plot_limit <= 0 ) { 
        plot_limit = 1 
        cat("WARNING: 'Top rank filter' cannot be 0 or less; its value was changed to 1.\n")
    }
        
    gsea_df <- gsea_df %>% dplyr::group_by(contrast, collection) %>% dplyr::mutate(p_rank = rank(pval, ties.method="min")) %>% filter(p_rank <= plot_limit) %>% dplyr::select(-p_rank) %>% dplyr::ungroup()

    if(n_input > nrow(gsea_df)) {
        cat(sprintf("WARNING: Preparing top-ranked plots from each contrast and collection (%g out of %g) based on max P-value rank of %g\n\t Change 'Top rank filter'parameter if you intend to generate more plots\n", nrow(gsea_df), n_input, plot_limit))
    } else {
        cat(sprintf("Preparing top-ranked plots from each contrast and collection (%g out of %g) based on max P-value rank of %g\n", nrow(gsea_df), n_input, plot_limit))
    }

} else {
    
    cat(sprintf("Preparing all available plots (%g)\n", nrow(gsea_df)))
}

gsea_list<- dplyr::group_split(gsea_df, contrast, pathway) 
names(gsea_list) = sapply(gsea_list, function(x) paste(x$pathway, x$contrast, sep="_"))

# auto removed: output <- new.output()
# auto removed: output_fs <- output$fileSystem()
conn <- file("GSEA-RunningES.csv", 'w')
header=c('contrast', 'pathway', 'geneRank', 'runningES', 'gene', 'leadingEdge', 'geneScore')
write.table(t(header), conn, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')

for ( i in 1:length(gsea_list) ) {

    gsea = gsea_list[[i]]
     txt <- sprintf( "%s, %s, ES=%g, NES=%g, pval=%g, padj=%g", gsea$collection, gsea$contrast,  round(gsea$ES,2), round(gsea$NES,2), signif(gsea$pval,2), signif(gsea$padj,2) )
    name = gsea$pathway

    fontsize_row = 6
    fontsize_col = 8

    counter1 = seq(0,length(gsea_list),100)
    counter2 = seq(0,length(gsea_list),25)
    if (i %in% counter1 | i == length(gsea_list)) { cat(i,"\n") } else if (i %in% counter2) { cat(i, " ") } else { cat(".")}

    plotES = gg.plotES(ranks=genescore_df, gsea=gsea, ntop=10, add_top=TRUE, cex_top=3, image_size='reduced', gset=gset, line_top=topline, ES_colo = 'ES sign')
    gsea_out <- plotES[[3]]
    write.table(gsea_out, conn, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',', append=TRUE)

       
    if (which_plot == "ES+RNK+LE") {

        # ES with heatbar and RNK barplot
        pRunes <- plotES[[1]]
        pRank <- plotES[[2]]        
        
        # gex transformation
        le = unlist(strsplit(gsea$leadingEdge, ","))
        gex_le = gex_df %>% dplyr::filter(Gene %in% le) %>% tibble::column_to_rownames("Gene") 
        metadata = prep.metadata(meta=metadata_df, gsea=gsea, samples=colnames(gex_le), Sample="Sample", Group="Genotype_Plating", transformation=gex_transformation, dropREF=drop_ref, reference=reference_phenotype, own_palette=c(), label_palette='Accent')
        gex_le = gex_le[, match(metadata$Sample, colnames(gex_le))]
        gex_trans = transform.gex(gex_le, meta=metadata, transformation = gex_transformation )
        if (drop_ref) {
            if(grepl("reference", gex_transformation)) {
                gex_trans <- gex_trans[, metadata$Condition != 'Reference']
                metadata <- metadata[metadata$Condition != 'Reference', ]
            }
        }
        
        # clustering        
        clustOrder <- prep.clustOrder( df=gex_trans, linkage='complete', distance='Euclidean', way='rows and columns' )
        rowv = clustOrder$rowv
        colv = clustOrder$colv

        # font size
        if (fontsize_row == 0) {
        fontsize_row = find.fontsize(max_size=8, min_size=3, max_n=c(500, gsea$size_leadingEdge)[which.max(c(500, gsea$size_leadingEdge))], min_n=25, stepdown=0.6, n=gsea$size_leadingEdge) }
        if (fontsize_col == 0) {
        fontsize_col = find.fontsize(max_size=4, min_size=3, max_n=c(500, nrow(metadata))[which.max(c(500, nrow(metadata)))], min_n=25, stepdown=0.6, n=nrow(metadata)) }
              
        # heatmap
        pLEdge = plot.heatmap( gex=gex_trans, meta=metadata, heat_colors=NULL, limit=NULL, rowv=rowv, colv=colv, transformation=gex_transformation, show_rownames = TRUE, show_colnames = FALSE, show_coldend = FALSE, show_rowdend = TRUE, row_size=fontsize_row, col_size=fontsize_col, heatmap_leg = TRUE, sample_leg = TRUE)

        # plot layout     
        lay <- c(
        area(t = 1, l = 1, b = 2, r = 2), 
        area(t = 3, l = 1, b = 4, r = 2),
        area(t = 1, l = 3, b = 4, r = 3))

        # plot
        patch = pRunes + pRank + patchwork::wrap_elements(pLEdge) + plot_annotation(title = name, subtitle = txt) + plot_layout( design=lay) + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(size = 14))

    } else if (which_plot == "ES+RNK") {

        # ES with heatbar and RNK barplot
        pRunes <- plotES[[1]]
        pRank <- plotES[[2]]
                
        # plot layout
        lay <- c(
        area(t = 1, l = 1, b = 2, r = 2), 
        area(t = 3, l = 1, b = 4, r = 2))
        
        # plot
        patch = pRunes + pRank + plot_annotation(title = name, subtitle = txt) + plot_layout( design=lay) + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(size = 14))
        
    } else if (which_plot == "ES") {

        # ES with heatbar and RNK barplot
        pRunes <- plotES[[1]]
                
        # plot
        patch = pRunes + plot_annotation(title = name, subtitle = txt)
       
    } else if (which_plot == "LE") {

        # gex transformation
        le = unlist(strsplit(gsea$leadingEdge, ","))
        gex_le = gex_df %>% dplyr::filter(Gene %in% le) %>% tibble::column_to_rownames("Gene") 
        metadata = prep.metadata(meta=metadata_df, gsea=gsea, samples=colnames(gex_le), Sample="Sample", Group="Genotype_Plating", transformation=gex_transformation, dropREF=drop_ref, reference=reference_phenotype, own_palette=c(), label_palette='Accent')
        gex_le = gex_le[, match(metadata$Sample, colnames(gex_le))]
        gex_trans = transform.gex(gex_le, meta=metadata, transformation = gex_transformation )
        if (drop_ref) {
            if(grepl("reference", gex_transformation)) {
                gex_trans <- gex_trans[, metadata$Condition != 'Reference']
                metadata <- metadata[metadata$Condition != 'Reference', ]
            }
        }
        
        # clustering        
        clustOrder <- prep.clustOrder( df=gex_trans, linkage='complete', distance='Euclidean', way='rows and columns' )
        rowv = clustOrder$rowv
        colv = clustOrder$colv
      
        # font size
        if (fontsize_row == 0) {
        fontsize_row = find.fontsize(max_size=8, min_size=3, max_n=c(500, gsea$size_leadingEdge)[which.max(c(500, gsea$size_leadingEdge))], min_n=25, stepdown=0.6, n=gsea$size_leadingEdge)
        }
        if (fontsize_col == 0) {
        fontsize_col = find.fontsize(max_size=4, min_size=3, max_n=c(500, nrow(metadata))[which.max(c(500, nrow(metadata)))], min_n=25, stepdown=0.6, n=nrow(metadata))
        }

        # heatmap
        pLEdge = plot.heatmap( gex=gex_trans, meta=metadata, heat_colors=NULL, limit=NULL, rowv=rowv, colv=colv, transformation=gex_transformation, show_rownames = TRUE, show_colnames = FALSE, show_coldend = FALSE, show_rowdend = TRUE, row_size=fontsize_row, col_size=fontsize_col, heatmap_leg = TRUE, sample_leg = TRUE)

        # plot
        patch = wrap_elements(pLEdge) + plot_annotation(title = name, subtitle = txt)
        
    }

    if (pdfType == 'common PDF') {
        
        gsea_list[[i]] = patch
        if (i == 1) { preview = patch }
    
    } else {
        
        fileName =  make.names(sprintf("%s.pdf", names(gsea_list)[i]))
        pdf(file(fileName, 'w'), height=pageHeight, width=pageWidth)
        print(patch)
        dev.off()
        if (i == 1) { preview = patch }
        
    }
}

if (pdfType == 'common PDF') {

    fileName=sprintf("GSEA-Plot_%s.pdf", make.names(which_plot))
    cat(sprintf("\n%s\n\n", fileName))          
    pdf(file(fileName, 'w'), height=pageHeight, width=pageWidth)
    lapply(gsea_list, function(patch) {
        print(patch)
        text_contrast = sapply(strsplit(patch$patches$annotation$subtitle, ", "), function(x) x[2])
        text_pathway = patch$patches$annotation$title
        cat(sprintf("%s (%s)\n", text_pathway, text_contrast))
    })
    dev.off()

} else {
    cat(sprintf("\n%s.pdf", make.names(names(gsea_list))))
}

png(filename="GSEA_Visualization_KOvWT_Activated.png", width=pageWidth*300, height=pageHeight*300, units="px", pointsize=4, bg="white", res=300, type="cairo") 
print(preview)

}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

gg.plotES <- function(ranks, gsea, ntop, add_top, image_size, gset, ES_colo, line_top, cex_top) {

    # get all ranks
    ranks <- ranks %>% dplyr::filter(contrast %in% gsea$contrast)
    rnk=ranks$genescores
    names(rnk) = ranks$geneid
    le = gsea %>% dplyr::select(leadingEdge) %>% tidyr::separate_rows(leadingEdge) 
        lep=is.element(names(rnk), le$leadingEdge) & rnk>0
    len=is.element(names(rnk), le$leadingEdge) & rnk<0
    zero = sum(rnk>0)
    yrange = diff(range(rnk))
    ymax = max(rnk)
    ymin = min(rnk)

    # get gene set ranks and running scores
    ES_sign = sign(gsea$NES)
    gset = gsea %>% dplyr::select(inPathway) %>% tidyr::separate_rows(inPathway)
    gs = intersect(paste(gset$inPathway),names(rnk))
    es.data <- plotEnrichment(pathway=gs, stats=rnk, gseaParam = 1)
    x <- es.data$data$x
    y <- es.data$data$y
    xranks = sort(unname(as.vector(na.omit(match(gs, names(rnk))))))
    
    gs_ranks = data.frame(xranks=xranks)
    gs_scores = data.frame(x = x, y = y)
    ys = round( y[ -c(1, length(y)) ], 10)
    if(ES_sign > 0) {
        keep=seq(2,length(ys), 2)
    } else {
        keep=seq(1,length(ys), 2)
    }    
    es_data=data.frame(rank=xranks,running_es=ys[keep])

    # set positive/negative params    
    
    if( ES_sign > 0 ){
        
        all_ranks = data.frame(Index=1:length(rnk), Rank=rnk, LE=ifelse(lep==TRUE,'LE','Outside'), order=ifelse(lep==TRUE,2,1))
        
        if(add_top==TRUE & ntop > 0){
            if(sum(lep) < ntop) {
                ntop = sum(lep)                
                warning(sprintf("Max number of leading edge genes available is %g", sum(lep)))
            }
            top = sort(rnk[lep], decreasing=TRUE)[1:ntop]; top = paste(names(top), collapse='\n')
            topx = 1; topy= ymin
            nx=which(all_ranks$LE=="LE")[1]; ny = 0-yrange/30
            v=1; h=0; hn=0
            colorMargin = color.margins()['up']
            angletop=90
        }
    
    } else if ( ES_sign < 0) {
        
        all_ranks = data.frame(Index=1:length(rnk), Rank=rnk, LE=ifelse(len==TRUE,'LE','Outside'), order=ifelse(len==TRUE,2,1))
        
        if(add_top==TRUE & ntop > 0){
            if(sum(len) < ntop) {
                ntop = sum(len)
                warning(sprintf("Max number of leading edge genes available is %g", sum(len)))
            }
            top = sort(rnk[len],decreasing=FALSE)[1:ntop]; top = paste((names(top)), collapse='\n') # rev(names(top)) if angle=90
            topx = length(rnk); topy= ymax
            nx=which(all_ranks$LE=="LE")[sum(len)]; ny=0+yrange/30
            v=1; h=0; hn=1
            colorMargin = color.margins()['dn']
            angletop=-90
        } 
    }
   
    # set miscelanous

    #.. zero arrow
    df_arrow <- data.frame(x1 = zero+0.5, x2 = zero+0.5, y1 = 0, y2 = 0+yrange/16)

    #.. keep only gene set ranks
    if(image_size == 'reduced') { all_ranks$Rank = ifelse(all_ranks$LE=='LE', all_ranks$Rank, NA) } 

    #.. color ranks    
    qua = quantile(abs(rnk), 0.95)
    newrnk = ifelse(rnk > qua, qua, rnk); newrnk[newrnk < -qua] = -qua
    all_ranks$Ranklimit = newrnk

    #.. running score line color
    if(ES_colo == "green") {
            line_colo = "green2"
        } else {
            line_colo = colorMargin }

    #.. top line type
    topL = ifelse(ES_sign==1 , max(y), min(y))
    which_topL = which(y==topL)
    df_topL = data.frame(x1 = x[which_topL], y1 = topL, x2 = x[which_topL], y2=0) 
    
    #.. base text size
    base = 12
    
    # generate rank subplot (barplot + heatbar)

    p <- ggplot( all_ranks, aes( x=Index, y=Rank ) ) +
            
            geom_bar(stat="identity", aes(fill=LE), width=10, order=order, color=NA, show.legend=FALSE) + 
        
            scale_fill_manual(values=c("Outside"="#D8D8D855","LE"=paste(colorMargin))) +  
            
            scale_y_continuous(limits=c( min(rnk), max(rnk) )) +        
            
            xlab("Rank") + ylab("Gene score") + 
            
            annotate(geom="text", x=topx, y=topy, label=top, angle=angletop, color=colorMargin, hjust=h, size=cex_top, vjust=v) + 
            
            theme_bw() + theme( text = element_text(size = base+3), axis.text = element_text(size = base)) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            
            annotate(geom="text", x=nx, y=ny, label=sprintf("N=%g",sum(all_ranks$LE=="LE")), color=colorMargin, hjust=hn, size=3) +
            
            geom_segment(data=df_arrow, aes(x = x1, y = y2, xend = x2, yend = y1), arrow =arrow(length = unit(0.2, "cm"), angle=25, type='closed'), size=0.3, color='black', inherit.aes=FALSE) +            
            annotate(geom="text", x=zero+0.5, y=0+yrange/11, label=paste('zero crossed at', zero+1), color="black", size=3)
                    
    if(image_size == 'reduced'){
        
            p = p +  geom_segment(aes(xend = 1, y = 0, x = length(rnk), yend = 0), col="black", size=0.1) +        
            annotate(geom='text',x=-Inf,y= -Inf, label="+", hjust=-0.4, vjust=-0.2, size=7) +
            annotate(geom='text',x=Inf, y=-Inf, label="_", hjust=1.7, vjust=-1, size=5.5, fontface='bold')
    }

    if(length(rnk) >= 1000){
        
        p = p + scale_x_continuous(position='bottom', limits=c(0,length(rnk)), labels = function(l) {trans=l/1000; paste0(trans, "K")})
    
    } else {
        
        p = p + scale_x_continuous(position='bottom', limits=c(0,length(rnk)))
    }
    
    pRank <- p +  theme(plot.margin = unit(c(l=0,r=0.03,t=0,b=0), "npc"))

    # generate running ES sublot
    
    q <- ggplot(gs_scores, aes(x=x, y=y)) + geom_line(color=line_colo) + 
        
        theme_bw() +  theme( text = element_text(size = base+3), axis.text = element_text(size = base)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +

        scale_y_continuous(expand=expand_scale(add = c(0.1, 0)) ) +

        #annotate(geom='text',-Inf, Inf, label="+", hjust=-0.1, vjust=0, size=7) +            
        #annotate(geom='text',Inf, Inf, label="\u2013", hjust=1.2, vjust=0, size=4, fontface='bold') +

        ylab("Enrichment Score (ES)") + xlab("Rank") +
    
        geom_hline(yintercept=0)

        if(image_size == 'reduced') {

            q = q + geom_rug(data=gs_ranks, aes(x=xranks), sides='b', show.legend=FALSE, length = unit(0.04, "npc"), inherit.aes=FALSE, size=0.3)

        } else {

            q = q + geom_rug(data=gs_ranks, aes(x=xranks), sides='b', show.legend=FALSE, length = unit(0.08, "npc"), inherit.aes=FALSE, size=0.3) +
            
            geom_rug(data=all_ranks, aes(x=Index, color=Ranklimit), sides='b', show.legend=FALSE, length = unit(0.05, "npc"), inherit.aes=FALSE) +         
            scale_color_gradient2(high=color.margins()['up'], low=color.margins()['dn'], mid=color.margins()['md'], midpoint=0, limits=c(-1,1)*qua )
        }

        if(line_top == 'horizontal'){
            
            q = q + geom_hline(yintercept=topL, linetype='dashed')
        
        } else if(line_top == "coordinate") {

            q = q + geom_segment(data=df_topL, aes(x = x1, y = y2, xend = x2, yend = y1), size=0.3, color=line_colo, linetype='dashed', inherit.aes=FALSE)
        }

        if(length(rnk) >= 1000){
            
            q = q + scale_x_continuous(position='bottom', limits=c(0,length(rnk)+1), labels = function(l) {trans=l/1000; paste0(trans, "K")})
        
        } else { 
    
            q = q + scale_x_continuous(position='bottom', limits=c(0,length(rnk)+1))
                     
        }
        pRunes <- q + theme(plot.margin = unit(c(l=0,r=0.03,t=0,b=0), "npc")) +        
        annotate(geom='text',x=-Inf,y= -Inf, label="+", hjust=-0.4, vjust=-0.2, size=7) +
        annotate(geom='text',x=Inf, y=-Inf, label="_", hjust=1.7, vjust=-1, size=5.5, fontface='bold')

    # grobs + output table
    
    all_ranks$Rank = round(ranks$genescores[match(rownames(all_ranks), ranks$geneid)],10)
    gsea_out = tibble::rownames_to_column(all_ranks, 'geneid') %>% dplyr::filter(Index %in% gs_ranks$xranks) %>% dplyr::mutate("leadingEdge"=ifelse(LE=='Outside',FALSE,TRUE)) %>% dplyr::inner_join(es_data, by=c('Index'='rank'))  %>% dplyr::rename('geneScore'='Rank', 'geneRank'="Index", "runningES"="running_es", "gene"='geneid') %>% dplyr::mutate(contrast=gsea$contrast, pathway=gsea$pathway) %>% dplyr::select(contrast,pathway, geneRank, runningES, gene, leadingEdge, geneScore) 
    
    return(list( pRunes = pRunes, pRank = pRank, gsea_out = gsea_out ))
}

## prep sample metadata
prep.metadata <- function(meta, gsea, samples, Sample, Group, transformation, dropREF, reference, own_palette, label_palette) {
    
    meta <- dplyr::select(meta, c(as.name(Sample),as.name(Group))) %>% dplyr::rename("Sample"=Sample, "Group"=Group)
    meta <- meta[ na.omit(match(samples, meta$Sample)), ]
    groups = unlist(strsplit(gsea$contrast, "-"))
    meta <- meta %>% dplyr::filter(Group %in% groups)
    if (reference) { reference_group = groups[2] }
    
    if (!reference) {

        if (grepl('reference', transformation)) {
            stop("\nERROR: Reference phenotype equels FALSE, but gene transformation with reference phenotype selected (Heatmap parameters)\n")
        
        } else {
            meta$Condition <- rep('Experiment', nrow(meta))
        }
    
    } else if (reference) {

        if (! grepl('reference', transformation)) {
            stop("\nERROR: Reference phenotype equels TRUE, but gene transformation with reference phenotype not selected (Heatmap parameters)\n")
        
        } else {

            if (reference_group %in% meta$Group) {            
                meta$Condition <- ifelse(meta$Group==reference_group, 'Reference', 'Experiment')
        
            } else if (! reference_group %in% meta$Group) { 
                stop(sprintf("\nERROR: Reference phenotype (%s) not found \nSUGGESTION: check if input datasets are correct that is sample metadata (Sample id and Group id), gene expression (column names), and (GSEA contrasts)\n", reference))
            }
        }
    }  
     
    color_label <- factor(meta$Group)
    if (label_palette == "Custom") {
        colors <- own_palette
        if( length(colors) < length(levels(color_label)) & dropREF == TRUE) {
            colors <- c(colors,rep('grey'),1)
        } else if(length(colors) != length(levels(color_label)) ){
            stop('ERROR: select number of colors in Custom palette equal to the number of phenotype labels')
        }
    } else {        
        n = ifelse(is.element(label_palette, c('Paired','Set3')), 12, 8)
        n_color <- length(levels(color_label))
        colors <- colors <- brewer.pal(n, label_palette)[1:n_color]
    }
    levels(color_label) <- colors
    meta$Group_color <- paste(color_label)
    
    return(meta)
}

## transfer gene expression by row with R (input r data.frame)
transform.gex <- function(df, transformation, meta) {

    if ( ! all(colnames(df) == meta$sample) ) { 
        
        stop("\nERROR: column names in gene expression dataset and sample metadata are not matched correctly: contact template maintainer at michaloa@mail.nih.gov")
    }
                 
    if (transformation == 'median centering') {
        mat <- t(apply(df, 1, function(y) y-median(y, na.rm=TRUE)) )

    } else if ( transformation == 'mean centering') {
        mat <- t(apply(df, 1, function(y) y-mean(y, na.rm=TRUE)) )
        
    } else if (transformation == 'z-score') {
        mat <- t(apply(df, 1, function(y) (y-mean(y, na.rm=TRUE))/sd(y, na.rm=TRUE)) )
        
    } else if (transformation == 'reference median centering' & !is.null(meta)) {
        mat <- df - apply(df[, meta$Condition == 'Reference'], 1, median, na.rm=TRUE)
        
    } else if (transformation == 'reference mean centering' & !is.null(meta) ) {
        mat <- df - apply(df[, meta$Condition == 'Reference'], 1, mean, na.rm=TRUE)
        
    } else {
        mat <- df
    } 

    return(mat)
}

## clustering order
prep.clustOrder <- function(df, linkage, distance, way) {

    rowv=FALSE
    colv=FALSE

    if( (way == 'rows'| way == 'rows and columns') ){

    if (distance=='1-Spearman') {    
        rowv = hclust(as.dist(1 - cor(t(df), method='spearman')), method = linkage)
    } else if(distance == "1-Pearson") {
        rowv = hclust(as.dist(1 - cor(t(df), method='pearson')), method = linkage)
    } else if (distance == 'Euclidean'){
        rowv = hclust(dist(df), method = linkage)
    } else { 
        rowv = FALSE
    }}

    if( (way == 'columns'| way == 'rows and columns') ){
    
    if (distance=='1-Spearman') {    
        colv = hclust(as.dist(1 - cor(df, method='spearman')), method = linkage)
    } else if(distance == "1-Pearson") {
        colv = hclust(as.dist(1 - cor(df, method='pearson')), method = linkage)
    } else if (distance == 'Euclidean'){
        colv = hclust(dist(t(df)), method = linkage)
    } else { 
        colv = FALSE
    }}

    return(list(rowv=rowv,colv=colv))
}

# plot heatmap
plot.heatmap <- function(gex, meta, heat_colors=NULL, limit, rowv, colv, transformation, show_rownames, show_colnames, show_coldend, show_rowdend, row_size, col_size, heatmap_leg, sample_leg) {
      
    meta_label = meta$Group
    meta_color = unlist(lapply(split(meta$Group_color, meta_label), unique))
    ha = HeatmapAnnotation(Class = meta_label, col=list(Class=meta_color), annotation_height = unit(rep(0.3,1), "cm"), annotation_legend_param = list(title_gp = gpar(fontsize = 7), grid_width=unit( c(0.3), "cm"), labels_gp = gpar(fontsize = 8)))

    pal = color.margins()
    if(is.null(limit)) {
        limit = quantile(abs(as.matrix(gex)), 0.95)
    } else { limit = limit }
    palette_function <- circlize::colorRamp2( c(-limit, 0, limit), c(pal['dn'],pal['md'],pal['up']), space="LAB")
    if(transformation=='z-score') { heat_legend ='SD' } else { heat_legend = 'log2' }
    
# set clustering matrix

    h1 = Heatmap(gex, col=palette_function, cluster_rows=rowv, cluster_columns=colv, name = heat_legend, column_title = "", column_title_gp = gpar(fontsize = 12), column_title_side = "top", show_column_names = show_colnames, column_names_gp = gpar(fontsize = col_size), row_names_gp = gpar(fontsize = row_size), show_row_names = show_rownames, column_dend_height = unit(1, "cm"), column_dend_reorder = T, row_dend_reorder=T, show_column_dend=show_coldend, show_row_dend=show_rowdend, top_annotation=ha,heatmap_legend_param = list(color_bar = "continuous", title_gp = gpar(fontsize = 7),labels_gp = gpar(fontsize = 6)))

    gb_heatmap = grid.grabExpr(draw(h1, heatmap_legend_side='right',  annotation_legend_side='right' , show_heatmap_legend = heatmap_leg, show_annotation_legend = sample_leg) )

   return(gb_heatmap)
}     

# color ranks
color.ranks <- function(val) {
    posVal=val[val>=0]
    up=sequential_hcl(length(posVal),h=0,c.=c(180,0),l=c(30,90),power=1.5,gamma=NULL,fixup=TRUE,alpha=1)[rev(rank(posVal))]
    downVal=val[val<0]
    down= sequential_hcl(length(downVal),h=260,c.=c(90,0),l=c(30,90),power=2,gamma=NULL,fixup=TRUE,alpha=1)[rank(downVal)]
    return(c(up,down))
}
    
# set marginal colors (default: red, blue, whitesmoke, no transparency)
color.margins <- function(dn_hew=260, up_hew=0, md_hew=94, dn_c=90, up_c=180, md_c=0, dn_l=30, up_l=30, md_l=97, a=1) {
    dn = sequential_hcl(1, h = dn_hew, c. = c(dn_c), l = c(dn_l), fixup = TRUE, alpha = a)
    up = sequential_hcl(1, h = up_hew, c. = c(up_c), l = c(up_l), fixup = TRUE, alpha = a)
    md = sequential_hcl(1, h = md_hew, c. = c(md_c), l = c(md_l), fixup = TRUE, alpha = a)
    return(c(dn=dn,md=md,up=up))
}

# set font size in LE heatmap
find.fontsize <- function(n, max_size, min_size, stepdown, max_n, min_n) {
            
            numbers = seq(min_n, max_n, by=25)
            font_sizes = seq(max_size, min_size, by=-stepdown)
            index = which.min(abs(numbers-n))
            if (index > length(font_sizes)) { index = length(font_sizes) }
            return(font_sizes[index])
        }

# adjust old GSEA table (released v67 or lower)
adjust.v67 <- function(input, gS, tC, sp, db, required_columns=c("contrast","geneScore","fdr_correction_mode","collection","pathway","pval","padj","ES","NES","nMoreExtreme","size","leadingEdge","size_leadingEdge","inPathway")) {

    missing_columns = setdiff(required_columns, colnames(input))

    if ( length(missing_columns) > 0 ) {

        cat("WARNING: Output from outdated 'Preranked GSEA [CCBR]' detected.\n")
        cat(sprintf("\nWARNING: Missing columns to be added (%s)\n\t 'inPathway' column will include all genes annotated to a gene set - run the latest released version of 'Preranked GSEA [CCBR]'\n\t to return only genes mapped in the dataset.\n", paste(missing_columns, collapse=", ")))

        if ( is.null(gS) ) {
            stop("\nERROR: 'Gene score' parameter not provided (Advanced Parameters)\n")        
        } else {
            if ( length(gS) > 1 ) {
                gS = gS[1]
                cat(sprintf("\nWARNING: too many values for 'Gene score' provided; only the first one used, '%s'\n", gS))
            }
            has.underscore = substring(gS, 1, 1) == "_" 
            if ( !has.underscore ) {
                gS = paste0("_",gS)
                cat(sprintf("\nWARNING: 'Gene score' (Advanced Parameters) should start with underscore to match column naming convention in DEG table used for GSEA run (e.g. '_tstat');\n\t column 'geneScore' is assigned the provided value with '_' added ('%s');\n\t this template will fail if this is not an adequate specification.\n", gS))
                
            }
        }

        if ( is.null(tC) ) {
            tC = "NULL"
            stop("\nERROR: 'Tested contrast' parameter not provided (Advanced Parameters)\n")        
        } else {
            if ( length(tC) > 1 ) {
                tC = tC[1]
                cat(sprintf("\nWARNING: too many values for 'Tested contrast' provided; only the first one used, '%s'\n",tC))
            }
            is.contrast = any(grepl("-", tC))
            if ( !is.contrast ) {
                stop(sprintf("\nERROR: 'Tested contrast' (Advanced Parameters) should be specified exactly the same way as when Preranked GSEA was run (e.g. treated-control);\n\t provided value is '%s'", tC))
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

        return(input)
    }
}

print("template_function_GSEA_Visualization_KOvWT_Activated.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Analysis<-readRDS(paste0(rds_output,"/var_DEG_Analysis.rds"))
var_DEG_Analysis<-as.data.frame(var_DEG_Analysis)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_GSEA_Filters_KOVWT_Activated<-readRDS(paste0(rds_output,"/var_GSEA_Filters_KOVWT_Activated.rds"))
var_GSEA_Filters_KOVWT_Activated<-as.data.frame(var_GSEA_Filters_KOVWT_Activated)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Normalized_Counts<-readRDS(paste0(rds_output,"/var_Normalized_Counts.rds"))
var_Normalized_Counts<-as.data.frame(var_Normalized_Counts)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr804_metadata_for_NIDAP<-readRDS(paste0(rds_output,"/var_ccbr804_metadata_for_NIDAP.rds"))
var_ccbr804_metadata_for_NIDAP<-as.data.frame(var_ccbr804_metadata_for_NIDAP)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_msigdb_v6_2_with_orthologs<-readRDS(paste0(rds_output,"/var_msigdb_v6_2_with_orthologs.rds"))
var_msigdb_v6_2_with_orthologs<-as.data.frame(var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
var_GSEA_Visualization_KOvWT_Activated<-GSEA_Visualization_KOvWT_Activated(var_DEG_Analysis,var_GSEA_Filters_KOVWT_Activated,var_Normalized_Counts,var_ccbr804_metadata_for_NIDAP,var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
saveRDS(var_GSEA_Visualization_KOvWT_Activated, paste0(rds_output,"/var_GSEA_Visualization_KOvWT_Activated.rds"))
