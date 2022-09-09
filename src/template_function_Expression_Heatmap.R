Expression_Heatmap <- function(Normalized_Counts,ccbr804_metadata_for_NIDAP) {
    #This function uses pheatmap to draw a heatmap, scaling first by rows
    #(with samples in columns and genes in rows)

    # image: png
    suppressMessages(library(colorspace))
    suppressMessages(library(dendsort))
    suppressMessages(library(pheatmap))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(dplyr))
    suppressMessages(library(grid))
    suppressMessages(library(gtable))
    suppressMessages(library(gridExtra))
    suppressMessages(library(gridGraphics))

    # Toggle sets output matrix (TRUE = Zscores, FALSE = expression).
    return_zscores_toggle <- FALSE

    pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
        if (n < 1L) 
            return(character(0L))
        h <- rep(h, length.out = 2L)
        c <- c[1L]
        l <- rep(l, length.out = 2L)
        power <- rep(power, length.out = 2L)
        rval <- seq(1, -1, length = n)
        rval <- hex(
            polarLUV(
                L = l[2L] - diff(l) * abs(rval)^power[2L], 
                C = c * abs(rval)^power[1L],
                H = ifelse(rval > 0, h[1L], h[2L])
            ),
            fixup=fixup, ...
        )
        if (!missing(alpha)) {
            alpha <- pmax(pmin(alpha, 1), 0)
            alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                width = 2L, upper.case = TRUE)
            rval <- paste(rval, alpha, sep = "")
        }
        return(rval)
    }
    #Color selections for heatmap:
    np0 = pal(100)
    np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1)  #Blue to Red
    np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2))  #Red to Vanilla
    np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) #Violet to Pink
    np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))
    np5 = colorRampPalette(c("steelblue","white", "red"))(100) #Steelblue to White to Red

    np = list(np0, np1, np2, np3, np4, np5)
    names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")

    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col) {
        require(pheatmap)
        require(dendsort)
        
        col.pal <- np[[col]]
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # define metrics for clustering
        drows1 <- "correlation"
        dcols1 <- "correlation"
        minx = min(dat)
        maxx = max(dat)

        if (TRUE) {
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            #absmax = ceiling(max(abs(c(minx, maxx))))
            #breaks = c(-1*absmax, seq(0, 1, length=98), absmax)
            #legbreaks = c(-1*absmax, 0, absmax)
            breaks = seq(-2, 2, length=100)
            legbreaks = seq(-2, 2, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)

        #Run cluster method using 
        hc = hclust(dist(t(dat)), method="average")
        hcrow = hclust(dist(dat), method="average")
        
        if (FALSE) {
            sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
        } else {
            sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        }
        if (clus) {
            colclus <- sort_hclust(hc)
        } else {
            colclus = FALSE
        }
        if (clus2) {
            rowclus <- sort_hclust(hcrow)
        } else {
            rowclus = FALSE
        }
        #print('sorted the clusters')
        if (TRUE) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        hm.parameters <- list(
            dat, 
            color=col.pal,
            legend_breaks=legbreaks,
            legend=FALSE,
            # cellwidth=15, 
            # cellheight=10, 
            scale="none",
            treeheight_col=treeheight,
            treeheight_row=treeheight,
            kmeans_k=NA,
            breaks=breaks,
            # height=80,
            fontsize=10,
            fontsize_row=8,
            fontsize_col=8,
            show_rownames=rn, 
            show_colnames=cn,
            cluster_rows=rowclus, 
            cluster_cols=clus,
            cutree_rows=1,
            clustering_distance_rows=drows1, 
            clustering_distance_cols=dcols1,
            annotation_col = annotation_col,
            annotation_colors = annot_col,
            labels_col = labels_col
        )

        mat = t(dat)
        # print('calculated mat')
        
        callback = function(hc, mat) {
            # print('inside the callback')
            dend=rev(dendsort(as.dendrogram(hc)))
            # print ('reversed the dendsorted hc')
            dend %>% dendextend::rotate(c(1:nobs(dend))) -> dend
            as.hclust(dend)
        }
        phm<-do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
        
        ##Do color spectrum legend
        mat_breaks = breaks
        mat_legend = legbreaks
        #mat_breaks=seq(1,100, length.out = 100)
        #mat_legend=seq(0,100, length.out = 11)
        names(mat_legend)<-mat_legend
        mat_color             = colorRampPalette(c("steelblue","white", "red"))(length(mat_breaks)-1)
        dl<-draw_legend(color=col.pal,breaks=mat_breaks,legend=mat_legend,fontsize=5)
        
        ###Do group plots
       # print(mat_legend)  
       # print(mat_breaks)
       # print(col.pal)     
        #print(unique(annotation_col)) 
       # print(annot_col)
        legd<- draw_annotation_legend(annotation = annot_col, annotation_colors = annot_col,border_color = "black", fontsize = 20)

        ###put them together

        tg <- phm[[4]]
        print(tg)

        tg<-gtable_filter(tg, "legend",invert=TRUE,trim=TRUE)

        print(tg)

        #######
        tg
        grid.newpage()
        grid.draw(tg)
        eg<-grid.grab()
        grid.newpage()
        grid.draw(dl)
        egdl<-grid.grab()
        grid.newpage()
        grid.draw(legd)
        egdlkey<-grid.grab()
        lay <- rbind(c(1,NA,2,3),
                    c(1,NA,2,3))
        grid.newpage()
        grid.arrange(grobs=list(eg,egdlkey,egdl),
                    widths=c(unit(2000, "bigpts"),unit(50, "bigpts"),unit(250, "bigpts"),unit(250, "bigpts")),
                            heights=c(unit(1000, "bigpts"),unit(1000, "bigpts")),layout_matrix = lay)
    }

    df1 <- Normalized_Counts
    if("Gene"!="Gene"){
        #drop original Gene column
        df1 = df1[,!(colnames(df1)%in% c("Gene")) ]
        #rename column to Gene
        colnames(df1)[which(colnames(df1) == "Gene")] <- 'Gene'

    }
    samples_to_include = c("WT_UT_1","WT_CD3_CD28_1","CD47KO_UT_1","WT_UT_2","WT_CD3_CD28_2","CD47KO_UT_2","CD47KO_CD3_CD28_2","WT_UT_3","WT_CD3_CD28_3","CD47KO_CD3_CD28_3","CD47KO_UT_3","WT_UT_4","WT_CD3_CD28_4","CD47KO_UT_4","CD47KO_CD3_CD28_4","CD47KO_CD3_CD28_1")
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

    genes_to_include_parsed = c()
    genes_to_include_string = "RP23-186O3.13 Il2 Cpeb1 Spats2 Nkain1 Tox2 Gm11224 E130218I03Rik Cd209a Prokr1 Fzd9 Gm37790 Spin4 Ccdc85c Frmd4b Lif Dlc1 Tbkbp1 1700001L05Rik Nkx1-1 Phactr2 Batf2 A430072P03Rik Gm37787 Atg9b Aldh1l2 Galnt3 Phlda3 Nrg1 2810429I04Rik Pcolce BC048644 Hivep3 Gm37422 Spire1 Adora2b Hoxb9 Eea1 P2ry2 Pros1 Vcan Abca9 Bag2 Irf4 Abhd15 Pdzk1 Ak7 Kctd12b Cited4"
    genes_to_include_parsed = strsplit(genes_to_include_string, " ")[[1]]
        
    if (FALSE == FALSE) {
       # spark_df %>% filter(`%in%`(spark_df$Gene, genes_to_include_parsed)) -> spark_df
        df1[ df1$Gene %in% genes_to_include_parsed,] -> df1
    }
    
    #spark_df %>% dplyr::select(append("Gene", samples_to_include)) -> spark_df
    df1 <- df1[,append("Gene", samples_to_include)]

    df.orig = df1 #dplyr::collect(spark_df, stringsAsFactors=FALSE)
    df.orig %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean)) -> df
    df.mat = df[ , (colnames(df) != "Gene" )] %>% as.data.frame
    df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
    df %>% dplyr::mutate(Gene = stringr::str_wrap(Gene,50)) -> df
    row.names(df.mat) <- df$Gene
    df.mat <- as.data.frame(df.mat)
    
    filter_by_variance = FALSE
    if (filter_by_variance) {
        df.mat = as.matrix(df.mat)
        var <- matrixStats::rowVars(df.mat)
        df <- as.data.frame(df.mat)
        rownames(df) <- rownames(df.mat)
        df.mat <- df
        df.mat$var <- var
        df.mat %>% rownames_to_column("Gene") -> df.mat 
        df.mat %>% dplyr::arrange(desc(var)) -> df.mat
        df.mat.extra.genes = dplyr::filter(df.mat, Gene %in% genes_to_include_parsed)
        df.mat = df.mat[1:500,]
        df.mat = df.mat[complete.cases(df.mat),]
        df.mat <- rbind(df.mat, df.mat.extra.genes)
        df.mat <- df.mat[!duplicated(df.mat),] 
        rownames(df.mat) <- df.mat$Gene
        df.mat$Gene <- NULL
        df.mat$var <- NULL
    }
    

    manually_replace_sample_names = FALSE
    if (manually_replace_sample_names) {
        replacements = c("")
        old <- c()
        new <- c()
        for (i in 1:length(replacements)) {
            old[i] <- strsplit(replacements[i], ": ?")[[1]][1]
            new[i] <- strsplit(replacements[i], ": ?")[[1]][2]
        }
        df.relabel <- as.data.frame(cbind(old, new), stringsAsFactors=FALSE)
        labels_col %>% replace(match(df.relabel$old, labels_col), df.relabel$new) -> labels_col
    }

    annot <- ccbr804_metadata_for_NIDAP
    annot %>% dplyr::filter(Sample %in% samples_to_include) -> annot

    if(TRUE){
      annot %>% dplyr::arrange_(.dots=c("Genotype_Plating")) -> annot
      df.mat <- df.mat[,match(annot$Sample,colnames(df.mat))] 
    }

    groups = c("Genotype_Plating")

    annot %>% dplyr::select(groups) -> annotation_col    
    annotation_col = as.data.frame(unclass(annotation_col))
    annotation_col[] <- lapply(annotation_col,factor)
    rownames(annotation_col) <- annot$Sample
    annot_col = list()
    colors <- c("greenyellow","darkviolet","darkorange","darkturquoise","darkblue","chocolate","deeppink","royalblue","darkgoldenrod","red")
    b=1
    i=1
    while (i <= length(groups)){
        nam <- groups[i]
        grp <- as.factor(annotation_col[,i])
        c <- b+length(levels(grp))-1
        col = colors[b:c]
        names(col) <- levels(grp)
        assign(nam,col)
        annot_col = append(annot_col,mget(nam))
        b = b+c
        i=i+1
    }

    print(paste0("The total number of genes in heatmap: ", nrow(df.mat)))

    labels_col <- colnames(df.mat)

    imageWidth = 3000
    imageHeight = 1500
    dpi = 300

    png(
      filename="Expression_Heatmap.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    # Optionally apply centering and rescaling (default TRUE).
    if (TRUE) {
            tmean.scale = t(scale(t(df.mat)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
    } else {
            tmean.scale=df.mat
        }

    doheatmap(dat=tmean.scale, clus=FALSE, clus2=TRUE, ht=50, rn=TRUE, cn=TRUE, col="Blue to Red")
    
    # If user sets toggle to TRUE, return Z-scores.
    # Else return input counts matrix by default (toggle FALSE).
    # Returned matrix includes only genes & samples used in heatmap.
    if(return_zscores_toggle){
        df.new <- data.frame(tmean.scale) # Convert to Z-scores.
        df.new %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    } else {
        df.mat %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    }
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

draw_annotation_legend = function(annotation, annotation_colors, border_color, fontsize,...){

  y = unit(1, "npc")
  text_height = unit(4, "grobheight", textGrob("FGH",gp=gpar(fontsize=10,...)))
  #text_height = min(unit(1, "npc"), unit(300, "bigpts"))
   fontsize=15
  res = gList()
  
  for(i in names(annotation)){
      print(i)
    res[[i]] = textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontsize=fontsize,fontface="bold"))
    
    y = y - 1.5 * text_height
    if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
      n = length(annotation_colors[[i]])
      yy = y - (1:n - 1) * 2 * text_height
      res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = 2 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]]))
      res[[paste(i, "t")]] = textGrob(names(annotation_colors[[i]]), x = text_height * 2.4, y = yy - text_height, hjust = 0, vjust = 0.5, gp = gpar(fontsize=fontsize))
      
      y = y - n * 2 * text_height
      
    }
    else{
      yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
      h = 8 * text_height * 0.25
      
      res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = 2 * text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
      res[[paste(i, "r2")]] = rectGrob(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = 8 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = NA))
      txt = rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
      yy = y - c(1, 7) * text_height
      res[[paste(i, "t")]]  = textGrob(txt, x = text_height * 2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(fontsize,...))
      y = y - 8 * text_height
    }
    y = y - 1.5 * text_height
  }
  
  res = gTree(children = res)
  
  return(res)
}

draw_legend = function(color, breaks, legend, ...){
  color = color[!is.infinite(breaks)]
  breaks = breaks[!is.infinite(breaks)]
  
  height = min(unit(1, "npc"), unit(300, "bigpts"))
  
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height*1.02)
  
  breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)
  
  h = breaks[-1] - breaks[-length(breaks)]
  #print(color)
  rect = rectGrob(x = 0, y = breaks[-length(breaks)]- height*.02, width = unit(40, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  text = textGrob(names(legend), x = unit(44, "bigpts"), y = legend_pos, hjust = 0,gp = gpar(fontsize = 10))
  
  res = grobTree(rect, text)
  
  return(res)
}

print("template_function_Expression_Heatmap.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Normalized_Counts<-readRDS(paste0(rds_output,"/var_Normalized_Counts.rds"))
var_Normalized_Counts<-as.data.frame(var_Normalized_Counts)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr804_metadata_for_NIDAP<-readRDS(paste0(rds_output,"/var_ccbr804_metadata_for_NIDAP.rds"))
var_ccbr804_metadata_for_NIDAP<-as.data.frame(var_ccbr804_metadata_for_NIDAP)
invisible(graphics.off())
var_Expression_Heatmap<-Expression_Heatmap(var_Normalized_Counts,var_ccbr804_metadata_for_NIDAP)
invisible(graphics.off())
saveRDS(var_Expression_Heatmap, paste0(rds_output,"/var_Expression_Heatmap.rds"))
