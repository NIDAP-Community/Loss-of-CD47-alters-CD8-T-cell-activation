Normalized_Counts <- function(Filtered_Counts, ccbr804_metadata_for_NIDAP) {
    #image: png
    imageType = "png"
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(colorspace))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(reshape))
    suppressMessages(library(gridExtra))

    df <- Filtered_Counts
    
    samples_for_deg_analysis = c("Gene","WT_UT_1","WT_CD3_CD28_1","CD47KO_UT_1","CD47KO_CD3_CD28_1","WT_UT_2","WT_CD3_CD28_2","CD47KO_UT_2","CD47KO_CD3_CD28_2","WT_UT_3","WT_CD3_CD28_3","CD47KO_UT_3","CD47KO_CD3_CD28_3","WT_UT_4","WT_CD3_CD28_4","CD47KO_UT_4","CD47KO_CD3_CD28_4")
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis !=            "Gene"]
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis !=            "GeneName"]
    df.m <- df[,samples_for_deg_analysis]
    gene_names <- NULL
    gene_names$GeneID <- df[,1]
    targetfile <- ccbr804_metadata_for_NIDAP
    Sampinfo <- targetfile
    targetfile <- targetfile[match(colnames(df.m),targetfile$Sample),]
    targetfile <- targetfile[rowSums(is.na(targetfile)) != ncol(targetfile), ]
    df.m <- df.m[,match(targetfile$Sample,colnames(df.m))]
    if(FALSE){
    x <- DGEList(counts=2^df.m, genes=gene_names)
    } else {
    x <- DGEList(counts=df.m, genes=gene_names)     
    }
    dm.formula <- as.formula(paste("~0 +", paste(c("Genotype_Plating"), sep="+", collapse="+")))
    design=model.matrix(dm.formula, targetfile)
    colnames(design) <- str_replace_all(colnames(design), c("Genotype_Plating"), "")
    v <- voom(x,design=design,normalize="quantile")
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    print(paste0("Total number of genes included: ", nrow(df.voom)))

    samples_to_include = c("Gene","WT_UT_1","WT_CD3_CD28_1","CD47KO_UT_1","CD47KO_CD3_CD28_1","WT_UT_2","WT_CD3_CD28_2","CD47KO_UT_2","CD47KO_CD3_CD28_2","WT_UT_3","WT_CD3_CD28_3","CD47KO_UT_3","CD47KO_CD3_CD28_3","WT_UT_4","WT_CD3_CD28_4","CD47KO_UT_4","CD47KO_CD3_CD28_4")
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

    #df %>% dplyr::select(append("Gene", samples_to_include)) -> spark.df
    df.voom -> df
    df -> edf.orig

    
    # cell <- Sampinfo$Sample
    Sampinfo = Sampinfo[, colSums(is.na(Sampinfo)) == 0]
    rownames(Sampinfo) <- Sampinfo$Sample
    Sampinfo <- Sampinfo[match(colnames(edf.orig[,-1]), Sampinfo$Sample), ]
    Sampinfo = na.omit(Sampinfo)
    print(paste0("Total number of genes included: ", nrow(edf.orig)))

    edf <- edf.orig[,match(Sampinfo$Sample,colnames(edf.orig))]
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2)
    pca.df$celltype <- Sampinfo$Genotype_Plating
    pca.df$sample <- Sampinfo$Genotype_Plating_Label
    
    # Manual changes to sample names
    replacements = c("")

    plotcolors <- c("aquamarine3","lightskyblue3","salmon1","darkolivegreen3","plum3","goldenrod1","burlywood2","gray70","firebrick2","steelblue","palegreen4","orchid4","darkorange1","yellow","sienna","palevioletred1","gray60","cyan4","darkorange3","mediumpurple3","violetred2","olivedrab","darkgoldenrod2","darkgoldenrod","gray40","palegreen3","thistle3","khaki1","deeppink2","paleturquoise3","chocolate3","wheat1","lightsteelblue","salmon","sandybrown","darkolivegreen2","thistle2","gray85","orchid3","darkseagreen1","lightgoldenrod1","lightskyblue2","dodgerblue3","darkseagreen3","forestgreen","lightpink2","lightpink1","mediumpurple4","thistle","navajowhite","lemonchiffon","bisque2","mistyrose","gray95","lightcyan3","peachpuff2","lightyellow2","lightsteelblue2","moccasin","antiquewhite2","gray80","lightgrey")
    if (length(unique(Sampinfo$Genotype_Plating)) > length(plotcolors)) {
        ## Original color-picking code.
        k=length(unique(Sampinfo$Genotype_Plating))-length(plotcolors)
       
        more_cols<- getourrandomcolors(k) 
        plotcolors <- c(plotcolors , more_cols)
    }

    if (!is.null(replacements)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0("PC1 ", perc.var[1],"%")
    pc.y.lab <- paste0("PC2 ", perc.var[2],"%")
    
    labelpos <- pca.df
    labelpos$mean_y <- pca.df$PC2+2
    labelpos$mean_x <- pca.df$PC1+2
    
    if (TRUE){
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=2) +
      geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
          label=sample, color=celltype, vjust="inward", hjust="inward"), size=3, show.legend=FALSE) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)
    } else {
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=2) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
    }
    par(mfrow = c(2,1))
   
    #samples_to_include = c("Gene","WT_UT_1","WT_CD3_CD28_1","CD47KO_UT_1","CD47KO_CD3_CD28_1","WT_UT_2","WT_CD3_CD28_2","CD47KO_UT_2","CD47KO_CD3_CD28_2","WT_UT_3","WT_CD3_CD28_3","CD47KO_UT_3","CD47KO_CD3_CD28_3","WT_UT_4","WT_CD3_CD28_4","CD47KO_UT_4","CD47KO_CD3_CD28_4")
   # samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    df.filt <- df %>% dplyr::select(samples_to_include)
   
        # cell <- Sampinfo$SampleName
        rownames(Sampinfo) <- Sampinfo$Sample
        Sampinfo = Sampinfo[complete.cases(Sampinfo[, "Sample"]),]
        print(paste0("Total number of samples included: ", nrow(Sampinfo)))
    df.filt <- df.filt[,match(Sampinfo$Sample,colnames(df.filt))]
    rownames(df.filt)=df[,1]
    df.filt %>% rownames_to_column("Gene") -> df.filt
    df.m <- melt(df.filt,id.vars=c("Gene"))
    if(FALSE){
    df.m %>% mutate(value = log2(value+0.5)) -> df.m
    }
    n <- 40
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    if (length(unique(Sampinfo$Sample)) > length(cols)) {
        cols <- c(cols, rep("black", length(unique(Sampinfo$Sample)) - length(cols)))
    }

    if(FALSE){
        xmin = -1
        xmax = 1
    } else {
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

 
 if(FALSE ){

 df.m %>% mutate(colgroup = Sampinfo[variable,"Genotype_Plating"]) -> df.m
 df.m = df.m[complete.cases(df.m[, "colgroup"]),]
 df.m$colgroup = gsub("\\s","_",df.m$colgroup)
 df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
 if("Genotype_Plating"=="Genotype_Plating"){
     print("Using same colors as PCA plot")
     altplotcolors<-plotcolors
 }else{
     print("Generating new colors")
     altplotcolors<- getourrandomcolors(length(unique(df.m$colgroup))) 
 }
 print(altplotcolors)
 print(unique(df.m$variable))
# plot Density 
g2 = ggplot(df.m, aes(x=value, group=variable)) + 
    geom_density(aes(colour = colgroup, linetype = colgroup)) +
    xlab("Filtered Counts") + ylab("Density") +
    theme_bw() +
    theme(legend.position='top',legend.text = element_text(size = 10)) + #scale_x_log10() + 
    ggtitle("Frequency Histogram") +
    xlim(xmin,xmax) +
    scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
    scale_colour_manual(values=altplotcolors)
 }
 else{
     df.m$variable = Sampinfo[df.m$variable,"Genotype_Plating_Label"]
     n=length(unique(df.m$variable))

     

   
    cols<- getourrandomcolors(n) 

   
  g2 = ggplot(df.m, aes(x=value, group=variable)) + 
    geom_density(aes(colour = variable, linetype = variable)) +
    xlab("Filtered Counts") + ylab("Density") +
    theme_bw() +
    theme(legend.position='top',legend.text = element_text(size = 10)) + #scale_x_log10() + 
    ggtitle("Frequency Histogram") +
    xlim(xmin,xmax) +
    scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
    scale_colour_manual(values=cols)#+
    guides(linetype = guide_legend(ncol = 6))
 }

    imageWidth = 3000
    imageHeight = 1500*2
    dpi = 300

    ## Choice of image output format: PNG or SVG.
    dev.off()
    if (imageType == 'png') {
    png(
      filename="Normalized_Counts.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
    } else {
        library(svglite)
        svglite::svglite(
        file="Normalized_Counts.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    require(gridExtra)
     gh<-make_heatmap(df.filt[,samples_to_include],targetfile,plotcolors)
 grid.arrange(g,g2,gh, nrow=2)
  dev.off()

    return(df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

#################################################
## Global imports and functions included below ##
#################################################

#################################################
## Global imports and functions included below ##
#################################################
make_heatmap <- function(df, metadata,plotcolors) {
  
  suppressMessages(library(amap))
  suppressMessages(library(lattice))
  suppressMessages(library(gplots))
  suppressMessages(library(gridGraphics))
  suppressMessages(library(dendsort))
  
  
  mat <- as.matrix(df) 
  tcounts=t(mat)
  d=Dist(tcounts,method="correlation",diag=TRUE)
  dend = rev(dendsort(as.dendrogram(hclust( d,method="average"))))
  m=as.matrix(d)
  targetfile <- metadata
  rownames(targetfile) = targetfile$Sample
  idx = as.factor(targetfile[rownames(m),"Genotype_Plating"])
  col = plotcolors
  cols <- col[idx]
  new.palette=colorRampPalette(c("blue","green","yellow"),space="rgb")
  #levelplot(m[1:ncol(m),ncol(m):1],col.regions=new.palette(20))

   

  mk<-function(){
        if(length(colnames(m))>20){
            par(mar=c(0,0,0,0))
    
            heatmap.2(m,labRow = NA, labCol = NA,col=new.palette(20),trace="none",colRow = col[idx], colCol = col[idx],rowDendrogram=dend,colDendrogram=dend,RowSideColors = col[idx],ColSideColors = col[idx],dendrogram = "row",cexRow=3,cexCol=3,margins=c(0,0),   lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), lhei=c( .2,4,2), lwid=c(1, .2,4 ),
              key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),key.xlab = "Correlation",key.ylab="Count",key.title=" ")       }else{
            heatmap.2(m,col=new.palette(20),trace="none",colRow = col[idx], colCol = col[idx],rowDendrogram=dend,colDendrogram=dend,RowSideColors = col[idx],ColSideColors = col[idx],dendrogram = "row",cexRow=3,cexCol=3,margins=c(4,1),  lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), lhei=c( .2,4,2), lwid=c(1, .2,4),
              key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),key.xlab = "Correlation",key.ylab="Count",key.title=" ")
        }
    #heatmap.2(m,col=new.palette(20),trace="none",colRow = col[idx], colCol = col[idx],rowDendrogram=dend,colDendrogram=dend,RowSideColors = col[idx],ColSideColors = col[idx],dendrogram = "row",cexRow=3,cexCol=3,margins=c(4,1),  lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ))
    
    #  heatmap.2(m,col=new.palette(20),trace="none",dendrogram = "row",margins=c(4,1),  lmat=rbind( c(0, 4, 3), c(2,1,0 ) ), lwid=c(1, 4, 2 ),lhei=c(1,4),cexRow=3,cexCol=3) 

  }
  
  tg<-mk()
  grid.echo(mk)
  gh1<-grid.grab()
  mklegend<-function(){
      plot.new()
          legend(x="top", legend=levels(idx), col=col[as.factor(levels(idx))],pch=15,x.intersp=3,bty ="n",cex=2)
  }
  grid.echo(    mklegend )
  gh2<-grid.grab()
  lay <- c(1,3)
grid.newpage()
grid.arrange(gh1,gh2,nrow=1,widths=c(unit(1000, "bigpts"),unit(300, "bigpts")))
    gh<-grid.grab()
  return(gh)
  
}
 getourrandomcolors<-function(k){
        seed=10
        n <- 2e3
        ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
        ourColorSpace <- as(ourColorSpace, "LAB")
        currentColorSpace <- ourColorSpace@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)

        km <- kmeans(currentColorSpace, k, iter.max=20)
        return( unname(hex(LAB(km$centers))))
    }

print("template_function_Normalized_Counts.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Filtered_Counts<-readRDS(paste0(rds_output,"/var_Filtered_Counts.rds"))
var_Filtered_Counts<-as.data.frame(var_Filtered_Counts)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr804_metadata_for_NIDAP<-readRDS(paste0(rds_output,"/var_ccbr804_metadata_for_NIDAP.rds"))
var_ccbr804_metadata_for_NIDAP<-as.data.frame(var_ccbr804_metadata_for_NIDAP)
invisible(graphics.off())
var_Normalized_Counts<-Normalized_Counts(var_Filtered_Counts,var_ccbr804_metadata_for_NIDAP)
invisible(graphics.off())
saveRDS(var_Normalized_Counts, paste0(rds_output,"/var_Normalized_Counts.rds"))
