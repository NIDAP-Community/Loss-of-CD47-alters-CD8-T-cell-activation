Filtered_Counts <- function(ccbr804_rawcounts_for_NIDAP_csv, ccbr804_metadata_for_NIDAP) {
   imageType = "png"
  spark.df = ccbr804_rawcounts_for_NIDAP_csv
    samples_to_include = c("WT_UT_1","WT_CD3_CD28_1","CD47KO_UT_1","CD47KO_CD3_CD28_1","WT_UT_2","WT_CD3_CD28_2","CD47KO_UT_2","CD47KO_CD3_CD28_2","WT_UT_3","WT_CD3_CD28_3","CD47KO_UT_3","CD47KO_CD3_CD28_3","WT_UT_4","WT_CD3_CD28_4","CD47KO_UT_4","CD47KO_CD3_CD28_4")
    gene_col="Gene"
    sample_col="Sample"
    metadata <- ccbr804_metadata_for_NIDAP
    group_col<- "Genotype_Plating"
    plot_labels_col<-"Genotype_Plating_Label"
    OPT_Minimum_Count_Value_to_be_Considered_Nonzero<-1
    OP_Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total<-1
    OPT_Minimum_Number_of_Samples_with_Nonzero_Counts_in_a_Group<-3
    outlier_samples_to_remove = c()
    OPT_Use_Group_Based_Filtering<-TRUE
    opt_use_cpm_counts_to_filter<-TRUE
    OPT_Principal_Component_on_X_axis<-1
    opt_pcx<-paste0("PC",OPT_Principal_Component_on_X_axis)
    OPT_Principal_Component_on_Y_axis<-2
    opt_pcy<-paste0("PC",OPT_Principal_Component_on_Y_axis)
    OPT_Add_Labels<-TRUE
    OPT_Label_Offset_Y_axis_<-2
    OPT_Label_Offset_X_axis_<-2
    OPT_Point_Size<-1
    OPT_Label_Font_Size<-3
    plotcolors <- c("aquamarine3","salmon1","lightskyblue3","plum3","darkolivegreen3","goldenrod1","burlywood2","gray70","firebrick2","steelblue","palegreen4","orchid4","darkorange1","yellow","sienna","palevioletred1","gray60","cyan4","darkorange3","mediumpurple3","violetred2","olivedrab","darkgoldenrod2","darkgoldenrod","gray40","palegreen3","thistle3","khaki1","deeppink2","chocolate3","paleturquoise3","wheat1","lightsteelblue","salmon","sandybrown","darkolivegreen2","thistle2","gray85","orchid3","darkseagreen1","lightgoldenrod1","lightskyblue2","dodgerblue3","darkseagreen3","forestgreen","lightpink2","mediumpurple4","lightpink1","thistle","navajowhite","lemonchiffon","bisque2","mistyrose","gray95","lightcyan3","peachpuff2","lightsteelblue2","lightyellow2","moccasin","antiquewhite2","gray80","lightgrey")
    opt_color_by_group <-FALSE 
    opt_group_to_color<-"Genotype_Plating"
    opt_legend_position<-"top"
    opt_legend_osition<-'top'
    opt_legend_font_size<-10
    Opt_Set_Min_Max_for_X_axis<-FALSE
    Opt_Maximum_for_X_axis<-1
    Opt_Minimum_for_X_axis<--1
    opt_perform_log_transform<-TRUE
    opt_num_legend_cols<- 6
    OPT_Number_of_Image_Rows<-2
    replacements = c("")
  
    # packages necessary to compute remove low count genes 
    suppressMessages(library(tibble))
    suppressMessages(library(magrittr))
    suppressMessages(library(dplyr))
    suppressMessages(library(reshape2))
    suppressMessages(library(edgeR))
    # packages necessary to display the PCA and Freq/Density plots
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(colorspace))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(gridExtra))
  suppressMessages(library(gridGraphics))
    # suppress warnings
    options(warn=-1)

# get the sample data
  

    # remove outliers
   
    samples_to_include <- samples_to_include[! samples_to_include %in% outlier_samples_to_remove]
    #
    samples_to_include <- samples_to_include[samples_to_include != gene_col]
    spark.df[,append(gene_col, samples_to_include)]-> spark.df

#### remove low count genes ########
    df <- spark.df[complete.cases(spark.df),]
    #df <- dplyr::collect(df)
    df %>% dplyr::group_by_at(gene_col) %>% summarize_all(sum) -> df
    print(paste0("Number of genes before filtering: ", nrow(df)))
    df$isexpr1 <- rowSums(edgeR::cpm(df[, -1]) > 1) >= 1
    df <- as.data.frame(df[df$isexpr1, ])
    df$isexpr1 <- NULL

###################################################
#This code block does input data validation

metadata <- metadata[metadata[,sample_col]%in%samples_to_include,]
print(nrow(metadata))
  Sampinfo <- metadata

   if(make.names(colnames(df))!=colnames(df)){
        print("Error: The following counts matrix column names are not valid:\n")
        print(colnames(df)[make.names(colnames(df))!=colnames(df)])

        print("Likely causes are columns starting with numbers or other special characters eg spaces.")
        stop("Bad column names.")
    }
    if(sum(grepl("-",metadata[,sample_col]))!=0){
        print("The sample names cannot contain dashes.")
        print(metadata[,sample_col][grepl("-",metadata[,sample_col])])
        stop("No dashes allowed in column names")
    }
   if(sum(duplicated(colnames(df)))!=0){
        print("Duplicate column names are not allowed, the following columns were duplicated.")
        colnames(df)[duplicated(colnames(df))]
        stop("Duplicated columns")
    }
    if(sum(duplicated(metadata[,sample_col]))!=0){
        print("Duplicate sample names are not allowed, the following samples were duplicated.")
        metadata[,sample_col][duplicated(metadata[,sample_col])]
        stop("Duplicated samples")
    }
    print(paste0("There are ",length(metadata[,sample_col])," samples"))
    if(sum(!(metadata[,sample_col] %in% colnames(df) ))!=0 ){
        print("All metadata sample names should match one column")
        print("the following was an additional column")
        print(paste0("Metadata sample column: ",sample_col))
        print(metadata[,sample_col][!(metadata[,sample_col] %in% colnames(df) )])
        stop("Metadata data mismatch")
    }
    
    #geneids<-df[,gene_col]
    #Original_IDs, Gene, Ensembl_ID
    if(sum(!grepl("ENSG\\d{11}", df[,gene_col]))==0){
        print("ENSEMBL ids found.")
        #ENSEMBL ids
        if(sum(grepl("\\|", df[,gene_col]))>0){
            print("ENSEMBL ids found containing bars.")
            #with bars
            #df[,"Orginal_IDS"]<-df[,gene_col]
            #df[,"Ensembl_ID"]<-gsub("(.*?)(\\|.*)", "\\1",df[,gene_col])
            print("Removing ENSEMBL prefix")
            df[,gene_col]<-gsub(".*\\|","",df[,gene_col])
            print("Aggregating the counts for the same gen in different chromosome locations.")
            df<-aggregate(as.formula(paste( ".~ ",gene_col)), df, sum)
            print("Done modifying gene ids.")
         }else if(sum(grepl("_", df[,gene_col]))>0){
            #with underscores
            print("ENSEMBL ids found containing underscores.")
            #df[,"Orginal_IDS"]<-df[,gene_col]
           # df[,"Ensembl_ID"]<-gsub("(.*?)(_.*)", "\\1",df[,gene_col])
            print("Removing ENSEMBL prefix")
            df[,gene_col]<-gsub(".*_","",df[,gene_col])
            print("Aggregating the counts for the same gen in different chromosome locations.")
             df<-aggregate(as.formula(paste( ".~ ",gene_col)), df, sum)
             print("Done modifying gene ids.")
         }else{
             stop("Ensembl IDs detected, but not a known format, such as ENSG00000288601.1|5S_rRNA or ENSG00000288601.1_5S_rRNA")
         }
    }else{
        print("Using orginal Gene IDs")
    }
#End inpute validataion
###################################################
    # not assume group based filtering; keep only the gens that have at least a certain number of 
     
    if (OPT_Use_Group_Based_Filtering) {
        rownames(df) <- df[,gene_col]
        df[,gene_col] <- NULL
        if(opt_use_cpm_counts_to_filter){
            counts <- edgeR::cpm(df) >  OPT_Minimum_Count_Value_to_be_Considered_Nonzero # boolean matrix
        }else{
            counts <- as.matrix(df) > OPT_Minimum_Count_Value_to_be_Considered_Nonzero
        }
        tcounts <- as.data.frame(t(counts))
        colnum <- dim(counts)[1] # number of genes
     

        rownames(metadata) <- metadata[,sample_col]
        tcounts <- merge(metadata[group_col], tcounts, by="row.names")
        tcounts$Row.names <- NULL
        melted <- melt(tcounts, id.vars=group_col)
        tcounts.tot <- dplyr::summarise(dplyr::group_by_at(melted, c(group_col, "variable")), sum=sum(value))
        tcounts.tot %>% tidyr::spread(variable, sum) -> tcounts.group
        colSums(tcounts.group[(1:colnum+1)]>=OPT_Minimum_Number_of_Samples_with_Nonzero_Counts_in_a_Group) >= 1 -> tcounts.keep ## change default value 3 here to adjust the minimum number of samples per group
        df <- df[tcounts.keep, ]
        df %>% rownames_to_column(gene_col) -> df

   } else {
        if (opt_use_cpm_counts_to_filter){
df$isexpr1 <- rowSums(edgeR::cpm(as.matrix(df[, -1])) > OPT_Minimum_Count_Value_to_be_Considered_Nonzero) >= OP_Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total} else {
df$isexpr1 <- rowSums(as.matrix(df[, -1]) > OPT_Minimum_Count_Value_to_be_Considered_Nonzero) >= OP_Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total}

    df <- as.data.frame(df[df$isexpr1, ])
    df$isexpr1 <- NULL
    }

    colnames(df)[colnames(df)==gene_col] <- "Gene"

    print(paste0("Number of genes after filtering: ", nrow(df)))

######## end of remove low count section ###############

df -> edf.orig

###### QC evaluation: PCA plot + Density plot ##############
# evaluate and display PCA plot

    rownames(Sampinfo) <- Sampinfo[,sample_col]
    Sampinfo <- Sampinfo[match(colnames(edf.orig[,-1]), Sampinfo[,sample_col]), ]
    Sampinfo = Sampinfo[complete.cases(Sampinfo[, sample_col]),]
    print(paste0("Total number of genes included in PCA evaluation: ", nrow(edf.orig)))

    edf <- edf.orig[,match(Sampinfo[,sample_col],colnames(edf.orig))]
    idx=!rowMeans(edf)==0
    edf=edgeR::cpm(edf[idx,])
    edf.orig = edf.orig[idx,]
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select_at(c(opt_pcx,opt_pcy))
    pca.df$celltype <- Sampinfo[,group_col]
    pca.df$sample <- Sampinfo[,plot_labels_col]

    
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
    pc.x.lab <- paste0("PC",OPT_Principal_Component_on_X_axis," ", perc.var[OPT_Principal_Component_on_X_axis],"%")
    pc.y.lab <- paste0("PC",OPT_Principal_Component_on_Y_axis," ", perc.var[OPT_Principal_Component_on_Y_axis],"%")
    
    labelpos <- pca.df
    labelpos$mean_y <- pca.df[,opt_pcy]+OPT_Label_Offset_Y_axis_
    labelpos$mean_x <- pca.df[,opt_pcx]+OPT_Label_Offset_X_axis_
    
    # set the colors to be used in the plot 
    
    if (length(unique(Sampinfo[,group_col])) > length(plotcolors)) {
        ## Original color-picking code.
        k=length(unique(Sampinfo[,group_col]))-length(plotcolors)
       
        more_cols<- getourrandomcolors(k) 
        plotcolors <- c(plotcolors , more_cols)
    }
    names(plotcolors)<-unique(Sampinfo[,group_col])
    print("Starting PCA")
# plot PCA  
    if (OPT_Add_Labels){
    g1 <- ggplot(pca.df, aes_string(x=opt_pcx, y=opt_pcy)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position=opt_legend_position) +
      geom_point(aes(color=pca.df$celltype), size=OPT_Point_Size) +
      geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
          label=sample, color=celltype, vjust="inward", hjust="inward"), size=OPT_Label_Font_Size, show.legend=FALSE) +
      ggtitle("PCA Plot") +    
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)
    }

    else{
    g1 <- ggplot(pca.df, aes_string(x=opt_pcx, y=opt_pcy)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position=opt_legend_position) +
      geom_point(aes(color=pca.df$celltype), size=OPT_Point_Size) +
      ggtitle("PCA plot") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
    }
if (TRUE){
    print(g1)
}

###### eval and display Density plot
df -> df.filt
  df.filt <- df.filt[,match(Sampinfo[,sample_col],colnames(df.filt))]
  print(paste0("Total number of samples included in the density plot: ", ncol(df.filt)))
  rownames(df.filt)=df[,1]
  df.filt %>% rownames_to_column(gene_col) -> df.filt
  df.m <- melt(df.filt,id.vars=c(gene_col))
  
  if(opt_perform_log_transform){
  df.m %>% mutate(value = log2(value+0.5)) -> df.m 
  }
  n <- 40
 

# set xlimits
if(Opt_Set_Min_Max_for_X_axis){
    xmin = Opt_Minimum_for_X_axis
    xmax = Opt_Maximum_for_X_axis
}else{
        xmin = min(df.m$value)
        xmax = max(df.m$value)
}

 if(opt_color_by_group ){

 df.m %>% mutate(colgroup = Sampinfo[variable,opt_group_to_color]) -> df.m
 df.m = df.m[complete.cases(df.m[, "colgroup"]),]
 df.m$colgroup = gsub("\\s","_",df.m$colgroup)
 df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
 if(opt_group_to_color==group_col){
     print("Using same colors as PCA plot")
     altplotcolors<-plotcolors
 }else{
     print("Generating new colors")
     altplotcolors<- getourrandomcolors(length(unique(df.m$colgroup))) 
 }
# plot Density 
g2 = ggplot(df.m, aes(x=value, group=variable)) + 
    geom_density(aes(colour = colgroup, linetype = colgroup)) +
    xlab("Filtered Counts") + ylab("Density") +
    theme_bw() +
    theme(legend.position=opt_legend_osition,legend.text = element_text(size = opt_legend_font_size)) + #scale_x_log10() + 
    ggtitle("Frequency Histogram") +
    xlim(xmin,xmax) +
    scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
    scale_colour_manual(values=altplotcolors)
 }
 else{
     df.m$variable = Sampinfo[df.m$variable,plot_labels_col]
     n=length(unique(df.m$variable))
     #m=ceiling(n/4)
     #n=m*4
     #color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
     

   
    cols<- getourrandomcolors(n) 

   
  g2 = ggplot(df.m, aes(x=value, group=variable)) + 
    geom_density(aes(colour = variable, linetype = variable)) +
    xlab("Filtered Counts") + ylab("Density") +
    theme_bw() +
    theme(legend.position=opt_legend_osition,legend.text = element_text(size = opt_legend_font_size)) + #scale_x_log10() + 
    ggtitle("Frequency Histogram") +
    xlim(xmin,xmax) +
    scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
    scale_colour_manual(values=cols)#+
    guides(linetype = guide_legend(ncol = opt_num_legend_cols))
 }
 
if (TRUE){ 
  dev.off()

# format image 
imageType = "png"
imageWidth = 3000
imageHeight = 1500*OPT_Number_of_Image_Rows
dpi = 300
if (imageType == 'png') {
    png(
      filename="Filtered_Counts.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
} else{
        library(svglite)
        svglite::svglite(
        file="Filtered_Counts.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }
 
  gh<-make_heatmap(df.filt[,samples_to_include],metadata,plotcolors,group_col,sample_col)
  grid.arrange(g1,g2,gh, nrow=OPT_Number_of_Image_Rows)
dev.off()

}
 #add heatmap
  
  
  
  
  

  
return(df)
}

#################################################
## Global imports and functions included below ##
#################################################
make_heatmap <- function(df, metadata,plotcolors,group_col,sample_col) {
  
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
  rownames(targetfile) = targetfile[,sample_col]
  idx = as.factor(targetfile[rownames(m),group_col])
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
# Functions defined here will be available to call in
# the code for any table.

#install_bioconductor_package <- function(pkg) {
#}

print("template_function_Filtered_Counts.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr804_rawcounts_for_NIDAP_csv<-readRDS(paste0(rds_output,"/var_ccbr804_rawcounts_for_NIDAP_csv.rds"))
var_ccbr804_rawcounts_for_NIDAP_csv<-as.data.frame(var_ccbr804_rawcounts_for_NIDAP_csv)
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr804_metadata_for_NIDAP<-readRDS(paste0(rds_output,"/var_ccbr804_metadata_for_NIDAP.rds"))
var_ccbr804_metadata_for_NIDAP<-as.data.frame(var_ccbr804_metadata_for_NIDAP)
invisible(graphics.off())
var_Filtered_Counts<-Filtered_Counts(var_ccbr804_rawcounts_for_NIDAP_csv,var_ccbr804_metadata_for_NIDAP)
invisible(graphics.off())
saveRDS(var_Filtered_Counts, paste0(rds_output,"/var_Filtered_Counts.rds"))
