set -e
Rscript template_function_Filtered_Counts.R
Rscript template_function_Normalized_Counts.R
Rscript template_function_DEG_Analysis.R
Rscript template_function_Expression_Heatmap.R
Rscript template_function_GSEA_Preranked_KOvWT_Activated.R
Rscript template_function_Volcano_Summary.R
Rscript template_function_GSEA_Filters_KOVWT_Activated.R
Rscript template_function_GSEA_Visualization_KOvWT_Activated.R
source("workbook_start_globals.R")
