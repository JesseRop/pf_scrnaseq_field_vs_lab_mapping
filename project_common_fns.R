# ## !!! NOTE - Function to run SCT within notebook instead of farm

HVG <- 500
UMAP_NM <- "umap"
RESOLN <- 0.5
CL_NM <- "seurat_clusters"
UMAP_REDCTN_NM <- "umap"
REDCTN <- "pca"
INTEGRTD_PC30 <- "no"
FINDPC_MTHD <- "perpendicular line"

sct_run_fn <- function(Seurat_object, hvg = HVG){
  
  ## Set assay for seurat object
  DefaultAssay(Seurat_object) <- "RNA"
  
  ## Runs SCT normalisation
  Seurat_object <- SCTransform(Seurat_object, vst.flavor = "v2", verbose = FALSE, variable.features.n = hvg, return.only.var.genes = FALSE) ## Added return.only.var.genes = FALSE on 30 March 2025 to ensure all genes are returned in data
  Seurat_object <- RunPCA(Seurat_object, verbose = F)
  
  sdev <- prcomp(t(Seurat_object@reductions$pca@cell.embeddings),scale. = T)$sdev[1:30]
  sdev <- sort(sdev, decreasing = TRUE)
  # npc <- findPC(sdev = sdev, number = 30, figure = T,method = 'all',aggregate = 'voting')
  npc <- findPC(sdev = sdev, number = 30, figure = T, method = FINDPC_MTHD)
  
  
  print(last_plot())
  
  npc <- ifelse(as.numeric(npc) < 5, 5, as.numeric(npc))
  npc <- ifelse(INTEGRTD_PC30 == "yes", 30, npc)
  
  print(npc)
  
  ## Run UMAP 
  Seurat_object <- RunUMAP(Seurat_object, reduction = REDCTN, dims = 1:npc, n.components = 3L, verbose=F, seed.use = 949, reduction.key = UMAP_NM, reduction.name = UMAP_REDCTN_NM)
  Seurat_object <- FindNeighbors(Seurat_object, reduction = REDCTN, dims = 1:npc, verbose = FALSE)
  Seurat_object <- FindClusters(Seurat_object, resolution = RESOLN, verbose = FALSE, cluster.name = CL_NM)
  
  DefaultAssay(Seurat_object) <- "RNA"
  
  return(Seurat_object)
  
}


## Common plotting functions 

stat_box_data <- function(y) {
  return( 
    data.frame(
      y=quantile(y,probs=1)*1.4,  #may need to modify this depending on your data
      label = paste(
        'n:', length(y), '\n',
        'mu:', round(mean(y), 1)
      )
    )
  )
}

umap_theme <- theme(legend.position = "right",
                    axis.title=element_text(size=12,face="bold",hjust = 0),
                    # plot.title = element_blank(),
                    legend.title=element_text(size=14,face="bold"),
                    legend.text=element_text(size=14),
                    legend.spacing.y = unit(0.25, "cm"),
                    axis.line = element_line(arrow = arrow(angle = 20, length = unit(0.15, "inches"), ends = "last", type = "closed"))
)

trunc_axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2.5, "cm")
)




methds_ln_plt <- function(tbl, x_var = "donor", y_var, col_var = "mthds_union") {
  ggplot(data = tbl, aes(x= !!rlang::sym(x_var), y = !!rlang::sym(y_var), colour = !!rlang::sym(col_var)))  +
    geom_point(size = 3) +
    geom_line(aes(group = !!rlang::sym(col_var)), size = 1) +
    geom_hline(yintercept = .50) +
    geom_hline(yintercept = .70) +
    theme_classic() + 
    theme(text = element_text(size =20),
          axis.text.x = element_text(size =20, angle = 45, hjust = 1),
          legend.position = "bottom")
}

## Function to make line plots to compare the different cell calling methods and the overlap between them

methds_ln_plt_nohline <- function(tbl, x_var = "donor", y_var, col_var = "mthds_union") {
  ggplot(data = tbl, aes(x= !!rlang::sym(x_var), y = !!rlang::sym(y_var), colour = !!rlang::sym(col_var)))  +
    geom_point(size = 4) +
    geom_line(aes(group = !!rlang::sym(col_var)), size = 1) +
    theme_classic() + 
    theme(text = element_text(size =20),
          axis.text.x = element_text(size =20, angle = 45, hjust = 1),
          legend.position = "bottom")
}

# Knee plot function for independent plots not facetting to check for the called cells
strn_knee_p_10_fn <- function(dset = seu_mdata, donor_nm,  ft_count = "nCount_RNA", ft_name = "Transcript", cell_method = "is_cell_cr") {
  dset %>% 
    mutate(rank = rank(desc(!!rlang::sym(ft_count)), ties.method = "average")) %>% 
    filter(!duplicated(rank)) %>%
    ggplot(., aes(y = !!rlang::sym(ft_count), x= rank, colour = !!rlang::sym(cell_method))) +
    ggrastr::rasterise(geom_point(size=0.5), dpi=100) + ## Reducing the value increases the speed but makes the plot more blurry
    geom_hline(yintercept = 100) +
    geom_hline(yintercept = 50) +
    theme_classic() +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    # scale_color_manual(name= 'Strain', values = strain_cols_n) +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"), 
          legend.text = element_text(size=20), 
          legend.title = element_text(size=20,face="bold"),
          legend.position = "none"
    )+ 
    annotation_logticks() +
    labs(x= 'Cell rank', y= paste0(ft_name, ' count'),
         title = donor_nm
    )
  
}

# Knee plot function faceting by a donor category and also by another level with each donor to check for the called cells

kneeplt_facet_fn <- function(ncount_tbl, cell_call_methd, donor_group, numbr_of_cols = 4, mthd_col = cell_call_mthd_col, rank_var = 'rank', cell_var = 'nCount_RNA', cut_offs = 50, thresh = "man_thresh", thresh_ln_size = 0.5, point_size = 0.5){ 
  
  ncount_tbl %>%
    filter(!is.na(!!rlang::sym(cell_call_methd))) %>%
    ggplot(., aes(y = !!rlang::sym(cell_var), x= !!rlang::sym(rank_var), colour = !!rlang::sym(cell_call_methd))) +
    # ggrastr::rasterise(geom_point(size=0.5), dpi=100) + ## Reducing the value increases the speed but makes the plot more blurry
    geom_point(size=point_size) + ## Reducing the value increases the speed but makes the plot more blurry
    # geom_hline(aes(yintercept = !!rlang::sym(thresh)), colour = "red") +
    geom_hline(aes(yintercept = !!rlang::sym(thresh)), size = thresh_ln_size, colour = "red") +
    # geom_hline(yintercept = 100) +
    geom_hline(yintercept = cut_offs) +
    theme_classic() +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                  labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                  labels = trans_format("log10", math_format(10^.x)))+
    # facet_nested_wrap(vars(!!rlang::sym(donor_group), donor), ncol = 4) +
    facet_nested_wrap(vars(!!!rlang::syms(unique(c(donor_group, "donor")))), ncol = numbr_of_cols) +
    scale_color_manual(name= 'CellCall', values = mthd_col) +
    theme(axis.text=element_text(size=15), 
          axis.title=element_text(size=20,face="bold"), 
          legend.text = element_text(size=20), 
          legend.title = element_text(size=15,face="bold"),
          legend.position = "bottom",
          strip.text.x = element_text(size = 10, face = "bold")
    )+ 
    # annotation_logticks() +
    labs(x= 'Cell rank', y= paste0("Transcript", ' count'),
         # title = donor_nm
    )
}

# function for making scatter plots (ncount VS nfeature) faceting by a donor category and also by each donor to check for the behavior of called cells

scatter_facet_fn <- function(ncount_tbl, cell_call_methd, donor_group, numbr_of_cols = 4, raster_resln = 50, use_raster = T, mthd_col = cell_call_mthd_col){ 
  
  ncount_tbl %>%
    ggplot(aes(x = nFeature_RNA, y= nCount_RNA, color = !!rlang::sym(cell_call_methd))) +
    (if (use_raster) ggrastr::rasterise(geom_point(size=0.1), dpi=raster_resln) else geom_point(size=0.1)) +
    geom_hline(yintercept = c(50, 100)) + 
    geom_vline(xintercept = c(50, 100)) +
    theme_classic() +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
    facet_nested_wrap(vars(!!rlang::sym(donor_group), donor), ncol = 4) +
    scale_color_manual(name= 'CellCall', values = mthd_col) +
    theme(axis.text=element_text(size=15), 
          axis.title=element_text(size=20,face="bold"), 
          legend.text = element_text(size=20), 
          legend.title = element_text(size=15,face="bold"),
          legend.position = "bottom",
          strip.text.x = element_text(size = 10, face = "bold")
    )+ 
    # annotation_logticks() +
    labs(x= 'Gene count', y= paste0("Transcript", ' count'),
         # title = donor_nm
    )+
    labs(title = paste0(str_remove(cell_call_methd, "is_cell_"), "_", donor_group))
}


# function for making boxplots/violinplots of either genes or transcripts faceting by a donor category and also by each donor to check for the behavior of called cells

umi_box_vln_plt_fn <- function(ncount_tbl, cell_call_methd, numbr_of_cols = 4, box_or_vln, donor_group = "frac_load", gn_umi = "nCount_RNA", gn_umi_nm = "Transcript count", scales_free = "free_y", mthd_col = cell_call_mthd_col){ 
  
  ncount_tbl %>%
    ggplot(aes( y= !!rlang::sym(gn_umi), x = !!rlang::sym(cell_call_methd), color = !!rlang::sym(cell_call_methd))) + 
    (if (box_or_vln == "boxplot") {geom_boxplot()} else if (box_or_vln == "violin") {geom_violin()}) +
    geom_hline(yintercept = 20000) + 
    geom_hline(yintercept = 100) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+ 
    theme_classic() +
    # facet_nested_wrap(vars(frac_load, donor), ncol = numbr_of_cols, scales = "free_y") +
    facet_nested_wrap(vars(!!!rlang::syms(unique(c(donor_group, "donor")))), ncol = numbr_of_cols, scales = scales_free) +
    scale_color_manual(name= 'CellCall', values = mthd_col) +
    theme(axis.text=element_text(size=15), 
          axis.text.x=element_text(angle=45, hjust = 1),
          axis.title=element_text(size=20,face="bold"), 
          legend.text = element_text(size=20), 
          legend.title = element_text(size=15,face="bold"),
          legend.position = "bottom",
          strip.text.x = element_text(size = 10, face = "bold")
    )+ 
    labs(x= 'Gene count', y= gn_umi_nm,
         title = paste0(str_remove(cell_call_methd, "is_cell_"), "_", donor_group)
    )
}


snk_pt_fn <- function(decod_tbl, strns_2_plt = c('Strain_2nd', 'n_strains'), labls = paste0("K=", opt_clusters_nms[[3]][1:2]), bp_face = c("bold", "plain"), ttl = "plot", strain_colours = strain_cols_n, labl_size = 3.6, bs_size = 9, strns_2_plt_lvls = sort(names(strain_cols_n)), labl_sep = ": ", rm_from_lab = "nothing") {
  decod_tbl %>% 
    mutate(across(contains("train"), ~str_replace_all(., c("Doublet" = "Dbt", "Negative" = "Neg")))) %>%
    make_long(all_of(strns_2_plt)) %>%
    # filter(node != "NA", !is.na(node)) %>% 
    mutate(node = droplevels(factor(node, strns_2_plt_lvls))) %>% 
    # arrange(node) %>% 
    add_count(x,node, "n") %>%
    ggplot(., aes(x = x, next_x = next_x,node = node,next_node = next_node,fill = node,label = paste0(str_remove(node, rm_from_lab), labl_sep,n))) +
    # ggplot(., aes(x = x, next_x = next_x,node = node,next_node = next_node,fill = factor(node),label = paste0(node, labl_sep,n))) +
    geom_sankey() +
    theme_void()+
    scale_fill_manual(values = strain_colours) +
    scale_x_discrete(labels=labls)+
    geom_sankey_label(size = labl_size, color = 1, fill = "white") +
    labs(title = ttl) +
    theme_sankey(base_size = bs_size) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 15, face = bp_face))
  
}



##Function to save tables in excel as sheets

createWorksheet <- function(x, y){
  wb <- createWorkbook()
  for(i in seq_along(x)){
    addWorksheet(wb, names(x)[i])
    writeData(wb, sheet = names(x)[i], x[[i]])
    saveWorkbook(wb, file = y, overwrite = TRUE)
  }
}

## Function for plotting the average expression of a group of cell

# Function 

## Function to calculate average expression score from Sunil
makeSetScore <- function(seur, gene.set,assay='RNA', slt = 'data') {
  # Get mean expression of genes of interest per cell
  gene.set <- gene.set[gene.set %in% rownames(seur) ]
  mean.exp <- colMeans(x = as.matrix(GetAssayData(seur, slot = slt, assay = assay))[gene.set, ], na.rm = TRUE)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = seur@meta.data))) {
    return(mean.exp)
  }
}


## Function to plot UMAPs and color clusters while labelling inside the plot
umap_col_plot_fn <- function(dset = asx_gns_mdata_umap, col_nm = "Stage", clstr_labl = "leiden", group_by_var = "leiden", dim1 = "umap_1", dim2 = "umap_2", pt_size = 0.01, clustr_cols = leiden_cols, is_continuous = "no", cluster_text_lab_size = 9) {
  
  cluster_centers <- dset %>%
    group_by(!!rlang::sym(clstr_labl)) %>%
    dplyr::summarize(umap_cntr_1 = median(!!rlang::sym(dim1)), umap_cntr_2 = median(!!rlang::sym(dim2)))
  
  dset %>%
    # group_by(!!!rlang::syms(grp_var)) %>%
    # arrange(!!rlang::sym(gene)) %>%
    ggplot(aes(x = !!rlang::sym(dim1), y = !!rlang::sym(dim2), colour = !!rlang::sym(group_by_var))) +
    geom_point(size = pt_size) +
    geom_text(data = cluster_centers,      # Add labels at centers
              aes(x = umap_cntr_1, y = umap_cntr_2, label = !!rlang::sym(clstr_labl)), 
              color = "black", 
              fontface = "bold",
              size = cluster_text_lab_size) +
    {if(is_continuous == "yes") {
      scale_color_distiller(direction = 1) 
    } else {
      scale_color_manual(values = clustr_cols)
    }} +
    # facet_nested_wrap(vars(!!!rlang::syms(grp_var)), nrow = row_num) +
    labs(x = str_replace(dim1, "umapharmonyPCman_", "UMAP "), y = str_replace(dim2, "umapharmonyPCman_", "UMAP ")) +
    theme_classic() +
    theme(text = element_text(size = 18),
          legend.position = "bottom",
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    guides(colour = guide_colorbar(title = col_nm, direction = "horizontal"))
}


# Theme to remove umalabels

umap_axis_rm_theme <- theme(axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.title = element_blank())


## ggally plotting functions
#### GGgally for pairwise comparisons across all clusters

lowerfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 0.7, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "blue")
}

# Add this custom function to your environment
cor_heatmap <- function(data, mapping, ...) {
  GGally::ggally_cor(data, mapping, ...) + 
    theme(panel.background = element_rect(fill = "white")) # Ensures background doesn't override colors
}

# Define a helper to ensure panel borders are visible
my_ggpairs_theme <- theme_bw() + # Starting with theme_bw provides the panel boxes
  theme(
    axis.text.x = element_text(angle = 90, size = 16),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    text = element_text(size = 16),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), # Grid surrounding facets
    legend.position = "none"
  )


### !!!!!!NOTE - BELOW FUNCTIONS ARE FOR FOR CELL CALLING PLOTS
## COPIED FROM ~/mali22_code_lnk/2_cell_calling/plotting_fns_n_variables.R - to have one script with all the key functions
## VARIABLES

rank_colours <- c("CellEd" = "red", "Non-Cell" = "black", "Below_Limit" = "grey60")


# FUNCTIONS
## Barcode ranks plotting function incorporating thresholds for the UMI cut off for the knee and inflection for the number of likely cells and total droplets calculated by droplet utils - Added as geom hlines

bcranks_cutoff_plot_fn <-  function(dset, hline_droplets = "umi_counts_cutoff_total_droplets_", hline_cells = "umi_counts_cutoff_cell_estimate_", hline_cells_lab = "umi_cutoff_cells:", hline_droplets_lab = "umi_cutoff_droplets:", line_part = "knee", ncols = 6) {
  label_data = distinct(dset, donor, species, jumpcode, !!rlang::sym(paste0(hline_cells, line_part)), !!rlang::sym(paste0(hline_droplets, line_part)))
  
  dset %>%
    ggplot(., aes(x = rank, y = total)) +
    geom_point(size = 0.3, alpha = 0.6) +
    # scale_colour_manual(values = rank_colours) +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(aes(yintercept = !!rlang::sym(paste0(hline_cells, line_part)))) +
    geom_hline(aes(yintercept = !!rlang::sym(paste0(hline_droplets, line_part)))) +
    geom_text(data = label_data,
              aes(x = 500, 
                  y = !!rlang::sym(paste0(hline_cells, line_part)), 
                  label = paste(hline_cells_lab, !!rlang::sym(paste0(hline_cells, line_part)))),
              vjust = -0.2) +
    geom_text(data = label_data,
              aes(x = 500, 
                  y = !!rlang::sym(paste0(hline_droplets, line_part)), 
                  label = paste(hline_droplets_lab, !!rlang::sym(paste0(hline_droplets, line_part)))),
              vjust = -0.2) +
    facet_wrap2(vars(donor, species, jumpcode), ncol = ncols) +
    labs(title = paste0(line_part)) +
    theme_classic()
  
  
}



## Barcode ranks plotting function incorporating thresholds for the number of likely cells and total droplets calculated by droplet utils - Added as geom vlines

# Generated in 'm22_EDrops_ambient_cutoff_est.Rmd'
bcranks_ncells_ndrops_plot_fn <-  function(dset, don_nm, vline_cells_lab = "cell_estimate: ", vline_droplets_lab = "total_droplets: ",  lab_size = 5, ncols = 6, c_pos=2, d_pos = 1500) {
  # allow datasets with or without `jumpcode` column
  has_jump <- "jumpcode" %in% colnames(dset)
  if (has_jump) {
    label_data <- distinct(dset, donor, species, jumpcode, n_cells_cell_estimate_estimated_ncells, n_cells_total_droplets_estimated_ndroplets, umi_counts_cutoff_cell_estimate_knee, umi_counts_cutoff_total_droplets_knee, umi_counts_cutoff_cell_estimate_inflection)
    
  } else {
    label_data <- distinct(dset, donor, species, n_cells_cell_estimate_estimated_ncells, n_cells_total_droplets_estimated_ndroplets, umi_counts_cutoff_cell_estimate_knee, umi_counts_cutoff_total_droplets_knee, umi_counts_cutoff_cell_estimate_inflection)
    
  }
  
  p <- dset %>%
    data.frame() %>%
    # filter(!(duplicated(rank))) %>%
    ggplot(., aes(x = rank, y = total)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                  labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_vline(aes(xintercept = n_cells_cell_estimate_estimated_ncells)) +
    geom_vline(aes(xintercept = n_cells_total_droplets_estimated_ndroplets), linetype="dashed") +
    # geom_hline(aes(yintercept = umi_counts_cutoff_cell_estimate_inflection)) +
    # labs(title = don_nm) +
    geom_text(data = label_data,
              aes(y = c_pos,
                  x = n_cells_cell_estimate_estimated_ncells,
                  label = paste0(vline_cells_lab, n_cells_cell_estimate_estimated_ncells)),
              vjust = -0.2, hjust =0.2, angle = 90, size = lab_size) +
    geom_text(data = label_data,
              aes(y = d_pos,
                  x = n_cells_total_droplets_estimated_ndroplets,
                  label = paste0(vline_droplets_lab, n_cells_total_droplets_estimated_ndroplets)),
              vjust = 1.2, hjust =0.2, angle = 90, size = lab_size) +
    theme_classic()
  
  if (has_jump) {
    # p <- p + facet_wrap2(vars(donor, species, jumpcode), ncol = ncols, scales = "free_x")
    p <- p + facet_wrap2(vars(donor, species, jumpcode), ncol = ncols)
  } else {
    p <- p + facet_wrap2(vars(donor, species), ncol = ncols)
  }
  
  p
  
}



## Function to plot how different cell calling methods including different iterations of cellbender overlap
call_methods_overlap_barplot_fn <- function(tbl = cr_cbender_mdata_slct_uni, mthd = "mthds_union", pf_umi_filt = 0, min_grp_size =0) {
  map(tbl, ~.x %>% 
        filter(pf_umi_count >= pf_umi_filt) %>%
        count(!!rlang::sym(mthd)) %>%
        filter(n > min_grp_size)
  ) %>%
    bind_rows(., .id = "donor") %>%
    # filter(n > 100) %>%
    ggplot(aes(y = n/1000, x = !!rlang::sym(mthd), fill = !!rlang::sym(mthd))) +
    geom_col() +
    # geom_text(aes(label=n), hjust=0.5, vjust= -0.3, position = position_dodge(width = .9), size = 5, show.legend = F)+
    geom_text(aes(label=n), angle=90, hjust=-0.1, vjust= 0.5, position = position_dodge(width = .9), size = 3.8, show.legend = F)+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.8))) +
    scale_fill_manual(values = c_combos_cols) +
    facet_wrap(. ~ donor, ncol = 5, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 11),
          legend.position = "bottom")
  
}


## Function to plot the scatter of the human versus plasmodium umi or gene count coloured by different cell calling method
call_methods_overlap_scatter_fn <- function(tbl = cr_cbender_uni_pf100_mrg, gn_or_umi = "umi",  mthd = "mthds_union") {
  
  ggplot(tbl, aes(x = log10(!!rlang::sym(paste0('pf_',gn_or_umi,'_count'))), y = log10(!!rlang::sym(paste0('hs_',gn_or_umi,'_count'))), colour = !!rlang::sym(mthd))) +
    geom_point(size = 0.01) +
    # geom_density_2d()  +
    scale_colour_manual(values = c_combos_cols) +
    facet_wrap(donor ~ .,  ncol = 5) +
    geom_hline(yintercept = log10(100)) +
    geom_vline(xintercept = log10(100)) +
    guides(color = guide_legend(title = 'Cell call', override.aes = list(size = 4), ncol = 1))
}




