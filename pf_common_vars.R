library(MetBrewer)
library(RColorBrewer)

## Common tables and variables 

##Strain cluster colours
# strain_cols_n <-c("#910000", "#ff1a8d", "#007c71", "#0ea300", "#00deca", "#9a96c7", "#5d57a4", "#a8875b","#708db3", "#c39a00", "#ffd94d", "#9770b3", "#C9E0F1", "#a40000" , "#00b7a7", "#ffcd12","lightgrey","#a40000" , "#00b7a7","#00b7a7", "#e7f1f9", "#D2C1DC", 'black', '#050301','grey','lightgrey','#ededed')
#, 'lavender', 'lightgray', 'grey'

strain_cols_n <-c("#910000", "#ff1a8d", "#007c71", "#0ea300", "#00deca", "#9a96c7", "#5d57a4", "#a8875b","#708db3", "#c39a00", "#ffd94d", "#9770b3", "#C9E0F1", "#a40000" , "#00b7a7", "#ffcd12","lightgrey","#a40000" , "darkgrey","darkgrey", "#e7f1f9", "#D2C1DC", 'black', '#050301','grey','lightgrey','#ededed',"#afe1ed","lightgrey","lightgrey")


strain_nms <- c('SC2', 'SC1','SC5',  'SC4', 'SC3',  'SC6', 'SC7', 'SC8', 'SC9', 'SC10', 'SC11', 'SC12', 'S3_b', 'S1', 'S2', 'Singlet', 'Negative', 'SC1_3i','Doublet','Dbt', 'Lab', 'Field', 'Less5', 'PoorQC', NA, 'Neg', 'Other', 'SC0', "Missing", "Mis")

names(strain_cols_n) <- strain_nms
strain_nms_lvls <- str_sort(strain_nms, numeric = T)
strain_nms_lvls <- c(strain_nms_lvls[str_detect(strain_nms_lvls, "SC")], strain_nms_lvls[!str_detect(strain_nms_lvls, "SC")])

## Add stage to strains and get colours
strain_stg_nms <- map(strain_nms, ~paste0(.x, c("_A", "_G"))) %>% unlist()

strain_stg_cols_n <- c("#0110b5", "#69e095", "#d65559", "#3de4f7", "#9751e2", "#0b66c1", 
  "#12eadc", "#b220cc", "#d34593", "#a5431c", "#e2b3f2", "#2538b2", 
  "#ef3b0e", "#17a597", "#88bc18", "#3ed80f", "#bc66d6", "#5149b2", 
  "#e2b34d", "#f9dbb3", "#f7eb65", "#0e8759", "#dd21dd", "#91f7cb", 
  "#d35d6e", "#f23242", "#2f7bd8", "#e8c97a", "#a5ff89", "#a8ffc8", 
  "#83fca1", "#b0f449", "#ccfc9c", "#b75b2d", "#ce7d44", "#7ff232", 
  "#ffcdba", "#95edce", "#2727aa", "#3e529b", "#f7f972", "#ff59a6", 
  "#1d0fe2", "#95fc7e", "#d32e3c", "#ed8e78", "#acc0ad", "#9b4788", "#afd0ad", "#9b4677", "#adf1ad", "#9f6687", "#adf1ed", "#9a6687", "#Bdf8ed", "#9a5087", "#0120b5", "#60e005", "#d66569", "#1af2e7") %>% set_names(strain_stg_nms)

dblt_cols <- c(brewer.pal(name = 'Set1', n = 8)[c(1:5,7:8)], 'lightblue')
names(dblt_cols) <- c('ScDfSt', 'ScDf', 'ScSt', 'DfSt', 'Sc', 'Df','St', 'Singlet')

cell_qc_filt_nm_lvls <- c("sc_neg", "lo_gn_cnt",  "hi_gn_cnt",  "hi_mtp", "ss_db",  "strn_db", "stg_db", "pass", "cell_prob.9", "cb_bg.5", "stg_unassignd")
cell_qc_cols <- c(brewer.pal(name = 'Set1', n = 8)[c(1:5,7:8)], 'lightblue', "#009E73", "#eebe99", "#00a7d9") %>% set_names(cell_qc_filt_nm_lvls)

dblt_ss_cols <- c(brewer.pal(name = 'Set1', n = 8)[c(1:3)], 'lightgrey')
names(dblt_ss_cols) <- c('Strain & stage', 'Strain', 'Stage', 'Singlet')

cluster_cols <- brewer.pal(8, name = 'Set2')[c(1:4,7,8)] %>% set_names(as.character(c(0:5)))

stg_refnd_lvls = c("Early ring", "early ring", "Late ring", "late ring", "Early trophozoite", "early trophozoite", "MidTrophozoite", "Mid trophozoite", "Late trophozoite", "late trophozoite", "Early schizont", "early schizont",  "Late schizont","late schizont", "Schizont",  "Asexual", "Gametocyte (committed)", "Gametocyte (developing)", "gametocyte (developing)", "Gametocyte (branching)", "Female (HE)", "Female", "Female (LE)", "Female LE", "Gametocyte (female)", "gametocyte (female)", "Gametocyte (early female)", "Gametocyte (late female)", "Gametocyte (male)", "Male (HE)", "Male", "Male (LE)", "Male LE", "Gametocyte (early male)", "Gametocyte (late male)", "gametocyte (male)", "Sexual", "Unassigned", "Failed qc", "Not assigned", "Atlas", "stgQC_fail", "Other", "Clean")

# stg_lvls = c("Early ring", "Late ring", "Early trophozoite","Late trophozoite" , "Early schizont", "Late schizont",  "Gametocyte (developing)", "Female (HE)", "Female", "Female (LE)",  "Gametocyte (female)", "Male (HE)", "Male", "Male (LE)", "Gametocyte (male)", "Asexual", "Unassigned", "Failed QC", "Not assigned")


# stg_refnd_lvls2 = c('Early ring','Late ring','Early trophozoite','Late trophozoite','Early schizont','Late schizont', "Gametocyte (committed)","Gametocyte (developing)","Gametocyte (branching)", "Female (HE)", "Female", "Female (LE)", "Gametocyte (early female)","Gametocyte (late female)", "Male (HE)", "Male", "Male (LE)", "Gametocyte (early male)","Gametocyte (late male)", "Asexual", "Unassigned", "Failed QC", "Not assigned") 

sp_colors <- c('Early ring'='#D1EC9F',
               'Late ring'='#78C679',
               'Early trophozoite'='#FEEEAA',
               'Late trophozoite'='#FEB24C',
               'Early schizont'='#C9E8F1',
               'Late schizont'='#85B1D3',
               'Gametocyte (developing)'='thistle',
               'Gametocyte (female)'='#551A8B',
               'Gametocyte (male)'='mediumpurple',
               'Unassigned'='red', 
               'Not assigned'='#D6CFC7',
               'Female (LE)'='#208cf7',
               'Female (HE)'='#551A8B',
               'Male (LE)'='#a894d1',
               'Male (HE)'='mediumpurple',
               'Failed QC' = 'grey',
               'Gametocyte (committed)' = '#C399A2', 
               'Gametocyte (developing)' = '#66444C', 
               'Gametocyte (branching)' = '#66444C', 
               'Gametocyte (early female)' = '#749E89', 
               'Gametocyte (early male)' = '#7D87B2', 
               'Gametocyte (late female)' = '#4E6D58', 
               'Gametocyte (late male)' = '#41507B',
               "Other" = '#DCDCDC',
               "Clean" = '#DCDCDC'
)


sp_fld_colors <- c('Late ring'= '#78C679',#met.brewer('Java')[2],
                   'Early trophozoite'='#FEEEAA',#met.brewer('Java')[3],
                   'Female LE'= met.brewer('Ingres')[5],
                   'Female (LE)'= met.brewer('Ingres')[5],
                   'Female'= '#551A8B', #met.brewer('Ingres')[7],
                   # 'Female (HE)'= '#551A8B', #met.brewer('Ingres')[7],
                   # 'Male (LE)'=met.brewer('Nizami')[6],
                   # 'Male (HE)'=met.brewer('Nizami')[8],
                   # 'Gametocyte (male)'=met.brewer('Austria')[5],
                   'Male LE'=met.brewer('Ingres')[7],#'#a2e19e',
                   'Male (LE)'=met.brewer('Ingres')[7],#'#a2e19e',
                   'Male'= 'mediumpurple',#'#039911',
                   ## 'Male (HE)'= 'mediumpurple',#'#039911',
                   'Gametocyte (male)'='#1a6c8b',
                   'Gametocyte (female)'=met.brewer('Ingres')[6],
                   'Lab' = '#D3D3D3',
                   'Unassigned'='#DCDCDC',
                   'Failed QC' = 'grey',
                   'Not assigned' = 'grey',
                   'Atlas' = 'grey',
                   'Early ring'='#D1EC9F',
                   'Late trophozoite'='#FEB24C',
                   'Early schizont'='#C9E8F1',
                   'Late schizont'='#85B1D3',
                   'Gametocyte (developing)'='thistle',
                   'Developing \ngametocyte'='thistle',
                   'early female' = '#749E89',
                   'late female'= '#4E6D58',
                   "Asexual" = "#FA8072",
                   NULL = 'lightgrey',
                   'Field unlabelled'='#D2C1DC',
                   "Schizont" = '#85B1D3',
                   # New assignments:
                   "fem14" = "#2A9D8F",              # muted teal
                   "fem8_activated_like" = "#E76F51",# warm coral
                   "fem_non8_asx_dblt" = "#6C757D",  # cool grey
                   "fem10_low2" = "#BC6C25",         # earthy amber
                   "fem9_low2" = "#457B9D",          # muted steel blue
                   "fem8_asx_dblt" = "#8AB17D",      # soft sage green
                   "fem15_egress" = "#F4A261",        # soft orange
                   "ActivFemLE" = "#d900a7",   # rich green
                   "ActivFem"   = "#77e0ff",
                   "Female4" = "#1a508b",   
                   "Female2"   = "#1a898b",
                   "Other" = '#DCDCDC',
                   "Clean" = '#DCDCDC',
                   " " = 'lightgrey'
)

sp_fld_colors2 <- c('Late ring'= '#78C679',#met.brewer('Java')[2],
                   'Early trophozoite'='#FEEEAA',#met.brewer('Java')[3],
                   'Female LE'= met.brewer('Ingres')[5],
                   'Female (LE)'= met.brewer('Ingres')[5],
                   'Female'= '#551A8B', #met.brewer('Ingres')[7],
                   # 'Female (HE)'= '#551A8B', #met.brewer('Ingres')[7],
                   # 'Male (LE)'=met.brewer('Nizami')[6],
                   # 'Male (HE)'=met.brewer('Nizami')[8],
                   # 'Gametocyte (male)'=met.brewer('Austria')[5],
                   'Male LE'=met.brewer('Ingres')[7],#'#a2e19e',
                   'Male (LE)'=met.brewer('Ingres')[7],#'#a2e19e',
                   'Male'= 'mediumpurple',#'#039911',
                   ## 'Male (HE)'= 'mediumpurple',#'#039911',
                   'Gametocyte (male)'='#1a6c8b',
                   'Gametocyte (female)'="#f5688e",
                   'Lab' = '#D3D3D3',
                   'Unassigned'='#DCDCDC',
                   'Failed QC' = 'grey',
                   'Not assigned' = 'grey',
                   'Atlas' = 'grey',
                   'Early ring'='#D1EC9F',
                   'Late trophozoite'='#FEB24C',
                   'Early schizont'='#C9E8F1',
                   'Late schizont'='#85B1D3',
                   'Gametocyte (developing)'='thistle',
                   'Developing \ngametocyte'='thistle',
                   'early female' = '#749E89',
                   'late female'= '#4E6D58',
                   "Asexual" = "#FA8072",
                   NULL = 'lightgrey',
                   'Field unlabelled'='#D2C1DC',
                   "Schizont" = '#85B1D3',
                   "DevelopingGametocyte"='thistle',
                   "EarlyRing"='#D1EC9F',
                   "EarlyTrophozoite"='#FEEEAA',
                   "FemaleGametocyte"= '#551A8B',
                   "LateRing"= '#78C679',
                   "LateTrophozoite"='#fe9911',
                   "MidTrophozoite" = "#fec373",
                   "Mid trophozoite" = "#fec373",
                   # New assignments:
                   "fem14" = "#8AB17D",              # muted teal
                   "fem8_activated_like" = "#E76F51",# warm coral
                   "fem_non8_asx_dblt" = "#d4d7da",  # cool grey
                   "fem10_low2" = "#BC6C25",         # earthy amber
                   "fem9_low2" = "#457B9D",          # muted steel blue
                   "fem8_asx_dblt" = "#e9eaec",      # soft sage green
                   "fem_asx_dblt" = "#e9eaec",
                   "fem15_egress" = "#f7f086",        
                   "ActivFemLE" = "#d900a7",   
                   "ActivFem"   = "#77e0ff",
                   "Female4" = "#1a508b",   
                   "Female2"   = "#1a898b",
                   "Other" = '#DCDCDC',
                   "Clean" = '#DCDCDC',
                   " " = 'lightgrey', "weak_score" = "#c9c9c9", "unclear" = "#dddddd",`discrepantFM` =  "#e5d0ff", "Missing" = '#DCDCDC'
)

sp_fld_colors3 <- c('Late ring'= '#78C679',#met.brewer('Java')[2],
                    'Early trophozoite'='#FA8072',#met.brewer('Java')[3],
                    'Female LE'= met.brewer('Ingres')[5],
                    'Female (LE)'= met.brewer('Ingres')[5],
                    'Female'= '#551A8B', #met.brewer('Ingres')[7],
                    # 'Female (HE)'= '#551A8B', #met.brewer('Ingres')[7],
                    # 'Male (LE)'=met.brewer('Nizami')[6],
                    # 'Male (HE)'=met.brewer('Nizami')[8],
                    # 'Gametocyte (male)'=met.brewer('Austria')[5],
                    'Male LE'=met.brewer('Ingres')[7],#'#a2e19e',
                    'Male (LE)'=met.brewer('Ingres')[7],#'#a2e19e',
                    'Male'= 'mediumpurple',#'#039911',
                    ## 'Male (HE)'= 'mediumpurple',#'#039911',
                    'Gametocyte (male)'='#1a6c8b',
                    'Gametocyte (female)'="#4169E1",
                    'Lab' = '#D3D3D3',
                    'Unassigned'='#DCDCDC',
                    'Failed QC' = 'grey',
                    'Not assigned' = 'grey',
                    'Atlas' = 'grey',
                    'Early ring'='#D1EC9F',
                    'Late trophozoite'='#DC143C',
                    'Early schizont'='#C9E8F1',
                    'Late schizont'='#85B1D3',
                    'Gametocyte (developing)'='thistle',
                    'early female' = '#749E89',
                    'late female'= '#4E6D58',
                    "Asexual" = "#FBC81D",
                    NULL = 'lightgrey',
                    'Field unlabelled'='#D2C1DC',
                    "Schizont" = '#85B1D3',
                    # New assignments:
                    "fem14" = "#8AB17D",              # muted teal
                    "fem8_activated_like" = "#E76F51",# warm coral
                    "fem_non8_asx_dblt" = "#6C757D",  # cool grey
                    "fem10_low2" = "#BC6C25",         # earthy amber
                    "fem9_low2" = "#457B9D",          # muted steel blue
                    "fem8_asx_dblt" = "#2A9D8F",      # soft sage green
                    "fem15_egress" = "#F4A261",        # soft orange
                    "ActivFemLE" = "#d900a7",   # rich green
                    "ActivFem"   = "#77e0ff",
                    "Female4" = "#1a508b",   
                    "Female2"   = "#1a898b",
                    "Other" = '#DCDCDC',
                    "Clean" = '#DCDCDC' ,
                    " " = 'lightgrey'
)
# "#DC143C"
# "#B22222"
# "#F08080"
# "#FA8072"
# "#4169E1"

cluster_cols <- brewer.pal(7, name = 'Set2') %>% set_names(as.character(c(0:6)))


## donor colors 
# donor_cols_m = c(met.brewer("Thomas")[c(1)],c('#F39EC8','#99ceff', '#b9e589')) %>% set_names(paste0(c("MSC1", "MSC3", "MSC13", "MSC14" ), '.mrna'))
# # donor_cols = met.brewer("Egypt")[c(2,4,1,3)] %>% set_names(paste0(c("MSC1", "MSC3", "MSC13", "MSC14" )))
# donor_cols = c(met.brewer("Thomas")[c(1)],c('#F39EC8','#99ceff', '#b9e589')) %>% set_names(paste0(c("MSC1", "MSC3", "MSC13", "MSC14" )))

c_stg_cols= c('Asexual' = "#D7B7C7", 'Male' = "#C4CAA4", 'Female'="#6289A9", 'Developing' = "#C9F6E1", 'Unassigned'="lightgrey")


##Metadata from all p.falciparum genes downloaded from plasmoDB
# NOTE: Using read.csv() instead of read_csv() for speed
if (!exists("Pf3D7_genes_plasmoDB")) {
  Pf3D7_genes_plasmoDB <- read.csv("/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pf3d7_Genes_Summary_270922.csv",
                                     stringsAsFactors = FALSE, check.names = F) %>%
    as_tibble() %>%
    dplyr::rename("gene_id" = `Gene ID`,
                  'Gene_name' = `Gene Name or Symbol`, 
                  'Ka_Ks' = `NonSyn/Syn SNP Ratio All Strains`, 
                  'TM_domains' = `# TM Domains`) %>%
    mutate(across(where(is.character), ~na_if(., "N/A"))) %>%
    add_count(gene_id, name = 'gn_iso') %>%
    mutate(gene_id = if_else(gn_iso > 1, source_id, gene_id),
           gene_id = str_replace(gene_id, "_", "-")) %>%
    mutate(Gene = coalesce(Gene_name, gene_id)) %>%
    add_count(Gene, name = 'n') %>%
    mutate(Gene = if_else(n > 1, paste0(Gene, '_', gene_id), Gene))
}


## Irods decode file plus metadata
irods_id_sk21_mdata_cln <- read.csv("/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/irods_id_sk21_mdata_cln.csv",
                                     stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  mutate(sample_nm = str_remove_all(sample_name, "_JC")) %>%
  arrange(sample_nm)

smpl_nms_all <-irods_id_sk21_mdata_cln$sample_nm

## Read in msc annotations 
# donor_cols = met.brewer(name = "Renoir")[1:length(sample_name_nms)] %>% set_names(str_remove(sample_name_nms, "SC"))
sample_name_nms_all1 <-sort(c("MSC33", "MSC48","MSC49", "MSC55", "MSC60", "MSC66",  "MSC40", "MSC50_MACS", "MSC53", "MSC57", "MSC24", "MSC37", "MSC45", "MSC50_SLO", "MSC54", "MSC67", "MSC68", "MSC70", "MSC41", "MSC51", "MSC25", "MSC3", "MSC13", "MSC14", "MSC1", "Other"))
donor_cols = c("#C1EAE5", "#A7B19D", "#73935B", "#669CE2", "#62E1D0", "#DB6F87", "#DC48A2", "#E593D6", "#E65840", "#B29BDB", "#66EA54", "#B8E9B5", "#AFE685", "#D6ACCC", "#76586A", "#D8A46B", "#736AD9", "#62919E", "#65E297", "#E3AE3B", "#DA4CE6", "#BF6BD1", "#6FCEE8", "#8537E1", "#E3E544", "#D4D0E8", "#B0E653", "#E1A797", "#DFDC7F", "#E7E3C4" , "#DAEC3F",'#ededed')[1:length(sample_name_nms_all1)] %>% set_names(str_remove(sample_name_nms_all1, "SC"))

strain_don_cols = c("#db7dcd", "#71e8da", "#fcc2d0", "#4366a8", "#09efef", "#7cc435", 
  "#7736c1", "#286ea0", "#d4f77b", "#49ddb8", "#7d88ed", "#f2ef63", 
  "#9be076", "#373ce5", "#d8a6ea", "#9cd64a", "#83e04a", "#166170", 
  "#f7afc0", "#8058c9", "#f9728b", "#3ec0d1", "#f9917c", "#ead87e", 
  "#ffcccd", "#daef99", "#fcd916", "#c9e88d", "#646bef", "#e018c9", 
  "#6af455", "#e031cb", "#016ba8", "#56d3c9", "#562abf", "#55f4d4", 
  "#d15780", "#244b93", "#8aea98", "#70f9e7", "#5294ce", "#fccdbd", 
  "#c68c1f", "#f41916", "#008080", "#dcbeff", "#f7b2c6", "#cdc4fc", 
  "#c68af7", "#c10583", "#82001a", "#1d3daf", "#b2ffbb", "#bc5c32", 
  "#d1b5ff", "#3ae859", "#57d690", "#267693", "#9baddd", "#42e2f7", 
  "#eae13a", "#53b0b5", "#ed253d", "#ef4a5a", "#1a0d91", "#edcb97", 
  "#354ae8", "#162f72", "#ddb66c", "#ba6a1f", "#0c3675", "#2fce94", 
  "#d8d634", "#5de2b8", "#b77217", "#f9e952", "#ffb560", "#e8ce81", 
  "#bd49fc", "#d13c95", "#75dd85", "#610b89", "#e9a295", "#767d85", "#a20b89") %>% 
  set_names(c("M1_SC2", "M1_SC3", "M1_SC4", "M1_SC1", "M1_SC5", "M3_SC2", 
                "M3_SC1", "M13_SC2", "M13_SC1", "M14_SC3", "M14_SC1", "M14_SC2", 
                "M14_SC5", "M14_SC7", "M14_SC8", "M14_SC6", "M14_SC4", "M24_SC0", 
                "M25_SC1", "M25_SC2", "M33_SC3", "M33_SC1", "M33_SC7", "M33_SC2", 
                "M33_SC6", "M33_SC5", "M33_SC4", "M33_SC8", "M40_SC6", "M40_SC4", 
                "M40_SC1", "M40_SC5", "M40_SC2", "M40_SC3", "M40_SC8", "M40_SC7", 
                "M41_SC1", "M41_SC2", "M41_SC3", "M45_SC0", "M48_SC3", "M48_SC2", 
                "M48_SC1", "M48_SC4", "M49_SC5", "M49_SC1", "M49_SC2", "M49_SC4", 
                "M49_SC6", "M49_SC3", "M50_MACS_SC0", "M50_SLO_SC2", "M50_SLO_SC1", 
                "M51_SC0", "M53_SC0", "M54_SC0", "M55_SC4", "M55_SC3", "M55_SC2", 
                "M55_SC6", "M55_SC1", "M57_SC0", "M60_SC4", "M60_SC2", "M60_SC3", 
                "M60_SC6", "M60_SC1", "M60_SC7", "M60_SC5", "M66_SC5", "M66_SC3", 
                "M66_SC2", "M66_SC8", "M66_SC6", "M66_SC4", "M66_SC7", "M66_SC1", 
                "M67_SC0", "M68_SC2", "M68_SC1", "M68_SC3", "M70_SC0", "M50_SC0", "M50_SC1", "M50_SC2"))

strain_don_cols <- rep(strain_don_cols, times = 2) %>% set_names(c(names(strain_don_cols), str_replace(names(strain_don_cols), "_SC", "C_SC")))


sample_name_nms_color <- c("MSC33", "MSC48", "MSC49", "MSC55", "MSC60", "MSC66", "MSC38", 
  "MSC40", "MSC50_MACS", "MSC53", "MSC57", "MSC24", "MSC37", "MSC39", 
  "MSC45", "MSC50_SLO", "MSC54", "MSC67", "MSC68", "MSC70", "MSC41", 
  "MSC51", "MSC25", "MSC3", "MSC13", "MSC14", "MSC1", "MSC50", "Other")

# donor_cols = c("#C1EAE5", "#A7B19D", "#73935B", "#669CE2", "#62E1D0", "#DB6F87", "#DC48A2", "#E593D6", "#E65840", "#B29BDB", "#66EA54", "#B8E9B5", "#AFE685", "#D6ACCC", "#76586A", "#D8A46B", "#736AD9", "#62919E", "#65E297", "#E3AE3B", "#DA4CE6", "#BF6BD1", "#6FCEE8", "#8537E1", "#E3E544", "#D4D0E8", "#B0E653", "#E1A797", "#DFDC7F", "#E7E3C4", "#DAEC3F", "#D8A46B",'#ededed' )[1:length(sample_name_nms_color)] %>% set_names(str_remove(sample_name_nms_color, "SC"))
donor_cols_msc = c("#C1EAE5", "#A7B19D", "#73935B", "#669CE2", "#62E1D0", "#DB6F87", "#DC48A2", "#E593D6", "#E65840", "#B29BDB", "#66EA54", "#B8E9B5", "#AFE685", "#D6ACCC", "#76586A", "#D8A46B", "#736AD9", "#62919E", "#65E297", "#E3AE3B", "#DA4CE6", "#BF6BD1", "#6FCEE8", "#8537E1", "#E3E544", "#D4D0E8", "#B0E653", "#E1A797", '#ededed' )[1:length(sample_name_nms_color)] %>% set_names(sample_name_nms_color)
donor_cols_msc <- rep(donor_cols_msc, times = 2) %>% set_names(c(names(donor_cols_msc), paste0(names(donor_cols_msc), "C")))

donor_cols = donor_cols_msc
names(donor_cols) = str_remove(names(donor_cols_msc), "SC")
# donor_cols <- rep(donor_cols, times = 2) %>% set_names(c(names(donor_cols), paste0(names(donor_cols), "C")))



## strain plotting factor levels
fld_strns = paste0('SC', 1:12)
strain_dset_lvls = c( fld_strns, 'Field', 'Lab', 'Doublet')

sp_col_new <- c(`Late ring` = "#78C679", 
                # `Early trophozoite` = "#FEEEAA",
                `Early trophozoite` = "#FA8072",
                `Female (LE)` = "#CBB282",
                `Early ring` = "#D1EC9F",
                # `Late trophozoite` = "#FEB24C", 
                # `Late trophozoite` = "#FEB24C", 
                `Late trophozoite` = "#ee2008", 
Female = "#551A8B", `Male (LE)` = "#7e5522", Male = "mediumpurple",
`Gametocyte (male)` = "#1a6c8b", `Gametocyte (female)` = "#d1b252",
Lab = "#D3D3D3", Unassigned = "#DCDCDC", `Failed QC` = "grey",
`Not assigned` = "grey", Atlas = "grey", `Early schizont` = "#C9E8F1",
`Late schizont` = "#85B1D3", 
# `Gametocyte (developing)` = "thistle",
# `Gametocyte (developing)` = "#ffb2e5",
`Gametocyte (developing)` = "#ed009d",
`early female` = "#749E89", `late female` = "#4E6D58", Asexual = "#FF6961",
`NULL` = "lightgrey", `Field unlabelled` = "#D2C1DC", "weak_score" = "#c9c9c9", "unclear" = "#dddddd",`discrepantFM` =  "#e5d0ff", "Schizont" = '#85B1D3')

# Create abbreviated color map for plotting (short names used in figures)
sp_col_new_abrv <- sp_fld_colors2
names(sp_col_new_abrv) <- str_replace_all(names(sp_fld_colors2), c("discrepantFM" = "dFM", "weak_score" = "WS", "Early " = "E", "Late " = "L", "ring" = "R", "trophozoite" = "T", "Female" = "F", "Male" = "M", " \\(LE\\)" = "LE", "Gametocyte \\(developing\\)" = "Gd", "Missing" = "Mis"))


## Colors for the different methods for calling cells
cell_call_mthd_col <- c("Cell" = "#FFA808", "CellCr" = "#BC60D5", "CellEd" = "#99E564", "CellCb.5"= "#93DACC", "CellCb.8" = "#DC7886", other_calls = "#FF474C", "nonCell" = 'grey', "empty_bc" = '#ededed')

## SEurat clusters for integrated asexuals
# seu_cols <- brewer.pal(12, name = "Set3") %>% set_names(c(0:11))
seu_cols_base <- c(met.brewer(name = "Juarez"), brewer.pal(5, name = "Set3"))
seu_cols <- setNames(rev(colorRampPalette(seu_cols_base)(19)), 0:18)
seu_cols

## Random colours for pseudotime chunks for downsampling along pseudotime
# rand_color(n = 21) %>% dput()
ptc_cols <- c("#068B4CFF", "#4FB0C0FF", "#74DC09FF", "#3377ff", "#E769ACFF", 
              "#8A1708FF", "#048461FF", "#3BDF3FFF", "#5559B5FF", "#557406FF", 
              "#2A158CFF", "#B48606FF", "#CA2F4DFF", "#E75FFDFF", "#E35655FF", 
              "#42CFC3FF", "#1BA0B9FF", "#D481E8FF", "#958FEBFF", "#4A7406FF", 
              "#540B9CFF") %>% set_names(paste0("PT.", 1:21))


## Random colours for pseudotime chunks for DE along pseudotime
# rand_color(n = 9, transparency = 0.2) %>% dput()
ptc_cols2 <- c("#8FBE68CC", "#85F3D6CC", "#92C70ACC", "#CC0D28CC", "#958FEBFF", 
               "#0C9592CC", "#CAE37DCC", "#116575CC") %>% set_names(paste0("PT.", 2:9))


## Colours for the 16 leiden clusters for the UMAP
green_shades <- c("#D9F0A3", "#ADDD8E", "#41AB5D", "#45AB8D", "#238443", "#005A32", "#009E73") %>% set_names(c(15, 10,9,4,2,8, 14))
medium_purple_shades <- c("#8A66D2", "#6A4DBA", "#8d97f4") %>% set_names(c(0, 6, 11))
deep_purple_shades <- c("#F5DEB3","#E3C58F", "#C6A15B", "#A67C00", "#8B5A2B", "#845422") %>% set_names(c(1, 5,7,12,13,16)) 
femLE <- c("#CC79A7") %>% set_names(c(3))
empty_labl <- c("lightgrey") %>% set_names(c(" "))
# femLE <- c("#7cf4ec") %>% set_names(c(3))

leiden_cols <- c(green_shades, medium_purple_shades, deep_purple_shades, femLE, empty_labl)

## Key samples chosen for good quality to include in the manuscript
key_samples <- c("MSC3", "MSC13", "MSC14", "MSC33C", "MSC40C", "MSC45", "MSC48C", "MSC49C", "MSC50_MACS", "MSC50_SLOC", "MSC53", "MSC54", "MSC57", "MSC60C", "MSC66C", "MSC67", "MSC68C", "MSC70C")
key_samples <- c(key_samples, str_remove(key_samples, "SC") )
don_order <- unique(str_remove(key_samples, "C$"))
don_order_m <- unique(str_remove_all(key_samples, "SC|C$"))

### Asymptomatic Symptomatic colours
sympto_cols <- c("Asymptomatic" = "#97CBEFFF",  "Symptomatic" = "#EA2B60FF")
# c("Asymptomatic" = "#80F097FF",  "Symptomatic" = "#B59FF1FF")

sample_cols = c("#630CA4FF","#37C4BEFF", "#D647CFFF", "#F7C0AAFF",  "#0BFEAAFF", "#F0F695FF", "#39F6D3FF", "#FC65CCFF", "#99F868FF", "#09C89FFF", "#C3EB72FF", "#007343FF", "#C058DDFF", "#D46382FF", "#F31061FF", "#BEACF8FF") %>% set_names(c("D0", "D10a", "D3", "D5", "M14", "M33C", "M40C", "M48C", "M49C", "M50_SLOC", "M54", "M60C", "M66C", "M67", "M68C", "M70C"))


## strain coloors
# strain_dset_lvls_m = paste(rep(str_remove(sample_id_lst[1:4], "SC"), each=length(fld_strns)), fld_strns, sep = '_')
# strain_dset_lvls_don = c(strain_dset_lvls_m, strain_dset_lvls)
# strain_dset_lvls_don_col = randomcoloR::distinctColorPalette(k = length(strain_dset_lvls_m)) %>% set_names(strain_dset_lvls_m) %>% c(., strain_cols_n)