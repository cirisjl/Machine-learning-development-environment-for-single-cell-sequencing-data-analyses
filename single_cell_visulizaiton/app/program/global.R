# ----------------------------------------
# --          PROGRAM global.R          --
# ----------------------------------------
# USE: Global variables and functions
#
# NOTEs:
#   - All variables/functions here are
#     globally scoped and will be available
#     to server, UI and session scopes
# ----------------------------------------


# -- Setup your Application --
g_app_version <- "0.1.0"

set_app_parameters(title = "Single Cell RNA-Seq Data Explorer (SCRSDE)",
                   titleinfo = HTML(sprintf(readLines("program/data/about.html"), g_app_version)),
                   loglevel = "DEBUG",
                   showlog = FALSE,
                   app_version = g_app_version)

rm(list=ls())

library(plotly)
require(Matrix)

# -- PROGRAM --
suppressPackageStartupMessages(library(canvasXpress))

source('program/fxn/module_heatmapDownloadableTable.R')
source('program/fxn/supporting_data.R')
source('program/fxn/supporting_plots.R')
source('program/fxn/supporting_misc.R')
source('program/fxn/diff_expression.R')

# Variables

g_differential_logfc_threshold <- 0.5
g_differential_gene_threshold  <- 1000
g_differential_pct_threshold   <- 1
g_differential_min_no_cells    <- 3

# ui
g_add_top_genes_options  <- c('Off' = 'off', 'Top 10' = 'top10', 'Top 30' = 'top30')

# File upload
options(shiny.maxRequestSize = 500*1024^2) #max 500 MB
g_error_field    <- "error"
g_missing_fields <- "missing_fields"

#load("/home/philge/singlecell/iPSC/cell_count_iPSC.RData")
load("/srv/shiny-server/program/data/cell_count_SFG_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_EC_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_EC_Exc.RData")
load("/srv/shiny-server/program/data/cell_count_EC_Inh.RData")
load("/srv/shiny-server/program/data/cell_count_SFG_Exc.RData")
load("/srv/shiny-server/program/data/cell_count_SFG_Inh.RData")
load("/srv/shiny-server/program/data/cell_count_MTG_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_V1_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_ACC_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_VISp_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_ALM_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_ACA_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_MOpC_eis4.RData")
load("/srv/shiny-server/program/data/cell_count_MOpN_eis4.RData")

# load("/home/philge/singlecell/allen/odds_ratio/aca_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/aca_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/acc1_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/acc1_alm_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/acc1_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/acc1_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/acc1_visp_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/alm_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/alm_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/alm_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mopc_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg1_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg1_alm_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg1_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg1_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg1_visp_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg_acc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/mtg_v1_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/v11_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/v11_alm_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/v11_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/v11_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/v11_visp_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/v1_acc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/visp_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/visp_alm_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/visp_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/visp_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/ec_exc_sfg_exc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/allen/odds_ratio/ec_exc_sfg_inh_odds_ratio_v4.RData")
# load("/home/philge/singlecell/allen/odds_ratio/ec_inh_sfg_exc_odds_ratio_v4.RData")
# load("/home/philge/singlecell/allen/odds_ratio/ec_inh_sfg_inh_odds_ratio_v4.RData")
# load("/home/philge/singlecell/allen/odds_ratio/sfg_exc_sfg_inh_odds_ratio_v4.RData")

# load("/home/philge/singlecell/leng_et_al/odds_ratio/v1_ec_inh_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/acc_ec_exc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/acc_ec_inh_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/v1_ec_exc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_exc1_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_exc1_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_inh1_visp_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_exc_ec_inh_odds_ratio_v4.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_exc1_alm_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_inh1_aca_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_inh1_alm_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_inh1_mopn_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_exc1_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_exc1_visp_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/ec_inh1_mopc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/mtg_ec_exc_odds_ratio_v3.RData")
# load("/home/philge/singlecell/leng_et_al/odds_ratio/mtg_ec_inh_odds_ratio_v3.RData")
