# ==================================================================================================================
# File name:      Catch_reconstruction_uncertainty.R      (ONE script for the whole study)
#
# TITLE:          How much does the catch-reconstruction method change the reconstructed catch series of
#                 Brazilian marine shrimp? A hide-and-refill test and its consequence for stock-status indicators
#
# AUTHOR:         Silva, M.L.S. (UFRPE/PGRPAq, CAMMARADA project) - matheus.ssilva@ufrpe.br
# CO-AUTHORS:     [to be added]
#
# CREATED:        2026-10-30
# LAST MODIFIED:  2026-11-12
# VERSION:        0.3 (written BEFORE any inferential analysis; see analysis_plan_OSF_EN.docx)
# LICENSE:        CC-BY 4.0
#
# PURPOSE
#   Quantify how much the CHOICE OF CATCH-RECONSTRUCTION METHOD changes the reconstructed catch of Brazilian marine
#   shrimp, how large the error of each method is when it fills gaps of known size and position, and whether the
#   differences are large enough to change a simple stock-status indicator (ZBRT biomass depletion, B/K).
#   Analysis plan written before the inferential analyses: analysis_plan_OSF_EN.docx (secondary data).
#
#   UNIT OF ANALYSIS = the TOTAL shrimp catch of each STATE (17 states): the sum of every bulletin category (generic
#   "camaroes" + rosa + branco + sete-barbas). State totals are summed to the 4 coastal regions (N, NE, SE, S) and to Brazil.
#   Why totals and not species stocks? Most of the observed catch is in the generic category, which the published
#   reconstructions split among species; the total is the only quantity that can be compared with the observed data
#   without that split.
#
#   Idea in plain words:
#     - Three published reconstructions exist (Freire et al. 2021, Silva & Andrade 2023, MPA 2026). They use different
#       methods: LOESS/interpolation (Freire), GLM (Silva & Andrade), neural networks (MPA).
#     - OBSERVED catches = national bulletins 1962-2007 + MPA collected landings 2023-2024. Everything else is estimated.
#     - We refill the gaps ourselves with four methods applied identically to every state series:
#         M0 = naive: the nearest observed value is carried over (the reference every real method must beat)
#         M1 = linear interpolation + LOESS at the edges     (the family of Freire et al. 2021)
#         M2 = gamma GLM, log link, year^2 and sin(year)      (Silva & Andrade 2023, chosen by AIC)
#         M3 = small neural network                          (the family of the MPA)
#       They are SIMPLIFIED versions: same input for all (only the series' own years), no environmental covariates, no
#       expert opinion, no allocation of the generic category. Differences therefore come from the method alone.
#
# THE THREE STEPS (and the question each one answers)
#   STEP 1  Compare the published reconstructions (descriptive, no hypothesis).
#           Do they agree with each other WHERE THEY HAD TO FILL A GAP (state-years without an observed value)? In the years with a
#           bulletin the three keep the reported value, so agreement there would be by construction. Metrics: ratio, typical
#           difference, trend agreement. By state (median across the 17), by region and for Brazil. Also a check of how well each
#           reproduces the bulletins, and the share of the observed total that is identified to species, by state.
#   STEP 2  Hide-and-refill test on the 17 state series (hypothesis: EDGE ERROR).
#           Known observed years are hidden in blocks and refilled by M0-M3; the error is |ln(predicted+1) - ln(observed+1)|.
#           Gap types: interior4 and interior10 (4- and 10-year blocks inside the bulletin era), interior_long (the last
#           5 bulletin years, which leaves a ~20-year hole bounded by the 2023-24 data), start (first 5 observed years)
#           and end (2023-24 hidden and predicted from the bulletin era alone).
#           Edge-error hypothesis: mean error when filling the START or END of a series (edge) is >= 1.5 x the mean
#           error when filling INSIDE it, and the 95% CI of the difference excludes 0.
#           Model: abs_e ~ method * gap_type + (1 | series_id/mask_id)   [lme4 + emmeans]
#   STEP 3  Rebuild every state series ONCE with M0-M3 (observed values kept, gaps filled), sum to regions and Brazil,
#           compare the methods with each other and with the published reconstructions summed the same way, and
#           apply ZBRT (Zhou et al. 2017, datalimited2::zbrt) to each series.
#           ZBRT is run on every state, region and Brazil, over the WHOLE series (B/K for every year).
#           Divergence hypothesis (ONE, decided on the 17 STATES): the spread of B/K among M0-M3 exceeds ZBRT's own 95%
#           interval width in more than half of the years, in more than half of the states (at least 9 of 17).
#           Regions and Brazil are sums of the states, reported for summary and figures, not used for the decision.
#           The share of diverging years in each era (bulletin era, long gap 2008-2022, recent 2023-24) is reported as a
#           descriptive result, not as a hypothesis.
#   (Steps 2 and 3 use the SAME four methods: Step 2 applies them to years whose true value we know, to measure their
#    error; Step 3 applies them to the real gaps.)
#
# INPUT (all in the working folder 'Code and results')
#   FINAL_RECONSTRUCTED_Brazil_1950_2015_CommercialEtapaII_04072021.xlsx   Freire et al. 2021 (row level)
#   Banco_reconstrucao_marinha_2026.xlsx   MPA: reconstruction 1950-2022 + collected landings 2023-2024
#   data_AL.csv ... data_SP.csv            Silva & Andrade 2023 GLM outputs by STATE (17 files, 1946-2011). The columns
#        c.rosa.recon, c.branco.recon and c.7barbas.recon are the GLM prediction for ALL years (also where a bulletin exists). Their
#        series here = observed state total where a bulletin exists + the sum of these three predictions elsewhere (section 1.8.2).
#   Camaroes_final-3.csv                   OBSERVED bulletin landings by state x year, 1962-2007, in commercial categories
#        (generic shrimp, white, pink, seabob). Years of all zeros (1946-1961, 1965, 2008-2011) = no bulletin. A zero is
#        treated as 'no data' (it cannot be told from a real zero).
#
# DESIGN DECISIONS (author, 31 Oct - 3 Nov 2026)
#   - observed = bulletins + MPA collected landings 2023-2024; ALL observed data are used even though the collection
#     methods differ between eras (stated limitation)
#   - the generic category is NOT allocated to species; the unit is the total
#   - CMSY++ and the rainfall test were removed from this study (limitations)
#
# HOW TO RUN
#   Set the switches in RUN (section 0). Section 1 builds the database from the raw files and writes it to data/;
#   later steps read what it wrote, so any step can be rerun alone.
#   Neural networks are the heavy part: each series runs in a fresh R process and is saved to models/ (restart-safe).
#   ZBRT needs the package datalimited2 (remotes::install_github('cfree14/datalimited2')) and JAGS.
#
# STRUCTURE
#   0.  Setup and configuration
#   1.  Data: ONE unified database, observed vs estimated, with source          -> data/
#   2.  STEP 1 - compare the published reconstructions (state totals)           -> results/Table_1-2, Fig_1-3, Fig_9
#   3.  Load the observed state series used in Steps 2-3
#   4.  The four reconstruction methods M0-M3 (functions)
#   5.  STEP 2 - hide-and-refill test                                           -> results/Table_3, Fig_4
#   6.  STEP 3 - rebuild, compare, ZBRT                                         -> results/Table_4, Table_S2, Fig_5-8
#   7.  Session info
#
# NOTE ON DEVIATIONS: any change to the analysis plan after it was written is listed here with date and reason.
#   2026-11-10  Step 3: the divergence hypothesis is decided on the 17 states (more than half of the states), not on the
#               4 regions (3 of 4), so that the whole study is at the scale of the state, as in Steps 1-2. The second
#               hypothesis (concentration in the long gap) became a descriptive result (share of diverging years by era).
#               Both versions were computed on the same results and give the same conclusion (not supported).
#   2026-11-12  Step 1 and the Silva & Andrade series. The first version used the full Silva & Andrade GLM prediction (all years) as their
#               series and compared the three reconstructions in the bulletin era. That was unfair: Freire et al. 2021 and MPA 2026
#               keep the reported values and fill only the gaps, so in the bulletin years they match the bulletins by construction,
#               while the full model fit does not. Now (a) the Silva & Andrade series is observed total where a bulletin exists and the GLM
#               prediction elsewhere, and (b) the reconstructions are compared with each other only in the FILLED years (state-years with no
#               observed value). The comparison with the bulletins stays as a check. Step 2 is unchanged (it does not use the published
#               series); in Step 3 only the published series (Silva & Andrade) change, so their ZBRT and the spread among the published
#               reconstructions are recomputed.
#   2026-11-12  Maranhao moved from the North to the Northeast (North = Amapa and Para). The study is at the scale of the state and is not
#               a stock analysis; regions are only sums of states used for presentation, so they follow the usual geographic division.
#               State-level results do not change; regional sums (Figs 2, 5, 7 and the regional tables) change for North and Northeast.
# ==================================================================================================================

# ------------------------------------------------------------------------------------------------------------------
# 0. SETUP AND CONFIGURATION
# ------------------------------------------------------------------------------------------------------------------
# Working folder = 'Code and results'. The environment variable lets the script run on another machine.
wd <- Sys.getenv('CRU_DIR', 'C:/Matheus/Universidade/Posdoc/Camarada/Catch reconstruction uncertainty/Code and results')
setwd(wd)

library(data.table)   # fast table handling
library(readxl)       # read the raw .xlsx inputs
library(openxlsx)     # write the organised .xlsx and the result tables
library(ggplot2)      # figures
library(lme4)         # mixed model of Step 2
library(emmeans)      # contrasts of Step 2
library(callr)        # run each series in a fresh R process (keeps neural-network memory under control)
library(lmerTest)     # same models as lme4 plus Satterthwaite degrees of freedom (df) in the tables; install.packages('lmerTest')
library(datalimited2) # ZBRT (zbrt) for Step 3; remotes::install_github('cfree14/datalimited2'); needs JAGS installed on the computer

dir.create('data',    showWarnings = FALSE)   # unified database and the files derived from it
dir.create('models',  showWarnings = FALSE)   # per-series checkpoints (.rds)
dir.create('results', showWarnings = FALSE)   # tables and figures

# ---- 0.1 Look of all figures (same style as the life-history paper: Shrimp_life_history.R) --------------------------------
# Minimal theme, light horizontal grid only, bold centred panel titles, legend at the bottom, a short note under the
# figure. Colours are fixed per meaning so a colour never changes meaning between figures:
#   one colour per published reconstruction and one per reconstruction method M0-M3; observed values are dark points.
# Files: PNG at 300 dpi (publication quality, white background) and a vector PDF with the same name.
ink <- '#262626'; ink2 <- '#52514e'; grid_col <- '#E5E5E2'
source_cols <- c('Freire et al. 2021' = '#2a78d6', 'Silva & Andrade 2023' = '#eb6834', 'MPA 2026' = '#1baf7a')
# The four methods use hues that do NOT repeat the three source colours above (blue/orange/green belong to the published
# reconstructions in every figure), so a method line is never mistaken for a publication.
method_cols  <- c(M0 = '#8c8b86', M1 = '#7d5ba6', M2 = '#d9a21b', M3 = '#d1457c')   # grey, purple, gold, magenta
obs_col     <- ink

theme_cru <- function() {
  theme_minimal(base_size = 11) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          strip.text = element_text(face = 'bold', hjust = 0.5, colour = ink),
          legend.position = 'bottom', legend.title = element_text(face = 'bold', size = 9),
          plot.title = element_text(face = 'bold', size = 12), plot.subtitle = element_text(colour = ink2, size = 9),
          plot.caption = element_text(colour = ink2, size = 8, hjust = 0), plot.caption.position = 'plot',
          panel.spacing = unit(1.1, 'lines'), plot.margin = margin(8, 12, 8, 8))
}
say <- function(...) if (isTRUE(CFG$verbose)) cat(format(Sys.time(), '%H:%M:%S'), '|', ..., '\n', sep = ' ')   # timestamped progress message
save_fig <- function(plot, name, width, height) {      # PNG 300 dpi + vector PDF, same name, in results/
  say('saving figure', name)
  ggsave(file.path('results', paste0(name, '.png')), plot, width = width, height = height, dpi = 300, bg = 'white')
  try(ggsave(file.path('results', paste0(name, '.pdf')), plot, width = width, height = height, device = cairo_pdf), silent = TRUE)
}

# ---- 0.2 Switches: which parts to run ------------------------------------------------------------------------------
RUN <- list(
  build_data = TRUE,    # section 1: build the unified database from the raw files
  step1      = TRUE,    # section 2: compare the published reconstructions
  step2      = TRUE,    # sections 3-5: hide-and-refill test
  step3      = TRUE    # sections 3-4 and 6: rebuild the series, compare, ZBRT
)

# ---- 0.3 Settings fixed in the analysis plan ------------------------------------------
CFG <- list(
  # --- raw files ---
  f_freire   = 'FINAL_RECONSTRUCTED_Brazil_1950_2015_CommercialEtapaII_04072021.xlsx',
  f_mpa      = 'Banco_reconstrucao_marinha_2026.xlsx',
  f_camaroes = 'Camaroes_final-3.csv',
  f_silva_state = 'data_%s.csv',     # Silva & Andrade state files; %s = two-letter state code

  # --- observed data and windows ---
  mpa_obs_from = 2023,   # MPA rows from this year on are collected landings (observed); before = reconstruction
  mpa_recon_last = 2022, # last year of the MPA reconstruction
  silva_last   = 2011,   # last year of the Silva & Andrade state files (data_XX.csv). The SAME state files are used in every figure and table
                         # (states, regions and Brazil = sum of states). In 2006-2011 the N and NE files hold GLM extrapolations that the
                         # authors later corrected in the published series; they are used as they are (flagged in the plan).
  zero_as_missing = TRUE,# a ZERO in the bulletin file cannot be told apart from 'not reported' -> treated as missing
  bulletin_last = 2007,  # last year with a national bulletin
  year_first   = 1962,   # first year of the state series used in Steps 2-3 (first bulletin year)
  year_last    = 2024,   # last year of the state series used in Steps 2-3 (MPA collected landings)
  src_obs        = 'Silva & Andrade 2023 (bulletins)',   # observed national bulletins (name inside the unified records)
  obs_source     = c('Silva & Andrade 2023 (bulletins)', 'MPA 2026'),   # sources that provide OBSERVED values
  src_recon      = c('Freire et al. 2021', 'Silva & Andrade 2023', 'MPA 2026'),   # the three published reconstructions

  # --- Step 2 (hide-and-refill) ---
  methods_to_run = c('M0', 'M1', 'M2', 'M3'),   # use c('M0', 'M1', 'M2') for a quick test without neural networks
  gap_interior  = c(4, 10), # interior block lengths hidden (years); 10 reflects the long gap 2008-2022
  n_place       = 5,        # random placements per interior length
  gap_edge      = 5,        # 'start' block length and the block hidden for 'interior_long' (years)
  min_obs_hide  = 15,       # a series is used in Step 2 only if it has this many observed years
  seed_sample   = 2026,     # seed of masks and networks
  nn_reps_A     = 2,        # network repeats in Step 2
  nn_reps_B     = 3,        # network repeats in Step 3
  edge_ratio    = 1.5,      # edge-error hypothesis: mean |e| edge / interior >= this value

  # --- messages ---
  verbose = TRUE,           # TRUE = print what is running (time, state, mask, method) while the script works; FALSE = silent

  # --- Step 3 ---
  min_share_states = 0.5,   # divergence hypothesis: supported if MORE than this share of the states meet the criterion
  diverge_share = 0.5,      # a unit 'diverges' when the methods differ by more than ZBRT's own uncertainty in MORE than this share of the years
  zbrt_states = TRUE        # also run ZBRT on each of the 17 states (about 120 more runs, ~8 min); FALSE = regions and Brazil only
)

# ---- 0.4 Names used everywhere ------------------------------------------------------------------------------------
# State -> region. The whole study works at the scale of the STATE (17 coastal states); this is NOT a stock analysis and no fishery
# system or stock is defined. Regions (and Brazil) are only SUMS OF STATES used to present results. Regions follow the usual
# geographic division of Brazil, so Maranhao belongs to the Northeast (North = Amapa and Para; Northeast = Maranhao to Bahia).
# (Until 2026-11-12 Maranhao was counted with the North, following the stock assessment of Silva et al. 2026.)
REGION <- c('Amapá' = 'N', 'Pará' = 'N',
            'Maranhão' = 'NE', 'Piauí' = 'NE', 'Ceará' = 'NE', 'Rio Grande do Norte' = 'NE', 'Paraíba' = 'NE', 'Pernambuco' = 'NE',
            'Alagoas' = 'NE', 'Sergipe' = 'NE', 'Bahia' = 'NE',
            'Espírito Santo' = 'SE', 'Rio de Janeiro' = 'SE', 'São Paulo' = 'SE',
            'Paraná' = 'S', 'Santa Catarina' = 'S', 'Rio Grande do Sul' = 'S')
# Two-letter code -> state name (bulletin file and Silva & Andrade state files).
UF  <- c(AP = 'Amapá', PA = 'Pará', MA = 'Maranhão', PI = 'Piauí', CE = 'Ceará', RN = 'Rio Grande do Norte', PB = 'Paraíba',
         PE = 'Pernambuco', AL = 'Alagoas', SE = 'Sergipe', BA = 'Bahia', ES = 'Espírito Santo', RJ = 'Rio de Janeiro',
         SP = 'São Paulo', PR = 'Paraná', SC = 'Santa Catarina', RS = 'Rio Grande do Sul')
REG_LAB <- c(N = 'North', NE = 'Northeast', SE = 'Southeast', S = 'South')
UNITS_REGION <- c('N', 'NE', 'SE', 'S')            # north to south
thousand <- function(x) format(x / 1000, big.mark = ',', trim = TRUE)    # axis in thousand tonnes

# ---- 0.5 Helper functions shared by Steps 1 and 3 ---------------------------------------------------------------------
# agree_fun: three agreement metrics between two catch series a and b over the SAME years (yr sorted)
#   ratio              total of a / total of b (1 = same total)
#   typical_difference exp(mean(|ln(a+1) - ln(b+1)|)) - 1  (0.25 = typically about 25% apart)
#   trend_agreement    correlation of year-to-year changes of ln(a+1), ln(b+1), CONSECUTIVE years only (+1 = move together)
agree_fun <- function(yr, a, b) {
  la <- log(a + 1); lb <- log(b + 1)
  i  <- which(diff(yr) == 1)                          # consecutive years only: gaps are skipped, not bridged
  ta <- NA_real_
  if (length(i) >= 5) { da <- la[i + 1] - la[i]; db <- lb[i + 1] - lb[i]; if (sd(da) > 0 && sd(db) > 0) ta <- cor(da, db) }
  data.table(n_years = length(yr), ratio = sum(a) / sum(b), typical_difference = exp(mean(abs(la - lb))) - 1, trend_agreement = ta)
}
# fmt3: "median (min to max)" of a vector, used in the summary tables
fmt3 <- function(x) if (all(is.na(x))) 'n/a' else sprintf('%.2f (%.2f to %.2f)', median(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
sum_na <- function(x) if (anyNA(x)) NA_real_ else sum(x)     # a sum that stays NA if any element is missing


# ------------------------------------------------------------------------------------------------------------------
# 1. DATA: ONE UNIFIED DATABASE, OBSERVED vs ESTIMATED, WITH SOURCE
# ------------------------------------------------------------------------------------------------------------------
# Reads the four raw files, standardises names, flags every record as observed or estimated and writes:
#   data/catch_unified_base.xlsx   the organised database (records, annual tables, coverage, QC)
#   data/catch_state_total.csv   TOTAL shrimp per source x state x year (observed + three reconstructions): input of Steps 1-3
# STATUS FLAG (the key column)
#   observed  = a landing reported in the national bulletins (Camaroes_final-3.csv, 1962-2007) or collected by the MPA
#               in 2023-2024. Freire's 'Landing (reported)' rows are kept for COMPARISON only: she attributes species to
#               categories that the bulletins do not separate (and 2008-2015 'reported' are local programmes in 5 states).
#   estimated = anything produced by a reconstruction (interpolation, trend, expert opinion, MPA 1950-2022 ...)
#   For Freire rows 'status_detail' keeps the finer classes: observed | observed_split | reported_modified | estimated.
if (RUN$build_data) {
say('SECTION 1: reading raw files and building the unified database (about 1-2 min)')

# ---- 1.0 The 12 stock units (group x region) of Silva et al. 2026 ------------------------------------------------
# By-product of the shared compiled database only (sheet 'annual_stock_unit'); NOT used by the analysis of this study, which works with
# state totals. The unit regions use the REGION map above (Maranhao in the Northeast).
STOCKS <- data.table(group  = c('pink', 'pink', 'brown', 'brown', rep('seabob', 4), rep('white', 4)),
                     region = c('S', 'SE', 'N', 'NE', rep(c('N', 'NE', 'S', 'SE'), 2)))
STOCKS[, stock := paste(group, region, sep = '_')]

# ------------------------------------------------------------------------------------------------------------------
# 1.1 READ THE THREE SPREADSHEETS
# ------------------------------------------------------------------------------------------------------------------
# Only the first 15 columns of Freire's sheet hold data; reading just those avoids harmless readxl warnings about
# stray text in the two last (almost empty) columns ('Expecting logical ... got Checar').
fr  <- as.data.table(read_excel(CFG$f_freire, sheet = 'Reconstruction', range = cell_cols(1:15), guess_max = 1e6,
                                .name_repair = 'unique_quiet'))
mpa <- as.data.table(read_excel(CFG$f_mpa))
cat('Freire rows:', nrow(fr), '| MPA rows:', nrow(mpa), '\n')

# ------------------------------------------------------------------------------------------------------------------
# 1.2 NAME STANDARDISATION
# ------------------------------------------------------------------------------------------------------------------
# 1.2.1 Taxon: the MPA file uses the genus Penaeus for everything; the life-history paper uses Litopenaeus and
#     Farfantepenaeus. We convert to the life-history names so the two files can be compared.
std_taxon <- function(x) {
  x <- trimws(x)
  x <- sub('\\s+spp\\.?$', ' spp', x)                       # 'spp.' and 'spp' -> 'spp'
  x[x == 'Penaeus schmitti']    <- 'Litopenaeus schmitti'
  x[x == 'Penaeus subtilis']    <- 'Farfantepenaeus subtilis'
  x[x == 'Penaeus paulensis']   <- 'Farfantepenaeus paulensis'
  x[x == 'Penaeus brasiliensis'] <- 'Farfantepenaeus brasiliensis'
  x
}

# 1.2.2 Stock group. Rules (in this order):
#   - X. kroyeri -> seabob; L. schmitti -> white
#   - PINK-SHRIMP COMPLEX (F. subtilis, F. paulensis, F. brasiliensis, and the genus-level 'Farfantepenaeus spp' /
#     'Penaeus spp' with common name 'rosa' or 'camarão'): the bulletins do not separate these species, so, as in
#     Freire et al. 2021 and Silva & Andrade 2023, the complex is assigned by REGION:
#         N and NE -> brown (F. subtilis);   S and SE -> pink (F. brasiliensis + F. paulensis, aggregated)
#   - Artemesia longinaris (barba-ruça) and Pleoticus muelleri (santana): kept, but out of the 12 stocks
#   - Penaeidae / 'Penaeus spp' / 'camarões' / 'estuarino' (no species, not rosa) -> 'generic' (the undifferentiated
#     shrimp category of the bulletins); it is tracked because the reconstructions use it to estimate species catches
#   - anything else (e.g. common name 'branco' with a pink taxon): unassigned
#   - non-penaeid and farmed/exotic shrimps are excluded (list written to the 'excluded' sheet)
assign_group <- function(taxon, common, region) {
  common <- tolower(common)
  g <- rep(NA_character_, length(taxon)); rule <- rep(NA_character_, length(taxon))
  sp <- function(t, grp) { i <- taxon == t; g[i] <<- grp; rule[i] <<- 'species' }
  sp('Xiphopenaeus kroyeri', 'seabob');  sp('Litopenaeus schmitti', 'white')
  sp('Artemesia longinaris', 'barba_ruca'); sp('Pleoticus muelleri', 'santana')
  rosa <- taxon %in% c('Farfantepenaeus subtilis', 'Farfantepenaeus paulensis', 'Farfantepenaeus brasiliensis') |
          (taxon %in% c('Farfantepenaeus spp', 'Penaeus spp') & (grepl('rosa|vermelho', common) | common %in% c('camarão', 'camarao')))
  g[rosa] <- ifelse(region[rosa] %in% c('N', 'NE'), 'brown', 'pink'); rule[rosa] <- 'regional_rosa'
  gen <- is.na(g) & (taxon == 'Penaeidae' | taxon == 'Penaeus spp')
  g[gen] <- 'generic'; rule[gen] <- 'generic_category'
  g[is.na(g)] <- 'unassigned'; rule[is.na(rule)] <- 'unassigned'
  list(group = g, rule = rule)
}
NON_TARGET <- 'Plesionika|Aristaeopsis|Aristaeomorpha|Aristeus|Sicyonia|Macrobrachium|Parapenaeus|Decapoda|vannamei|monodon'
IS_PENAEID <- 'Xiphopenaeus|Litopenaeus|Farfantepenaeus|Penaeus|Penaeidae|Artemesia|Pleoticus'

# ------------------------------------------------------------------------------------------------------------------
# 1.3 FREIRE ET AL. 2021  (row-level: landings reported and estimated, 1950-2015)
# ------------------------------------------------------------------------------------------------------------------
fr[, taxon_orig := TaxonName]
fz <- fr[grepl(IS_PENAEID, TaxonName) & !grepl(NON_TARGET, TaxonName)]       # shrimp rows
fz[, `:=`(source = 'Freire et al. 2021', state = OtherArea, year = Year, catch_t = CatchAmount_t,
          common_orig = PortugueseCommonName, sector = fifelse(tolower(Sector) == 'artisanal (ss, c)', 'Artisanal', 'Industrial'),
          origin = Input, reference = Reference, notes = fifelse(is.na(Notes), '', as.character(Notes)))]
fz[, taxon := std_taxon(TaxonName)]
fz[, region := REGION[state]]

# 1.3.1 STATUS. 'Landing (reported)' is NOT automatically observed: the Notes of many reported rows say the value was
#     replaced by a linear trend, or that a reported total was split among species/sectors. We read the notes.
RE_MODIFIED <- 'linear trend|replaced|interpolat|spline|Kalman|Akima|running average|trend from previous|monetary|Missing bulletin|ZERO|proportion in relation'
RE_SPLIT    <- 'split|values estimated based on fishery registered pattern|assumed to be artisanal'
fz[, status_detail := fcase(
  CatchType == 'Landing (estimated)',                        'estimated',
  CatchType == 'Landing (reported)' & grepl(RE_MODIFIED, notes, ignore.case = TRUE), 'reported_modified',
  CatchType == 'Landing (reported)' & grepl(RE_SPLIT,    notes, ignore.case = TRUE), 'observed_split',
  CatchType == 'Landing (reported)',                         'observed',
  default = NA_character_)]
fz[, status := fifelse(status_detail == 'observed', 'observed', 'estimated')]   # strict rule
fz[, record_note := notes]

# ------------------------------------------------------------------------------------------------------------------
# 1.3b SILVA & ANDRADE 2023 - OBSERVED BULLETIN LANDINGS (Camaroes_final-3.csv, commercial categories)
# ------------------------------------------------------------------------------------------------------------------
# The file is semicolon-separated with decimal comma. Columns used: state, year and the category totals
# c.marinho (generic marine shrimp), c.branco (white), c.rosa (pink) and c.7barbas (seabob).
cam <- fread(CFG$f_camaroes, sep = ';', dec = ',', encoding = 'UTF-8')
cm <- melt(cam[, .(fonte, estado, ano, c.marinho, c.branco, c.rosa, c.7barbas)], id.vars = c('fonte', 'estado', 'ano'),
           variable.name = 'cat', value.name = 'catch_t')
cm[, cat := as.character(cat)]
# 1.3b.1 Years with a bulletin: the national total of the year is > 0 (1946-1961, 1965 and 2008-2011 are all zeros = no data)
yr_ok <- cm[, .(tot = sum(catch_t)), by = ano][tot > 0, ano]
# 1.3b.2 A category exists from the first year in which its national total is > 0 (before: not reported, not zero)
cat_start <- cm[, .(first = min(ano[catch_t > 0 & ano %in% yr_ok])), by = cat]
cm <- merge(cm, cat_start, by = 'cat')[ano %in% yr_ok & ano >= first]
cz <- cm[, .(source = 'Silva & Andrade 2023 (bulletins)', state = UF[estado], year = ano, catch_t = catch_t,
             common_orig = fcase(cat == 'c.marinho', 'generic marine shrimp', cat == 'c.branco', 'white', cat == 'c.rosa', 'pink (rosa)', default = 'seabob (sete-barbas)'),
             taxon_orig = NA_character_, taxon = NA_character_, sector = 'all',
             status = fifelse(catch_t > 0 | !CFG$zero_as_missing, 'observed', 'estimated'),
             status_detail = fifelse(catch_t > 0 | !CFG$zero_as_missing, 'observed', 'zero_not_reported'),
             origin = paste('bulletin', fonte), reference = paste('Bulletins of', fonte, '(Camaroes_final-3.csv)'), record_note = '',
             rule = fifelse(cat == 'c.rosa', 'regional_rosa', 'bulletin_category'))]
cz[, region := REGION[state]]
cz[, group := fcase(common_orig == 'generic marine shrimp', 'generic', common_orig == 'white', 'white',
                    common_orig == 'seabob (sete-barbas)', 'seabob', common_orig == 'pink (rosa)', fifelse(region %in% c('N', 'NE'), 'brown', 'pink'))]

# ------------------------------------------------------------------------------------------------------------------
# 1.4 MPA 2026  (reconstruction 1950-2022 + collected landings 2023-2024)
# ------------------------------------------------------------------------------------------------------------------
mz <- mpa[grepl(IS_PENAEID, Táxon) & !grepl(NON_TARGET, Táxon) & !is.na(Estado)]
mz[, `:=`(source = 'MPA 2026', state = Estado, year = Ano, catch_t = `Captura (t)`, common_orig = `Nome Comum`,
          taxon_orig = Táxon, sector = Setor, origin = fifelse(Ano >= CFG$mpa_obs_from, 'MPA collected landings', 'MPA reconstruction'),
          reference = fifelse(Ano >= CFG$mpa_obs_from, 'MPA Boletim 2023-2024', 'MPA Boletim da Reconstrução 1950-2022'),
          record_note = fifelse(is.na(Notas), '', as.character(Notas)))]
mz[, taxon := std_taxon(Táxon)]
mz[, region := REGION[state]]
mz[, status_detail := fifelse(year >= CFG$mpa_obs_from, 'observed', 'estimated')]
mz[, status := status_detail]

# ------------------------------------------------------------------------------------------------------------------
# 1.6 STACK, AGGREGATE, COVERAGE
# ------------------------------------------------------------------------------------------------------------------
keep <- c('source', 'state', 'region', 'year', 'common_orig', 'taxon_orig', 'taxon', 'sector', 'catch_t', 'status', 'status_detail',
          'origin', 'reference', 'record_note')
rec <- rbind(fz[, ..keep], mz[, ..keep], fill = TRUE)
a   <- assign_group(rec$taxon, rec$common_orig, rec$region)
rec[, `:=`(group = a$group, group_rule = a$rule)]
rec <- rbind(rec, cz[, c(keep, 'group', 'rule'), with = FALSE][, setnames(.SD, 'rule', 'group_rule')], fill = TRUE)
rec[, in_scope := group %in% c('seabob', 'white', 'brown', 'pink')]
rec[, tracked := in_scope | group == 'generic']      # 'generic' is tracked in the annual tables but is not a stock
rec <- merge(rec, STOCKS[, .(group, region, stock)], by = c('group', 'region'), all.x = TRUE)   # NA = outside the 12 stocks
rec[, record_id := .I]
setcolorder(rec, c('record_id', 'source', 'status', 'status_detail', 'state', 'region', 'year', 'group', 'stock', 'in_scope',
                   'common_orig', 'taxon_orig', 'taxon', 'group_rule', 'sector', 'catch_t', 'origin', 'reference', 'record_note'))

# ---- 1.6.1 State x group x year (Freire and MPA only; Silva has no state) --------------------------------------------
# A state-group-year is 'observed' ONLY if every record that composes it is observed. If part of it is estimated
# (e.g. one species observed and a related one interpolated), the whole cell is 'estimated', so that an observed
# value is never a partial sum.
ann_state <- rec[source != 'Silva & Andrade 2023' & tracked,
                 .(catch_t = sum(catch_t), n_records = .N,
                   n_observed = sum(status == 'observed'),
                   catch_observed_t = sum(catch_t[status == 'observed']),
                   catch_estimated_t = sum(catch_t[status != 'observed'])),
                 by = .(source, state, region, group, year)]
ann_state[, status := fifelse(n_observed == n_records, 'observed', 'estimated')]
ann_state[, status_cell := fifelse(n_observed == n_records, 'all observed', fifelse(n_observed == 0, 'all estimated', 'mixed'))]

# ---- 1.6.2 Stock unit x year, by source -----------------------------------------------------------------------------
n_states <- table(REGION)
unit_f <- merge(ann_state, STOCKS[, .(group, region, stock)], by = c('group', 'region'))
unit <- unit_f[, .(catch_total_t = sum(catch_t), catch_observed_t = sum(catch_observed_t),
                   n_states_with_data = uniqueN(state), n_states_all_observed = sum(status == 'observed')),
               by = .(source, stock, group, region, year)]
unit[, n_states_region := as.integer(n_states[region])]
unit[, share_states_observed := n_states_all_observed / n_states_region]
setorder(unit, source, stock, year)

# ---- 1.6.3 Coverage of OBSERVED values per source x state x group ---------------------------------------------------
longest_run <- function(y) {
  y <- sort(unique(y)); if (!length(y)) return(0L)
  brk <- c(0, which(diff(y) != 1), length(y)); max(diff(brk))
}
cover <- ann_state[, {
  o <- year[status == 'observed']
  .(n_years_observed = as.integer(length(o)), first_obs = if (length(o)) as.integer(min(o)) else NA_integer_, last_obs = if (length(o)) as.integer(max(o)) else NA_integer_,
    longest_run_observed = as.integer(longest_run(o)), n_obs_2008_2015 = as.integer(sum(o >= 2008 & o <= 2015)),
    n_obs_2016_2022 = as.integer(sum(o >= 2016 & o <= 2022)), n_obs_2023_2024 = as.integer(sum(o >= 2023 & o <= 2024)))
}, by = .(source, state, region, group)]
setorder(cover, source, group, region, state)

# ---- 1.6.4 Name map (what each original name became) --------------------------------------------------------------
name_map <- rec[source != 'Silva & Andrade 2023', .(n_records = .N, catch_t_total = round(sum(catch_t), 1),
                years = paste(range(year), collapse = '-')),
                by = .(source, common_orig, taxon_orig, taxon, group, group_rule, in_scope)]
setorder(name_map, source, -catch_t_total)
excluded <- rbind(fr[!(grepl(IS_PENAEID, TaxonName) & !grepl(NON_TARGET, TaxonName)) & grepl('amar', PortugueseCommonName),
                     .(source = 'Freire et al. 2021', common_orig = PortugueseCommonName, taxon_orig = TaxonName, n_records = .N), by = .(PortugueseCommonName, TaxonName)][, .(source, common_orig, taxon_orig, n_records)],
                  mpa[!(grepl(IS_PENAEID, Táxon) & !grepl(NON_TARGET, Táxon)) & grepl('amar', `Nome Comum`),
                      .(source = 'MPA 2026', common_orig = `Nome Comum`, taxon_orig = Táxon, n_records = .N), by = .(`Nome Comum`, Táxon)][, .(source, common_orig, taxon_orig, n_records)])
excluded[, reason := 'non-penaeid, deep-water, freshwater or farmed/exotic shrimp (or fish named "bacamarte")']

# ------------------------------------------------------------------------------------------------------------------
# 1.7 CROSS-CHECKS AND QC LOG
# ------------------------------------------------------------------------------------------------------------------
qc <- list()
add_qc <- function(item, n, detail) qc[[length(qc) + 1]] <<- data.table(item = item, n = n, detail = detail)

add_qc('Freire rows labelled reported', fz[CatchType == 'Landing (reported)', .N], 'CatchType = Landing (reported)')
for (s in c('observed', 'observed_split', 'reported_modified', 'estimated'))
  add_qc(paste('Freire status_detail =', s), fz[status_detail == s, .N], 'rows after reading the Notes column')
add_qc('Freire rows reported but modified by a note', fz[status_detail == 'reported_modified', .N],
       'Notes mention trend/replacement/interpolation; treated as estimated under the strict rule')
add_qc('Freire estimated rows with note "Original data"', fz[CatchType == 'Landing (estimated)' & grepl('^Original data', notes), .N],
       'labelled estimated but note says original data; kept as estimated, to be checked')
add_qc('Freire: common name x taxon conflict (branco with Farfantepenaeus; sete barbas with Artemesia)',
       fz[(grepl('branco', tolower(common_orig)) & grepl('Farfantepenaeus', taxon)) | (grepl('sete barbas', tolower(common_orig)) & grepl('Artemesia', taxon)), .N],
       'taxon column used when species-level; genus-level with common name branco -> unassigned')
add_qc('MPA rows 2023-2024 (observed)', mz[status == 'observed', .N], 'status set by year, as informed by the author; confirm in Boletim 2023-2024')
add_qc('MPA rows 1950-2022 (estimated)', mz[status == 'estimated', .N], 'reconstruction; no observed/estimated column in the file')
add_qc('Records in scope but outside the 12 stock units', rec[in_scope & is.na(stock) & source != 'Silva & Andrade 2023', .N],
       paste('t:', round(rec[in_scope & is.na(stock) & source != 'Silva & Andrade 2023', sum(catch_t)]),
             '- e.g. pink in N/NE, brown in S/SE; excluded from stock units, kept in records'))
add_qc('Pink-shrimp complex records assigned by region (N/NE brown, S/SE pink)', rec[group_rule == 'regional_rosa', .N],
       paste('t:', round(rec[group_rule == 'regional_rosa', sum(catch_t)]), '- includes named F. subtilis/paulensis/brasiliensis and generic rosa, as in Freire 2021 / Silva & Andrade 2023'))
add_qc('Bulletin file: years without bulletin (all zeros) removed', length(setdiff(min(cam$ano):max(cam$ano), yr_ok)),
       paste(paste(setdiff(min(cam$ano):max(cam$ano), yr_ok), collapse = ','), '- zeros there mean no data, not zero catch'))
add_qc('Bulletin file: category start years', nrow(cat_start), paste(cat_start$cat, cat_start$first, collapse = '; '))
add_qc('Bulletin records with catch = 0 after the category start', cz[catch_t == 0, .N], 'CFG$zero_as_missing = TRUE: treated as missing (zero = not reported cannot be told apart from zero catch)')
add_qc('Records in generic category (Penaeidae, camarões, ...)', rec[group == 'generic', .N], paste('t:', round(rec[group == 'generic', sum(catch_t)])))
add_qc('Records unassigned', rec[group == 'unassigned', .N], 'e.g. common name branco with a pink taxon')
add_qc('Records out of scope: barba-ruça / santana', rec[group %in% c('barba_ruca', 'santana'), .N], 'kept in records, not in annual tables')
add_qc('Duplicated keys in Freire (state, year, name, taxon, sector, type)', fz[, sum(duplicated(.SD)), .SDcols = c('state', 'year', 'common_orig', 'taxon_orig', 'sector', 'CatchType')], 'should be 0')

# ---- 1.7.2 Overlap of MPA and Freire (both reconstructions 1950-2015) --------------------------------------------------
chk2 <- merge(unit[source == 'Freire et al. 2021', .(stock, year, freire = catch_total_t)],
              unit[source == 'MPA 2026', .(stock, year, mpa = catch_total_t)], by = c('stock', 'year'))
chk2_sum <- chk2[, .(n_years = .N, median_ratio_mpa_over_freire = round(median(mpa / freire, na.rm = TRUE), 2),
                     cor = round(cor(mpa, freire), 3)), by = stock]
# ---- 1.7.3 Freire 'reported' vs the bulletin file (same bulletins; do the totals agree?) ------------------------------
fr_rep <- rec[source == 'Freire et al. 2021' & status_detail %in% c('observed', 'observed_split', 'reported_modified') & tracked,
              .(freire_reported_t = sum(catch_t)), by = .(state, year)]
bl_tot <- rec[source == 'Silva & Andrade 2023 (bulletins)', .(bulletin_t = sum(catch_t)), by = .(state, year)]
chk3 <- merge(fr_rep, bl_tot, by = c('state', 'year'))
chk3[, period := cut(year, c(1961, 1977, 1989, 1999, 2007), labels = c('1962-77', '1978-89', '1990-99', '2000-07'))]
chk3_sum <- chk3[, .(n_state_years = .N, median_ratio_freire_over_bulletin = round(median(freire_reported_t / bulletin_t, na.rm = TRUE), 2),
                     share_within_10pct = round(mean(abs(freire_reported_t / bulletin_t - 1) < 0.10, na.rm = TRUE), 2)), by = period][order(period)]

# ---- 1.7.4 How much of the OBSERVED shrimp is generic (undifferentiated) ---------------------------------------------
gsh <- rec[source == 'Silva & Andrade 2023 (bulletins)', .(t = sum(catch_t)),
           by = .(region, period = cut(year, c(1961, 1977, 1989, 1999, 2007), labels = c('1962-77', '1978-89', '1990-99', '2000-07')),
                  generic = group == 'generic')]
gsh <- dcast(gsh, region + period ~ generic, value.var = 't', fill = 0)
setnames(gsh, c('FALSE', 'TRUE'), c('species_categories_t', 'generic_t'))
gsh[, generic_share := round(generic_t / (generic_t + species_categories_t), 2)]
qc <- rbindlist(qc)

# ------------------------------------------------------------------------------------------------------------------
# 1.8 WRITE THE FILES
# ------------------------------------------------------------------------------------------------------------------
# 1.8.2 STATE TOTALS: the unit of analysis of Steps 1-3
#     ONE row per source x state x year with the TOTAL shrimp catch of the state = sum of every tracked record
#     (generic + white + seabob + pink/brown; barba-ruca, santana and other non-listed penaeids are outside, in ALL sources).
#     source 'Observed'  = bulletin totals (1962-2007) + MPA collected landings (2023-2024). A state-year whose total is zero
#                          is NOT observed (cannot be told from 'no bulletin'), so it has no row.
#     the three reconstructions = their state totals, as published.
#       Freire et al. 2021 and MPA 2026: sum of the unified records above.
#       Silva & Andrade 2023: a MERGED series, built like the other two (reported values kept, model only where there is no data):
#         - state-years with an observed bulletin total (generic + species categories, the same total used for Freire and MPA):
#           the OBSERVED total;
#         - all other state-years (1950-2011): the sum of their three GLM predictions (c.rosa.recon + c.branco.recon + c.7barbas.recon).
#       Their GLM models the three species categories and uses the general shrimp category as a covariate (Silva & Andrade 2023); the
#       generic category is never added to the prediction. The full GLM prediction of ALL years (including the bulletin years) is kept as
#       a separate source, 'Silva & Andrade 2023 (model fit)', only to check how well the model reproduces the bulletins in Step 1.
#       (An earlier version of this script used the full model fit as the Silva & Andrade series; corrected on 2026-11-12, see the note
#       on deviations at the top.)
#     'unassigned' penaeid records (Freire: 'Camarao branco' with a Farfantepenaeus name, 9.6 thousand t, NE) ARE in the total: for a
#     total only 'is it a marine penaeid shrimp of the three bulletin groups?' matters. Barba-ruca and santana (S, SE) stay OUT in every
#     source: the bulletins do not report them and including them makes the S/SE reconstructions disagree MORE with the bulletins
#     (typical difference S: Freire 0.17 -> 0.34, MPA 0.11 -> 0.40; checked 1978-2007).
#     A state-year absent from a reconstruction inside the years that source covers is a ZERO catch (filled with 0, counted in qc).
obs_bul <- rec[source == CFG$src_obs & tracked, .(catch_t = sum(catch_t)), by = .(state, region, year)][catch_t > 0]
obs_bul[, origin := 'bulletin']
obs_mpa <- rec[source == 'MPA 2026' & year >= CFG$mpa_obs_from & tracked, .(catch_t = sum(catch_t)), by = .(state, region, year)][catch_t > 0]
obs_mpa[, origin := 'MPA collected landings']
st_obs <- rbind(obs_bul, obs_mpa)[, `:=`(source = 'Observed', status = 'observed')]
# A second observed series, bulletins only, WITHOUT the generic category (white + seabob + pink/brown). Why: the Silva & Andrade
# state files reconstruct only the three species categories (rosa, branco, sete-barbas); their reconstruction does not include the
# generic category (e.g. Ceara: generic = 99% of the observed total, S&A reconstruction ~30 t against ~2,800 t observed). Compared
# with this series the S&A reconstruction is judged on what it actually models. Used in Step 1 only.
obs_sp <- rec[source == CFG$src_obs & tracked & group != 'generic', .(catch_t = sum(catch_t)), by = .(state, region, year)][catch_t > 0]
obs_sp[, `:=`(origin = 'bulletin, species categories only', source = 'Observed species only', status = 'observed')]
st_obs <- rbind(st_obs, obs_sp)

st_rec <- rec[source %in% c('Freire et al. 2021', 'MPA 2026') & (tracked | group == 'unassigned') & !is.na(state),
              .(catch_t = sum(catch_t)), by = .(source, state, region, year)]
sa <- rbindlist(lapply(names(UF), function(ab) {
  f <- sprintf(CFG$f_silva_state, ab)
  if (!file.exists(f)) { warning('missing file ', f); return(NULL) }
  x <- fread(f, select = c('ano', 'c.rosa.recon', 'c.branco.recon', 'c.7barbas.recon'))
  x[, .(source = 'Silva & Andrade 2023', state = UF[[ab]], region = REGION[[UF[[ab]]]], year = ano,
        catch_t = c.rosa.recon + c.branco.recon + c.7barbas.recon)]
}))
sa <- sa[year >= 1950 & year <= CFG$silva_last]
sa_fit <- copy(sa)[, source := 'Silva & Andrade 2023 (model fit)']          # full GLM prediction, all years: diagnostic only
# merged series: observed bulletin total where it exists, GLM prediction where it does not
sa <- merge(sa, obs_bul[, .(state, year, obs_total_t = catch_t)], by = c('state', 'year'), all.x = TRUE)
n_sa_obs <- sa[!is.na(obs_total_t), .N]; n_sa_model <- sa[is.na(obs_total_t), .N]
sa[, catch_t := fifelse(!is.na(obs_total_t), obs_total_t, catch_t)][, obs_total_t := NULL]
st_rec <- rbind(st_rec, sa, sa_fit, use.names = TRUE)
# fill zeros inside the years each source covers
rng <- st_rec[, .(y0 = min(year), y1 = max(year)), by = source]
grid <- rng[, CJ(state = names(REGION), year = y0:y1), by = source]
st_rec <- merge(grid, st_rec[, .(source, state, year, catch_t)], by = c('source', 'state', 'year'), all.x = TRUE)
n_filled <- st_rec[is.na(catch_t), .N, by = source]
st_rec[is.na(catch_t), catch_t := 0]
st_rec[, `:=`(region = REGION[state], status = 'estimated', origin = 'reconstruction')]
st_all <- rbind(st_obs, st_rec, fill = TRUE)[, .(source, state, region, year, catch_t, status, origin)]
setorder(st_all, source, state, year)
fwrite(st_all, 'data/catch_state_total.csv')
fwrite(rec, 'data/catch_records.csv')          # every shrimp record (also in the workbook), read by Step 1 when section 1 is not rerun
qc <- rbind(qc, data.table(item = c('State-years filled with zero (absent from a reconstruction inside its own years)', 'Silva & Andrade state files read',
                                    'Silva & Andrade merged series: state-years with the observed total', 'Silva & Andrade merged series: state-years with the GLM prediction'),
                           n = c(sum(n_filled$N), sa[, uniqueN(state)], n_sa_obs, n_sa_model),
                           detail = c(paste(n_filled$source, n_filled$N, collapse = '; '), paste0('data_XX.csv, 1950-', CFG$silva_last), 'bulletin total kept', 'c.rosa.recon + c.branco.recon + c.7barbas.recon')), fill = TRUE)
cat('State totals written. Observed state-years:', nrow(st_obs), '\n')


# 1.8.3 The organised workbook
readme <- data.table(sheet = c('records', 'annual_state_group', 'annual_stock_unit', 'coverage_observed', 'name_map', 'excluded',
                               'qc_log', 'overlap_mpa_freire', 'freire_vs_bulletins', 'generic_share', 'state_total'),
  description = c('every shrimp record: source, status (observed/estimated), state, year, names (original and standard), group, stock, catch, reference, note',
                  'Freire and MPA summed to state x group x year; status = observed only if every record of the cell is observed',
                  'stock unit (12 stocks) x year by source; catch_observed_t and share of states observed',
                  'for each source, state and group: how many years are observed, first/last, longest run, observed years in 2008-2015, 2016-2022, 2023-2024',
                  'how each original common name/taxon was standardised and which group it received',
                  'names removed (non-penaeid, deep-water, freshwater or farmed shrimps)',
                  'checks and counts that need your eyes',
                  'MPA vs Freire annual totals per stock in the years both cover',
                  'Freire reported totals vs the bulletin file, per state-year, by period',
                  'share of observed shrimp in the generic (undifferentiated) category by region and period',
                  'TOTAL shrimp per source x state x year: Observed (bulletins + MPA 2023-24) and the three reconstructions (unit of Steps 1-3)'))
wb <- createWorkbook()
put <- function(nm, x) { addWorksheet(wb, nm); writeData(wb, nm, x); freezePane(wb, nm, firstRow = TRUE); setColWidths(wb, nm, cols = 1:ncol(x), widths = 'auto') }
put('README', readme); put('records', rec); put('annual_state_group', ann_state); put('annual_stock_unit', unit)
put('coverage_observed', cover); put('name_map', name_map); put('excluded', excluded); put('qc_log', qc)
put('overlap_mpa_freire', chk2_sum); put('freire_vs_bulletins', chk3_sum); put('generic_share', gsh); put('state_total', st_all)
saveWorkbook(wb, 'data/catch_unified_base.xlsx', overwrite = TRUE)

cat('
Data ready. Records:', nrow(rec), '| state-group-year cells:', nrow(ann_state), '| state-total rows:', nrow(st_all), '
')
print(qc[, .(item, n)])
print(chk2_sum); print(chk3_sum); print(gsh)

say('SECTION 1 done: data/catch_state_total.csv and data/catch_unified_base.xlsx written')
}  # end of section 1 (RUN$build_data)

# ------------------------------------------------------------------------------------------------------------------
# 2. STEP 1 - COMPARE THE PUBLISHED RECONSTRUCTIONS (state totals)
# ------------------------------------------------------------------------------------------------------------------
# Nothing is refitted here. Two plain questions about the reconstructions that already exist, at the level of the
# TOTAL shrimp catch of each state (and of the regions and Brazil, obtained by summing states):
#   QUESTION 1 (main)   Where the three reconstructions HAD TO FILL A GAP, do Freire et al. 2021, Silva & Andrade 2023 and MPA 2026
#                       agree with EACH OTHER?
#   QUESTION 2 (check)  Where a bulletin exists, do they reproduce it?
# There is no hypothesis in this step.
#
# WHY ONLY THE FILLED YEARS ARE COMPARED. In the years with a bulletin the three reconstructions are anchored on the reported values
# (Silva & Andrade 2023 kept the observed value and used the GLM only where there was no data; Freire et al. 2021 and MPA 2026
# reproduce the bulletin total exactly in about 55% of the state-years and within 10% in about 78%, see 'reproduces_observed' in
# Table 1), so in those years they agree with each other and with the bulletins largely BY CONSTRUCTION, and the agreement says
# nothing about the methods. The methods differ only where there is no observed value. A FILLED YEAR is therefore a state-year
# with no observed value (no bulletin and, from 2023, no MPA collected landing): 1950-1961, 1965, 2008-2022 and the state-years
# inside the bulletin era in which a state did not report. The same set of filled state-years is used for all three
# reconstructions (it is defined by the OBSERVED data, not by a reconstruction).
#
# The Silva & Andrade series used here is the one built in section 1.8.2: the observed state total (with the generic category) in
# years with a bulletin, and the sum of their three GLM predictions (c.rosa.recon + c.branco.recon + c.7barbas.recon) in the other
# years (1950-2011). Their full model fit (all years) is kept as 'Silva & Andrade 2023 (model fit)' ONLY to show how well the model
# reproduces the bulletins (QUESTION 2); it is not one of the three reconstructions compared.
#
# The three metrics (a = first series of the pair, b = second):
#   ratio              total of a / total of b over the common filled years (1 = same total)
#   typical_difference exp(mean(|ln(a+1) - ln(b+1)|)) - 1: typical yearly difference as a fraction (0.25 = about 25%)
#   trend_agreement    correlation of the year-to-year changes of ln(a+1) and ln(b+1), CONSECUTIVE filled years only
#                      (+1 = up and down together, 0 = unrelated); NA with fewer than 5 pairs of consecutive years
#                      (a state whose filled years are scattered, e.g. only 1965 and 2008-2011, can have no trend value)
# For regions and Brazil the value is the SUM of the reconstruction over the states of the unit, and it exists only in years in which EVERY
# state of the unit is filled (a year in which part of the region was observed would mix observed and reconstructed states and make the
# reconstructions agree artificially). For Brazil these are 1950-1961, 1965 and 2008-2022.
if (RUN$step1) {
say('STEP 1: comparing the published reconstructions in the years they had to fill')

st <- fread('data/catch_state_total.csv')
if (!exists('rec')) rec <- fread('data/catch_records.csv')

# 2.1 Wide tables. SA_FIT is the full Silva & Andrade model fit (diagnostic only).
SA_FIT <- 'Silva & Andrade 2023 (model fit)'
SRC <- c('Observed', 'Observed species only', CFG$src_recon, SA_FIT)
wide_state <- dcast(st[year <= CFG$mpa_recon_last], state + region + year ~ source, value.var = 'catch_t')
wide_state <- wide_state[, c('state', 'region', 'year', SRC), with = FALSE]
# a year with an observed total but no species record is a year with ZERO identified species (only the generic category was reported)
wide_state[!is.na(Observed) & is.na(`Observed species only`), `Observed species only` := 0]

# 2.1a Filled years: the reconstructions are blanked where an observed value exists
wide_fill <- copy(wide_state)
wide_fill[!is.na(Observed), (CFG$src_recon) := NA_real_]
sum_gap <- function(x) if (any(is.na(x))) NA_real_ else sum(x)     # blanked (= observed) or uncovered state-years are NA, so the sum exists only if every state is filled
FCOLS <- c('unit', 'unit_type', 'year', CFG$src_recon)
wide_fill_unit <- rbind(
  copy(wide_fill)[, `:=`(unit = state, unit_type = 'state')][, ..FCOLS],
  wide_fill[, lapply(.SD, sum_gap), by = .(unit = region, year), .SDcols = CFG$src_recon][, unit_type := 'region'][, ..FCOLS],
  wide_fill[, lapply(.SD, sum_gap), by = .(year), .SDcols = CFG$src_recon][, `:=`(unit = 'Brazil', unit_type = 'Brazil')][, ..FCOLS])
setorder(wide_fill_unit, unit_type, unit, year)

# 2.1b Observed-era table (state, region, Brazil): used only for the check against the bulletins. Region and Brazil values are SUMS OF STATES
#      and exist only in years where every state of the unit has a value (a partial sum would be a different quantity).
wide_region <- wide_state[, lapply(.SD, sum_na), by = .(unit = region, year), .SDcols = SRC]
wide_brazil <- wide_state[, lapply(.SD, sum_na), by = .(year), .SDcols = SRC][, unit := 'Brazil']
UCOLS <- c('unit', 'unit_type', 'year', SRC)
wide_unit <- rbind(copy(wide_state)[, `:=`(unit = state, unit_type = 'state')][, ..UCOLS],
                   wide_region[, unit_type := 'region'][, ..UCOLS],
                   wide_brazil[, unit_type := 'Brazil'][, ..UCOLS])
setorder(wide_unit, unit_type, unit, year)

# 2.2 QUESTION 1: the three pairs of reconstructions in the filled years (all, before 1962, 2008 onwards)
WINDOWS_F <- list('filled years, all' = c(1950, CFG$mpa_recon_last),
                  'filled years before 1962' = c(1950, CFG$year_first - 1),
                  'filled years 2008 onwards' = c(CFG$bulletin_last + 1, CFG$mpa_recon_last))
PAIRS <- combn(CFG$src_recon, 2, simplify = FALSE)      # Freire-Silva, Freire-MPA, Silva-MPA
agree <- rbindlist(lapply(PAIRS, function(p) rbindlist(lapply(names(WINDOWS_F), function(wn) {
  w <- wide_fill_unit[year >= WINDOWS_F[[wn]][1] & year <= WINDOWS_F[[wn]][2]]
  w <- w[!is.na(get(p[1])) & !is.na(get(p[2]))]
  w[, c(list(pair = paste(p[1], 'vs', p[2]), window = wn), agree_fun(year, get(p[1]), get(p[2]))), by = .(unit_type, unit)]
}))))
agree[, `:=`(ratio = round(ratio, 2), typical_difference = round(typical_difference, 2), trend_agreement = round(trend_agreement, 2))]

# 2.3 Summary: one line per pair and window, median (minimum to maximum) ACROSS THE 17 STATES; Brazil printed apart.
#     n_states_trend = states with enough consecutive filled years to compute the trend agreement.
agree_summary <- agree[unit_type == 'state', .(n_states = .N, n_states_trend = sum(!is.na(trend_agreement)), ratio = fmt3(ratio),
                                              typical_difference = fmt3(typical_difference), trend_agreement = fmt3(trend_agreement)), by = .(pair, window)]
agree_brazil <- agree[unit_type == 'Brazil', .(pair, window, n_years, ratio, typical_difference, trend_agreement)]
setorder(agree_summary, window, pair); setorder(agree_brazil, window, pair)
cat('\nStep 1 - agreement in the FILLED years, median (min to max) across the 17 states:\n'); print(agree_summary)
cat('\nStep 1 - agreement in the filled years for Brazil (sum of the 17 states, years when all are filled):\n'); print(agree_brazil)
# How many state-years are actually filled (the comparison rests on these)
# (zero_* = filled state-years in which that reconstruction has NO RECORD for the state and the catch was set to zero (section 1.8). The original
#  spreadsheets were checked: no shrimp record exists under any scientific or common name; Freire et al. start the series of a state with its first
#  record, e.g. 1962 in AL, BA, PB, PE, RN and SE. A zero makes ln(catch + 1) = 0 and inflates the typical difference.)
filled_by_state <- wide_state[, .(years_filled_1950_2011 = sum(is.na(Observed) & year <= CFG$silva_last),
                                  years_filled_1950_2022 = sum(is.na(Observed)),
                                  zero_Freire = sum(is.na(Observed) & `Freire et al. 2021` == 0, na.rm = TRUE),
                                  zero_Silva_Andrade = sum(is.na(Observed) & year <= CFG$silva_last & `Silva & Andrade 2023` == 0, na.rm = TRUE),
                                  zero_MPA = sum(is.na(Observed) & `MPA 2026` == 0, na.rm = TRUE)), by = .(state, region)]
setorder(filled_by_state, region, state)
cat('\nFilled state-years per state:\n'); print(filled_by_state)

# 2.3b QUESTION 2 (check): do the reconstructions reproduce the bulletins (1962-2007, state-years with an observed total)?
#      Also the diagnostic of the Silva & Andrade model: what share of the observed total do the three species categories make up, and
#      how large is the model fit relative to the observed total and relative to the species categories alone?
CHK_SRC  <- c(CFG$src_recon, SA_FIT)
PAIRS_OBS <- c(lapply(CHK_SRC, function(s) c(s, 'Observed')), lapply(CHK_SRC, function(s) c(s, 'Observed species only')))
agree_obs <- rbindlist(lapply(PAIRS_OBS, function(p) {
  w <- wide_unit[year >= CFG$year_first & year <= CFG$bulletin_last]
  w <- w[!is.na(get(p[1])) & !is.na(get(p[2]))]
  w[, c(list(pair = paste(p[1], 'vs', p[2]), window = 'bulletin era 1962-2007'), agree_fun(year, get(p[1]), get(p[2]))), by = .(unit_type, unit)]
}))
agree_obs[, `:=`(ratio = round(ratio, 2), typical_difference = round(typical_difference, 2), trend_agreement = round(trend_agreement, 2))]
agree_obs_summary <- agree_obs[unit_type == 'state', .(n_states = .N, ratio = fmt3(ratio), typical_difference = fmt3(typical_difference),
                                                      trend_agreement = fmt3(trend_agreement)), by = .(pair, window)]
cat('\nStep 1 - check against the bulletins (1962-2007), median (min to max) across states:\n'); print(agree_obs_summary)
bul <- wide_state[year >= CFG$year_first & year <= CFG$bulletin_last & !is.na(Observed)]
repro <- rbindlist(lapply(CHK_SRC, function(s) bul[, .(source = s, n_state_years = .N,
          exactly_equal = round(mean(abs(get(s) - Observed) < 0.5), 2), within_10pct = round(mean(abs(get(s) / Observed - 1) <= 0.1), 2))]))
cat('\nShare of the observed state-years that each series reproduces (|difference| < 0.5 t; within 10%):\n'); print(repro)
fit_vs_obs <- bul[, .(n_years = .N, observed_total_t = sum(Observed), species_categories_t = sum(`Observed species only`),
                      model_fit_t = sum(get(SA_FIT))), by = .(state, region)]
fit_vs_obs[, `:=`(species_share = species_categories_t / observed_total_t, model_fit_over_total = model_fit_t / observed_total_t,
                  model_fit_over_species = model_fit_t / species_categories_t)]
setorder(fit_vs_obs, region, state)
cat('\nShare of the observed total made of the three species categories, and Silva & Andrade model fit / observed total (bulletin years):\n')
print(fit_vs_obs[, .(state, region, n_years, species_share = round(species_share, 2), model_fit_over_total = round(model_fit_over_total, 2),
                     model_fit_over_species = round(model_fit_over_species, 2))])
cat(sprintf('Median across states: species share %.2f; model fit / observed total %.2f. All states together: species share %.2f.\n',
            median(fit_vs_obs$species_share), median(fit_vs_obs$model_fit_over_total), sum(fit_vs_obs$species_categories_t) / sum(fit_vs_obs$observed_total_t)))

# 2.4 FIGURES
state_levels <- c(names(REGION), 'Brazil')
wide_unit[, unit_f := factor(unit, levels = state_levels)]
long_rc <- melt(wide_unit[unit %in% state_levels & !is.na(unit_f)], id.vars = c('unit', 'unit_f', 'unit_type', 'year'),
                measure.vars = CFG$src_recon, variable.name = 'source', value.name = 'catch_t')
long_rc <- long_rc[!is.na(catch_t)]
# the full range of each reconstruction is drawn (years after the MPA reconstruction included), states in all of 1950-2015
rc_full <- rbind(
  st[source %in% CFG$src_recon, .(unit = state, year, source, catch_t)],
  st[source %in% CFG$src_recon, .(catch_t = sum_na(catch_t), n = .N), by = .(source, year)][n == 17][, .(unit = 'Brazil', year, source, catch_t)])
rc_full[, unit_f := factor(unit, levels = state_levels)]
obs_pts <- wide_unit[unit %in% state_levels & !is.na(Observed), .(unit_f = factor(unit, levels = state_levels), year, observed_t = Observed)]
obs_pts <- rbind(obs_pts, st[source == 'Observed' & year > CFG$mpa_recon_last, .(unit_f = factor(state, levels = state_levels), year, observed_t = catch_t)])

# Fig 1: the three reconstructions (lines) and the observed total (points), one panel per state + Brazil
p1 <- ggplot(rc_full, aes(year, catch_t, colour = source)) +
  geom_line(linewidth = 0.5, alpha = 0.95) +
  geom_point(data = obs_pts, aes(year, observed_t), inherit.aes = FALSE, colour = obs_col, size = 0.8) +
  facet_wrap(~ unit_f, scales = 'free_y', ncol = 6) +
  scale_colour_manual(values = source_cols, name = 'Reconstruction') +
  scale_x_continuous(breaks = seq(1950, 2020, 30)) +
  scale_y_continuous(labels = thousand, limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = 'Total shrimp catch (thousand t)',
       caption = paste('Lines: published reconstructions, total shrimp per state (Silva & Andrade: observed bulletin total where reported, GLM prediction elsewhere). Points: observed total (bulletins 1962-2007, MPA collected landings 2023-2024).',
                       'Brazil = sum of the 17 states; its observed total is shown only in years when all 17 states were observed.', sep = '\n')) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2))) + theme_cru()
save_fig(p1, 'Fig_1_state_totals', width = 15, height = 8.5)

# Fig 2: what the observed total is made of, by region (stacked bars) + the reconstructed totals (lines)
comp_levels <- c('generic (not identified)', 'seabob', 'white shrimp', 'pink/brown shrimp ("rosa")')
comp_fills  <- c('generic (not identified)' = '#D5D4CF', 'seabob' = '#F1D79B', 'white shrimp' = '#CBBEE0', 'pink/brown shrimp ("rosa")' = '#F0B5C0')
obs_comp <- rec[source == CFG$src_obs & !is.na(region) & !is.na(year) & catch_t > 0, .(catch_t = sum(catch_t)), by = .(region, year, group)]
obs_comp[, category := factor(fcase(group == 'generic', comp_levels[1], group == 'seabob', comp_levels[2],
                                    group == 'white', comp_levels[3], group %in% c('brown', 'pink'), comp_levels[4]), levels = comp_levels)]
obs_comp <- obs_comp[!is.na(category), .(catch_t = sum(catch_t)), by = .(region, year, category)]
obs_comp[, region_f := factor(region, levels = UNITS_REGION, labels = REG_LAB)]
tot_reg <- st[source %in% CFG$src_recon & year <= CFG$mpa_recon_last, .(catch_t = sum(catch_t)), by = .(source, region, year)]
tot_reg[, region_f := factor(region, levels = UNITS_REGION, labels = REG_LAB)]
generic_share <- obs_comp[year <= CFG$bulletin_last, .(generic_share = round(sum(catch_t[category == comp_levels[1]]) / sum(catch_t), 2)), by = region_f]
p2 <- ggplot() +
  geom_col(data = obs_comp, aes(year, catch_t, fill = category), width = 0.9) +
  geom_line(data = tot_reg, aes(year, catch_t, colour = source), linewidth = 0.6, alpha = 0.95) +
  facet_wrap(~ region_f, scales = 'free_y', ncol = 2) +
  scale_fill_manual(values = comp_fills, name = 'Observed bulletin catch') +
  scale_colour_manual(values = source_cols, name = 'Reconstruction (total)') +
  scale_x_continuous(breaks = seq(1950, 2020, 10)) +
  scale_y_continuous(labels = thousand, limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = 'Total shrimp catch (thousand t)',
       caption = paste('Bars: observed catch by bulletin category (bar height = observed total of the states that reported; years without a bulletin are empty).',
                       'Lines: totals of the published reconstructions (sum of the states of the region; Silva & Andrade: observed bulletin total where reported, GLM prediction elsewhere).', sep = '\n')) +
  guides(colour = guide_legend(order = 1, override.aes = list(linewidth = 1.2)), fill = guide_legend(order = 2)) +
  theme_cru() + theme(legend.box = 'vertical', legend.spacing.y = unit(0.1, 'lines'))
save_fig(p2, 'Fig_2_total_shrimp_by_region', width = 12, height = 7.5)
cat('\nShare of the generic category in the observed total, bulletin era, by region:\n'); print(generic_share)

# Fig 3: the reconstructions compared with each other in the FILLED years only. Each point is one filled state-year:
#   d = ln(reconstruction + 1) - ln(MPA 2026 + 1)   (reference = MPA 2026, the series that covers most years; the reference is
#   a yardstick, NOT the truth: no reconstruction is validated against the other). d = 0 -> same value; +0.69 -> twice MPA (+100%);
#   -0.69 -> half (-50%). Axis shown as a percentage.  Each panel prints, per reconstruction, the typical difference and trend agreement r.
d_long <- rbindlist(lapply(c('Freire et al. 2021', 'Silva & Andrade 2023'), function(s) {
  w <- wide_fill_unit[unit %in% state_levels & !is.na(get(s)) & !is.na(`MPA 2026`)]
  w[, .(unit_f = factor(unit, levels = state_levels), year, source = s, d = log(get(s) + 1) - log(`MPA 2026` + 1), rc = get(s), ref = `MPA 2026`)]
}))
setorder(d_long, source, unit_f, year)
d_long[, run := cumsum(c(1, diff(year) != 1)), by = .(source, unit_f)]          # lines break where years are missing
d_lab <- agree[window == 'filled years, all' & pair %in% paste(c('Freire et al. 2021', 'Silva & Andrade 2023'), 'vs MPA 2026') & unit %in% state_levels,
               .(unit_f = factor(unit, levels = state_levels), source = sub(' vs MPA 2026', '', pair), typical_difference, trend_agreement)]
d_lab[, `:=`(line = match(source, c('Freire et al. 2021', 'Silva & Andrade 2023')),
             text = sprintf('%s: %s%%, r = %s', c('Freire', 'Silva & Andrade')[match(source, c('Freire et al. 2021', 'Silva & Andrade 2023'))],
                            round(100 * typical_difference), ifelse(is.na(trend_agreement), 'n/a', sprintf('%.2f', trend_agreement))))]
pct_axis   <- function(b) ifelse(abs(b) < 1e-9, '0%', sprintf('%+.0f%%', (exp(b) - 1) * 100))
pct_breaks <- function(lim) {            # percentage ticks inside the panel's range, at most 5 so the labels never overlap
  b <- log(c(0.02, 0.1, 0.25, 0.5, 1, 2, 4, 10, 25, 100)); b <- b[b >= lim[1] & b <= lim[2]]
  if (length(b) > 5) b[unique(round(seq(1, length(b), length.out = 5)))] else b
}
p3 <- ggplot(d_long, aes(year, d, colour = source, group = interaction(source, run))) +
  geom_hline(yintercept = 0, colour = ink, linewidth = 0.4) +
  geom_line(linewidth = 0.45, alpha = 0.9) + geom_point(size = 0.7) +
  geom_text(data = d_lab, aes(x = -Inf, y = Inf, label = text, vjust = 1.5 + 1.6 * (line - 1)), inherit.aes = FALSE,
            colour = source_cols[d_lab$source], hjust = -0.03, size = 2.3) +
  facet_wrap(~ unit_f, scales = 'free_y', ncol = 4) +
  scale_colour_manual(values = source_cols, name = 'Reconstruction') +
  scale_x_continuous(breaks = seq(1950, 2020, 20)) +
  scale_y_continuous(breaks = pct_breaks, labels = pct_axis, expand = expansion(mult = c(0.05, 0.45))) +
  labs(x = NULL, y = 'Reconstruction relative to MPA 2026',
       caption = paste('Each point: ln(reconstruction + 1) - ln(MPA 2026 + 1) for one FILLED state-year (no observed value), shown as a percentage (0% = equal to MPA 2026).',
                       'In each panel: typical difference and trend agreement r (correlation of year-to-year changes) over the filled years. Brazil: sum of the 17 states, only in years when all are filled.', sep = '\n')) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2))) + theme_cru()
save_fig(p3, 'Fig_3_difference_filled_years', width = 15, height = 13)

# Fig 9: what the observed bulletin total is made of, by state: share of the three identified species categories (pink/brown, white,
# seabob) and share of the generic category (accumulated over the years with a reported total, 1962-2007).
sp_levels <- fit_vs_obs[order(match(region, c('N', 'NE', 'SE', 'S')), species_share)]$state       # lowest share at the bottom of each panel
sp_plot <- rbind(fit_vs_obs[, .(state, region, part = 'species identified (pink/brown, white, seabob)', pct = 100 * species_share)],
                 fit_vs_obs[, .(state, region, part = 'generic (not identified)', pct = 100 * (1 - species_share))])
sp_plot[, `:=`(state_f = factor(state, levels = sp_levels), region_f = factor(region, levels = UNITS_REGION, labels = REG_LAB),
               part = factor(part, levels = c('generic (not identified)', 'species identified (pink/brown, white, seabob)')))]
sp_lab <- fit_vs_obs[, .(state_f = factor(state, levels = sp_levels), region_f = factor(region, levels = UNITS_REGION, labels = REG_LAB),
                         pct = 100 * species_share, label = sprintf('%.0f%%', 100 * species_share))]
p9 <- ggplot(sp_plot, aes(pct, state_f, fill = part)) +
  geom_col(width = 0.7) +
  geom_text(data = sp_lab, aes(x = pct + 1.5, y = state_f, label = label), inherit.aes = FALSE, hjust = 0, size = 3, colour = ink) +
  facet_grid(region_f ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_manual(values = c('generic (not identified)' = '#D5D4CF', 'species identified (pink/brown, white, seabob)' = '#3D5A80'), name = NULL) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = 'Share of the observed shrimp total (%)', y = NULL,
       caption = 'Accumulated bulletin catch of each state over the years with a reported total (1962-2007). Bars sorted within region.') +
  theme_cru() + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.3))
save_fig(p9, 'Fig_9_species_share_by_state', width = 8.5, height = 9)

write.xlsx(list(filled_years_summary = agree_summary, filled_years_brazil = agree_brazil, filled_years_by_state = filled_by_state,
                agreement_all_units = agree, check_bulletins_summary = agree_obs_summary, check_bulletins_all_units = agree_obs,
                reproduces_observed = repro),
           'results/Table_1_agreement_state_totals.xlsx')
write.xlsx(list(log_difference_filled_years = d_long[, .(source, unit = unit_f, year, reconstruction_t = rc, mpa_t = ref, log_difference = round(d, 4))],
                species_share_and_model_fit = fit_vs_obs[, .(state, region, n_years, observed_total_t = round(observed_total_t), species_categories_t = round(species_categories_t),
                                                           model_fit_t = round(model_fit_t), species_share = round(species_share, 3),
                                                           model_fit_over_total = round(model_fit_over_total, 3), model_fit_over_species = round(model_fit_over_species, 3))],
                generic_share_by_region_period = generic_share),
           'results/Table_2_difference_vs_observed_by_year.xlsx')
cat('\nStep 1 done. Tables and figures are in results/.\n')

}  # end of section 2 (RUN$step1)



if (RUN$step2 || RUN$step3) {

# ------------------------------------------------------------------------------------------------------------------
# 3. LOAD THE OBSERVED STATE SERIES USED IN STEPS 2-3
# ------------------------------------------------------------------------------------------------------------------
# One series per state: years year_first..year_last, catch = OBSERVED total (NA = not observed, i.e. to be filled).
st <- fread('data/catch_state_total.csv')
yrs_all <- CFG$year_first:CFG$year_last
catch <- merge(CJ(state = names(REGION), year = yrs_all), st[source == 'Observed', .(state, year, catch = catch_t)],
               by = c('state', 'year'), all.x = TRUE)
catch[, region := REGION[state]]
longest_run <- function(y) { y <- sort(unique(y)); if (!length(y)) return(0L); brk <- c(0, which(diff(y) != 1), length(y)); max(diff(brk)) }
coverage <- catch[, {
  o <- year[!is.na(catch)]
  .(region = region[1], n_obs = length(o), first_obs = min(o), last_obs = max(o), longest_run = longest_run(o),
    n_obs_bulletin_era = sum(o <= CFG$bulletin_last), n_obs_2023_24 = sum(o >= CFG$mpa_obs_from),
    usable_step2 = length(o) >= CFG$min_obs_hide)
}, by = state]
setorder(coverage, region, state)
fwrite(coverage, 'data/state_series_coverage.csv')
cat('\nObserved years per state series:\n'); print(coverage)

# ------------------------------------------------------------------------------------------------------------------
# 4. THE FOUR RECONSTRUCTION METHODS (they use ONLY the series' own observed values and its time index)
# ------------------------------------------------------------------------------------------------------------------
# Each function takes a data.frame with columns year and catch (NA = to be filled) and returns a numeric vector with a
# prediction for EVERY year (the caller decides which years to use). Fitting is on ln(C + 1) so that zeros are allowed.
# Every method caps its prediction at 1e6 t (written as a literal inside each function, because each function runs in a
# fresh R process): it only prevents numeric overflow if a GLM extrapolates absurdly and is far above any Brazilian state catch.

# ---- 4.0 M0: naive. Each year receives the value of the nearest observed year (the earlier one if two are equally near). ----
# Why: a method that cannot beat 'repeat the nearest known value' adds nothing. It is the reference of the whole test.
fit_M0 <- function(d, ...) {
  obs <- !is.na(d$catch); yo <- d$year[obs]; co <- d$catch[obs]
  vapply(d$year, function(y) co[which.min(abs(yo - y))], numeric(1))
}

# ---- 4.1 M1: linear interpolation inside gaps, LOESS extrapolation at the edges (Freire et al. 2021 logic) ---------
fit_M1 <- function(d, ...) {
  y   <- log(d$catch + 1); obs <- !is.na(y)
  out <- rep(NA_real_, nrow(d))
  inside <- d$year >= min(d$year[obs]) & d$year <= max(d$year[obs])
  out[inside] <- approx(d$year[obs], y[obs], xout = d$year[inside])$y     # interior: straight line between neighbours
  dd  <- data.frame(year = d$year[obs], y = y[obs])
  lo  <- loess(y ~ year, data = dd, span = 0.75, degree = 2,              # edges: smooth curve (span 0.75 = 75% of the points) ...
               control = loess.control(surface = 'direct'))               # ... 'direct' allows extrapolation outside the data range
  out[!inside] <- predict(lo, newdata = data.frame(year = d$year[!inside]))
  pmin(pmax(exp(out) - 1, 0), 1e6)
}

# ---- 4.2 M2: gamma GLM, as in Silva & Andrade 2023 (without their environmental covariates) -----------------------
# Silva & Andrade (2023, Bol. Inst. Pesca 49:e747) fitted generalized linear models with a gamma distribution and log link.
# Time entered as a quadratic term (year^2) and a sinusoidal term (sin of year); terms were chosen by AIC, starting from the
# intercept-only model and adding one term at a time while the AIC decreased. Their other candidate terms (general shrimp
# category, sea-surface temperature, rainfall) are NOT used: all methods see only the series' own years.
#   - response: catch + 1 (a gamma needs values > 0; the same +1 used by M1 and M3 on the log scale)
#   - year^2 is scaled by 1e-6 only for numerical stability (this does not change the fit)
#   - sin(year): year in radians, as named in their Table 1 (exact definition to be confirmed with the authors)
fit_M2 <- function(d, ...) {
  obs <- !is.na(d$catch)
  dd  <- data.frame(y = d$catch[obs] + 1, year2 = d$year[obs]^2 / 1e6, sinyear = sin(d$year[obs]))
  nd  <- data.frame(year2 = d$year^2 / 1e6, sinyear = sin(d$year))
  fam <- Gamma(link = 'log')
  fit_terms <- function(tt) try(glm(reformulate(if (length(tt)) tt else '1', response = 'y'), family = fam, data = dd,
                                   control = glm.control(maxit = 100)), silent = TRUE)
  chosen <- character(0); cur <- fit_terms(chosen)
  if (inherits(cur, 'try-error')) return(rep(mean(d$catch[obs]), nrow(d)))
  repeat {                                                                             # forward selection by AIC
    cand <- setdiff(c('year2', 'sinyear'), chosen)
    if (!length(cand)) break
    fits <- lapply(cand, function(v) fit_terms(c(chosen, v)))
    ok   <- !vapply(fits, inherits, logical(1), what = 'try-error')
    if (!any(ok)) break
    aic  <- vapply(fits[ok], AIC, numeric(1))
    if (min(aic) < AIC(cur) - 1e-8) { chosen <- c(chosen, cand[ok][which.min(aic)]); cur <- fits[ok][[which.min(aic)]] } else break
  }
  pmin(pmax(as.numeric(predict(cur, newdata = nd, type = 'response')) - 1, 0), 1e6)
}

# ---- 4.3 M3: small feedforward neural network, year as the only input (MPA approach) ------------------------------
# keras is loaded INSIDE the function because each call runs in a fresh R process (see 5.2).
#   3 hidden layers x 32 tanh units, dropout 0.2, Adam, mean squared error on the standardised ln(C+1);
#   early stopping (patience 30) on a RANDOM 15% of the observed years. The prediction is the mean of `reps` networks.
fit_M3 <- function(d, reps = 2, seed = 1, ...) {
  # Silence the TensorFlow/keras chatter (C++ info lines, optimizer and function-call tips, Python deprecation warnings).
  # The environment variables must be set BEFORE TensorFlow is loaded, so they come first; errors are still shown.
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = '3',      # 3 = hide TensorFlow's C++ INFO / WARNING / ERROR lines
             TF_ENABLE_ONEDNN_OPTS = '0',     # hide the oneDNN "custom operations" notice
             PYTHONWARNINGS = 'ignore')       # hide Python warnings (input_shape tip, optimizer tips, etc.)
  suppressPackageStartupMessages(library(keras))
  try({ reticulate::py_run_string("import warnings, logging; warnings.filterwarnings('ignore'); logging.getLogger('tensorflow').setLevel(logging.ERROR)")
        tensorflow::tf$get_logger()$setLevel('ERROR') }, silent = TRUE)
  obs <- !is.na(d$catch)
  x   <- (d$year - mean(d$year[obs])) / sd(d$year[obs])                   # standardised year
  y   <- log(d$catch + 1)
  mu  <- mean(y[obs]); s <- sd(y[obs]); if (!is.finite(s) || s == 0) s <- 1
  ys  <- (y - mu) / s
  preds <- matrix(NA_real_, nrow(d), reps)
  for (r in seq_len(reps)) {
    set.seed(seed + r); tensorflow::set_random_seed(seed + r)
    m <- suppressWarnings(suppressMessages(keras_model_sequential() |>
      layer_dense(32, activation = 'tanh', input_shape = 1) |> layer_dropout(0.2) |>
      layer_dense(32, activation = 'tanh') |> layer_dropout(0.2) |>
      layer_dense(32, activation = 'tanh') |> layer_dropout(0.2) |>
      layer_dense(1)))
    suppressWarnings(suppressMessages(m |> compile(optimizer = 'adam', loss = 'mse')))
    # Early-stopping hold-out: a RANDOM 15% of the observed years. (keras' validation_split would take the LAST 15% of the
    # rows, i.e. the most recent years, and never train on them: that would leave the end of every series out of the fit
    # and bias the edge-error test.) The held-out years are used only to decide when to stop.
    xo <- x[obs]; yo <- ys[obs]
    iv <- sample(length(xo), max(2, round(0.15 * length(xo))))
    suppressWarnings(suppressMessages(
      m |> fit(xo[-iv], yo[-iv], epochs = 500, batch_size = 16, validation_data = list(xo[iv], yo[iv]), verbose = 0,
               callbacks = list(callback_early_stopping(patience = 30, restore_best_weights = TRUE)))))
    preds[, r] <- suppressWarnings(as.numeric(predict(m, x, verbose = 0))) * s + mu
    keras::k_clear_session(); rm(m); gc()                                 # release memory after EVERY network
  }
  pmin(pmax(exp(rowMeans(preds)) - 1, 0), 1e6)
}

METHODS <- list(M0 = fit_M0, M1 = fit_M1, M2 = fit_M2, M3 = fit_M3)
# The neural network needs keras/tensorflow. If they are not installed the script stops here instead of failing later.
if ('M3' %in% CFG$methods_to_run && !requireNamespace('keras', quietly = TRUE))
  stop("M3 (neural network) needs the 'keras' package. Install keras/tensorflow, or set CFG$methods_to_run <- c('M0', 'M1', 'M2') for a test run.")
# tag written in the checkpoint file names, so results of different method sets are never mixed
TAG <- paste(CFG$methods_to_run, collapse = '')

}  # end of sections 3-4

# ------------------------------------------------------------------------------------------------------------------
# 5. STEP 2 - HIDE-AND-REFILL TEST (edge-error hypothesis)
# ------------------------------------------------------------------------------------------------------------------
# For each state series, observed years are hidden in blocks (a 'mask'), the series is refilled with M0-M3 using ONLY the
# years that remain, and the prediction is compared with the hidden truth. No hidden value ever enters a fit.
#   interior4, interior10  5 random placements each of a 4- and a 10-year block, strictly inside the bulletin era
#   interior_long          the last 5 bulletin years (2003-2007): leaves a ~20-year hole (2003-2022) bounded by the 2023-24 data.
#                          Only for states with 2023-24 data (otherwise it would be an end block, not an interior gap).
#   start                  the first 5 observed years (predicted from the later ones)
#   end                    2023-24 hidden and predicted from the bulletin era alone (only for states with 2023-24 data)
# Error of one hidden year:  e = ln(predicted + 1) - ln(observed + 1)   (0 = exact; +0.69 = twice too high; -0.69 = half)
if (RUN$step2) {
make_masks <- function(d, cfg, seed) {
  set.seed(seed)
  yrs <- d$year[!is.na(d$catch)]                       # observed years of this state
  bul <- yrs[yrs <= cfg$bulletin_last]                 # observed years of the bulletin era
  late <- yrs[yrs >= cfg$mpa_obs_from]                 # observed years of the MPA collected landings (2023-24)
  masks <- list()
  for (L in cfg$gap_interior) {                        # interior blocks (4 and 10 years), strictly inside the bulletin era
    ok <- bul[vapply(bul, function(s) all(s:(s + L - 1) %in% bul), logical(1))]
    ok <- ok[ok > min(bul) & (ok + L - 1) < max(bul)]
    if (length(ok) > 0)
      for (s in ok[sample.int(length(ok), min(cfg$n_place, length(ok)))])
        masks[[length(masks) + 1]] <- list(type = paste0('interior', L), years = s:(s + L - 1))
  }
  if (length(late) > 0 && length(bul) > cfg$gap_edge + 5) {
    masks[[length(masks) + 1]] <- list(type = 'interior_long', years = tail(bul, cfg$gap_edge))
    masks[[length(masks) + 1]] <- list(type = 'end', years = late)
  }
  masks[[length(masks) + 1]] <- list(type = 'start', years = yrs[1:cfg$gap_edge])
  masks
}

# ---- 5.2 One series = one fresh R process (memory is returned to the system when the process ends) -------------
run_series_A <- function(sid, d, cfg, seed) {
  f <- file.path('models', paste0('stateA_', TAG, '_', gsub('[^A-Za-z0-9]', '_', sid), '.rds'))
  if (file.exists(f)) { say('   checkpoint found, skipping', sid); return(f) }   # checkpoint: skip finished series
  res <- callr::r(function(d, cfg, seed, fit_M0, fit_M1, fit_M2, fit_M3, make_masks) {
    out <- list(); mask_id <- 0
    masks <- make_masks(d, cfg, seed)
    for (mk in masks) {
      mask_id <- mask_id + 1
      if (isTRUE(cfg$verbose)) cat(sprintf('      mask %d/%d  %-13s  hides %d-%d (%d yr)\n', mask_id, length(masks), mk$type, min(mk$years), max(mk$years), length(mk$years)))
      dm <- d; truth <- dm$catch[dm$year %in% mk$years]; dm$catch[dm$year %in% mk$years] <- NA
      for (mn in cfg$methods_to_run) {
        fit  <- get(paste0('fit_', mn))
        pred <- if (mn == 'M3') fit(dm, reps = cfg$nn_reps_A, seed = seed) else fit(dm)
        pr   <- pred[d$year %in% mk$years]
        # SAVE the prediction itself (not only the error): metrics can be recomputed later without refitting anything
        out[[length(out) + 1]] <- data.frame(method = mn, gap_type = mk$type, mask_id = mask_id, year = mk$years,
                                             observed = truth, predicted = pr, e = log(pr + 1) - log(truth + 1))
      }
    }
    do.call(rbind, out)
  }, args = list(d, cfg, seed, fit_M0, fit_M1, fit_M2, fit_M3, make_masks), show = isTRUE(cfg$verbose),
    env = c(callr::rcmd_safe_env(), TF_CPP_MIN_LOG_LEVEL = '3', TF_ENABLE_ONEDNN_OPTS = '0', PYTHONWARNINGS = 'ignore'))   # quiet TensorFlow in the child process   # show = print the subprocess messages here
  res$series_id <- sid
  saveRDS(res, f); f
}
states_A <- coverage[usable_step2 == TRUE, state]
say('STEP 2: hide-and-refill on', length(states_A), 'states; methods:', paste(CFG$methods_to_run, collapse = ', '),
    '| each state runs in its own R process; finished states are skipped (checkpoints in models/)')
for (i in seq_along(states_A)) {
  sid <- states_A[i]; t0 <- Sys.time()
  say(sprintf('STEP 2 [%d/%d] %s: start', i, length(states_A), sid))
  run_series_A(sid, as.data.frame(catch[state == sid][order(year), .(year, catch)]), CFG, CFG$seed_sample)
  say(sprintf('STEP 2 [%d/%d] %s: done in %.1f min', i, length(states_A), sid, as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
}

errA <- rbindlist(lapply(list.files('models', paste0('^stateA_', TAG, '_'), full.names = TRUE), readRDS))
errA[, abs_e := abs(e)]
errA[, edge := fifelse(gap_type %in% c('start', 'end'), 'edge', 'interior')]     # interior_long is interior (bounded on both sides)
errA[, method := factor(method, levels = CFG$methods_to_run)]
errA[, gap_type := factor(gap_type, levels = c('interior4', 'interior10', 'interior_long', 'start', 'end'))]
write.xlsx(errA, 'results/Table_S1_step2_hidden_year_predictions.xlsx')   # every hidden year: observed, predicted, error

# ---- 5.3 Edge-error hypothesis: mixed model and decision rule --------------------------------------------------------
#   abs_e ~ method * gap_type + (1 | series_id/mask_id)   error by method and gap type (descriptive, Fig_4)
#   abs_e ~ method * edge     + (1 | series_id/mask_id)   edge (start + end) versus interior
#   The years of one hidden block (mask) share the same gap and are not independent, so the random effects are the state
#   AND the mask nested in the state; with only the state, the confidence intervals would be too narrow.
#   SUPPORTED if mean |e| edge / mean |e| interior >= CFG$edge_ratio AND the 95% CI of the difference excludes 0.
# Degrees of freedom (df): Satterthwaite approximation when the package lmerTest is installed (df is then a number in the tables).
# Without lmerTest, intervals are asymptotic (normal) and df is reported as Inf. (lmerTest.limit raised: emmeans switches the
# approximation off above 3000 observations by default, and here there are ~5,600 hidden years.)
say('STEP 2: fitting the mixed models (a minute or two)')
use_sat  <- requireNamespace('lmerTest', quietly = TRUE)
fit_lmer <- if (use_sat) lmerTest::lmer else lme4::lmer
emm_options(lmer.df = if (use_sat) 'satterthwaite' else 'asymptotic', lmerTest.limit = 1e5, pbkrtest.limit = 1e5)
fix_df <- function(x) { if ('df' %in% names(x)) x$df <- ifelse(is.finite(x$df), as.character(round(x$df, 1)), 'Inf (asymptotic)'); x }   # keep df readable in the xlsx
mA  <- fit_lmer(abs_e ~ method * gap_type + (1 | series_id/mask_id), data = errA)
mA2 <- fit_lmer(abs_e ~ method * edge + (1 | series_id/mask_id), data = errA)
saveRDS(mA, 'models/mA_abs_error.rds')
em_edge   <- emmeans(mA2, ~ edge)
edge_mean <- fix_df(as.data.frame(em_edge)); edge_mean$edge <- as.character(edge_mean$edge)
edge_diff_raw <- as.data.frame(summary(pairs(em_edge), infer = c(TRUE, TRUE)))   # edge - interior, with 95% CI
edge_diff <- edge_diff_raw
names(edge_diff)[grep('LCL|lower', names(edge_diff))] <- 'lower.CL'; names(edge_diff)[grep('UCL|upper', names(edge_diff))] <- 'upper.CL'
ratio_edge <- edge_mean$emmean[edge_mean$edge == 'edge'] / edge_mean$emmean[edge_mean$edge == 'interior']
edge_supported <- ratio_edge >= CFG$edge_ratio && edge_diff$lower.CL > 0
cat('\nEdge-error hypothesis: mean |e| edge / interior =', round(ratio_edge, 2),
    '| difference', round(edge_diff$estimate, 3), '[', round(edge_diff$lower.CL, 3), ',', round(edge_diff$upper.CL, 3), ']  ->',
    ifelse(edge_supported, 'SUPPORTED', 'NOT supported'), '\n')
em_gap <- as.data.frame(emmeans(mA, ~ method * gap_type))
names(em_gap)[grep('LCL', names(em_gap))] <- 'lower.CL'; names(em_gap)[grep('UCL', names(em_gap))] <- 'upper.CL'
em_met <- as.data.frame(emmeans(mA2, ~ method * edge))
# Edge / interior ratio of each METHOD separately (the pooled ratio above averages the four methods equally)
ratio_by_method <- dcast(as.data.table(em_met)[, .(method, edge = as.character(edge), emmean)], method ~ edge, value.var = 'emmean')
ratio_by_method[, `edge / interior` := round(edge / interior, 2)]
# Which method has the lowest mean |e| in each hidden block (one block = one state x one mask): a simple win count
best_by_block <- errA[, .(mae = mean(abs_e)), by = .(series_id, mask_id, gap_type, method)][, .SD[which.min(mae)], by = .(series_id, mask_id, gap_type)][, .N, by = .(gap_type, method)]
best_by_block <- dcast(best_by_block, gap_type ~ method, value.var = 'N', fill = 0)
em_gap <- fix_df(em_gap); em_met <- fix_df(em_met)
raw_sum <- errA[, .(n_hidden_years = .N, mean_abs_e = round(mean(abs_e), 3), median_abs_e = round(median(abs_e), 3),
                    typical_error = round(exp(mean(abs_e)) - 1, 3),
                    n_prediction_at_cap = sum(predicted >= 1e6 - 1)), by = .(method, gap_type)]   # cap = 1e6 t: a prediction at the cap is a failed extrapolation
# Statistical tests reported in the paper:
#  (a) F tests of the fixed effects (method, block type or edge/interior, and their interaction), Satterthwaite df;
#  (b) standard deviations of the random effects (mask within state, state) and of the residual, on the scale of |e|;
#  (c) pairwise differences among the four methods inside edge blocks and inside interior blocks (Tukey-adjusted, 95% CI);
#  (d) the edge - interior difference of each method separately.
tests_fixed <- if (use_sat) rbind(data.frame(model = 'method x gap_type', term = rownames(anova(mA)), as.data.frame(anova(mA)), check.names = FALSE),
                                  data.frame(model = 'method x edge', term = rownames(anova(mA2)), as.data.frame(anova(mA2)), check.names = FALSE)) else data.frame(note = 'lmerTest not installed: no F tests')
vc <- function(m, nm) { v <- as.data.frame(VarCorr(m)); data.frame(model = nm, component = ifelse(is.na(v$var1), 'residual', paste(v$grp, 'intercept')), sd = round(v$sdcor, 3)) }
var_components   <- rbind(vc(mA, 'method x gap_type'), vc(mA2, 'method x edge'))
method_contrasts <- as.data.frame(summary(pairs(emmeans(mA2, ~ method | edge), adjust = 'tukey'), infer = c(TRUE, TRUE)))
edge_by_method   <- as.data.frame(summary(contrast(emmeans(mA2, ~ edge | method), 'pairwise'), infer = c(TRUE, TRUE)))
write.xlsx(list(edge_hypothesis = data.frame(ratio_edge_over_interior = ratio_edge, diff = edge_diff$estimate, df = fix_df(edge_diff)$df, lower = edge_diff$lower.CL,
                                             upper = edge_diff$upper.CL, supported = edge_supported),
                model_means_by_method_gap = em_gap, model_means_by_method_edge = em_met, edge_ratio_by_method = ratio_by_method,
                best_method_per_block = best_by_block, raw_summary = raw_sum,
                tests_fixed_effects = tests_fixed, random_effect_sd = var_components,
                method_contrasts = method_contrasts, edge_minus_interior_by_method = edge_by_method),
           'results/Table_3_step2_hide_and_refill.xlsx')

# Fig 4: error by method and gap type (model means with 95% CI) and the raw errors as light points
# The y axis is ZOOMED to 0-fig4_ymax so the model means and CIs are readable; a few extreme hidden years (|e| up to ~8, i.e.
# a prediction e^8 = ~3000 times off) fall outside the window. coord_cartesian only zooms: the points are NOT removed and the
# means/CIs are unchanged. The number of points outside the window is printed in the caption.
fig4_ymax  <- 3.5
n_outside  <- errA[abs_e > fig4_ymax, .N]
p4 <- ggplot(em_gap, aes(gap_type, emmean, colour = method)) +
  geom_jitter(data = errA, aes(gap_type, abs_e, colour = method), inherit.aes = FALSE, alpha = 0.12, size = 0.6,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.6), size = 0.5) +
  scale_colour_manual(values = method_cols, name = 'Method') +
  coord_cartesian(ylim = c(0, fig4_ymax), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, fig4_ymax, 0.5)) +
  labs(x = 'Gap type (hidden block)', y = NULL,                                  # y title removed on purpose (explained in the caption)
       caption = paste('Y axis: absolute error on the log scale (0 = perfect). Large points: model mean and 95% CI (mixed model, random intercepts for state and block). Small points: every hidden year;',
                       sprintf('axis zoomed to 0-%.1f, so %d of %d extreme hidden years (%.1f%%) lie above the window (means and CIs use all of them).', fig4_ymax, n_outside, nrow(errA), 100 * n_outside / nrow(errA)),
                       'interior4/interior10: blocks inside the bulletin era; interior_long: 2003-2007 hidden (hole 2003-2022); start: first 5 observed years; end: 2023-24 hidden.', sep = '\n')) +
  theme_cru()
save_fig(p4, 'Fig_4_hide_and_refill_error', width = 11, height = 6)
}  # end of section 5 (RUN$step2)

# ------------------------------------------------------------------------------------------------------------------
# 6. STEP 3 - REBUILD ALL STATE SERIES, SUM TO REGIONS AND BRAZIL, COMPARE, ZBRT
# ------------------------------------------------------------------------------------------------------------------
# Each state series is refitted ONCE with every method on all its observed years; observed values are kept and the method
# fills only the gaps. The 17 filled series are summed to the 4 regions and to Brazil (years 1962-2024). The published
# reconstructions are summed the same way (and cut to 1962 on, the same start).
#
# ZBRT (Zhou et al. 2017, Fish and Fisheries 18(4): datalimited2::zbrt) estimates relative biomass B/K (saturation) from the
# catch history alone, with a boosted regression tree trained on simulated stocks; it returns B/K with a 95% interval.
# Here it is a RULER: it turns each catch series into the quantity a manager reads (how depleted?) so that differences among
# methods can be judged against ZBRT's own uncertainty. It is applied to regional/national TOTALS, which are not stocks, and the
# series start in 1962, not at the start of the fishery: the absolute B/K is therefore NOT an assessment. Only differences among
# series (same unit, same ruler) are interpreted.
if (RUN$step3) {
run_series_B <- function(sid, d, cfg, seed) {
  f <- file.path('models', paste0('stateB_', TAG, '_', gsub('[^A-Za-z0-9]', '_', sid), '.rds'))
  if (file.exists(f)) { say('   checkpoint found, skipping', sid); return(f) }
  res <- callr::r(function(d, cfg, seed, fit_M0, fit_M1, fit_M2, fit_M3) {
    out <- data.frame(year = d$year)
    for (mn in cfg$methods_to_run) {                                       # one column per method
      if (isTRUE(cfg$verbose)) cat('      fitting', mn, '\n')
      out[[mn]] <- if (mn == 'M3') fit_M3(d, reps = cfg$nn_reps_B, seed = seed) else get(paste0('fit_', mn))(d)
    }
    out
  }, args = list(d, cfg, seed, fit_M0, fit_M1, fit_M2, fit_M3), show = isTRUE(cfg$verbose),
    env = c(callr::rcmd_safe_env(), TF_CPP_MIN_LOG_LEVEL = '3', TF_ENABLE_ONEDNN_OPTS = '0', PYTHONWARNINGS = 'ignore'))
  res$state <- sid; saveRDS(res, f); f
}
say('STEP 3: refitting each state series ONCE per method on all its observed years (', length(REGION), 'states)')
for (i in seq_along(REGION)) {
  sid <- names(REGION)[i]; t0 <- Sys.time()
  say(sprintf('STEP 3 [%d/%d] %s: start', i, length(REGION), sid))
  run_series_B(sid, as.data.frame(catch[state == sid][order(year), .(year, catch)]), CFG, CFG$seed_sample)
  say(sprintf('STEP 3 [%d/%d] %s: done in %.1f min', i, length(REGION), sid, as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
}
say('STEP 3: summing states to regions and Brazil, comparing series, then ZBRT')

fit_all <- rbindlist(lapply(list.files('models', paste0('^stateB_', TAG, '_'), full.names = TRUE), readRDS))
fit_all <- melt(fit_all, id.vars = c('state', 'year'), variable.name = 'series', value.name = 'catch_hat')
fit_all <- merge(fit_all, catch[, .(state, year, obs = catch)], by = c('state', 'year'), all.x = TRUE)
fit_all[, `:=`(catch_t = fifelse(!is.na(obs), obs, catch_hat), observed = !is.na(obs))]     # observed kept, method fills gaps
fit_all[, series := as.character(series)]
state_series <- rbind(fit_all[, .(series, state, year, catch_t, observed)],
                      st[source %in% CFG$src_recon & year >= CFG$year_first & year <= CFG$year_last,
                         .(series = source, state, year, catch_t, observed = FALSE)])
state_series[, region := REGION[state]]
fwrite(state_series, 'data/state_series_rebuilt.csv')       # the rebuilt database: one row per series x state x year

# ---- 6.1 Sum to regions and Brazil -----------------------------------------------------------------------------------------
# A published series that does not cover a year (Silva & Andrade ends in 2011, Freire in 2015) has no value there (NA). The Silva & Andrade
# series is the merged one (observed total where a bulletin exists, GLM prediction elsewhere; section 1.8.2).
n_state <- c(table(REGION)); n_state['Brazil'] <- 17L
reg_series <- state_series[, .(catch_t = sum(catch_t), n = .N), by = .(series, unit = region, year)]
reg_series <- reg_series[n == n_state[unit]][, n := NULL]
bra_series <- state_series[, .(catch_t = sum(catch_t), n = .N), by = .(series, year)][n == 17][, `:=`(unit = 'Brazil', n = NULL)]
unit_series <- rbind(reg_series, bra_series, use.names = TRUE)
state_u <- state_series[, .(series, unit = state, year, catch_t)]
all_series <- rbind(unit_series, state_u, use.names = TRUE)
SERIES_M <- CFG$methods_to_run; SERIES_P <- CFG$src_recon

# ---- 6.2 Do the methods agree with each other and with the published reconstructions? --------------------------------------
# Same three metrics as Step 1 (ratio, typical difference, trend agreement), over the years both series cover.
# All years both series cover up to 2022 (years observed in a state are identical in M0-M3 by construction).
wide_s <- dcast(all_series, unit + year ~ series, value.var = 'catch_t')
pairs_s <- combn(c(SERIES_M, SERIES_P), 2, simplify = FALSE)
cmp <- rbindlist(lapply(pairs_s, function(p) {
  w <- wide_s[!is.na(get(p[1])) & !is.na(get(p[2])) & year <= CFG$mpa_recon_last]
  w[, c(list(pair = paste(p[1], 'vs', p[2])), agree_fun(year, get(p[1]), get(p[2]))), by = unit]
}))
cmp[, `:=`(unit_type = fifelse(unit %in% names(REGION), 'state', fifelse(unit == 'Brazil', 'Brazil', 'region')),
           ratio = round(ratio, 2), typical_difference = round(typical_difference, 2), trend_agreement = round(trend_agreement, 2))]
cmp_summary <- cmp[unit_type == 'state', .(n_states = .N, ratio = fmt3(ratio), typical_difference = fmt3(typical_difference),
                                          trend_agreement = fmt3(trend_agreement)), by = pair]
cat('\nStep 3 - agreement among methods and published reconstructions, median (min to max) across states:\n'); print(cmp_summary)

# ---- 6.3 Supporting descriptors: peak, year of the peak, recent catch over peak ---------------------------------------------
desc <- unit_series[unit %in% c(UNITS_REGION, 'Brazil'), .(peak_t = round(max(catch_t)), year_of_peak = year[which.max(catch_t)],
                                                           recent_over_peak = round(mean(tail(catch_t[order(year)], 5)) / max(catch_t), 2),
                                                           last_year = max(year)), by = .(unit, series)]

# ---- 6.4 ZBRT B/K for every state, region and Brazil series -------------------------------------------------------------
# zbrt_fun(year, catch) must return a data.frame with columns year, s (B/K), s_lo and s_hi (95% interval).
# The study unit is the STATE, so ZBRT is run on each state series too (CFG$zbrt_states); regions and Brazil are sums of the states.
# datalimited2 prints two harmless warnings that clutter the console: "Breakpoint estimate(s) outdistanced..." (its segmented
# regression on the last years of catch has no clear breakpoint) and "essentially perfect fit" (the last years of a series are an exact
# straight line or constant, which is exactly what a filled gap produces). They do not change the B/K, so they are silenced here;
# an actual failure of a run is NOT silenced: it is reported as 'ZBRT failed' and that series is skipped.
if (!exists('zbrt_fun')) {
  if (!requireNamespace('datalimited2', quietly = TRUE))
    stop("ZBRT needs the package datalimited2: remotes::install_github('cfree14/datalimited2') (it also needs JAGS).")
  zbrt_fun <- function(year, catch) { z <- datalimited2::zbrt(year = year, catch = catch)$ts; as.data.frame(z) }
}
zb_units <- c(if (isTRUE(CFG$zbrt_states)) names(REGION), UNITS_REGION, 'Brazil')
zb_input <- all_series[unit %in% zb_units & !is.na(catch_t)]
zb_list  <- split(zb_input, by = c('unit', 'series'))
say('ZBRT:', length(zb_list), 'series (unit x series), about 4 s each')
zb_i <- 0
zb <- rbindlist(lapply(zb_list, function(x) {
  x <- x[order(year)]
  zb_i <<- zb_i + 1; say(sprintf('   ZBRT %d/%d  %s  %s', zb_i, length(zb_list), x$unit[1], x$series[1]))
  if (nrow(x) < 15 || sum(x$catch_t) <= 0) { say('      ZBRT skipped (fewer than 15 years or no catch)'); return(NULL) }
  z <- tryCatch(suppressWarnings(as.data.table(zbrt_fun(x$year, x$catch_t))),
                error = function(e) { say('      ZBRT failed:', conditionMessage(e)); NULL })
  if (is.null(z)) return(NULL)
  z[, .(unit = x$unit[1], series = x$series[1], year, bk = s, bk_lo = s_lo, bk_hi = s_hi)]
}))
zb[, unit_type := fifelse(unit %in% names(REGION), 'state', fifelse(unit == 'Brazil', 'Brazil', 'region'))]
write.xlsx(list(zbrt_all_years = zb), 'results/Table_S2_zbrt_all_years.xlsx')

# ---- 6.5 Divergence hypothesis (decided on the 17 states) ---------------------------------------------------------------------
# ZBRT returns B/K for EVERY year (each year uses the catch history up to that year; the first estimate is for year 10 of the
# series), so the whole depletion trajectory of every series is compared, not only its last year.
# For a unit and a year, the four methods DIVERGE if the spread of B/K among them (max - min) is larger than the mean width of
# ZBRT's own 95% interval for those series that year (the difference would change the conclusion by more than ZBRT's uncertainty).
# Each year belongs to one of three ERAS of the data: the bulletin era (up to 2007), the long gap with no observed catch (2008-2022)
# and the recent MPA years (2023-24). The share of years that diverge is computed overall and in each era (descriptive).
#   Hypothesis - divergence along the series: the methods diverge in MORE than CFG$diverge_share of all years in MORE than
#                CFG$min_share_states of the 17 states. The state is the unit of the whole study; regions and Brazil (sums of
#                the states) and the published reconstructions get the same criterion as a descriptive reference.
era_of <- function(y) fcase(y <= CFG$bulletin_last, 'bulletin era (to 2007)', y <= CFG$mpa_recon_last, 'long gap (2008-2022)', default = 'recent MPA years (2023-24)')
# (ZBRT gives no B/K for the first years of a series: those years are NA and are left out here.)
spread <- function(x, series_set, label) {
  x[series %in% series_set & !is.na(bk) & !is.na(bk_lo) & !is.na(bk_hi), .(spread_bk = max(bk) - min(bk), mean_ci_width = mean(bk_hi - bk_lo), n_series = .N), by = .(unit, unit_type, year)][
    n_series == length(series_set)][, `:=`(group = label, diverge = spread_bk > mean_ci_width, era = era_of(year))]
}
spr <- rbind(spread(zb, SERIES_M, 'our methods (M0-M3)'), spread(zb, SERIES_P, 'published reconstructions'))
share_tab <- spr[, .(years = .N, share_all = round(mean(diverge), 2),
                     share_bulletin = round(mean(diverge[era == 'bulletin era (to 2007)']), 2),
                     share_long_gap = round(mean(diverge[era == 'long gap (2008-2022)']), 2),
                     share_recent = round(mean(diverge[era == 'recent MPA years (2023-24)']), 2),
                     spread_over_ci_all = round(mean(spread_bk / mean_ci_width), 2),
                     spread_over_ci_long_gap = round(mean((spread_bk / mean_ci_width)[era == 'long gap (2008-2022)']), 2)),
                 by = .(group, unit_type, unit)]
share_tab[, diverges_more_than_half_of_years := share_all > CFG$diverge_share]
share_tab[is.na(diverges_more_than_half_of_years), diverges_more_than_half_of_years := FALSE]
mine <- share_tab[group == 'our methods (M0-M3)']
n_states <- mine[unit_type == 'state', sum(diverges_more_than_half_of_years)]; of_states <- mine[unit_type == 'state', .N]
hyp_tab <- data.table(hypothesis = 'divergence along the series (more than half of the years in more than half of the states)',
                      states_meeting_criterion = n_states, of_states = of_states, needed = floor(CFG$min_share_states * of_states) + 1L,
                      supported = n_states > CFG$min_share_states * of_states,
                      regions_meeting_criterion_descriptive = mine[unit_type == 'region', sum(diverges_more_than_half_of_years)], of_regions = 4L,
                      brazil_meets_criterion_descriptive = mine[unit_type == 'Brazil', diverges_more_than_half_of_years])
cat('\nDivergence hypothesis (our methods M0-M3 vs ZBRT uncertainty; criterion: share of years > ', CFG$diverge_share, '):\n', sep = ''); print(hyp_tab)
cat('\nShare of years in which the methods diverge, by unit and era:\n'); print(mine[order(factor(unit_type, c('Brazil', 'region', 'state')), unit)])
zb_last <- zb[, .SD[year == max(year)], by = .(unit, series)]
write.xlsx(list(hypotheses = hyp_tab, divergence_by_unit_and_era = share_tab, spread_by_year = spr, bk_last_year_by_series = zb_last, descriptors = desc,
                agreement_summary = cmp_summary, agreement_all_units = cmp),
           'results/Table_4_step3_rebuilt_series_zbrt.xlsx')

# ---- 6.6 Figures -----------------------------------------------------------------------------------------------------------------
# Same look as Fig 1: one panel per unit, free y axis, thin grid, legend at the bottom, explanatory caption.
series_cols <- c(source_cols, method_cols)
unit_lab <- c(REG_LAB, Brazil = 'Brazil')
state_lvls <- c(names(REGION), 'Brazil')                     # panel order of Fig 1 (north to south, then Brazil)
reg_lvls   <- c(UNITS_REGION, 'Brazil')

# The rebuilt series themselves are NOT plotted: with the generic category in the totals, most years are observed and the methods
# differ only in a few gap years, so the published reconstructions (Fig 1) are the only catch-series figure. The consequence of the
# method is shown where it matters: on the ZBRT B/K (Figs 5-6) and on the disagreement among methods (Figs 7-8).

# Figs 5 (regions + Brazil) and 6 (the 17 states + Brazil): ZBRT B/K through the whole series, every series. B/K = 0.5 (dotted) is
# the B/K at maximum sustainable yield of a Schaefer model, drawn as a reference only; B/K here is not an assessment (see the caption).
bk_plot <- function(units, lvls, labs, ncol) {
  zr <- zb[unit %in% units]
  zr[, `:=`(unit_f = factor(unit, levels = lvls, labels = labs), kind = fifelse(series %in% SERIES_M, 'Our methods', 'Published reconstructions'))]
  ggplot(zr, aes(year, bk, colour = series)) +
    geom_hline(yintercept = 0.5, linetype = 'dotted', colour = ink2, linewidth = 0.4) +
    geom_line(aes(linewidth = kind), alpha = 0.95) +
    facet_wrap(~ unit_f, ncol = ncol) +
    scale_colour_manual(values = series_cols, breaks = names(series_cols), name = NULL) +
    scale_linewidth_manual(values = c('Our methods' = 0.8, 'Published reconstructions' = 0.4), guide = 'none') +
    scale_x_continuous(breaks = seq(1970, 2020, 30)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = 'ZBRT B/K (relative biomass)',
         caption = paste('B/K estimated year by year from the catch history up to that year (first estimate: year 10; series start in 1962). Dotted line: B/K = 0.5 (Schaefer reference).',
                         'Only differences among series are interpreted, not the absolute level (series start in 1962 and the units are totals, not stocks).', sep = '\n')) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2))) + theme_cru()
}
save_fig(bk_plot(reg_lvls, reg_lvls, unit_lab[reg_lvls], 5), 'Fig_5_zbrt_bk_regions', width = 15, height = 4.6)
if (isTRUE(CFG$zbrt_states)) save_fig(bk_plot(state_lvls, state_lvls, state_lvls, 6), 'Fig_6_zbrt_bk_states', width = 15, height = 8.5)

# Figs 7 (regions + Brazil) and 8 (states + Brazil): how much the series disagree in each year. The y axis is a DIFFERENCE between series (max B/K - min B/K), not a B/K:
# purple line = spread among our four methods; grey line = spread among the three published reconstructions; grey band = mean
# width of ZBRT's own 95% interval for our four series. Where the purple line leaves the band, the methods differ by more than
# ZBRT's uncertainty (= 'diverge'). The shaded stripe is the long gap with no observed catch (2008-2022).
spread_plot <- function(units, lvls, labs, ncol) {
  d <- spr[unit %in% units][, unit_f := factor(unit, levels = lvls, labels = labs)]
  band <- d[group == 'our methods (M0-M3)']
  ggplot(d, aes(year, spread_bk)) +
    annotate('rect', xmin = CFG$bulletin_last + 0.5, xmax = CFG$mpa_recon_last + 0.5, ymin = -Inf, ymax = Inf, fill = '#e8c98a', alpha = 0.3) +   # warm tint so the stripe is not confused with the grey ZBRT band
    geom_ribbon(data = band, aes(ymin = 0, ymax = mean_ci_width, fill = 'ZBRT 95% interval width'), alpha = 0.9) +
    geom_line(data = d[group == 'published reconstructions'], aes(colour = 'Among the 3 published reconstructions'), linewidth = 0.4) +
    geom_line(data = band, aes(colour = 'Among our 4 methods (M0-M3)'), linewidth = 0.8) +
    facet_wrap(~ unit_f, ncol = ncol) +
    scale_colour_manual(values = c('Among our 4 methods (M0-M3)' = unname(method_cols['M1']), 'Among the 3 published reconstructions' = ink2), name = 'Spread of B/K') +
    scale_fill_manual(values = c('ZBRT 95% interval width' = grid_col), name = NULL) +
    scale_x_continuous(breaks = seq(1970, 2020, 30)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = 'Spread of ZBRT B/K (max - min)',
         caption = paste('Spread = highest minus lowest B/K among the series in each year (a difference, not a B/K). Grey band = mean width of ZBRT\'s 95% interval: a line above the band means the series',
                         'differ by more than ZBRT\'s own uncertainty. Shaded stripe = long gap with no observed catch (2008-2022); to its left the bulletin era, to its right the MPA years.', sep = '\n')) +
    guides(colour = guide_legend(order = 1, override.aes = list(linewidth = 1.2)), fill = guide_legend(order = 2)) + theme_cru()
}
save_fig(spread_plot(reg_lvls, reg_lvls, unit_lab[reg_lvls], 5), 'Fig_7_zbrt_spread_regions', width = 15, height = 4.6)
if (isTRUE(CFG$zbrt_states)) save_fig(spread_plot(state_lvls, state_lvls, state_lvls, 6), 'Fig_8_zbrt_spread_states', width = 15, height = 8.5)
cat('\nStep 3 done.\n')
}  # end of section 6 (RUN$step3)

# ------------------------------------------------------------------------------------------------------------------
# 7. SESSION INFO
# ------------------------------------------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), 'results/session_info.txt')
