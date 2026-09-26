# ==============================================================================
# DIET COMPOSITION AND ENVIRONMENTAL CORRELATES OF THE GREY WOLF
# IN SLOVENIA AND CROATIA - SCRIPT UNICO DI ANALISI
#
# Alfio Tomarchio - Universita' di Bologna
# Ultima revisione: 26 settembre 2026
#
# COME SI USA
#   1. Metti questo script nella cartella "File definitivi", insieme ai file
#      elencati sotto, e aprila come cartella di lavoro. In RStudio:
#      Session > Set Working Directory > Choose Directory..., oppure
#        setwd("C:/Users/OEM/Desktop/TESI MAGISTRALE LJUBLJANA ERASMUS/File definitivi")
#   2. Esegui:  source("Script_Tesi_Lupo_COMPLETO.R", encoding = "UTF-8")
#      In RStudio i grafici compaiono anche nel pannello Plots.
#   3. I grafici (.png) e i CSV vengono riscritti nella cartella con gli stessi
#      nomi di prima; numeri_risultati.tex (i numeri della tesi) va caricato su
#      Overleaf al posto di quello vecchio, insieme ai grafici .png (se su
#      Overleaf ci sono anche .pdf con lo stesso nome, la tesi usa quelli:
#      sostituiscili con quelli di output_R/figure_pdf/ oppure cancellali).
#   4. In output_R/ trovi le versioni PDF dei grafici, le tabelle, i dati
#      derivati, il log (log_esecuzione.txt) e i confronti con la bozza.
#
# FILE NELLA CARTELLA
#   Stessi nomi della versione precedente:
#     Taxonomic_coverage_Slovenia.csv, Taxonomic_coverage_Croatia.csv
#         tabelle delle varianti (dati grezzi): da qui si ricostruiscono le matrici
#     Slovenia wolves FILE DEFINITIVO.csv, Croatian wolves FILE DEFINITIVO.csv
#         matrici finali: ora servono solo per il confronto con le matrici
#         ricostruite (le differenze vanno in output_R/dati_derivati/)
#     Wolf_NDVI_Dynamic_Dinaric.csv, individui_e_branchi_103_campioni.csv,
#     Punti_103_per_QGIS.csv (separato da ";" oppure da ",")
#   Nuovi (26/09/2026):
#     assegnazioni_varianti.csv   assegnazione finale e ruolo di ogni variante
#     gruppi_geografici_103.csv   aree di studio, gruppi, branchi e regioni
#     specifica_macro.csv         elenco dei numeri citati nella tesi
#   Facoltativi (se mancano, i relativi controlli vengono saltati):
#     Punti_130_campionamento.csv, risultati_chiave_python.csv,
#     numeri_risultati_provvisori.tex
#
# COSA E' CAMBIATO RISPETTO ALLA VERSIONE PRECEDENTE
#   - Le matrici si ricostruiscono dalle tabelle delle varianti, con filtri e
#     soglia applicati dallo script, invece di partire dai FILE DEFINITIVO.
#   - MSV136: ripristinati 45 reads di capriolo; 86 reads restano Cervinae non
#     risolti (nuova categoria "Cervinae", solo in questo campione).
#   - Numero minimo di prede distinte accanto al numero di categorie.
#   - Analisi per individuo, errori standard robusti, modelli misti, analisi di
#     sensibilita'; tutti i numeri della tesi in numeri_risultati.tex.
#   - Nei file delle torte: colonne CERVINAE e NPREY_MIN in piu'.
#
# SCELTE DOCUMENTATE (modificabili nel blocco di configurazione)
#   - filtro di profondita' minima: 2.000 reads totali per campione nella
#     tabella delle varianti (tutte le varianti); l'inclusione e' identica con
#     il denominatore delle sole prede (verificato e salvato).
#   - esclusione dai denominatori: predatore (Canis lupus), carnivori non
#     target (Vulpes vulpes), assegnazioni sopra la famiglia, varianti non
#     validate (variante Bovidae di 87 reads, esclusa).
#   - soglia dell'1%: quota di ogni categoria sul totale delle prede valide
#     del campione prima della soglia; si azzerano le quote strettamente < 1%;
#     un solo denominatore, nessuna iterazione; poi RRA ricalcolata. Analisi di
#     sensibilita' con una soglia comune del 5%, perche' la tabella croata era
#     gia' filtrata al 5% per campione e quella slovena no.
#   - categorie annidate (Cervinae > Cervus elaphus; Caprinae > Rupicapra,
#     Ovis, Capra; Ovis > Ovis aries): restano distinte nelle matrici; il
#     numero minimo di prede distinte non conta la categoria superiore quando
#     e' presente una sua categoria inclusa.
#   - dipendenza fra campioni dello stesso individuo: analisi a livello di
#     individuo (profili medi) per PERMANOVA, betadisper e RDA; errori
#     standard robusti per individuo (sandwich::vcovCL) per i modelli
#     logistici; GLMM con intercetta casuale per individuo (lme4::glmer) per
#     il modello sulla dieta mista.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURAZIONE
# ==============================================================================
rm(list = ls())

PROJECT_DIR <- getwd()          # la cartella "File definitivi"

# --- file di ingresso: stessi nomi della versione precedente ------------------
FILE_SLO       <- "Slovenia wolves FILE DEFINITIVO.csv"   # solo confronto
FILE_CRO       <- "Croatian wolves FILE DEFINITIVO.csv"   # solo confronto
FILE_AMBIENTE  <- "Wolf_NDVI_Dynamic_Dinaric.csv"
FILE_INDIVIDUI <- "individui_e_branchi_103_campioni.csv"
FILE_PUNTI     <- "Punti_103_per_QGIS.csv"
FILE_COV_SLO   <- "Taxonomic_coverage_Slovenia.csv"       # tabella delle varianti
FILE_COV_CRO   <- "Taxonomic_coverage_Croatia.csv"        # tabella delle varianti
# --- file nuovi (26/09/2026) --------------------------------------------------
FILE_ASSEGN    <- "assegnazioni_varianti.csv"
FILE_GRUPPI    <- "gruppi_geografici_103.csv"
FILE_SPEC      <- "specifica_macro.csv"
# --- facoltativi: se mancano, i controlli relativi vengono saltati -------------
FILE_P130      <- "Punti_130_campionamento.csv"
FILE_PROV      <- "risultati_chiave_python.csv"
FILE_TEX_PROV  <- "numeri_risultati_provvisori.tex"

# --- uscite con i nomi della versione precedente (nella cartella) --------------
OUT_FOO_RRA    <- "FOO_RRA_Slovenia_and_Croatia.csv"
OUT_COMMUNITY  <- "Community_Matrix_SloCro.csv"
OUT_TEX        <- "numeri_risultati.tex"
# --- uscite nuove (in output_R/) ------------------------------------------------
DIR_OUT  <- file.path(PROJECT_DIR, "output_R")
DIR_FIG  <- file.path(DIR_OUT, "figure_pdf")
DIR_TAB  <- file.path(DIR_OUT, "tabelle")
DIR_DAT  <- file.path(DIR_OUT, "dati_derivati")

MIN_READS   <- 2000        # filtro di profondita' minima (reads totali)
SOGLIA      <- 0.01        # soglia campione-specifica dell'1%
N_PERM      <- 9999        # permutazioni per tutti i test
SEED        <- 20260926    # seme unico; ogni test lo reimposta prima di partire
VARIANTE_BOVIDAE_87 <- "3833:11173"   # variante non validata (sensibilita')
VARIANTE_CERVINAE_CRO <- "21331:1117" # variante croata risolta a Cervus elaphus (sensibilita')

pacchetti <- c("dplyr", "tidyr", "readr", "ggplot2", "vegan", "permute",
               "MASS", "sandwich", "lme4", "patchwork", "ggrepel")
mancanti <- pacchetti[!vapply(pacchetti, requireNamespace, logical(1), quietly = TRUE)]
if (length(mancanti) > 0) {
  stop("Pacchetti mancanti: ", paste(mancanti, collapse = ", "),
       "\nInstallali con: install.packages(c(\"", paste(mancanti, collapse = "\", \""), "\"))")
}
# MASS prima di dplyr: dplyr::select non viene mascherato; le chiamate usano
# comunque il prefisso del pacchetto dove c'e' ambiguita'.
suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(vegan); library(permute); library(sandwich); library(lme4)
  library(patchwork); library(ggrepel)
})
if (getRversion() < "4.4.0") {
  message("Nota: con R < 4.4.0 confint() sui glm usa il profilo di MASS; i valori coincidono.")
}

for (d in c(DIR_OUT, DIR_FIG, DIR_TAB, DIR_DAT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
for (f in c(FILE_COV_SLO, FILE_COV_CRO, FILE_SLO, FILE_CRO, FILE_INDIVIDUI, FILE_PUNTI, FILE_AMBIENTE,
            FILE_ASSEGN, FILE_GRUPPI, FILE_SPEC)) {
  if (!file.exists(f)) stop("File mancante nella cartella di lavoro (", PROJECT_DIR, "): ", f)
}

# --- log di esecuzione --------------------------------------------------------
F_LOG <- file.path(DIR_OUT, "log_esecuzione.txt")
cat("Log di esecuzione - Script_Tesi_Lupo_COMPLETO.R\n", file = F_LOG)
LOG <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = F_LOG, append = TRUE)
}
AVVISI <- character(0)
avviso <- function(w) {                 # per withCallingHandlers(): registra e silenzia
  AVVISI <<- c(AVVISI, conditionMessage(w))
  LOG("AVVISO: ", conditionMessage(w))
  invokeRestart("muffleWarning")
}
segnala <- function(...) {              # avvisi dello script: nel log e in avvisi.txt
  msg <- paste0(...)
  AVVISI <<- c(AVVISI, msg)
  LOG("AVVISO: ", msg)
}
versione <- function(pk) utils::packageDescription(pk, fields = "Version")   # es. "2.7-5"
LOG("Cartella del progetto: ", PROJECT_DIR)
LOG("R ", as.character(getRversion()), "; seme ", SEED, "; permutazioni ", N_PERM)
# sessionInfo subito, cosi' resta anche se l'esecuzione si interrompe
writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "sessionInfo_inizio.txt"))

# --- archivio dei risultati chiave (stessi nomi della versione Python) --------
K <- list()
put <- function(key, value) { K[[key]] <<- value; invisible(value) }

CATS <- c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa", "Caprinae",
          "Rupicapra rupicapra", "Bos", "Capra", "Ovis", "Ovis aries", "Lepus", "Cervinae")
CODE <- c("Cervus elaphus" = "cerela", "Capreolus capreolus" = "capcap", "Sus scrofa" = "susscr",
          "Caprinae" = "caprinae", "Rupicapra rupicapra" = "ruprup", "Bos" = "bos", "Capra" = "capra",
          "Ovis" = "ovis", "Ovis aries" = "ovisaries", "Lepus" = "lepus", "Cervinae" = "cervinae",
          "Domestic" = "domestic", "Bovidae" = "bovidae")
DOMESTIC <- c("Bos", "Capra", "Ovis", "Ovis aries")
RDA_CATS <- c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa", "Caprinae",
              "Rupicapra rupicapra", "Domestic", "Lepus", "Cervinae")
AREE <- c("Slovenia", "Croatia")
ab <- function(a) c(Slovenia = "slo", Croatia = "cro", Total = "tot")[[a]]


# ==============================================================================
# 1. TABELLE DELLE VARIANTI: LETTURA PER NOME DI COLONNA
# ==============================================================================
short_id <- function(m) {
  x <- regmatches(m, regexpr("[0-9]+:[0-9]+(?=_SUB)", m, perl = TRUE))
  if (length(x) != length(m)) stop("Identificativo MOTU non riconosciuto")
  x
}
come_intero <- function(x, cosa) {
  y <- suppressWarnings(as.numeric(x))
  if (any(is.na(y))) stop("Valori mancanti o non numerici in ", cosa, ": nessuna sostituzione con zero")
  if (any(abs(y - round(y)) > 1e-9)) stop("Valori non interi in ", cosa)
  as.integer(round(y))
}

# lettura come testo, senza intestazione: le righe di intestazione si leggono a mano
leggi_grezzo <- function(f) {
  x <- read_delim(f, delim = ";", col_names = FALSE, col_types = cols(.default = col_character()),
                  locale = locale(encoding = "UTF-8"), trim_ws = TRUE, na = character(),
                  progress = FALSE, show_col_types = FALSE)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[] <- lapply(x, function(v) sub("^\ufeff", "", v))
  x
}
assegn <- read_csv(FILE_ASSEGN, show_col_types = FALSE, col_types = cols(.default = col_character()))
stopifnot(!anyDuplicated(paste(assegn$run, assegn$motu_id)))

# --- 1.1 Slovenia: Taxonomic_coverage_Slovenia.csv -------------------------------
# Stessa struttura del foglio READ COUNT del composition report: la riga 2 e'
# l'intestazione principale (66 campioni); il blocco di destra (MSV1M2 ripetuto,
# MSV1MU, sequence) ha l'intestazione nella riga 1 e contiene valori per taxon
# accorpato, non per variante. I conteggi di MSV1MU si ricostruiscono per
# differenza dal totale di ogni variante.
rc_raw <- leggi_grezzo(FILE_COV_SLO)
h1 <- as.character(unlist(rc_raw[1, ])); h2 <- as.character(unlist(rc_raw[2, ]))
col_sn <- which(h2 == "scientific_name"); col_ct <- which(h2 == "count_total")
stopifnot(length(col_sn) == 1, length(col_ct) == 1)
col_camp <- which(grepl("^(M1|M2|MSV)[0-9A-Z]+$", h2))
camp_slo <- h2[col_camp]
if (anyDuplicated(camp_slo)) stop("Nomi di campione duplicati nella riga 2 di ", FILE_COV_SLO)
stopifnot("MSV1MU" %in% h1, !("MSV1MU" %in% h2))
LOG("SLO: ", length(camp_slo), " colonne-campione con intestazione nella riga 2; MSV1MU disallineato")
rc <- rc_raw[-(1:2), , drop = FALSE]
rc <- rc[rc[[col_sn]] != "", , drop = FALSE]
var_slo <- tibble(scientific_name = rc[[col_sn]],
                  count_total = come_intero(rc[[col_ct]], "count_total SLO"))
for (j in seq_along(col_camp)) {
  var_slo[[camp_slo[j]]] <- come_intero(rc[[col_camp[j]]], paste(FILE_COV_SLO, camp_slo[j]))
}
# MSV1MU: conteggi per variante = totale della variante - somma degli altri 66 campioni
var_slo$MSV1MU <- var_slo$count_total - rowSums(var_slo[, camp_slo])
if (any(var_slo$MSV1MU < 0)) stop("Residui negativi nella ricostruzione di MSV1MU")
camp_slo <- c(camp_slo, "MSV1MU")
if (sum(var_slo$count_total) != 839437) stop("Totale sloveno diverso da 839,437 reads: file cambiato?")
if (anyDuplicated(var_slo[, c("scientific_name", "count_total")])) stop("Varianti slovene non distinguibili")
# la tabella slovena non ha l'identificativo del MOTU: ogni variante si abbina alla
# sua assegnazione per nome tassonomico e totale di reads, che sono unici
a_slo <- assegn %>% filter(run == "SLO") %>%
  transmute(motu_id, short_id, scientific_name = pipeline_name, count_total = as.integer(reads_total_run))
var_slo <- inner_join(var_slo, a_slo, by = c("scientific_name", "count_total"), relationship = "one-to-one",
                      unmatched = "error")   # ogni variante deve avere la sua assegnazione, e viceversa
if (nrow(var_slo) != nrow(a_slo)) stop("Abbinamento varianti slovene - assegnazioni incompleto")
var_slo$run <- "SLO"
LOG("SLO: ", nrow(var_slo), " varianti abbinate alle assegnazioni per (scientific_name, count_total)")
# controllo di MSV1MU sui valori per taxon del blocco di destra
col_mu <- which(h1 == "MSV1MU")
mu_blocco <- suppressWarnings(as.numeric(unlist(rc_raw[-1, col_mu])))
mu_blocco <- sort(mu_blocco[!is.na(mu_blocco) & mu_blocco > 0])
LOG("MSV1MU, valori per taxon nel blocco di destra: ", paste(mu_blocco, collapse = ", "))

# --- 1.2 Croazia: Taxonomic_coverage_Croatia.csv ---------------------------------
rcc_raw <- leggi_grezzo(FILE_COV_CRO)
hc <- as.character(unlist(rcc_raw[1, ]))
cc <- function(n) { j <- which(hc == n); if (length(j) != 1) stop("Colonna ", n, " non trovata in ", FILE_COV_CRO); j }
rcc <- rcc_raw[-1, , drop = FALSE]
rcc <- rcc[rcc[[cc("motus")]] != "", , drop = FALSE]
camp_cro <- hc[grepl("^HRV", hc)]
if (anyDuplicated(camp_cro)) stop("Campioni croati duplicati")
var_cro <- tibble(motu_id = rcc[[cc("motus")]], count_total = come_intero(rcc[[cc("count_total")]], "count_total CRO"))
for (x in camp_cro) var_cro[[x]] <- come_intero(rcc[[cc(x)]], paste(FILE_COV_CRO, x))
stopifnot(all(rowSums(var_cro[, camp_cro]) == var_cro$count_total))
if (sum(var_cro$count_total) != 3980866) stop("Totale croato diverso da 3,980,866 reads: file cambiato?")
a_cro <- assegn %>% filter(run == "CRO") %>% dplyr::select(motu_id, short_id)
var_cro <- inner_join(var_cro, a_cro, by = "motu_id", relationship = "one-to-one", unmatched = "error")
if (nrow(var_cro) != nrow(a_cro)) stop("Abbinamento varianti croate - assegnazioni incompleto")
var_cro$run <- "CRO"
LOG("CRO: ", nrow(var_cro), " varianti, ", length(camp_cro), " campioni; somme verificate")


# ==============================================================================
# 2. ASSEGNAZIONI CURATE E PIPELINE DEI CAMPIONI
# ==============================================================================
chk <- function(v, run) {
  a <- assegn[assegn$run == run, ]
  if (!setequal(a$motu_id, v$motu_id)) stop("Le assegnazioni curate non corrispondono alle varianti ", run)
}
chk(var_slo, "SLO"); chk(var_cro, "CRO")
# categorie annidate: categoria -> categoria che la include
PARENT <- c("Cervus elaphus" = "Cervinae", "Rupicapra rupicapra" = "Caprinae", "Ovis aries" = "Ovis",
            "Ovis" = "Caprinae", "Capra" = "Caprinae")
LABEL <- c("Cervinae" = "Unresolved Cervinae", "Caprinae" = "Unresolved Caprinae",
           "Bovidae" = "Unresolved Bovidae")

formato_lungo <- function(v, camp, run) {
  v %>% dplyr::select(motu_id, short_id, all_of(camp)) %>%
    pivot_longer(all_of(camp), names_to = "Sample", values_to = "reads") %>%
    mutate(run = run)
}
LUNGO <- bind_rows(formato_lungo(var_slo, camp_slo, "SLO"), formato_lungo(var_cro, camp_cro, "CRO")) %>%
  left_join(assegn %>% dplyr::select(run, motu_id, final_name, role), by = c("run", "motu_id")) %>%
  mutate(Area = if_else(run == "SLO", "Slovenia", "Croatia"))
stopifnot(!any(is.na(LUNGO$final_name)))

pipeline <- function(lungo, prede_extra = character(0), cervinae_cro_irrisolti = FALSE, soglia = SOGLIA) {
  l <- lungo
  l$role[l$short_id %in% prede_extra] <- "prey"
  if (cervinae_cro_irrisolti) l$final_name[l$run == "CRO" & l$short_id == VARIANTE_CERVINAE_CRO] <- "Cervinae"
  ruoli <- c("prey", "predator", "non_target_carnivore", "not_attributable", "not_validated")
  tot <- l %>% group_by(Area, Sample, role) %>% summarise(r = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = role, values_from = r, values_fill = 0)
  for (r in ruoli) if (!r %in% names(tot)) tot[[r]] <- 0L
  tot <- tot %>% mutate(
    total_reads = prey + predator + non_target_carnivore + not_attributable + not_validated,
    pass_depth_total = total_reads >= MIN_READS,
    pass_depth_prey = prey >= MIN_READS,
    pass_depth_noCanis = (total_reads - predator) >= MIN_READS,
    pass_depth_noCanisVulpes = (total_reads - predator - non_target_carnivore) >= MIN_READS,
    status = case_when(!pass_depth_total ~ "Below 2,000 total reads",
                       prey == 0 ~ "No valid prey reads",
                       TRUE ~ "Analysed"))
  inclusi <- tot$Sample[tot$status == "Analysed"]
  c0 <- l %>% filter(role == "prey", Sample %in% inclusi) %>%
    group_by(Area, Sample, final_name) %>% summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = final_name, values_from = reads, values_fill = 0) %>%
    arrange(factor(Area, levels = AREE), Sample)
  categorie <- setdiff(names(c0), c("Area", "Sample"))
  m0 <- as.matrix(c0[, categorie]); rownames(m0) <- c0$Sample
  den <- rowSums(m0)
  quota <- m0 / den
  rimossi <- which(m0 > 0 & quota < soglia, arr.ind = TRUE)
  tab_rimossi <- tibble(Sample = rownames(m0)[rimossi[, 1]], category = categorie[rimossi[, 2]],
                        reads = m0[rimossi], valid_prey_reads_pre_threshold = den[rimossi[, 1]],
                        share_pct = 100 * quota[rimossi])
  m1 <- m0; m1[quota < soglia] <- 0L
  den1 <- rowSums(m1)
  if (any(den1 == 0)) stop("Campione senza reads dopo la soglia")
  rra <- 100 * m1 / den1
  tenere <- colSums(m1) > 0
  list(totals = tot, counts_pre = m0, counts = m1[, tenere, drop = FALSE], rra = rra[, tenere, drop = FALSE],
       area = setNames(c0$Area, c0$Sample), den_pre = den, den_post = den1, removed = tab_rimossi)
}
RES <- pipeline(LUNGO)
cat_trovate <- colnames(RES$counts)
LOG("Categorie analitiche dopo la soglia: ", length(cat_trovate), " (", paste(cat_trovate, collapse = ", "), ")")
if (!setequal(cat_trovate, CATS)) {
  stop("Le categorie trovate non coincidono con quelle previste nella tesi: controllare le tabelle LaTeX")
}
riorder <- function(m, cats) { out <- matrix(0, nrow(m), length(cats), dimnames = list(rownames(m), cats))
  comuni <- intersect(cats, colnames(m)); out[, comuni] <- m[, comuni]; out }
RRA <- riorder(RES$rra, CATS); CNT <- riorder(RES$counts, CATS); CNT0 <- riorder(RES$counts_pre, CATS)
AREA_S <- RES$area


# ==============================================================================
# 3. FLUSSO DEI CAMPIONI, CONFRONTI, ESPORTAZIONI
# ==============================================================================
# lettura di un CSV separato da ";" o da ",": tutto come testo
leggi_auto <- function(f) {
  prima <- readLines(f, n = 1, warn = FALSE, encoding = "UTF-8")
  car <- strsplit(prima, "")[[1]]
  sep <- if (sum(car == ";") > sum(car == ",")) ";" else ","
  x <- read_delim(f, delim = sep, col_types = cols(.default = col_character()), locale = locale(encoding = "UTF-8"),
                  trim_ws = TRUE, progress = FALSE, show_col_types = FALSE)
  names(x) <- sub("^\ufeff", "", names(x))
  x
}
numero <- function(x) as.numeric(gsub(",", ".", x, fixed = TRUE))   # accetta 0.97 e 0,97
if (file.exists(FILE_P130)) {
  p130 <- leggi_auto(FILE_P130)
  flusso <- p130 %>% dplyr::select(Sample, Area, Group, Status) %>%
    left_join(RES$totals %>% dplyr::select(-Area), by = "Sample") %>%
    mutate(status = if_else(is.na(status), "No sequence reads", status))
  if (!all(RES$totals$Sample %in% p130$Sample)) stop("Campioni con reads assenti dall'elenco dei 130 selezionati")
} else {
  LOG("File ", FILE_P130, " assente: flusso dei campioni dalle sole tabelle delle varianti")
  flusso <- RES$totals %>% mutate(Group = NA_character_, Status = status)
}
flusso$flag <- ""
flusso$flag[flusso$status == "Analysed" & flusso$non_target_carnivore > 0] <-
  "Non-target carnivore reads present (excluded from denominator)"
flusso$flag[flusso$Sample == "HRV018"] <-
  "Abundant Vulpes vulpes reads; retained in the main analysis, excluded in a sensitivity analysis"
write_csv(flusso, file.path(DIR_DAT, "flusso_campioni.csv"))
riep_flusso <- flusso %>% count(Area, status) %>% pivot_wider(names_from = status, values_from = n, values_fill = 0)
print(riep_flusso); write_csv(riep_flusso, file.path(DIR_DAT, "flusso_campioni_riepilogo.csv"))
if (any(flusso$status == "Analysed" & flusso$Status != "Analysed", na.rm = TRUE) ||
    any(flusso$status != "Analysed" & flusso$Status == "Analysed", na.rm = TRUE)) {
  segnala("Lo stato di alcuni campioni differisce dal file dei 130 punti: controllare")
}
den_chk <- flusso %>% filter(status != "No sequence reads") %>%
  mutate(incl_total = pass_depth_total & prey > 0, incl_prey = pass_depth_prey,
         incl_noCanisVulpes = pass_depth_noCanisVulpes & prey > 0)
n_diff_den <- sum(den_chk$incl_total != den_chk$incl_prey | den_chk$incl_total != den_chk$incl_noCanisVulpes)
LOG("Campioni la cui inclusione cambia con il denominatore del filtro di 2.000 reads: ", n_diff_den)
put("depth.differences", n_diff_den)
write_csv(den_chk, file.path(DIR_DAT, "verifica_denominatore_filtro_2000.csv"))

# confronto con le matrici finali precedenti (nessuna modifica agli originali)
confronta <- function(file, sep, area) {
  prev <- read_delim(file, delim = sep, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  prev <- as.data.frame(prev); rownames(prev) <- prev$scientific_name
  samp <- setdiff(names(prev), c("scientific_name", "rank"))
  nuovi <- names(AREA_S)[AREA_S == area]
  if (!setequal(samp, nuovi)) segnala("Insieme di campioni diverso dalla matrice precedente: ", area)
  diff <- list()
  for (t in union(rownames(prev), CATS)) for (s in samp) {
    a <- if (t %in% rownames(prev)) as.integer(prev[t, s]) else 0L
    b <- if (t %in% CATS) as.integer(CNT[s, t]) else 0L
    if (a != b) diff[[length(diff) + 1]] <- tibble(Area = area, Sample = s, category = t, previous = a, new = b)
  }
  bind_rows(diff)
}
diff_prev <- bind_rows(confronta(FILE_SLO, ";", "Slovenia"), confronta(FILE_CRO, ";", "Croatia"))
if (nrow(diff_prev) > 0) LOG("Differenze dai FILE DEFINITIVO: ", nrow(diff_prev),
                             " celle (elenco in output_R/dati_derivati/differenze_da_matrici_precedenti.csv)")
print(diff_prev); write_csv(diff_prev, file.path(DIR_DAT, "differenze_da_matrici_precedenti.csv"))

esporta <- function(m, nome) {
  out <- tibble(Sample = rownames(m), Area = AREA_S[rownames(m)]) %>% bind_cols(as_tibble(m))
  write_csv(out, file.path(DIR_DAT, nome))
}
esporta(CNT0, "matrice_conteggi_pre_soglia.csv"); esporta(CNT, "matrice_conteggi_finale.csv")
esporta(RRA, "matrice_RRA_percentuale.csv"); esporta((CNT > 0) * 1L, "matrice_presenza_assenza.csv")
write_csv(RES$removed, file.path(DIR_DAT, "rilevamenti_rimossi_soglia_1pct.csv"))

# riepilogo delle assegnazioni tassonomiche
riep_ass <- LUNGO %>% group_by(Area, final_name, role) %>%
  summarise(n_variants = n_distinct(motu_id), reads_run_all_samples = sum(reads), .groups = "drop")
riep_ass$reads_final_matrix <- vapply(seq_len(nrow(riep_ass)), function(i) {
  t <- riep_ass$final_name[i]; a <- riep_ass$Area[i]
  if (t %in% CATS) sum(CNT[AREA_S == a, t]) else 0 }, numeric(1))
riep_ass$n_samples_final <- vapply(seq_len(nrow(riep_ass)), function(i) {
  t <- riep_ass$final_name[i]; a <- riep_ass$Area[i]
  if (t %in% CATS) sum(CNT[AREA_S == a, t] > 0) else 0 }, numeric(1))
write_csv(riep_ass, file.path(DIR_DAT, "riepilogo_assegnazioni_tassonomiche.csv"))

# numeri del flusso e dei reads
for (a in AREE) {
  t <- RES$totals[RES$totals$Area == a, ]
  put(paste0("flow.withreads.", ab(a)), nrow(t))
  put(paste0("flow.below2000.", ab(a)), sum(t$status == "Below 2,000 total reads"))
  put(paste0("flow.noprey.", ab(a)), sum(t$status == "No valid prey reads"))
  put(paste0("reads.run.", ab(a)), sum(t$total_reads))
  put(paste0("reads.predator.", ab(a)), sum(t$predator))
  put(paste0("reads.fox.", ab(a)), sum(t$non_target_carnivore))
  s <- names(AREA_S)[AREA_S == a]
  rf <- rowSums(CNT[s, , drop = FALSE])
  put(paste0("reads.final.", ab(a)), sum(rf)); put(paste0("reads.pre.", ab(a)), sum(RES$den_pre[s]))
  put(paste0("reads.median.", ab(a)), median(rf)); put(paste0("reads.min.", ab(a)), min(rf))
  put(paste0("reads.max.", ab(a)), max(rf)); put(paste0("reads.below10k.pct.", ab(a)), 100 * mean(rf < 10000))
}
put("reads.final.tot", sum(CNT)); put("thr.removed.n", nrow(RES$removed))
quota_ret <- CNT / RES$den_pre[rownames(CNT)]
for (a in AREE) {
  q <- quota_ret[AREA_S == a, , drop = FALSE]; q[q == 0] <- NA
  w <- which(q == min(q, na.rm = TRUE), arr.ind = TRUE)[1, ]
  put(paste0("thr.minretained.pct.", ab(a)), 100 * min(q, na.rm = TRUE))
  put(paste0("thr.minretained.where.", ab(a)), paste(rownames(q)[w[1]], colnames(q)[w[2]]))
}
# pavimento di rilevamento delle due tabelle: quota di ogni occorrenza di una variante
# di preda sul totale delle prede valide del campione (campioni analizzati), e
# rilevamenti trattenuti dopo la soglia dell'1% con quota inferiore al 5%
occ <- LUNGO %>% filter(reads > 0, role == "prey") %>%
  inner_join(RES$totals %>% filter(status == "Analysed") %>% dplyr::select(Area, Sample, prey),
             by = c("Area", "Sample")) %>%
  mutate(quota = reads / prey)
for (a in AREE) {
  o <- occ[occ$Area == a, ]
  put(paste0("floor.", ab(a), ".occ"), nrow(o)); put(paste0("floor.", ab(a), ".below5"), sum(o$quota < 0.05))
  put(paste0("floor.", ab(a), ".minpct"), 100 * min(o$quota))
  q <- quota_ret[AREA_S == a, , drop = FALSE]
  put(paste0("floor.", ab(a), ".ret1to5"), sum(q > 0 & q < 0.05))
}
put("reads.table.min.slo", min(RES$totals$total_reads[RES$totals$Area == "Slovenia" & RES$totals$total_reads > 0]))


# ==============================================================================
# 4. METADATI, COVARIATE E VARIABILI DERIVATE (join per identificativo)
# ==============================================================================
ind  <- leggi_auto(FILE_INDIVIDUI)
pts  <- leggi_auto(FILE_PUNTI) %>%
  mutate(Year = as.integer(Year), Latitude = numero(Latitude), Longitude = numero(Longitude))
if (any(is.na(pts$Latitude) | is.na(pts$Longitude))) stop("Coordinate non leggibili in ", FILE_PUNTI)
env  <- read_delim(FILE_AMBIENTE, delim = ";", locale = locale(decimal_mark = ","), show_col_types = FALSE)
grp  <- read_csv(FILE_GRUPPI, show_col_types = FALSE)

DF <- tibble(Sample = rownames(RRA), Area = unname(AREA_S[rownames(RRA)])) %>% bind_cols(as_tibble(RRA))
DF$reads_final <- rowSums(CNT[DF$Sample, ]); DF$reads_prethreshold <- RES$den_pre[DF$Sample]
unisci <- function(x, y, nome) {
  mancano <- setdiff(x$Sample, y$Sample)
  if (length(mancano) > 0) stop("Campioni assenti da ", nome, ": ", paste(mancano, collapse = ", "))
  if (anyDuplicated(y$Sample)) stop("Campioni duplicati in ", nome)
  left_join(x, y, by = "Sample")
}
DF <- unisci(DF, ind %>% dplyr::select(Sample, Individual, Pack, GeneticSex, GenotypeQI), "individui")
DF <- unisci(DF, pts %>% dplyr::select(Sample, Date, Year, Latitude, Longitude), "punti")
DF <- unisci(DF, env %>% dplyr::select(Sample, ELEV_mean, NDVI_mean, NDVI_sd), "ambiente")
DF <- unisci(DF, grp %>% dplyr::select(Sample, Study_area, Group, Region_cluster_comparison), "gruppi")
stopifnot(nrow(DF) == nrow(RRA), !any(is.na(DF$ELEV_mean)), !any(is.na(DF$NDVI_sd)))
DF <- DF %>% mutate(
  Area = factor(Area, levels = AREE),
  cluster_id = if_else(is.na(Individual), paste0("unassigned_", Sample), Individual),
  ELEV_mean_z = as.numeric(scale(ELEV_mean)), NDVI_mean_z = as.numeric(scale(NDVI_mean)),
  NDVI_sd_z = as.numeric(scale(NDVI_sd)), Croatia = as.integer(Area == "Croatia")) %>%
  arrange(Area, Sample)

antenati <- function(c) { a <- character(0); while (c %in% names(PARENT)) { c <- PARENT[[c]]; a <- c(a, c) }; a }
n_prey_min <- function(det) { anc <- unique(unlist(lapply(det, antenati))); sum(!det %in% anc) }
lca <- function(det) {
  linea <- function(c) c(c, antenati(c))
  comuni <- Reduce(intersect, lapply(det, linea))
  if (length(comuni) == 0) return(NA_character_)
  l0 <- linea(det[1]); l0[l0 %in% comuni][1]
}
det_list <- lapply(seq_len(nrow(DF)), function(i) CATS[as.numeric(unlist(DF[i, CATS])) > 0])
DF$n_assign <- vapply(det_list, length, integer(1))
DF$n_prey_min <- vapply(det_list, n_prey_min, integer(1))
DF$mixed_prey <- as.integer(DF$n_prey_min >= 2)
DF$mixed_assign <- as.integer(DF$n_assign >= 2)
DF$diet_type <- vapply(seq_along(det_list), function(i) {
  if (DF$n_prey_min[i] >= 2) return("Mixed (more than one prey taxon)")
  top <- lca(det_list[[i]])
  if (top %in% DOMESTIC) "Domestic livestock only" else paste(top, "only")
}, character(1))
DF$dominant <- CATS[max.col(as.matrix(DF[, CATS]), ties.method = "first")]
write_csv(DF, file.path(DIR_DAT, "campioni_analitici_103.csv"))

nA <- table(DF$Area)
put("n.slo", nA[["Slovenia"]]); put("n.cro", nA[["Croatia"]]); put("n.tot", nrow(DF))
put("n.ind.slo", n_distinct(DF$Individual[DF$Area == "Slovenia"], na.rm = TRUE))
put("n.ind.cro", n_distinct(DF$Individual[DF$Area == "Croatia"], na.rm = TRUE))
put("n.units.ind", n_distinct(DF$cluster_id))
# campioni per individuo (descrizione del dataset e della dipendenza tra campioni)
cnt_ind <- table(DF$Individual[!is.na(DF$Individual)])
stopifnot(all(tapply(as.character(DF$Area[!is.na(DF$Individual)]), DF$Individual[!is.na(DF$Individual)],
                     function(x) length(unique(x))) == 1))            # ogni individuo in una sola area
put("n.ind.tot", length(cnt_ind)); put("n.ind.none", sum(is.na(DF$Individual)))
put("ind.single", sum(cnt_ind == 1)); put("ind.maxsamples", max(cnt_ind))
for (a in AREE) {
  cnt_a <- table(DF$Individual[DF$Area == a & !is.na(DF$Individual)])
  put(paste0("n.assigned.", ab(a)), sum(cnt_a)); put(paste0("ind.single.", ab(a)), sum(cnt_a == 1))
  put(paste0("ind.max.", ab(a)), max(cnt_a))
}
# composizione per anno e mese di raccolta; le date sono nel formato gg-mmm-aa con
# abbreviazioni italiane dei mesi (es. 22-giu-19), interpretate senza dipendere dal locale
MESI_IT <- c(gen = 1, feb = 2, mar = 3, apr = 4, mag = 5, giu = 6, lug = 7, ago = 8, set = 9, ott = 10, nov = 11, dic = 12)
DF$Month <- unname(MESI_IT[tolower(sub("^[0-9]+-([[:alpha:]]+)-[0-9]+$", "\\1", DF$Date))])
if (any(is.na(DF$Month))) stop("Date non interpretate: ", paste(DF$Sample[is.na(DF$Month)], collapse = ", "))
for (y in sort(unique(DF$Year))) for (a in AREE) put(sprintf("year.%s.%d", ab(a), as.integer(y)), sum(DF$Area == a & DF$Year == y))
put("cro.marapr", sum(DF$Area == "Croatia" & DF$Month %in% 3:4))
put("slo.sepjan", sum(DF$Area == "Slovenia" & DF$Month %in% c(9:12, 1)))
put("cro.north45", sum(DF$Area == "Croatia" & DF$Latitude >= 45))
LOG("Campioni analizzati: ", nrow(DF), " (", nA[["Slovenia"]], " + ", nA[["Croatia"]], "); unita' individuali: ", n_distinct(DF$cluster_id))


# ==============================================================================
# 5. DESCRITTORI: FOO, RRA, RICCHEZZA, TIPI DI DIETA
#    FOO = % di campioni positivi dopo i filtri; RRA di area = media aritmetica
#    delle proporzioni dei singoli campioni, zeri inclusi (non reads cumulati).
# ==============================================================================
foo_rra <- list()
for (a in AREE) {
  s <- DF[DF$Area == a, ]
  for (c in CATS) {
    det <- sum(s[[c]] > 0)
    put(paste0("det.", ab(a), ".", CODE[[c]]), det)
    put(paste0("foo.", ab(a), ".", CODE[[c]]), 100 * det / nrow(s))
    put(paste0("rra.", ab(a), ".", CODE[[c]]), mean(s[[c]]))
    foo_rra[[length(foo_rra) + 1]] <- tibble(Area = a, category = c, detections = det,
                                             FOO = 100 * det / nrow(s), RRA = mean(s[[c]]))
  }
  wild <- rowSums(s[, c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa", "Rupicapra rupicapra")])
  put(paste0("rra.", ab(a), ".wildresolved"), mean(wild))
  put(paste0("rra.", ab(a), ".domestic"), mean(rowSums(s[, DOMESTIC])))
  put(paste0("rra.", ab(a), ".domestic_plus_caprinae"), mean(rowSums(s[, c(DOMESTIC, "Caprinae")])))
  put(paste0("rra.", ab(a), ".wild_plus_caprinae"), mean(wild + s$Caprinae))
  put(paste0("foo.", ab(a), ".domestic"), 100 * mean(rowSums(s[, DOMESTIC]) > 0))
  put(paste0("det.", ab(a), ".domestic"), sum(rowSums(s[, DOMESTIC]) > 0))
  put(paste0("rra.", ab(a), ".three"), mean(rowSums(s[, c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa")])))
  for (c in CATS) put(paste0("only.", ab(a), ".", CODE[[c]]), sum(s[[c]] > 99.999999))
  cnt_a <- CNT[AREA_S == a, CATS, drop = FALSE]
  pooled <- 100 * colSums(cnt_a) / sum(cnt_a)
  for (c in CATS) put(paste0("pooled.", ab(a), ".", CODE[[c]]), pooled[[c]])
}
foo_rra <- bind_rows(foo_rra)
write_csv(foo_rra, file.path(DIR_TAB, "FOO_RRA_per_area.csv"))
put("n.categories", length(CATS))
put("n.categories.slo", sum(colSums(DF[DF$Area == "Slovenia", CATS]) > 0))
put("n.categories.cro", sum(colSums(DF[DF$Area == "Croatia", CATS]) > 0))
RANK <- c("Cervus elaphus" = "species", "Capreolus capreolus" = "species", "Sus scrofa" = "species",
          "Rupicapra rupicapra" = "species", "Ovis aries" = "species", "Caprinae" = "subfamily",
          "Cervinae" = "subfamily", "Bos" = "genus", "Capra" = "genus", "Ovis" = "genus", "Lepus" = "genus")
for (a in AREE) {
  pres <- CATS[colSums(DF[DF$Area == a, CATS]) > 0]
  for (r in c("species", "genus", "subfamily")) put(paste0("res.", ab(a), ".", r), sum(RANK[pres] == r))
}

for (a in c(AREE, "Total")) {
  s <- if (a == "Total") DF else DF[DF$Area == a, ]; k <- ab(a)
  for (v in 1:3) {
    put(paste0("nassign.", v, ".", k), sum(s$n_assign == v))
    put(paste0("nprey.", v, ".", k), sum(s$n_prey_min == v))
  }
  put(paste0("mix.prey.", k), sum(s$mixed_prey)); put(paste0("mix.assign.", k), sum(s$mixed_assign))
  put(paste0("single.prey.", k), sum(s$n_prey_min == 1)); put(paste0("single.assign.", k), sum(s$n_assign == 1))
  put(paste0("single.prey.pct.", k), 100 * mean(s$n_prey_min == 1))
  put(paste0("nassign.max.", k), max(s$n_assign)); put(paste0("nprey.max.", k), max(s$n_prey_min))
  put(paste0("nassign.mean.", k), mean(s$n_assign)); put(paste0("nprey.mean.", k), mean(s$n_prey_min))
}
put("mix.differ.samples", paste(DF$Sample[DF$mixed_prey != DF$mixed_assign], collapse = ", "))
annid <- DF[DF$n_assign != DF$n_prey_min, ]
put("nested.n", nrow(annid))
for (c in CATS) put(paste0("only.tot.", CODE[[c]]), sum(DF[[c]] > 99.999999))
put("nested.samples", paste(vapply(seq_len(nrow(annid)), function(i)
  paste0(annid$Sample[i], " (", paste(CATS[as.numeric(unlist(annid[i, CATS])) > 0], collapse = ", "), ")"),
  character(1)), collapse = "; "))
tipi <- DF %>% count(Area, diet_type) %>% pivot_wider(names_from = Area, values_from = n, values_fill = 0)
write_csv(tipi, file.path(DIR_TAB, "tipi_di_dieta.csv"))
for (i in seq_len(nrow(tipi))) for (a in AREE) {
  tk <- tolower(strsplit(tipi$diet_type[i], " ")[[1]][1])
  put(paste0("type.", ab(a), ".", tk), tipi[[a]][i])
  put(paste0("typepct.", ab(a), ".", tk), 100 * tipi[[a]][i] / sum(DF$Area == a))
}
foglie <- function(det) { anc <- unique(unlist(lapply(det, antenati))); det[!det %in% anc] }
mx <- DF[DF$mixed_prey == 1, ]
minore <- vapply(seq_len(nrow(mx)), function(i) {
  det <- CATS[as.numeric(unlist(mx[i, CATS])) > 0]; min(as.numeric(unlist(mx[i, foglie(det)]))) }, numeric(1))
put("mix.minor.median", median(minore)); put("mix.minor.min", min(minore))
put("mix.minor.max", max(minore)); put("mix.minor.below10", sum(minore < 10))
for (a in AREE) {
  s <- DF[DF$Area == a, ]
  put(paste0("reads.median.single.", ab(a)), median(s$reads_final[s$mixed_prey == 0]))
  put(paste0("reads.median.mixed.", ab(a)), median(s$reads_final[s$mixed_prey == 1]))
}


# ==============================================================================
# 6. AMPIEZZA E SOVRAPPOSIZIONE DI NICCHIA (categorie di assegnazione)
# ==============================================================================
levins <- function(p, n) { p <- p / sum(p); B <- 1 / sum(p^2); c(B = B, BA = (B - 1) / (n - 1)) }
pianka <- function(p, q) { p <- p / sum(p); q <- q / sum(q); sum(p * q) / sqrt(sum(p^2) * sum(q^2)) }
Pm <- DF %>% group_by(Area) %>% summarise(across(all_of(CATS), mean), .groups = "drop")
Pmat <- as.matrix(Pm[, CATS]); rownames(Pmat) <- as.character(Pm$Area)
for (a in AREE) {
  l <- levins(Pmat[a, ], length(CATS)); npres <- sum(Pmat[a, ] > 0)
  put(paste0("levins.B.", ab(a)), l[["B"]]); put(paste0("levins.BA.", ab(a)), l[["BA"]])
  put(paste0("levins.npres.", ab(a)), npres); put(paste0("levins.BAarea.", ab(a)), (l[["B"]] - 1) / (npres - 1))
}
put("pianka", pianka(Pmat["Slovenia", ], Pmat["Croatia", ]))
radice <- function(c) { while (c %in% names(PARENT)) c <- PARENT[[c]]; c }
ROOTS <- intersect(CATS, unique(vapply(CATS, radice, character(1))))
collassa <- function(M) { out <- sapply(ROOTS, function(r) rowSums(M[, CATS[vapply(CATS, radice, character(1)) == r], drop = FALSE])); out }
Pc <- collassa(Pmat)
put("collapsed.n", length(ROOTS)); put("collapsed.groups", paste(ROOTS, collapse = ", "))
for (a in AREE) {
  l <- levins(Pc[a, ], length(ROOTS))
  put(paste0("levins.coll.B.", ab(a)), l[["B"]]); put(paste0("levins.coll.BA.", ab(a)), l[["BA"]])
}
put("pianka.coll", pianka(Pc["Slovenia", ], Pc["Croatia", ]))


# ==============================================================================
# 7. ANALISI MULTIVARIATE (matrice delle 11 categorie di assegnazione)
#    Hellinger sulle proporzioni RRA, dissimilarita' di Bray-Curtis.
#    Permutazioni libere dei campioni: presuppongono campioni indipendenti;
#    la dipendenza fra campioni dello stesso individuo e' valutata al punto 7.5.
# ==============================================================================
H  <- decostand(as.matrix(DF[, CATS]), method = "hellinger")
DM <- vegdist(H, method = "bray")
set.seed(SEED)
perm <- adonis2(DM ~ Area, data = DF, permutations = how(nperm = N_PERM), by = "terms")
print(perm)
put("perm.R2", perm$R2[1]); put("perm.F", perm$F[1]); put("perm.p", perm$`Pr(>F)`[1])
put("perm.SS", perm$SumOfSqs[1]); put("perm.SSres", perm$SumOfSqs[2]); put("perm.SST", sum(perm$SumOfSqs[1:2]))
put("perm.df1", perm$Df[1]); put("perm.df2", perm$Df[2])

neg <- FALSE
bd <- withCallingHandlers(betadisper(DM, group = DF$Area, type = "median"),
                          warning = function(w) { if (grepl("negative", conditionMessage(w))) neg <<- TRUE
                                                  invokeRestart("muffleWarning") })
set.seed(SEED)
bdt <- permutest(bd, permutations = how(nperm = N_PERM))
print(bdt)
put("bd.F", bdt$tab$F[1]); put("bd.p", bdt$tab$`Pr(>F)`[1]); put("bd.df2", bdt$tab$Df[2])
put("bd.dist.slo", bd$group.distances[["Slovenia"]]); put("bd.dist.cro", bd$group.distances[["Croatia"]])
put("bd.negwarning", neg)

set.seed(SEED)
sim <- simper(H, group = as.character(DF$Area), permutations = how(nperm = N_PERM))
sims <- summary(sim)[[1]]
sims$category <- rownames(sims); sims$contribution_pct <- 100 * sims$average / sum(sims$average)
print(sims); write_csv(sims, file.path(DIR_TAB, "SIMPER.csv"))
put("simper.overall", sum(sims$average))
for (i in seq_len(nrow(sims))) {
  put(paste0("simper.", CODE[[sims$category[i]]], ".pct"), sims$contribution_pct[i])
  put(paste0("simper.", CODE[[sims$category[i]]], ".p"), sims$p[i])
}
put("simper.top4", paste(sims$category[1:4], collapse = ", ")); put("simper.top4.pct", sum(sims$contribution_pct[1:4]))
put("simper.rest.maxpct", max(sims$contribution_pct[-(1:4)])); put("simper.rest.pmin", min(sims$p[-(1:4)]))
put("simper.rest.pmax", max(sims$p[-(1:4)])); put("simper.rest.n", nrow(sims) - 4)

pc <- cmdscale(DM, k = 2, eig = TRUE)
eigpos <- pc$eig[pc$eig > 1e-10]
put("pcoa.ax1", 100 * pc$eig[1] / sum(eigpos)); put("pcoa.ax2", 100 * pc$eig[2] / sum(eigpos))
put("pcoa.plane", 100 * sum(pc$eig[1:2]) / sum(eigpos)); put("pcoa.rest", 100 - 100 * sum(pc$eig[1:2]) / sum(eigpos))
put("pcoa.npos", length(eigpos))
put("pcoa.neg", sum(pc$eig < -1e-10))
put("pcoa.profiles", nrow(unique(round(as.matrix(DF[, CATS]), 8))))
xy <- round(pc$points, 8)
put("pcoa.positions", nrow(unique(xy)))
pos_tab <- as.data.frame(xy) %>% count(V1, V2, name = "n")
put("pcoa.maxstack", max(pos_tab$n))
pc_all <- cmdscale(DM, k = length(eigpos), eig = TRUE)
put("pcoa.positions.allaxes", nrow(unique(round(pc_all$points, 8))))

# --- 7.5 dipendenza: analisi a livello di individuo (profilo medio per individuo)
IND <- DF %>% group_by(cluster_id) %>%
  summarise(Area = first(Area), across(all_of(CATS), mean),
            ELEV_mean = mean(ELEV_mean), NDVI_mean = mean(NDVI_mean), NDVI_sd = mean(NDVI_sd), .groups = "drop")
put("ind.n", nrow(IND)); put("ind.n.slo", sum(IND$Area == "Slovenia")); put("ind.n.cro", sum(IND$Area == "Croatia"))
DMi <- vegdist(decostand(as.matrix(IND[, CATS]), "hellinger"), "bray")
set.seed(SEED)
permi <- adonis2(DMi ~ Area, data = IND, permutations = how(nperm = N_PERM), by = "terms")
print(permi)
put("ind.perm.R2", permi$R2[1]); put("ind.perm.F", permi$F[1]); put("ind.perm.p", permi$`Pr(>F)`[1]); put("ind.perm.df2", permi$Df[2])
bdi <- suppressWarnings(betadisper(DMi, group = IND$Area, type = "median"))
set.seed(SEED)
bdit <- permutest(bdi, permutations = how(nperm = N_PERM))
put("ind.bd.F", bdit$tab$F[1]); put("ind.bd.p", bdit$tab$`Pr(>F)`[1])
put("ind.bd.dist.slo", bdi$group.distances[["Slovenia"]]); put("ind.bd.dist.cro", bdi$group.distances[["Croatia"]])


# ==============================================================================
# 8. MODELLI LOGISTICI
#    Coefficienti, IC 95% da profilo di verosimiglianza, p di Wald (summary).
#    Dipendenza: errori robusti per individuo (vcovCL, HC0, correzione G/(G-1))
#    e, per la dieta mista, GLMM con intercetta casuale per individuo.
# ==============================================================================
Z975 <- qnorm(0.975)
blocco_glm <- function(prefix, formula, dati, termine, cluster = NULL, glmm = FALSE) {
  m <- glm(formula, data = dati, family = binomial(link = "logit"))
  co <- summary(m)$coefficients
  ci <- suppressMessages(confint(m))
  for (t in rownames(co)[-1]) {
    pre <- if (t == termine) prefix else paste0(prefix, ".", t)
    put(paste0(pre, ".OR"), exp(co[t, "Estimate"])); put(paste0(pre, ".OR_lo"), exp(ci[t, 1]))
    put(paste0(pre, ".OR_hi"), exp(ci[t, 2])); put(paste0(pre, ".p"), co[t, "Pr(>|z|)"])
    if (t == termine) {
      put(paste0(pre, ".estimate"), co[t, "Estimate"]); put(paste0(pre, ".se"), co[t, "Std. Error"])
      put(paste0(pre, ".z"), co[t, "z value"])
    }
  }
  put(paste0(prefix, ".n"), nrow(dati)); put(paste0(prefix, ".events"), sum(m$y))
  if (!is.null(cluster)) {
    V <- sandwich::vcovCL(m, cluster = cluster, type = "HC0")
    se <- sqrt(diag(V))[termine]; b <- coef(m)[termine]
    put(paste0(prefix, ".cr.se"), se); put(paste0(prefix, ".cr.p"), 2 * pnorm(-abs(b / se)))
    put(paste0(prefix, ".cr.OR_lo"), exp(b - Z975 * se)); put(paste0(prefix, ".cr.OR_hi"), exp(b + Z975 * se))
    put(paste0(prefix, ".cr.G"), length(unique(cluster)))
    g <- tapply(m$y, cluster, function(v) c(length(v), sum(v)))
    g <- do.call(rbind, g); rep <- g[g[, 1] > 1, , drop = FALSE]
    put(paste0(prefix, ".rep.ind"), nrow(rep)); put(paste0(prefix, ".rep.discordant"), sum(rep[, 2] > 0 & rep[, 2] < rep[, 1]))
  }
  if (glmm) {
    dati$cluster_id_glmm <- cluster
    f2 <- update(formula, . ~ . + (1 | cluster_id_glmm))
    mm <- tryCatch(withCallingHandlers(
      glmer(f2, data = dati, family = binomial, control = glmerControl(optimizer = "bobyqa")),
      warning = function(w) { LOG("glmer ", prefix, ": ", conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) { LOG("glmer ", prefix, " non stimabile: ", conditionMessage(e)); NULL })
    if (!is.null(mm)) {
      b <- fixef(mm)[termine]; se <- sqrt(diag(as.matrix(vcov(mm))))[which(names(fixef(mm)) == termine)]
      put(paste0(prefix, ".glmm.OR"), exp(b)); put(paste0(prefix, ".glmm.OR_lo"), exp(b - Z975 * se))
      put(paste0(prefix, ".glmm.OR_hi"), exp(b + Z975 * se)); put(paste0(prefix, ".glmm.p"), 2 * pnorm(-abs(b / se)))
      put(paste0(prefix, ".glmm.sigma"), attr(VarCorr(mm)$cluster_id_glmm, "stddev")[[1]])
      put(paste0(prefix, ".glmm.singular"), isSingular(mm)); put(paste0(prefix, ".glmm.nclusters"), length(unique(cluster)))
    }
  }
  m
}

m_hhh <- blocco_glm("hhh", mixed_prey ~ NDVI_sd_z + ELEV_mean_z + Croatia, DF, "NDVI_sd_z",
                    cluster = DF$cluster_id, glmm = TRUE)
print(summary(m_hhh))
# modello con interazione stimato esplicitamente (update() rivaluterebbe la
# chiamata fuori dalla funzione blocco_glm)
m_hhh0 <- glm(mixed_prey ~ NDVI_sd_z + ELEV_mean_z + Croatia, data = DF, family = binomial)
m_int  <- glm(mixed_prey ~ NDVI_sd_z + ELEV_mean_z + Croatia + NDVI_sd_z:Croatia, data = DF, family = binomial)
lr <- anova(m_hhh0, m_int, test = "Chisq")
put("hhh.int.chi", lr$Deviance[2]); put("hhh.int.p", lr$`Pr(>Chi)`[2])
m_hhha <- blocco_glm("hhha", mixed_assign ~ NDVI_sd_z + ELEV_mean_z + Croatia, DF, "NDVI_sd_z",
                     cluster = DF$cluster_id, glmm = TRUE)

TRE <- c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa")
DOM <- DF[rowSums(DF[, TRE] > 0) > 0, ]
DOM$dom3 <- TRE[max.col(as.matrix(DOM[, TRE]), ties.method = "first")]
DOM$boar <- as.integer(DOM$dom3 == "Sus scrofa")
put("dom.n", nrow(DOM)); for (c in TRE) put(paste0("dom.", CODE[[c]]), sum(DOM$dom3 == c))
m_l1 <- blocco_glm("l1", boar ~ ELEV_mean_z + NDVI_mean_z, DOM, "ELEV_mean_z", cluster = DOM$cluster_id, glmm = TRUE)
print(summary(m_l1))
blocco_glm("l1.no2", boar ~ ELEV_mean_z + NDVI_mean_z, DOM[!DOM$Sample %in% c("HRV027", "HRV05X"), ], "ELEV_mean_z")
blocco_glm("l1.cro", boar ~ ELEV_mean_z + NDVI_mean_z, DOM[DOM$Area == "Croatia", ], "ELEV_mean_z")
CERV <- DOM[DOM$dom3 != "Sus scrofa", ]; CERV$red <- as.integer(CERV$dom3 == "Cervus elaphus")
m_l2 <- blocco_glm("l2", red ~ ELEV_mean_z + NDVI_mean_z, CERV, "ELEV_mean_z", cluster = CERV$cluster_id, glmm = TRUE)
print(summary(m_l2))
low <- DOM[order(DOM$ELEV_mean), ][1:3, ]
put("l1.lowest", paste(sprintf("%s %.0f m %s", low$Sample, low$ELEV_mean, low$dom3), collapse = "; "))
put("l1.min.cervid.elev", min(DOM$ELEV_mean[DOM$boar == 0]))


# ==============================================================================
# 9. RDA (matrice aggregata: 8 categorie; nessuna categoria esclusa, nessuna
#    rinormalizzazione: i taxa domestici sono sommati, Lepus e Cervinae restano)
# ==============================================================================
RD <- as.data.frame(DF[, CATS]); RD$Domestic <- rowSums(RD[, DOMESTIC])
stopifnot(all(abs(rowSums(RD[, RDA_CATS]) - 100) < 1e-8))
Y_rda <- decostand(as.matrix(RD[, RDA_CATS]), "hellinger")
X_rda <- as.data.frame(DF[, c("ELEV_mean", "NDVI_mean", "NDVI_sd")])
m_rda <- rda(Y_rda ~ ELEV_mean + NDVI_mean + NDVI_sd, data = X_rda)
print(m_rda)
set.seed(SEED); a_rda <- anova(m_rda, permutations = how(nperm = N_PERM))
set.seed(SEED); a_mar <- anova(m_rda, by = "margin", permutations = how(nperm = N_PERM))
print(a_rda); print(a_mar)
r2 <- RsquareAdj(m_rda)
put("rda.R2", r2$r.squared); put("rda.R2adj", r2$adj.r.squared)
put("rda.F", a_rda$F[1]); put("rda.p", a_rda$`Pr(>F)`[1]); put("rda.df2", a_rda$Df[2])
put("rda.ax1", 100 * m_rda$CCA$eig[1] / m_rda$tot.chi); put("rda.ax2", 100 * m_rda$CCA$eig[2] / m_rda$tot.chi)
for (t in c("ELEV_mean", "NDVI_mean", "NDVI_sd")) {
  put(paste0("rda.m.", t, ".F"), a_mar[t, "F"]); put(paste0("rda.m.", t, ".p"), a_mar[t, "Pr(>F)"])
}
v <- vif.cca(m_rda); for (t in names(v)) put(paste0("rda.vif.", t), v[[t]])
put("rda.ncat", length(RDA_CATS))
# livello di individuo: profilo e covariate medi per individuo
RDi <- as.data.frame(IND[, CATS]); RDi$Domestic <- rowSums(RDi[, DOMESTIC])
Yi <- decostand(as.matrix(RDi[, RDA_CATS]), "hellinger")
Xi <- as.data.frame(IND[, c("ELEV_mean", "NDVI_mean", "NDVI_sd")])
m_rdai <- rda(Yi ~ ELEV_mean + NDVI_mean + NDVI_sd, data = Xi)
set.seed(SEED); a_rdai <- anova(m_rdai, permutations = how(nperm = N_PERM))
set.seed(SEED); a_mari <- anova(m_rdai, by = "margin", permutations = how(nperm = N_PERM))
r2i <- RsquareAdj(m_rdai)
put("ind.rda.R2", r2i$r.squared); put("ind.rda.R2adj", r2i$adj.r.squared)
put("ind.rda.F", a_rdai$F[1]); put("ind.rda.p", a_rdai$`Pr(>F)`[1]); put("ind.rda.df2", a_rdai$Df[2])
for (t in c("ELEV_mean", "NDVI_mean", "NDVI_sd")) {
  put(paste0("ind.rda.m.", t, ".F"), a_mari[t, "F"]); put(paste0("ind.rda.m.", t, ".p"), a_mari[t, "Pr(>F)"])
}


# ==============================================================================
# 10. GRUPPI GEOGRAFICI E BRANCHI (descrittivi; PERMANOVA esplorative)
# ==============================================================================
ORD_G <- c("Jelovica pack", "Pokljuka pack", "Cerkljansko pack", "Other packs",
           "Zumberak group", "Gorski Kotar", "Southern Croatia")
TW <- DF %>% mutate(dom_any = rowSums(across(all_of(DOMESTIC))) > 0) %>%
  group_by(Group) %>%
  summarise(Area = first(Area), n = n(), individuals = n_distinct(Individual, na.rm = TRUE),
            no_individual = sum(is.na(Individual)),
            roe = sum(`Capreolus capreolus` > 0), red = sum(`Cervus elaphus` > 0),
            boar = sum(`Sus scrofa` > 0), chamois = sum(`Rupicapra rupicapra` > 0),
            caprinae = sum(Caprinae > 0), domestic = sum(dom_any), lepus = sum(Lepus > 0),
            cervinae = sum(Cervinae > 0), mixed_prey = sum(mixed_prey), mixed_assign = sum(mixed_assign),
            NDVI_sd_mean = mean(NDVI_sd), .groups = "drop") %>%
  mutate(Group = factor(Group, levels = ORD_G)) %>% arrange(Group)
print(TW, width = Inf); write_csv(TW, file.path(DIR_TAB, "tabella_gruppi_branchi.csv"))
put("mix.prey.hetgroups", sum(TW$mixed_prey[TW$Group %in% c("Jelovica pack", "Pokljuka pack", "Cerkljansko pack", "Southern Croatia")]))
for (i in seq_len(nrow(TW))) {
  g <- tolower(strsplit(as.character(TW$Group[i]), " ")[[1]][1])
  for (kk in c("n", "individuals", "roe", "red", "boar", "chamois", "caprinae", "domestic", "lepus",
               "cervinae", "mixed_prey", "mixed_assign", "NDVI_sd_mean")) put(paste0("grp.", g, ".", kk), TW[[kk]][i])
}
domg <- table(DF$Group, DF$dominant)
write.csv(as.data.frame.matrix(domg), file.path(DIR_TAB, "preda_dominante_per_gruppo.csv"))
for (g in ORD_G) for (c in colnames(domg)) {
  put(paste0("domgrp.", tolower(strsplit(g, " ")[[1]][1]), ".", CODE[[c]]), domg[g, c])
}
perm_gruppi <- function(dati) {
  d <- vegdist(decostand(as.matrix(dati[, CATS]), "hellinger"), "bray")
  set.seed(SEED); adonis2(d ~ Group, data = dati, permutations = how(nperm = N_PERM))
}
pgc <- perm_gruppi(DF[DF$Area == "Croatia", ]); print(pgc)
put("grp.cro.perm.R2", pgc$R2[1]); put("grp.cro.perm.F", pgc$F[1]); put("grp.cro.perm.p", pgc$`Pr(>F)`[1])
put("grp.cro.perm.df1", pgc$Df[1]); put("grp.cro.perm.df2", pgc$Df[2])
SLP <- DF[DF$Group %in% c("Jelovica pack", "Pokljuka pack", "Cerkljansko pack"), ]
pgs <- perm_gruppi(SLP); print(pgs)
put("grp.slo.perm.R2", pgs$R2[1]); put("grp.slo.perm.F", pgs$F[1]); put("grp.slo.perm.p", pgs$`Pr(>F)`[1])
put("grp.slo.perm.df1", pgs$Df[1]); put("grp.slo.perm.df2", pgs$Df[2]); put("grp.slo.perm.n", nrow(SLP))


# ==============================================================================
# 11. ANALISI DI SENSIBILITA'
#   main        analisi principale (riferimento)
#   hrv018      HRV018 escluso (campione con abbondante DNA di volpe)
#   crocervinae variante croata Cervinae lasciata non risolta
#   bovidae87   variante Bovidae di 87 reads trattenuta come Bovidae non risolto
#   collapsed   categorie annidate unite nella categoria superiore (gruppi esclusivi)
#   thr5        soglia comune del 5% al posto dell'1% (pavimento di rilevamento
#               uguale nelle due tabelle)
# ==============================================================================
SENS <- list()
sens <- function(tag, dati, cats, parent = PARENT) {
  dm <- vegdist(decostand(as.matrix(dati[, cats]), "hellinger"), "bray")
  set.seed(SEED); p <- adonis2(dm ~ Area, data = dati, permutations = how(nperm = N_PERM))
  P_ <- dati %>% group_by(Area) %>% summarise(across(all_of(cats), mean), .groups = "drop")
  Pm_ <- as.matrix(P_[, cats]); rownames(Pm_) <- as.character(P_$Area)
  ant <- function(c) { a <- character(0); while (c %in% names(parent)) { c <- parent[[c]]; a <- c(a, c) }; a }
  npm <- vapply(seq_len(nrow(dati)), function(i) {
    det <- cats[as.numeric(unlist(dati[i, cats])) > 0]; anc <- unique(unlist(lapply(det, ant))); sum(!det %in% anc) }, integer(1))
  dati$mixed_s <- as.integer(npm >= 2)
  m <- glm(mixed_s ~ NDVI_sd_z + ELEV_mean_z + Croatia, data = dati, family = binomial)
  ci <- suppressMessages(confint(m))
  cro <- dati[dati$Area == "Croatia", ]
  r <- list(n_slo = sum(dati$Area == "Slovenia"), n_cro = sum(dati$Area == "Croatia"), n_cat = length(cats),
            perm_R2 = p$R2[1], perm_F = p$F[1], perm_p = p$`Pr(>F)`[1],
            BA_slo = levins(Pm_["Slovenia", ], length(cats))[["BA"]],
            BA_cro = levins(Pm_["Croatia", ], length(cats))[["BA"]],
            pianka = pianka(Pm_["Slovenia", ], Pm_["Croatia", ]), mixed = sum(dati$mixed_s),
            hhh_OR = exp(coef(m)[["NDVI_sd_z"]]), hhh_lo = exp(ci["NDVI_sd_z", 1]), hhh_hi = exp(ci["NDVI_sd_z", 2]),
            hhh_p = summary(m)$coefficients["NDVI_sd_z", "Pr(>|z|)"],
            foo_cro_cerela = if ("Cervus elaphus" %in% cats) 100 * mean(cro$`Cervus elaphus` > 0) else NA_real_,
            rra_cro_cerela = if ("Cervus elaphus" %in% cats) mean(cro$`Cervus elaphus`) else NA_real_)
  for (k in names(r)) put(paste0("sens.", tag, ".", k), r[[k]])
  SENS[[tag]] <<- c(tag = tag, r)
}
covar <- DF %>% dplyr::select(Sample, Area, NDVI_sd_z, ELEV_mean_z, Croatia)
sens("main", DF, CATS)
DFh <- DF[DF$Sample != "HRV018", ] %>%
  mutate(ELEV_mean_z = as.numeric(scale(ELEV_mean)), NDVI_sd_z = as.numeric(scale(NDVI_sd)),
         NDVI_mean_z = as.numeric(scale(NDVI_mean)))
sens("hrv018", DFh, CATS)
RES_b <- pipeline(LUNGO, cervinae_cro_irrisolti = TRUE)
DFb <- covar %>% left_join(tibble(Sample = rownames(RES_b$rra)) %>% bind_cols(as_tibble(riorder(RES_b$rra, CATS))), by = "Sample")
sens("crocervinae", DFb, CATS)
elenco_e <- function(x) if (length(x) <= 1) paste(x, collapse = "") else
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
put("sens.crocervinae.samples", elenco_e(DFb$Sample[DFb$Area == "Croatia" & DFb$Cervinae > 0]))
RES_c <- pipeline(LUNGO, prede_extra = VARIANTE_BOVIDAE_87)
cats_c <- c(CATS, "Bovidae")
DFc <- covar %>% left_join(tibble(Sample = rownames(RES_c$rra)) %>% bind_cols(as_tibble(riorder(RES_c$rra, cats_c))), by = "Sample")
sens("bovidae87", DFc, cats_c, parent = c(PARENT, "Caprinae" = "Bovidae", "Bos" = "Bovidae"))
put("sens.bovidae87.msv17c", DFc$Bovidae[DFc$Sample == "MSV17C"])
put("sens.bovidae87.msv0l8", 100 * RES_c$counts_pre["MSV0L8", "Bovidae"] / sum(RES_c$counts_pre["MSV0L8", ]))
DFd <- covar %>% bind_cols(as_tibble(collassa(as.matrix(DF[, CATS]))))
sens("collapsed", DFd, ROOTS, parent = character(0))
RES_5 <- pipeline(LUNGO, soglia = 0.05)
stopifnot(setequal(rownames(RES_5$rra), DF$Sample))
cats_5 <- intersect(CATS, colnames(RES_5$rra))
DF5 <- covar %>% left_join(tibble(Sample = rownames(RES_5$rra)) %>% bind_cols(as_tibble(riorder(RES_5$rra, CATS))), by = "Sample")
sens("thr5", DF5, cats_5)
for (a in AREE) {
  put(paste0("sens.thr5.removed.", ab(a)), sum(RES_5$removed$Sample %in% names(RES_5$area)[RES_5$area == a]))
  npm5 <- vapply(which(DF5$Area == a), function(i) n_prey_min(CATS[as.numeric(unlist(DF5[i, CATS])) > 0]), integer(1))
  put(paste0("sens.thr5.mixed.", ab(a)), sum(npm5 >= 2))
}
put("sens.assign.mixed", sum(DF$mixed_assign))
SENS_T <- bind_rows(lapply(SENS, as_tibble))
print(SENS_T, width = Inf); write_csv(SENS_T, file.path(DIR_TAB, "analisi_sensibilita.csv"))

# dettagli usati nel testo
cap <- DF[DF$Caprinae > 0, ]
put("caprinae.samples", paste(sprintf("%s %s %s %s", cap$Sample, cap$Group, cap$Individual, cap$Date), collapse = "; "))
cham <- DF[DF$`Rupicapra rupicapra` > 0, ]
put("chamois.samples", paste(sprintf("%s %s", cham$Sample, cham$Group), collapse = "; "))
r136 <- unlist(DF[DF$Sample == "MSV136", CATS])
put("msv136", paste(sprintf("%s %.2f", names(r136)[r136 > 0], r136[r136 > 0]), collapse = "; "))
s4 <- DF[DF$Area == "Croatia", ] %>% arrange(Latitude) %>% head(4)
put("south4", paste(sprintf("%s %.3f", s4$Sample, s4$Latitude), collapse = "; "))
put("cro.south45", sum(DF$Area == "Croatia" & DF$Latitude < 45))


# ==============================================================================
# 12. FIGURE (PDF vettoriale + PNG 300 dpi, larghezza 16 cm come il testo)
#     Generi e specie in corsivo; Caprinae e Cervinae in tondo.
# ==============================================================================
COL_AREA <- c("Slovenia" = "#1F7FB5", "Croatia" = "#C73E1D")
GRP <- tibble(
  key   = c("ROE", "RED", "CERVINAE", "BOAR", "CHAMOIS", "CAPRINAE", "DOMESTIC", "LEPUS"),
  label = c("Roe deer", "Red deer", "Cervinae (unresolved)", "Wild boar", "Northern chamois",
            "Caprinae (unresolved)", "Taxa assigned to domestic livestock", "Lepus sp."),
  fill  = c("#eda100", "#8c510a", "#f3ebe1", "#4a3aa7", "#00a5a8", "#c51b7d", "#7570b3", "#4d9221"),
  edge  = c("white", "white", "#8c510a", "white", "white", "white", "white", "white"))
GRP_CATS <- list(ROE = "Capreolus capreolus", RED = "Cervus elaphus", CERVINAE = "Cervinae",
                 BOAR = "Sus scrofa", CHAMOIS = "Rupicapra rupicapra", CAPRINAE = "Caprinae",
                 DOMESTIC = DOMESTIC, LEPUS = "Lepus")
lab_tax <- function(x) {
  ifelse(x %in% c("Caprinae", "Cervinae"), paste0("plain('", x, "')"),
         ifelse(x == "Domestic", "plain('Domestic taxa')", paste0("italic('", x, "')")))
}
lab_tax_expr <- function(x) parse(text = lab_tax(x))
tema <- theme_minimal(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.line = element_line(linewidth = 0.3, colour = "black"), axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 9), legend.position = "bottom",
        legend.title = element_blank(), strip.text = element_text(face = "bold", size = 8.2))
salva <- function(p, nome, h_cm) {
  dev_pdf <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(DIR_FIG, paste0(nome, ".pdf")), p, width = 16, height = h_cm, units = "cm", device = dev_pdf)
  ggsave(file.path(PROJECT_DIR, paste0(nome, ".png")), p, width = 16, height = h_cm, units = "cm", dpi = 300)
  if (interactive()) print(p)                       # in RStudio: pannello Plots
  LOG("Figura salvata: ", nome, ".png (e .pdf in output_R/figure_pdf)")
}
lab_area_n <- c(Slovenia = paste0("Slovenia (n = ", nA[["Slovenia"]], ")"),
                Croatia = paste0("Croatia (n = ", nA[["Croatia"]], ")"))

# --- 12.1 composizione tassonomica dell'output di sequenziamento (per famiglia)
FAM <- c("Cervus elaphus" = "Cervidae", "Capreolus capreolus" = "Cervidae", "Cervinae" = "Cervidae",
         "Caprinae" = "Bovidae", "Rupicapra rupicapra" = "Bovidae", "Ovis aries" = "Bovidae", "Bos taurus" = "Bovidae",
         "Bovidae" = "Bovidae", "Bos" = "Bovidae", "Capra" = "Bovidae", "Ovis" = "Bovidae", "Sus scrofa" = "Suidae",
         "Vulpes vulpes" = "Canidae", "Gallus gallus" = "Phasianidae", "Lepus" = "Leporidae",
         "Laurasiatheria" = "Not assigned to a family", "Eukaryota" = "Not assigned to a family")
ORD_FAM <- c("Cervidae", "Bovidae", "Suidae", "Canidae", "Leporidae", "Phasianidae", "Not assigned to a family")
COL_FAM <- c("Cervidae" = "#8c510a", "Bovidae" = "#c51b7d", "Suidae" = "#4a3aa7", "Canidae" = "#9a9a9a",
             "Leporidae" = "#4d9221", "Phasianidae" = "#eda100", "Not assigned to a family" = "#d9d9d9")
fam_tab <- LUNGO %>% filter(role != "predator") %>% mutate(family = unname(FAM[final_name])) %>%
  group_by(Area, family) %>% summarise(reads = sum(reads), .groups = "drop") %>%
  group_by(Area) %>% mutate(pct = 100 * reads / sum(reads), tot = sum(reads)) %>% ungroup() %>%
  mutate(family = factor(family, levels = rev(ORD_FAM)),
         Area_lab = paste0(Area, "\n(", format(tot, big.mark = ","), " reads)"))
stopifnot(!any(is.na(fam_tab$family)))
write_csv(fam_tab, file.path(DIR_TAB, "copertura_tassonomica_famiglie.csv"))
fam_tab$Area_lab <- factor(fam_tab$Area_lab, levels = rev(unique(fam_tab$Area_lab[order(match(fam_tab$Area, AREE))])))
p <- ggplot(fam_tab, aes(x = pct, y = Area_lab, fill = family)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct >= 5, sprintf("%.1f%%", pct), "")), position = position_stack(vjust = 0.5),
            size = 2.6, colour = ifelse(fam_tab$family %in% c("Canidae", "Phasianidae", "Not assigned to a family"), "black", "white")) +
  scale_fill_manual(values = COL_FAM, breaks = ORD_FAM,
                    labels = c(ORD_FAM[1:3], expression(Canidae~(italic("Vulpes vulpes"))), ORD_FAM[5:7])) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  labs(x = "Share of the reads of the run (%)", y = NULL) + tema +
  guides(fill = guide_legend(nrow = 2))
salva(p, "Figura_TaxonomicCoverage", 5.6)

# --- 12.2 output del metabarcoding
barre <- function(dati, col, livelli, titolo) {
  # position_dodge con la stessa larghezza per barre ed etichette; l'ordine dei
  # livelli (Croatia, Slovenia) mette la Slovenia sopra, la legenda resta Slovenia-Croatia
  d <- dati %>% mutate(v = factor(.data[[col]], levels = livelli)) %>% count(Area, v, .drop = FALSE) %>%
    group_by(Area) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    mutate(Area = factor(Area, levels = rev(AREE)), v = factor(v, levels = rev(livelli)))
  ggplot(d, aes(x = pct, y = v, fill = Area)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.75), hjust = -0.15, size = 2.4) +
    scale_fill_manual(values = COL_AREA, breaks = AREE, labels = lab_area_n) +
    scale_x_continuous(limits = c(0, 108), expand = c(0, 0)) +
    labs(title = titolo, x = "Samples (%)", y = NULL) + tema
}
DFo <- DF %>% mutate(classe = cut(reads_final, c(0, 10000, 50000, 100000, Inf), right = FALSE,
                                  labels = c("< 10,000", "10,000-50,000", "50,000-100,000", "> 100,000")))
p <- barre(DFo, "n_assign", 1:3, "A  Assignment categories\n    per sample") +
  barre(DFo, "n_prey_min", 1:3, "B  Minimum number of distinct\n    prey taxa per sample") +
  barre(DFo, "classe", levels(DFo$classe), "C  Prey reads per sample\n    after filtering") +
  plot_layout(guides = "collect", widths = c(1, 1, 1.35)) & theme(legend.position = "bottom")
salva(p, "Figura_OutputMetabarcoding", 6.2)

# --- 12.3 FOO e RRA a confronto
fr <- foo_rra %>% pivot_longer(c(FOO, RRA), names_to = "metric", values_to = "value") %>%
  mutate(category = factor(category, levels = rev(CATS)), Area = factor(Area, levels = rev(AREE)))
pannello <- function(m, titolo, asse) {
  ggplot(fr[fr$metric == m, ], aes(x = value, y = category, fill = Area)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = ifelse(value > 0, sprintf("%.1f", value), "0")),
              position = position_dodge(width = 0.75), hjust = -0.2, size = 2.3) +
    scale_fill_manual(values = COL_AREA, breaks = AREE, labels = lab_area_n) +
    scale_y_discrete(labels = lab_tax_expr) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(fr$value[fr$metric == m]) * 1.18)) +
    labs(title = titolo, x = asse, y = NULL) + tema
}
p <- pannello("FOO", "A  Frequency of occurrence", "Samples (%)") +
  (pannello("RRA", "B  Mean relative read abundance", "Mean RRA (%)") + theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
salva(p, "Figura_Confronto_FOO_RRA", 8.6)

# --- 12.4 composizione dei singoli campioni
comp <- DF %>% dplyr::select(Sample, Study_area, all_of(CATS))
for (k in GRP$key) comp[[k]] <- rowSums(as.matrix(comp[, GRP_CATS[[k]], drop = FALSE]))
comp$dom <- GRP$key[max.col(as.matrix(comp[, GRP$key]), ties.method = "first")]
comp$domv <- apply(as.matrix(comp[, GRP$key]), 1, max)
comp <- comp %>% mutate(o = match(dom, GRP$key)) %>% arrange(Study_area, o, desc(domv))
comp$Sample <- factor(comp$Sample, levels = comp$Sample)
comp_l <- comp %>% dplyr::select(Sample, Study_area, all_of(GRP$key)) %>%
  pivot_longer(all_of(GRP$key), names_to = "key", values_to = "rra") %>%
  mutate(key = factor(key, levels = rev(GRP$key)))
n_sa <- table(comp$Study_area)
comp_l$facet <- factor(paste0(comp_l$Study_area, "\n(n = ", n_sa[comp_l$Study_area], ")"),
                       levels = paste0(c("Slovenia", "Northern Croatia", "Southern Croatia"), "\n(n = ",
                                       n_sa[c("Slovenia", "Northern Croatia", "Southern Croatia")], ")"))
p <- ggplot(comp_l, aes(x = Sample, y = rra, fill = key, colour = key)) +
  geom_col(width = 0.86, linewidth = 0.15) +
  facet_grid(. ~ facet, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = setNames(GRP$fill, GRP$key), breaks = GRP$key,
                    labels = c(GRP$label[1:7], expression(italic("Lepus")~"sp."))) +
  scale_colour_manual(values = setNames(GRP$edge, GRP$key), guide = "none") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  labs(x = "Samples, ordered by dominant prey category and its share", y = "Relative read abundance (%)") +
  tema + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
               panel.spacing = unit(0.3, "cm")) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(colour = GRP$edge)))
salva(p, "Figura_ComposizioneCampioni", 7.4)

# --- 12.5 PCoA: campioni, profili identici e posizioni
sc <- as.data.frame(pc$points); names(sc) <- c("x", "y")
sc$Area <- DF$Area
pos_pts <- sc %>% mutate(x = round(x, 8), y = round(y, 8)) %>% count(x, y, Area, name = "n") %>%
  group_by(x, y) %>% mutate(condivisa = n_distinct(Area) > 1) %>% ungroup() %>%
  mutate(xp = x + if_else(condivisa, if_else(Area == "Slovenia", -0.018, 0.018), 0))
cat_sc <- as.data.frame(wascores(pc$points, as.matrix(DF[, CATS]))); names(cat_sc) <- c("x", "y")
cat_sc$lab <- lab_tax(rownames(cat_sc))
p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_point(data = pos_pts, aes(x = xp, y = y, colour = Area, size = n), alpha = 0.62) +
  geom_point(data = cat_sc, aes(x = x, y = y), shape = 3, size = 2, stroke = 0.6) +
  ggrepel::geom_text_repel(data = cat_sc, aes(x = x, y = y, label = lab), parse = TRUE, size = 2.6,
                           seed = SEED, min.segment.length = 0, segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf) +
  scale_colour_manual(values = COL_AREA) +
  scale_size_area(max_size = 9, breaks = c(1, 5, 10, 20), name = "Samples at the same position") +
  labs(x = sprintf("PCoA axis 1 (%.2f%% of positive eigenvalues)", K[["pcoa.ax1"]]),
       y = sprintf("PCoA axis 2 (%.2f%%)", K[["pcoa.ax2"]])) +
  tema + theme(legend.box = "vertical", legend.title = element_text(size = 7.5), panel.grid.major.y = element_line(colour = NA))
salva(p, "Grafico_PCoA_final", 12.5)

# --- 12.6 tipi di dieta
ORD_T <- c("Cervus elaphus only", "Capreolus capreolus only", "Sus scrofa only", "Caprinae only",
           "Domestic livestock only", "Mixed (more than one prey taxon)")
if (!all(DF$diet_type %in% ORD_T)) stop("Tipo di dieta non previsto nella figura: ", paste(setdiff(DF$diet_type, ORD_T), collapse = ", "))
LAB_T <- c("italic('Cervus elaphus')~plain('only')", "italic('Capreolus capreolus')~plain('only')",
           "italic('Sus scrofa')~plain('only')", "plain('Caprinae only')",
           "plain('Taxa assigned to domestic livestock only')", "plain('Mixed (more than one prey taxon)')")
td <- DF %>% mutate(diet_type = factor(diet_type, levels = ORD_T)) %>% count(Area, diet_type, .drop = FALSE) %>%
  group_by(Area) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(diet_type = factor(diet_type, levels = rev(ORD_T)), Area = factor(Area, levels = rev(AREE)))
p <- ggplot(td, aes(x = pct, y = diet_type, fill = Area)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%% (%d)", pct, n)), position = position_dodge(width = 0.75), hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = COL_AREA, breaks = AREE, labels = lab_area_n) +
  scale_y_discrete(labels = function(b) parse(text = LAB_T[match(b, ORD_T)])) +
  scale_x_continuous(limits = c(0, 52), expand = c(0, 0)) +
  labs(x = "Samples in the geographic area (%)", y = NULL) + tema
salva(p, "Figura_TipiDieta", 7)

# --- 12.7 SIMPER
sp_df <- sims %>% mutate(category = factor(category, levels = rev(category)),
                         lab = sprintf("%.1f%%   (p %s)", contribution_pct,
                                       ifelse(p < 0.0001, "< 0.0001", ifelse(p < 0.01, sprintf("= %.4f", p), sprintf("= %.3f", p)))))
p <- ggplot(sp_df, aes(x = contribution_pct, y = category)) +
  geom_col(width = 0.62, fill = "#6b6b6b") +
  geom_text(aes(label = lab), hjust = -0.08, size = 2.4) +
  scale_y_discrete(labels = lab_tax_expr) +
  scale_x_continuous(limits = c(0, max(sp_df$contribution_pct) * 1.45), expand = c(0, 0)) +
  labs(x = "Contribution to the average Slovenia-Croatia dissimilarity (%)", y = NULL) + tema
salva(p, "Grafico_SIMPER_final", 7.2)

# --- 12.8 modello sulla dieta mista (curve del modello riportato, quota media)
g_hhh <- expand_grid(NDVI_sd = seq(min(DF$NDVI_sd), max(DF$NDVI_sd), length.out = 200), Croatia = 0:1) %>%
  mutate(NDVI_sd_z = (NDVI_sd - mean(DF$NDVI_sd)) / sd(DF$NDVI_sd), ELEV_mean_z = 0,
         Area = factor(if_else(Croatia == 1, "Croatia", "Slovenia"), levels = AREE))
pr <- predict(m_hhh0, newdata = g_hhh, type = "link", se.fit = TRUE)
g_hhh <- g_hhh %>% mutate(p = plogis(pr$fit), lo = plogis(pr$fit - Z975 * pr$se.fit), hi = plogis(pr$fit + Z975 * pr$se.fit))
set.seed(SEED)
p <- ggplot() +
  geom_ribbon(data = g_hhh, aes(x = NDVI_sd, ymin = lo, ymax = hi, fill = Area), alpha = 0.15) +
  geom_line(data = g_hhh, aes(x = NDVI_sd, y = p, colour = Area), linewidth = 0.8) +
  geom_point(data = DF, aes(x = NDVI_sd, y = mixed_prey, colour = Area), size = 1.3, alpha = 0.55,
             position = position_jitter(width = 0, height = 0.025, seed = SEED)) +
  scale_colour_manual(values = COL_AREA, labels = lab_area_n) + scale_fill_manual(values = COL_AREA, labels = lab_area_n) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0\n(single prey)", "0.25", "0.50", "0.75", "1\n(mixed)")) +
  labs(x = "NDVI standard deviation within 2 km of the sampling location", y = "Probability of a mixed-prey sample") +
  tema + theme(panel.grid.major.y = element_line(colour = "grey92"))
salva(p, "HHH_Logistic_Plot_final", 8.2)

# --- 12.9 modelli di dominanza (costanti di standardizzazione dei 103 campioni)
ELEV_MU <- mean(DF$ELEV_mean); ELEV_SD <- sd(DF$ELEV_mean)
curva <- function(m, dati) {
  g <- tibble(ELEV_mean = seq(min(dati$ELEV_mean), max(dati$ELEV_mean), length.out = 200)) %>%
    mutate(ELEV_mean_z = (ELEV_mean - ELEV_MU) / ELEV_SD, NDVI_mean_z = 0)
  pr <- predict(m, newdata = g, type = "link", se.fit = TRUE)
  g %>% mutate(p = plogis(pr$fit), lo = plogis(pr$fit - Z975 * pr$se.fit), hi = plogis(pr$fit + Z975 * pr$se.fit))
}
m_l1f <- glm(boar ~ ELEV_mean_z + NDVI_mean_z, data = DOM, family = binomial)
m_l2f <- glm(red ~ ELEV_mean_z + NDVI_mean_z, data = CERV, family = binomial)
fmtp <- function(p) ifelse(p < 0.0001, "< 0.0001", ifelse(p < 0.01, sprintf("= %.4f", p), sprintf("= %.3f", p)))
pan_dom <- function(m, dati, risp, titolo, etich, key) {
  set.seed(SEED)
  ggplot() +
    geom_ribbon(data = curva(m, dati), aes(x = ELEV_mean, ymin = lo, ymax = hi), fill = "grey70", alpha = 0.4) +
    geom_line(data = curva(m, dati), aes(x = ELEV_mean, y = p), linewidth = 0.7) +
    geom_point(data = dati, aes(x = ELEV_mean, y = .data[[risp]], colour = Area), size = 1.2, alpha = 0.6,
               position = position_jitter(width = 0, height = 0.03, seed = SEED)) +
    annotate("text", x = max(dati$ELEV_mean), y = 0.5, hjust = 1, size = 2.4,
             label = sprintf("n = %d; elevation OR %.3f\np %s", nrow(dati), K[[paste0(key, ".OR")]], fmtp(K[[paste0(key, ".p")]]))) +
    scale_colour_manual(values = COL_AREA) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(etich[1], "0.5", etich[2])) +
    labs(title = titolo, x = "Mean elevation within 2 km (m)", y = "Predicted probability") +
    tema + theme(panel.grid.major.y = element_line(colour = "grey92"))
}
p <- pan_dom(m_l1f, DOM, "boar", "A  Wild boar vs cervids", c("Cervid-\ndominated", "Wild boar-\ndominated"), "l1") +
  pan_dom(m_l2f, CERV, "red", "B  Red deer vs roe deer", c("Roe deer-\ndominated", "Red deer-\ndominated"), "l2") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
salva(p, "Grafico_Gerarchico_Cinghiale_Cervidi_final", 8.4)

# --- 12.10 RDA (scaling 2; frecce riscalate per la visualizzazione)
s_site <- as.data.frame(scores(m_rda, display = "sites", scaling = 2, choices = 1:2))
s_spec <- as.data.frame(scores(m_rda, display = "species", scaling = 2, choices = 1:2))
s_bp   <- as.data.frame(scores(m_rda, display = "bp", scaling = 2, choices = 1:2))
names(s_site) <- names(s_spec) <- names(s_bp) <- c("x", "y")
s_site$Area <- DF$Area
site_pts <- s_site %>% mutate(x = round(x, 8), y = round(y, 8)) %>% count(x, y, Area, name = "n") %>%
  group_by(x, y) %>% mutate(condivisa = n_distinct(Area) > 1) %>% ungroup() %>%
  mutate(xp = x + if_else(condivisa, if_else(Area == "Slovenia", -0.03, 0.03), 0))
r_site <- max(abs(c(s_site$x, s_site$y)))
mult_sp <- 0.75 * r_site / max(sqrt(s_spec$x^2 + s_spec$y^2))
mult_bp <- 0.75 * r_site / max(sqrt(s_bp$x^2 + s_bp$y^2))
put("rda.fig.mult.species", mult_sp); put("rda.fig.mult.env", mult_bp)
s_spec <- s_spec %>% mutate(x = x * mult_sp, y = y * mult_sp, lab = lab_tax(rownames(s_spec)))
s_bp <- s_bp %>% mutate(x = x * mult_bp, y = y * mult_bp,
                        lab = c(ELEV_mean = "Elevation", NDVI_mean = "NDVI mean", NDVI_sd = "NDVI SD")[rownames(s_bp)])
p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_point(data = site_pts, aes(x = xp, y = y, colour = Area, size = n), alpha = 0.5) +
  geom_segment(data = s_spec, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.35) +
  ggrepel::geom_text_repel(data = s_spec, aes(x = x, y = y, label = lab), parse = TRUE, size = 2.6, seed = SEED,
                           min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf) +
  geom_segment(data = s_bp, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 0.55, colour = "#8b1a1a") +
  ggrepel::geom_text_repel(data = s_bp, aes(x = x, y = y, label = lab), size = 2.6, fontface = "bold", colour = "#8b1a1a",
                           seed = SEED, min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf) +
  scale_colour_manual(values = COL_AREA) + scale_size_area(max_size = 7, guide = "none") +
  labs(x = sprintf("RDA1 (%.2f%% of total variance)", K[["rda.ax1"]]), y = sprintf("RDA2 (%.2f%%)", K[["rda.ax2"]])) +
  tema + theme(panel.grid.major.y = element_line(colour = NA))
salva(p, "RDA_triplot_final", 11.5)


# ==============================================================================
# 12b. FILE CON I NOMI DELLA VERSIONE PRECEDENTE (nella cartella di lavoro)
#      Stessi nomi e stesso formato dello script precedente, cosi' le mappe di
#      QGIS e gli altri file che li usano restano validi. Novita': la categoria
#      Cervinae (solo MSV136) e, nei file delle torte, le colonne CERVINAE e
#      NPREY_MIN (numero minimo di prede distinte).
# ==============================================================================
POP_IT <- c(Slovenia = "Slovenia", Croatia = "Croazia")
ricodifica <- function(x, mappa) { y <- x; k <- x %in% names(mappa); y[k] <- unname(mappa[x[k]]); y }
DFo <- DF %>% mutate(Area = as.character(Area)) %>% arrange(match(Area, AREE), Sample)

# --- FOO e RRA per area (FOO_RRA_Slovenia_and_Croatia.csv) ----------------------
foo_old <- bind_rows(lapply(AREE, function(a) {
  s_a <- DFo[DFo$Area == a, ]
  pres <- vapply(CATS, function(cc) sum(s_a[[cc]] > 0), numeric(1))
  tibble(Popolazione = POP_IT[[a]], Preda = CATS, Presenza_Assoluta = as.integer(pres),
         FOO_percentuale = round(100 * pres / nrow(s_a), 2),
         RRA_percentuale = round(vapply(CATS, function(cc) mean(s_a[[cc]]), numeric(1)), 2))
})) %>% filter(Presenza_Assoluta > 0) %>% arrange(Popolazione, desc(RRA_percentuale))
write.csv2(foo_old, file.path(PROJECT_DIR, OUT_FOO_RRA), row.names = FALSE)

# --- matrice RRA per campione (Community_Matrix_SloCro.csv) ---------------------
comm_old <- DFo %>% transmute(Sample_ID = Sample, Popolazione = unname(POP_IT[Area]), across(all_of(CATS)))
write.csv2(comm_old, file.path(PROJECT_DIR, OUT_COMMUNITY), row.names = FALSE)

# --- torte per QGIS: campioni e aree di studio -----------------------------------
torte_camp <- DFo %>% transmute(
  Sample, Area, Study = Study_area,
  Region = ricodifica(Region_cluster_comparison, c(Others = "Other packs", Banovina = "HRV027")),
  Latitude, Longitude,
  ROE = round(`Capreolus capreolus`, 3), RED = round(`Cervus elaphus`, 3), BOAR = round(`Sus scrofa`, 3),
  CHAMOIS = round(`Rupicapra rupicapra`, 3), CAPRINAE = round(Caprinae, 3),
  DOMESTIC = round(Bos + Capra + Ovis + `Ovis aries`, 3), LEPUS = round(Lepus, 3),
  CERVINAE = round(Cervinae, 3), NTAXA = n_assign, NPREY_MIN = n_prey_min)
write_csv(torte_camp, file.path(PROJECT_DIR, "Torte_campioni_103.csv"))
torte_aree <- DFo %>% group_by(Study = Study_area) %>%
  summarise(N = n(), N_IND = n_distinct(Individual[!is.na(Individual)]),
            ROE = round(mean(`Capreolus capreolus`), 2), RED = round(mean(`Cervus elaphus`), 2),
            BOAR = round(mean(`Sus scrofa`), 2), CHAMOIS = round(mean(`Rupicapra rupicapra`), 2),
            CAPRINAE = round(mean(Caprinae), 2), DOMESTIC = round(mean(Bos + Capra + Ovis + `Ovis aries`), 2),
            LEPUS = round(mean(Lepus), 2), CERVINAE = round(mean(Cervinae), 2),
            Lat_media = round(mean(Latitude), 5), Lon_media = round(mean(Longitude), 5), .groups = "drop") %>%
  transmute(Study, N, N_IND, Latitude = Lat_media, Longitude = Lon_media,
            ROE, RED, BOAR, CHAMOIS, CAPRINAE, DOMESTIC, LEPUS, CERVINAE) %>%
  arrange(match(Study, c("Slovenia", "Northern Croatia", "Southern Croatia")))
write_csv(torte_aree, file.path(PROJECT_DIR, "Torte_aree_3.csv"))

# --- FOO e RRA per area di studio (Tabella_FOO_RRA_3aree.csv) --------------------
AREE3 <- c("Slovenia", "Northern Croatia", "Southern Croatia")
tab3 <- tibble(Taxon = CATS)
for (a in AREE3) tab3[[paste0("FOO_", a)]] <-
  round(vapply(CATS, function(cc) 100 * mean(DFo[[cc]][DFo$Study_area == a] > 0), numeric(1)), 2)
for (a in AREE3) tab3[[paste0("RRA_", a)]] <-
  round(vapply(CATS, function(cc) mean(DFo[[cc]][DFo$Study_area == a]), numeric(1)), 2)
write.csv2(tab3, file.path(PROJECT_DIR, "Tabella_FOO_RRA_3aree.csv"), row.names = FALSE)

# --- copertura per famiglia (Tabella_coverage_famiglie.csv) ----------------------
ORD_FAM_OLD <- c("Cervidae", "Bovidae", "Suidae", "Canidae", "Leporidae", "Phasianidae", "Not assigned to family")
cov_old <- LUNGO %>% filter(role != "predator") %>%
  mutate(Famiglia = ricodifica(unname(FAM[final_name]), c("Not assigned to a family" = "Not assigned to family"))) %>%
  group_by(Area, Famiglia) %>%
  summarise(MOTU = n_distinct(motu_id), reads = sum(reads), .groups = "drop") %>%
  filter(reads > 0) %>%
  group_by(Area) %>% mutate(pct = round(100 * reads / sum(reads), 2)) %>% ungroup() %>%
  arrange(match(Area, AREE), match(Famiglia, ORD_FAM_OLD)) %>%
  dplyr::select(Area, Famiglia, reads, MOTU, pct)
write.csv2(cov_old, file.path(PROJECT_DIR, "Tabella_coverage_famiglie.csv"), row.names = FALSE)

# --- covariate ambientali per area (Tabella_covariate_ambientali.csv) ------------
cov_amb <- DFo %>% group_by(Area) %>%
  summarise(n = n(),
            ELEV_mean_min = min(ELEV_mean), ELEV_mean_mediana = median(ELEV_mean), ELEV_mean_max = max(ELEV_mean),
            NDVI_mean_min = min(NDVI_mean), NDVI_mean_mediana = median(NDVI_mean), NDVI_mean_max = max(NDVI_mean),
            NDVI_sd_min = min(NDVI_sd), NDVI_sd_mediana = median(NDVI_sd), NDVI_sd_max = max(NDVI_sd),
            .groups = "drop") %>%
  arrange(Area)
write.csv2(cov_amb, file.path(PROJECT_DIR, "Tabella_covariate_ambientali.csv"), row.names = FALSE)
LOG("CSV con i nomi della versione precedente riscritti nella cartella di lavoro")

# --- grafici a due pannelli per area (Grafico_2Pannelli_*_final.png) -------------
COL_CATEGORIA <- c("wild ungulates" = "#004C6D", "domestic livestock" = "#D95F59",
                   "unresolved Caprinae" = "#B84D9B", "unresolved Cervinae" = "#8C8C8C",
                   "other taxa" = "#E69F00")
categoria_trofica <- function(p) dplyr::case_when(
  p %in% c("Capreolus capreolus", "Cervus elaphus", "Sus scrofa", "Rupicapra rupicapra") ~ "wild ungulates",
  p %in% DOMESTIC ~ "domestic livestock",
  p == "Caprinae" ~ "unresolved Caprinae",
  p == "Cervinae" ~ "unresolved Cervinae",
  TRUE ~ "other taxa")
grafico_area <- function(a) {
  d <- foo_old %>% filter(Popolazione == POP_IT[[a]]) %>%
    mutate(Categoria = factor(categoria_trofica(Preda), levels = names(COL_CATEGORIA)))
  d$Preda <- factor(d$Preda, levels = d$Preda[order(d$RRA_percentuale)])
  dl <- d %>% pivot_longer(c(FOO_percentuale, RRA_percentuale), names_to = "Metrica", values_to = "Valore") %>%
    mutate(Metrica = factor(ifelse(Metrica == "FOO_percentuale", "A  Frequency of occurrence",
                                   "B  Relative read abundance"),
                            levels = c("A  Frequency of occurrence", "B  Relative read abundance")))
  sottotitolo <- if ("Cervinae" %in% d$Preda) {
    "; Caprinae and Cervinae are unresolved subfamilies, not attributed to wild or domestic taxa"
  } else "; Caprinae is not attributed to wild or domestic taxa"
  ggplot(dl, aes(x = Valore, y = Preda, fill = Categoria)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = sprintf("%.1f", Valore)), hjust = -0.2, size = 3.5) +
    facet_wrap(~ Metrica, scales = "free_x") +
    scale_fill_manual(values = COL_CATEGORIA, drop = TRUE) +
    scale_y_discrete(labels = function(x) parse(text = lab_tax(x))) +
    scale_x_continuous(limits = c(0, max(dl$Valore) * 1.15), expand = c(0, 0)) +
    labs(title = paste("Diet composition:", a),
         subtitle = paste0("n = ", sum(DFo$Area == a), " samples", sottotitolo),
         x = "Percentage (%)", y = NULL, fill = NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 11, colour = "black"),
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 4)),
          plot.subtitle = element_text(colour = "grey35", hjust = 0.5, margin = margin(b = 12)),
          strip.text = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")
}
for (a in AREE) {
  p2 <- grafico_area(a)
  nome2 <- c(Slovenia = "Grafico_2Pannelli_Slovenia_final.png", Croatia = "Grafico_2Pannelli_Croazia_final.png")[[a]]
  ggsave(file.path(PROJECT_DIR, nome2), p2, width = 12, height = 6, dpi = 300)
  if (interactive()) print(p2)
  LOG("Figura salvata: ", nome2)
}

# --- dispersione multivariata (Betadisper_Plot_final.png) ------------------------
# Figura descrittiva, non usata nella tesi: le ellissi non si usano per conclusioni.
COL_AREA_V <- unname(COL_AREA[levels(DF$Area)])
disegna_betadisper <- function() {
  plot(bd, hull = FALSE, ellipse = TRUE, main = "Multivariate dispersion of diet composition",
       sub = "Bray-Curtis dissimilarity on Hellinger-transformed RRA proportions",
       col = COL_AREA_V, lwd = 2, seg.col = "grey80", seg.lwd = 0.5)
  legend("topleft", legend = levels(DF$Area), col = COL_AREA_V, pch = 16, bty = "n", cex = 1.1)
}
png(file.path(PROJECT_DIR, "Betadisper_Plot_final.png"), width = 2000, height = 1600, res = 300)
disegna_betadisper()
invisible(dev.off())
if (interactive()) disegna_betadisper()
LOG("Figura salvata: Betadisper_Plot_final.png")


# ==============================================================================
# 13. RISULTATI CHIAVE, MACRO LATEX E CONFRONTO CON I VALORI PROVVISORI
# ==============================================================================
val_chr <- function(v) {
  if (length(v) == 0 || is.na(v)) return("NA")
  if (is.logical(v)) return(if (v) "True" else "False")
  if (is.numeric(v)) return(format(v, digits = 15))
  as.character(v)
}
chiavi <- tibble(key = names(K), value = vapply(K, function(v) val_chr(v[1]), character(1)))
write_csv(chiavi, file.path(DIR_OUT, "risultati_chiave_R.csv"))

spec <- read_csv(FILE_SPEC, show_col_types = FALSE, col_types = cols(.default = col_character()))
migliaia <- function(x) formatC(round(x), format = "d", big.mark = ",")
fmt_p <- function(p) ifelse(p < 0.0001, "<0.0001", ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.3f", p)))
fmt_pe <- function(p) ifelse(p < 0.0001, "< 0.0001", paste("=", fmt_p(p)))
evidenza <- function(p) ifelse(p < 0.001, "very strong", ifelse(p < 0.01, "strong", ifelse(p < 0.05, "moderate",
                                ifelse(p < 0.1, "weak", "little or no"))))
# numeri da zero a nove in lettere nel testo (wd; Wd a inizio frase), cifre da 10 in su
PAROLE <- c("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine")
in_lettere <- function(x, maiuscola = FALSE) {
  n <- round(x)
  s <- if (n >= 0 && n <= 9) PAROLE[n + 1] else migliaia(n)
  if (maiuscola) paste0(toupper(substr(s, 1, 1)), substring(s, 2)) else s
}
formatta <- function(v, f) {
  if (f == "str") return(as.character(v))
  x <- as.numeric(v)
  switch(f,
         int = migliaia(x), wd = in_lettere(x), Wd = in_lettere(x, TRUE),
         f0 = sprintf("%.0f", x), f1 = sprintf("%.1f", x), f2 = sprintf("%.2f", x),
         f3 = sprintf("%.3f", x), f4 = sprintf("%.4f", x), f5 = sprintf("%.5f", x),
         x100f1 = sprintf("%.1f", 100 * x), x100f2 = sprintf("%.2f", 100 * x),
         p = fmt_p(x), pe = fmt_pe(x), ev = evidenza(x),
         stop("Formato sconosciuto: ", f))
}
righe_tex <- c("% Numeri della tesi generati da Script_Tesi_Lupo_COMPLETO.R",
               paste0("% ", format(Sys.time(), "%Y-%m-%d %H:%M"), "; R ", as.character(getRversion()),
                      "; vegan ", versione("vegan"), "; permutazioni ", N_PERM, "; seme ", SEED))
mancano <- character(0)
for (i in seq_len(nrow(spec))) {
  rk <- spec$raw_key[i]
  if (!rk %in% names(K)) { mancano <- c(mancano, rk); next }
  righe_tex <- c(righe_tex, sprintf("\\defres{%s}{%s}", spec$macro[i], formatta(K[[rk]], spec$format[i])))
}
# versioni dei pacchetti per la sezione Software (con uno spazio iniziale, perche'
# nel testo la macro segue direttamente il nome del pacchetto)
for (pk in c("vegan", "permute", "MASS", "sandwich", "lme4", "ggplot2", "readr", "dplyr", "tidyr", "patchwork", "ggrepel")) {
  righe_tex <- c(righe_tex, sprintf("\\defres{ver.%s:str}{ %s}", pk, versione(pk)))
}
righe_tex <- c(righe_tex, sprintf("\\defres{ver.R:str}{%s}", as.character(getRversion())))
writeLines(righe_tex, file.path(PROJECT_DIR, OUT_TEX), useBytes = TRUE)
LOG("Numeri della tesi scritti in ", OUT_TEX, " (da caricare su Overleaf)")
if (length(mancano) > 0) {
  segnala("Chiavi della specifica non calcolate (nella tesi comparirebbero come [?chiave?]): ", paste(mancano, collapse = ", "))
  writeLines(mancano, file.path(DIR_OUT, "chiavi_mancanti.txt"))
}
LOG("Macro LaTeX scritte: ", sum(grepl("^\\\\defres", righe_tex)))

# confronto con i valori provvisori (Python) usati nella bozza
if (file.exists(FILE_PROV)) {
  prov <- read_csv(FILE_PROV, show_col_types = FALSE, col_types = cols(.default = col_character()))
  conf <- inner_join(prov, chiavi, by = "key", suffix = c("_provvisorio", "_R"))
  norm_na <- function(x) ifelse(tolower(x) %in% c("nan", "na", "none"), "NA", x)
  conf$value_provvisorio <- norm_na(conf$value_provvisorio); conf$value_R <- norm_na(conf$value_R)
  num <- suppressWarnings(cbind(as.numeric(conf$value_provvisorio), as.numeric(conf$value_R)))
  conf$differenza <- num[, 2] - num[, 1]
  conf$rel <- abs(conf$differenza) / pmax(abs(num[, 1]), 1e-12)
  # p-value per permutazione (errore Monte Carlo); tutto il resto e' deterministico
  mc <- "(^|\\.)perm\\.p$|^(bd|ind\\.bd|rda|ind\\.rda)\\.p$|^(ind\\.)?rda\\.m\\.[A-Za-z_]+\\.p$|^simper\\.[a-z]+\\.p$|^simper\\.rest\\.p(min|max)$|^sens\\.[a-z0-9]+\\.perm_p$"
  conf$tipo <- ifelse(grepl(mc, conf$key), "p per permutazione", "deterministico")
  conf$da_controllare <- (is.na(num[, 1]) & conf$value_provvisorio != conf$value_R) |
    (!is.na(conf$rel) & conf$rel > 1e-6 & conf$tipo == "deterministico")
  # cambi di categoria dell'evidenza (scala di Muff et al. 2022) per i p-value
  pk <- grepl("(^|\\.)p$|_p$", conf$key) & !is.na(num[, 1]) & !is.na(num[, 2])
  conf$cambio_evidenza <- FALSE
  conf$cambio_evidenza[pk] <- evidenza(num[pk, 1]) != evidenza(num[pk, 2])
  write_csv(conf, file.path(DIR_OUT, "confronto_R_vs_provvisori.csv"))
  LOG("Confronto con i valori provvisori: ", sum(conf$da_controllare, na.rm = TRUE), " valori deterministici diversi; ",
      sum(conf$cambio_evidenza, na.rm = TRUE), " p-value con categoria di evidenza diversa")
  if (any(conf$cambio_evidenza, na.rm = TRUE)) {
    LOG("ATTENZIONE: rileggere nel testo le frasi con questi p-value: ",
        paste(conf$key[conf$cambio_evidenza %in% TRUE], collapse = ", "))
  }
  solo_prov <- setdiff(prov$key, chiavi$key)
  if (length(solo_prov) > 0) LOG("Chiavi presenti solo nei valori provvisori: ", paste(solo_prov, collapse = ", "))
}
# confronto del testo stampato macro per macro con il file provvisorio della bozza:
# sono le differenze che cambiano la tesi
if (file.exists(FILE_TEX_PROV)) {
  leggi_macro <- function(righe) {
    righe <- righe[grepl("^\\\\defres\\{", righe)]
    tibble(macro = sub("^\\\\defres\\{([^}]*)\\}.*$", "\\1", righe),
           testo = sub("^\\\\defres\\{[^}]*\\}\\{(.*)\\}$", "\\1", righe))
  }
  m_prov <- leggi_macro(readLines(FILE_TEX_PROV, encoding = "UTF-8"))
  m_R <- leggi_macro(righe_tex)
  m_conf <- full_join(m_prov, m_R, by = "macro", suffix = c("_provvisorio", "_R")) %>%
    filter(is.na(testo_provvisorio) | is.na(testo_R) | testo_provvisorio != testo_R)
  write_csv(m_conf, file.path(DIR_OUT, "macro_cambiate_rispetto_ai_provvisori.csv"))
  LOG("Macro il cui testo stampato cambia rispetto alla bozza: ", nrow(m_conf),
      " (elenco in macro_cambiate_rispetto_ai_provvisori.csv; le versioni dei pacchetti cambiano sempre)")
}


# ==============================================================================
# 14. AMBIENTE DI ESECUZIONE
# ==============================================================================
si <- capture.output(sessionInfo())
writeLines(si, file.path(DIR_OUT, "sessionInfo.txt"))
cat(si, sep = "\n")
if (length(AVVISI) > 0) writeLines(AVVISI, file.path(DIR_OUT, "avvisi.txt"))
LOG("Fine. Output in ", DIR_OUT)
