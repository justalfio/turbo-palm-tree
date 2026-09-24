# ==============================================================================
# DIET COMPOSITION AND ENVIRONMENTAL CORRELATES OF THE GREY WOLF
# IN SLOVENIA AND CROATIA - SCRIPT UNICO DI ANALISI
#
# Alfio Tomarchio - Universita' di Bologna
# Ultima revisione: 23 settembre 2026
#
# COSA E' CAMBIATO RISPETTO ALLA VERSIONE PRECEDENTE (elenco completo)
#
#  1. library(MASS) caricato PRIMA di tidyverse. Serve per gli intervalli di
#     confidenza da profilo di verosimiglianza di confint() sui GLM; se MASS
#     non e' caricato, confint() restituisce intervalli di Wald, che sono
#     diversi da quelli riportati in tesi. Caricato prima perche' MASS::select
#     mascherebbe dplyr::select.
#  2. PCoA ricalcolata sui 10 TAXA, non sulle 6 categorie aggregate. La
#     versione a 6 categorie da' 52.70% + 27.06%; la figura della tesi riporta
#     52.70% + 26.72%, che e' la versione a 10 taxa. Ora lo script riproduce
#     esattamente la figura della tesi.
#  3. Grafico PCoA ridisegnato: i campioni con lo stesso profilo di preda
#     cadono nello stesso punto, quindi il jitter casuale nascondeva la
#     struttura del dato. Ora la dimensione del punto e' il numero di campioni
#     che occupano quella posizione (fino a 35 su un solo punto). ATTENZIONE
#     alla distinzione: i 100 campioni danno 25 PROFILI dietetici distinti ma
#     solo 24 POSIZIONI distinte nel piano, perche' il campione con solo Ovis e
#     quello con solo Capra hanno coordinate identiche sui primi cinque assi e
#     si separano solo sul sesto (1.51% degli autovalori).
#  4. Nuova Figura 3.1: output del metabarcoding (taxa per campione e reads
#     per campione). Risponde al commento "manca tutto il capitolo di
#     descrizione dei dati del metabarcoding".
#  5. Nuova Figura 3.3: FOO e RRA delle due aree sugli stessi assi. Risponde al
#     commento "utilizza qualche barplot per una migliore comprensione".
#  6. Nuova Figura 3.6: distribuzione dei campioni per tipo di dieta. E' la
#     figura che rende leggibile cio' che nella PCoA resta nascosto.
#  7. Grafico HHH corretto: la curva disegnata da geom_smooth() proveniva da un
#     modello univariato su NDVI_sd, NON dal modello riportato in tesi. Ora le
#     curve sono predette dal modello multivariato, separate per area e a
#     quota media, esattamente come dice la didascalia.
#  8. Aggiunto il test del rapporto di verosimiglianza sull'interazione
#     NDVI_sd x area (chi2 = 0.184, p = 0.668), riportato in tesi ma assente
#     dallo script.
#  9. Betadisper: titolo ed etichette in inglese.
# 10. Palette unica per tutte le figure. Il blu #2E86AB non supera il controllo
#     di croma (0.099: a stampa e in scala di grigi vira al grigio); sostituito
#     con #1F7FB5, visivamente quasi identico e con separazione adeguata anche
#     per daltonismo (Delta E 22.0 protanopia).
# 11. Percorsi dei file raccolti in un unico blocco in testa.
# 12. Controlli espliciti (stopifnot) sui numeri attesi, cosi' se un file viene
#     sostituito lo script si ferma invece di produrre numeri diversi in
#     silenzio.
# 13. sessionInfo() in coda, per la sezione Software and Reproducibility.
# 14. (22/09/2026) Nuova sezione 11: analisi esplorativa per regione dentro la
#     Croazia, per la Discussion (Sezione 4.2). Su indicazione della relatrice
#     l'attribuzione ai cluster genetici e il confronto fra regioni stanno nella
#     Discussion; lo script li rende riproducibili.
# 15. (23/09/2026) Recuperati 3 campioni sloveni (MSV134, MSV15H, MSV15M): nel
#     foglio "Merged scientific_name" del composition report le loro reads di
#     capriolo erano state azzerate, mentre nella tabella dei MOTU (READ COUNT)
#     e nella matrice LECA 2022 superano la soglia di 2.000 reads. Rispettano
#     tutte le regole di inclusione della tesi, quindi rientrano: Slovenia 52,
#     totale 103. Corretta anche la cella MSV12L (3850 -> 3580, trasposizione).
#     Le numerosita' attese sono ora in N_SLO / N_CRO / N_TOT; i valori "atteso"
#     che dipendono dai modelli ambientali vanno riletti dopo il rilancio.
# 16. (23/09/2026) Nuovo file ambientale Wolf_NDVI_Dynamic_Dinaric.csv con 103
#     righe, tutte ricalcolate con GEE_NDVI_quota_103_campioni.js (stesso
#     metodo per tutti i campioni). ELEV_mean = media nel buffer di 2 km (i
#     vecchi valori coincidono entro 1 m). I sottotitoli della figura
#     gerarchica leggono ora il p dal modello invece di un numero scritto a mano.
# 17. (24/09/2026) Nuove sezioni 13-17: tre aree di studio per le mappe e CSV
#     per QGIS; tabella FOO/RRA per area di studio; Figura G (composizione dei
#     singoli campioni); torte del taxonomic coverage per famiglia (dati grezzi);
#     tabella delle covariate per area; figura NDVI dei tre buffer di esempio
#     (facoltativa: serve terra + sf). La sessionInfo diventa la sezione 18.
#     In RStudio i grafici compaiono nel pannello Plots e le tabelle principali
#     nel visualizzatore dati (funzioni mostra() e mostra_tabella()).
#
# NOTA SULLA SOGLIA DELL'1%: i due file FILE DEFINITIVO hanno gia' la soglia
# applicata a monte. Lo script lo verifica invece di darlo per scontato.
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

# MASS PRIMA di tidyverse: confint() sui GLM usa il metodo di MASS (profilo di
# verosimiglianza). Invertendo l'ordine, MASS::select mascherebbe dplyr::select.
library(MASS)
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)

set.seed(123)

# --- percorsi dei file (unico punto da modificare se rinomini qualcosa) -------
FILE_SLO       <- "Slovenia wolves FILE DEFINITIVO.csv"
FILE_CRO       <- "Croatian wolves FILE DEFINITIVO.csv"
FILE_AMBIENTE  <- "Wolf_NDVI_Dynamic_Dinaric.csv"
OUT_FOO_RRA    <- "FOO_RRA_Slovenia_and_Croatia.csv"
OUT_COMMUNITY  <- "Community_Matrix_SloCro.csv"
FILE_INDIVIDUI <- "individui_e_branchi_103_campioni.csv"  # sezioni 11 e 13
FILE_PUNTI     <- "Punti_103_per_QGIS.csv"                # sezioni 13 e 17
FILE_COV_SLO   <- "Taxonomic_coverage_Slovenia.csv"       # sezione 15 (dati grezzi)
FILE_COV_CRO   <- "Taxonomic_coverage_Croatia.csv"        # sezione 15 (dati grezzi)

# --- controllo preliminare: tutti i file di input devono essere nella cartella -
# Se Windows nasconde le estensioni, un file che in Esplora risorse si vede come
# "Taxonomic_coverage_Slovenia.csv" in realta' si chiama "...csv.csv". Il
# controllo lo riconosce e dice quale file rinominare, prima di iniziare.
file_necessari <- c(FILE_SLO, FILE_CRO, FILE_AMBIENTE, FILE_INDIVIDUI,
                    FILE_PUNTI, FILE_COV_SLO, FILE_COV_CRO)
mancanti <- file_necessari[!file.exists(file_necessari)]
for (f in mancanti) {
  doppio <- paste0(f, ".csv")
  if (file.exists(doppio)) {
    message("Il file '", doppio, "' ha l'estensione doppia. Rinominalo con:\n",
            "  file.rename(\"", doppio, "\", \"", f, "\")")
  } else {
    message("Manca il file '", f, "' nella cartella ", getwd())
  }
}
if (length(mancanti) > 0) stop("File di input mancanti: lo script si ferma prima di iniziare.")

# --- numerosita' attese (23/09/2026: +3 campioni sloveni recuperati) ----------
N_SLO <- 52
N_CRO <- 51
N_TOT <- N_SLO + N_CRO     # 103

# --- palette unica per tutte le figure della tesi -----------------------------
# Slovenia #1F7FB5, Croazia #C73E1D. Controllate: banda di luminosita',
# soglia di croma, separazione per protanopia/deuteranopia/tritanopia,
# separazione a visione normale, contrasto sul fondo chiaro.
COL_AREA <- c("Slovenia" = "#1F7FB5", "Croazia" = "#C73E1D")
LAB_AREA <- c("Slovenia" = "Slovenia", "Croazia" = "Croatia")

# Palette per categoria trofica (grafici a due pannelli per area)
COL_CATEGORIA <- c(
  "wild ungulates"      = "#004C6D",
  "domestic livestock"  = "#D95F59",
  "unresolved Caprinae" = "#B84D9B",
  "other taxa"          = "#E69F00"
)

# Generi e specie in corsivo; Caprinae (sottofamiglia) in tondo
etichette_taxa_plotmath <- function(x) {
  parse(text = ifelse(x == "Caprinae",
                      "plain('Caprinae')",
                      sprintf("italic('%s')", x)))
}

# --- visualizzazione durante l'esecuzione -------------------------------------
# In RStudio ogni grafico compare nel pannello Plots (freccette per sfogliarli)
# e le tabelle principali si aprono in una scheda del visualizzatore dati.
# Lanciato con Rscript (non interattivo) lo script salva soltanto i file.
mostra <- function(grafico) {
  if (interactive()) print(grafico)
  invisible(grafico)
}
mostra_tabella <- function(tabella, titolo) {
  if (interactive()) utils::View(tabella, title = titolo)
  invisible(tabella)
}

tema_tesi <- theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey35"),
    panel.grid.minor = element_blank(),
    legend.title  = element_blank()
  )


# ==============================================================================
# 1. CARICAMENTO, FORMATO LUNGO E VERIFICA DELLA SOGLIA DELL'1%
# ==============================================================================

slo <- read_delim(FILE_SLO, delim = ";", locale = locale(encoding = "UTF-8"),
                  show_col_types = FALSE)
cro <- read_delim(FILE_CRO, delim = ";", locale = locale(encoding = "UTF-8"),
                  show_col_types = FALSE)

slo_long <- slo %>%
  dplyr::select(-rank) %>%
  pivot_longer(-scientific_name, names_to = "Sample_ID", values_to = "reads") %>%
  mutate(Popolazione = "Slovenia")

cro_long <- cro %>%
  dplyr::select(-rank) %>%
  pivot_longer(-scientific_name, names_to = "Sample_ID", values_to = "reads") %>%
  mutate(Popolazione = "Croazia")

grezzi <- bind_rows(slo_long, cro_long)

# Controlli: 103 campioni, 52 + 51
stopifnot(n_distinct(grezzi$Sample_ID) == N_TOT)
stopifnot(n_distinct(grezzi$Sample_ID[grezzi$Popolazione == "Slovenia"]) == N_SLO)
stopifnot(n_distinct(grezzi$Sample_ID[grezzi$Popolazione == "Croazia"])  == N_CRO)

grezzi <- grezzi %>%
  group_by(Sample_ID) %>%
  mutate(totale_grezzo = sum(reads)) %>%
  ungroup() %>%
  mutate(pct_grezza = if_else(totale_grezzo > 0, reads / totale_grezzo * 100, 0))

# Verifica della soglia dell'1%: nei FILE DEFINITIVO la soglia e' gia' stata
# applicata a monte, quindi qui non deve restare nessuna detection sotto l'1%.
sotto_soglia <- grezzi %>% filter(reads > 0 & pct_grezza < 1)
cat("\nDetection residue sotto la soglia dell'1%:", nrow(sotto_soglia), "\n")
cat("Detection piu' piccola trattenuta, per area (%):\n")
print(grezzi %>% filter(reads > 0) %>%
        group_by(Popolazione) %>%
        summarise(minima_pct = round(min(pct_grezza), 2), .groups = "drop"))
# Atteso: 0 detection sotto soglia; minima 1.60% in Slovenia, 5.63% in Croazia.
stopifnot(nrow(sotto_soglia) == 0)

# La soglia resta scritta nel codice (idempotente) per rendere esplicito il
# passaggio metodologico anche a chi rilegge solo lo script.
grezzi <- grezzi %>%
  mutate(reads_filtrate = if_else(pct_grezza < 1, 0, reads)) %>%
  group_by(Sample_ID) %>%
  mutate(totale_filtrato = sum(reads_filtrate)) %>%
  ungroup() %>%
  mutate(RRA_percentuale = if_else(totale_filtrato > 0,
                                   reads_filtrate / totale_filtrato * 100, 0))

# Ogni campione deve sommare a 100
controllo_somme <- grezzi %>%
  group_by(Sample_ID) %>%
  summarise(somma_RRA = sum(RRA_percentuale), .groups = "drop")
stopifnot(all(abs(controllo_somme$somma_RRA - 100) < 1e-6))


# ==============================================================================
# 2. DESCRITTORI DIETETICI: FOO, RRA, PRESENZE ASSOLUTE
# ==============================================================================

rra_pop <- grezzi %>%
  group_by(Popolazione, scientific_name) %>%
  summarise(RRA_percentuale = mean(RRA_percentuale), .groups = "drop")

foo_pop <- grezzi %>%
  group_by(Popolazione, scientific_name) %>%
  summarise(FOO_percentuale = mean(RRA_percentuale > 0) * 100, .groups = "drop")

presenza_pop <- grezzi %>%
  group_by(Popolazione, scientific_name) %>%
  summarise(Presenza_Assoluta = sum(RRA_percentuale > 0), .groups = "drop")

foo_rra_finale <- presenza_pop %>%
  left_join(foo_pop, by = c("Popolazione", "scientific_name")) %>%
  left_join(rra_pop, by = c("Popolazione", "scientific_name")) %>%
  filter(Presenza_Assoluta > 0) %>%
  rename(Preda = scientific_name) %>%
  mutate(FOO_percentuale = round(FOO_percentuale, 2),
         RRA_percentuale = round(RRA_percentuale, 2)) %>%
  arrange(Popolazione, desc(RRA_percentuale))

print(foo_rra_finale, n = Inf)
mostra_tabella(foo_rra_finale, "FOO e RRA per area")
write.csv2(foo_rra_finale, OUT_FOO_RRA, row.names = FALSE)

# Valori attesi in tesi (Tabella 3.2):
#   Slovenia (52): Capreolus 57.69 / 48.48 - Cervus 48.08 / 36.29 - Caprinae 9.62 / 8.87
#   Croazia (51, invariata): Cervus 41.18 / 40.39 - Capreolus 31.37 / 30.20 - Sus 19.61 / 18.84


# ==============================================================================
# 3. FIGURE 3.1 e 3.2 - COMPOSIZIONE DELLA DIETA PER AREA (due pannelli)
# ==============================================================================

df_indici <- read_csv2(OUT_FOO_RRA, show_col_types = FALSE)
stopifnot(all(c("Popolazione", "Preda", "Presenza_Assoluta",
                "FOO_percentuale", "RRA_percentuale") %in% names(df_indici)))
df_indici[is.na(df_indici)] <- 0

# Caprinae resta una categoria a se': non e' mai assegnato ne' al bestiame
# domestico ne' agli ungulati selvatici (Metodi, Sezione 2.5).
df_indici <- df_indici %>%
  mutate(Categoria = case_when(
    Preda %in% c("Capreolus capreolus", "Cervus elaphus",
                 "Sus scrofa", "Rupicapra rupicapra") ~ "wild ungulates",
    Preda %in% c("Ovis aries", "Ovis", "Capra", "Bos")  ~ "domestic livestock",
    Preda == "Caprinae"                                 ~ "unresolved Caprinae",
    TRUE                                                ~ "other taxa"
  ))
df_indici$Categoria <- factor(df_indici$Categoria,
  levels = c("wild ungulates", "domestic livestock",
             "unresolved Caprinae", "other taxa"))

crea_grafico_nazione <- function(dati, nazione) {

  df_nazione <- dati %>%
    filter(Popolazione == nazione,
           FOO_percentuale > 0 | RRA_percentuale > 0)

  ordine_prede <- df_nazione %>% arrange(RRA_percentuale) %>% pull(Preda)
  df_nazione$Preda <- factor(df_nazione$Preda, levels = ordine_prede)

  df_long <- df_nazione %>%
    pivot_longer(cols = c(FOO_percentuale, RRA_percentuale),
                 names_to = "Metrica", values_to = "Valore") %>%
    mutate(Metrica = recode(Metrica,
             "FOO_percentuale" = "A  Frequency of occurrence",
             "RRA_percentuale" = "B  Relative read abundance"))
  df_long$Metrica <- factor(df_long$Metrica,
    levels = c("A  Frequency of occurrence", "B  Relative read abundance"))

  ggplot(df_long, aes(x = Valore, y = Preda, fill = Categoria)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = sprintf("%.1f", Valore)), hjust = -0.2, size = 3.5) +
    facet_wrap(~ Metrica, scales = "free_x") +
    scale_fill_manual(values = COL_CATEGORIA) +
    scale_y_discrete(labels = etichette_taxa_plotmath) +
    scale_x_continuous(limits = c(0, max(df_long$Valore) * 1.15),
                       expand = c(0, 0)) +
    labs(title = paste("Diet composition:", LAB_AREA[[nazione]]),
         subtitle = paste0("n = ", n_distinct(grezzi$Sample_ID[grezzi$Popolazione == nazione]),
                           " samples; Caprinae is not attributed to wild or domestic taxa"),
         x = "Percentage (%)", y = NULL, fill = NULL) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 11, colour = "black"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey35", hjust = 0.5,
                                   margin = margin(b = 12)),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "white", colour = "black",
                                      linewidth = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

plot_slo <- crea_grafico_nazione(df_indici, "Slovenia")
plot_cro <- crea_grafico_nazione(df_indici, "Croazia")
mostra(plot_slo)
ggsave("Grafico_2Pannelli_Slovenia_final.png", plot_slo,
       width = 12, height = 6, dpi = 300)
mostra(plot_cro)
ggsave("Grafico_2Pannelli_Croazia_final.png", plot_cro,
       width = 12, height = 6, dpi = 300)


# ==============================================================================
# 4. FIGURA 3.3 (NUOVA) - FOO E RRA DELLE DUE AREE SUGLI STESSI ASSI
# Risponde al commento "utilizza qualche barplot per una migliore comprensione":
# i due grafici precedenti sono uno per area e non permettono il confronto.
# ==============================================================================

ordine_taxa <- c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa",
                 "Caprinae", "Rupicapra rupicapra", "Bos", "Capra",
                 "Ovis", "Ovis aries", "Lepus")

confronto <- expand_grid(Popolazione = c("Slovenia", "Croazia"),
                         Preda = ordine_taxa) %>%
  left_join(df_indici %>% dplyr::select(Popolazione, Preda,
                                        FOO_percentuale, RRA_percentuale),
            by = c("Popolazione", "Preda")) %>%
  mutate(across(c(FOO_percentuale, RRA_percentuale), ~ replace_na(.x, 0)),
         Preda = factor(Preda, levels = rev(ordine_taxa)),
         Popolazione = factor(Popolazione, levels = c("Slovenia", "Croazia")))

pannello_confronto <- function(dati, colonna, titolo, asse) {
  ggplot(dati, aes(x = .data[[colonna]], y = Preda, fill = Popolazione)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = ifelse(.data[[colonna]] > 0,
                                 format(round(.data[[colonna]], 1), nsmall = 1),
                                 "0")),
              position = position_dodge(width = 0.75),
              hjust = -0.15, size = 3.1) +
    scale_fill_manual(values = COL_AREA, labels = LAB_AREA) +
    scale_y_discrete(labels = etichette_taxa_plotmath) +
    scale_x_continuous(limits = c(0, max(dati[[colonna]]) * 1.2),
                       expand = c(0, 0)) +
    labs(title = titolo, x = asse, y = NULL) +
    tema_tesi
}

fig_confronto <-
  pannello_confronto(confronto, "FOO_percentuale",
                     "A  Frequency of occurrence",
                     "Frequency of occurrence (%)") +
  pannello_confronto(confronto, "RRA_percentuale",
                     "B  Mean relative read abundance",
                     "Mean relative read abundance (%)") +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Diet composition compared between geographic areas",
    subtitle = paste("Caprinae denotes assignments retained at subfamily level",
                     "and is not attributed to wild or domestic taxa"),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) &
  theme(legend.position = "bottom")

mostra(fig_confronto)
ggsave("Figura_Confronto_FOO_RRA.png", fig_confronto,
       width = 12, height = 5.5, dpi = 300)


# ==============================================================================
# 5. COMMUNITY MATRIX (RRA% per singolo campione)
# ==============================================================================

community_matrix <- grezzi %>%
  dplyr::select(Sample_ID, Popolazione, scientific_name, RRA_percentuale) %>%
  pivot_wider(names_from = scientific_name,
              values_from = RRA_percentuale,
              values_fill = 0) %>%
  dplyr::select(Sample_ID, Popolazione, all_of(ordine_taxa))

stopifnot(nrow(community_matrix) == N_TOT)
write.csv2(community_matrix, OUT_COMMUNITY, row.names = FALSE)
cat("\nCommunity matrix:", nrow(community_matrix), "campioni,",
    length(ordine_taxa), "taxa\n")


# ==============================================================================
# 6. FIGURA 3.4 (NUOVA) - OUTPUT DEL METABARCODING
# Risponde al commento "manca tutto il capitolo di descrizione dei dati del
# metabarcoding ottenuto dalle analisi".
# ==============================================================================

per_campione <- grezzi %>%
  group_by(Sample_ID, Popolazione) %>%
  summarise(n_taxa = sum(RRA_percentuale > 0),
            reads  = first(totale_filtrato), .groups = "drop") %>%
  mutate(Popolazione = factor(Popolazione, levels = c("Slovenia", "Croazia")),
         n_taxa = factor(n_taxa, levels = 1:3,
                         labels = c("1 taxon", "2 taxa", "3 taxa")),
         classe_reads = cut(reads,
            breaks = c(0, 10000, 50000, 100000, Inf),
            labels = c("< 10,000", "10,000-50,000",
                       "50,000-100,000", "> 100,000"),
            right = FALSE))

cat("\nTaxa per campione:\n");  print(table(per_campione$Popolazione,
                                            per_campione$n_taxa))
cat("\nReads per campione:\n"); print(table(per_campione$Popolazione,
                                            per_campione$classe_reads))
cat("\nMediana e intervallo dei reads per campione:\n")
print(per_campione %>% group_by(Popolazione) %>%
        summarise(mediana = median(reads), minimo = min(reads),
                  massimo = max(reads), .groups = "drop"))
# Atteso: 83 campioni su 103 con un solo taxon (38 Slovenia, 45 Croazia);
# mediana 8,931 reads in Slovenia (2,110-110,813) e 58,579 in Croazia
# (3,895-189,136).

barre_percentuali <- function(dati, colonna, titolo) {
  dati %>%
    mutate(categoria = .data[[colonna]]) %>%
    count(Popolazione, categoria, name = "n") %>%
    group_by(Popolazione) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = pct, y = fct_rev(factor(categoria)), fill = Popolazione)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.75),
              hjust = -0.15, size = 3.1) +
    scale_fill_manual(values = COL_AREA, labels = LAB_AREA) +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = titolo, x = "Samples in the geographic area (%)", y = NULL) +
    tema_tesi
}

fig_output <-
  barre_percentuali(per_campione, "n_taxa",       "A  Prey taxa detected per sample") +
  barre_percentuali(per_campione, "classe_reads", "B  Prey reads retained per sample") +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Metabarcoding output of the analytical dataset",
    subtitle = paste0("Percentages of the ", N_SLO, " Slovenian and ", N_CRO,
                      " Croatian samples retained after filtering"),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) &
  theme(legend.position = "bottom")

mostra(fig_output)
ggsave("Figura_OutputMetabarcoding.png", fig_output,
       width = 12, height = 4.8, dpi = 300)


# ==============================================================================
# 7. AMPIEZZA DI NICCHIA (LEVINS) E SOVRAPPOSIZIONE (PIANKA)
# ==============================================================================

rra_matrice <- rra_pop %>%
  pivot_wider(names_from = scientific_name,
              values_from = RRA_percentuale, values_fill = 0)

n_categorie <- ncol(rra_matrice) - 1          # 10, non piu' scritto a mano
stopifnot(n_categorie == 10)

levins_finale <- rra_matrice %>%
  rowwise() %>%
  mutate(B   = 1 / sum(c_across(-Popolazione)^2 / 10000),
         B_A = (B - 1) / (n_categorie - 1)) %>%
  ungroup() %>%
  dplyr::select(Popolazione, B_A)
print(levins_finale)      # atteso: Slovenia 0.184, Croazia 0.268
mostra_tabella(levins_finale, "Levins B_A")

p_slo <- as.numeric(rra_matrice[rra_matrice$Popolazione == "Slovenia", -1]) / 100
p_cro <- as.numeric(rra_matrice[rra_matrice$Popolazione == "Croazia",  -1]) / 100
pianka_finale <- sum(p_slo * p_cro) / sqrt(sum(p_slo^2) * sum(p_cro^2))
cat("\nIndice di Pianka:", round(pianka_finale, 3), "\n")   # atteso 0.915


# ==============================================================================
# 8. STATISTICA MULTIVARIATA: HELLINGER, PERMANOVA, BETADISPER, SIMPER
# Tutte e tre le analisi lavorano sulla STESSA matrice trasformata, cosi' il
# test di dispersione valida davvero le assunzioni del PERMANOVA.
# ==============================================================================

df_mv <- read_csv2(OUT_COMMUNITY, show_col_types = FALSE)
metadati       <- df_mv %>% dplyr::select(Sample_ID, Popolazione)
matrice_prede  <- df_mv %>% dplyr::select(all_of(ordine_taxa)) %>%
                   mutate(across(everything(), as.numeric)) %>% as.data.frame()

matrice_hellinger <- decostand(matrice_prede, method = "hellinger")
dist_matrix       <- vegdist(matrice_hellinger, method = "bray")

# --- PERMANOVA ---------------------------------------------------------------
set.seed(123)
permanova_risultato <- adonis2(dist_matrix ~ Popolazione, data = metadati,
                               permutations = 9999, by = "terms")
cat("\n=== PERMANOVA ===\n"); print(permanova_risultato)
mostra_tabella(as.data.frame(permanova_risultato), "PERMANOVA area")
# atteso: R2 = 0.0281, F = 2.92; il p va letto qui (una stima con 9999
# permutazioni in Python da' circa 0.04)

# --- BETADISPER (in inglese) --------------------------------------------------
# Le etichette dei centroidi nel grafico di betadisper sono i livelli del
# fattore: vanno rinominate qui, altrimenti dentro la figura resta "Croazia".
gruppo <- factor(metadati$Popolazione, levels = c("Slovenia", "Croazia"),
                 labels = c("Slovenia", "Croatia"))
dispersion_mod <- betadisper(dist_matrix, group = gruppo)

set.seed(123)
cat("\n=== BETADISPER ===\n")
print(permutest(dispersion_mod, permutations = 999))
print(dispersion_mod$group.distances)
# valori da leggere qui (quelli del 22/09, su 100 campioni, erano F = 3.744,
# p = 0.063; distanze medie 0.486 Slovenia, 0.571 Croazia).
# Nota: vegan avvisa che alcune distanze al quadrato sono negative e le porta a
# zero. E' il motivo per cui questi valori non si riproducono con un calcolo
# fatto a mano che tratti diversamente gli autovalori negativi della matrice di
# Bray-Curtis. I numeri da usare in tesi sono quelli stampati qui.

disegna_betadisper <- function() {
  plot(dispersion_mod, hull = FALSE, ellipse = TRUE,
       main = "Multivariate dispersion of diet composition",
       sub  = "Bray-Curtis dissimilarity on Hellinger-transformed RRA proportions",
       col = unname(COL_AREA),
       lwd = 2, seg.col = "grey80", seg.lwd = 0.5)
  legend("topleft", legend = levels(gruppo),
         col = unname(COL_AREA), pch = 16, bty = "n", cex = 1.1)
}
if (interactive()) disegna_betadisper()          # a schermo
png("Betadisper_Plot_final.png", width = 2000, height = 1600, res = 300)
disegna_betadisper()                              # su file
dev.off()

# --- SIMPER -------------------------------------------------------------------
set.seed(123)
simper_risultato <- simper(matrice_hellinger, group = metadati$Popolazione,
                           permutations = 999)
cat("\n=== SIMPER ===\n"); print(summary(simper_risultato))

simper_summary <- summary(simper_risultato)[[1]]
simper_summary$Taxon <- rownames(simper_summary)
simper_summary <- simper_summary %>%
  mutate(Taxon = gsub("\\.", " ", Taxon),
         contributo_pct = round(average / sum(average) * 100, 1),
         sig_label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                               p < 0.05  ~ "*",   TRUE     ~ ""),
         significativo = p < 0.05) %>%
  arrange(desc(contributo_pct))

plot_simper <- ggplot(simper_summary,
                      aes(x = reorder(Taxon, contributo_pct),
                          y = contributo_pct, fill = significativo)) +
  geom_col() +
  geom_text(aes(label = paste0(contributo_pct, "% ", sig_label)),
            hjust = -0.1, fontface = "bold", size = 3.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#C73E1D", "FALSE" = "#1F7FB5"),
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s.")) +
  expand_limits(y = max(simper_summary$contributo_pct) * 1.15) +
  labs(title = "SIMPER: taxon contributions to between-area dissimilarity",
       subtitle = paste("Bray-Curtis dissimilarity on Hellinger-transformed",
                        "RRA proportions (* p < 0.05, ** p < 0.01, *** p < 0.001)"),
       x = NULL, y = "Contribution to overall dissimilarity (%)") +
  tema_tesi + theme(legend.position = "top")

mostra(plot_simper)
ggsave("Grafico_SIMPER_final.png", plot_simper, width = 10, height = 7, dpi = 300)


# ==============================================================================
# 9. PCoA SUI 10 TAXA (Figura 3.5) E TIPI DI DIETA (Figura 3.6)
#
# PUNTO CENTRALE: la versione precedente calcolava la PCoA sulle 6 categorie
# aggregate e otteneva 52.70% + 27.06%. La figura della tesi riporta
# 52.70% + 26.72%, che e' la versione sui 10 taxa. Qui si usano i 10 taxa.
#
# Secondo punto: 81 campioni su 100 contengono un solo taxon, quindi campioni
# con lo stesso profilo cadono ESATTAMENTE nello stesso punto. I 100 campioni
# occupano 25 posizioni distinte, di cui una con 35 campioni e una con 33. Il
# jitter casuale della versione precedente nascondeva questo fatto; qui la
# dimensione del punto e' il numero di campioni sovrapposti.
# ==============================================================================

pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

eig_positivi <- pcoa_result$eig[pcoa_result$eig > 0]
var_asse1 <- round(100 * pcoa_result$eig[1] / sum(eig_positivi), 2)
var_asse2 <- round(100 * pcoa_result$eig[2] / sum(eig_positivi), 2)
cat("\nPCoA - varianza spiegata: asse 1 =", var_asse1,
    "%, asse 2 =", var_asse2, "%\n")
stopifnot(abs(var_asse1 - 53.35) < 0.01, abs(var_asse2 - 26.40) < 0.01)

site_scores <- as.data.frame(pcoa_result$points)
colnames(site_scores) <- c("PCoA1", "PCoA2")
site_scores$Popolazione <- factor(metadati$Popolazione,
                                  levels = c("Slovenia", "Croazia"))

# Quanti campioni occupano ciascuna posizione
site_points <- site_scores %>%
  mutate(x = round(PCoA1, 8), y = round(PCoA2, 8)) %>%
  count(x, y, Popolazione, name = "n_campioni") %>%
  # scostamento fisso (non casuale) per separare le due aree quando
  # condividono la stessa posizione: dichiarato in didascalia
  mutate(x_plot = x + if_else(Popolazione == "Slovenia", -0.012, 0.012))

# distinct() sulle coordinate non arrotondate conta 69 posizioni: e' rumore in
# virgola mobile, non struttura del dato. Va arrotondato, come per site_points.
cat("\nProfili dietetici distinti:",
    nrow(distinct(round(as.data.frame(matrice_hellinger), 10))), "\n")
cat("Posizioni distinte nel piano dei primi due assi:",
    nrow(distinct(site_points, x, y)), "\n")
cat("Campioni per posizione (prime righe):\n")
print(site_points %>% arrange(desc(n_campioni)) %>% head(8))
# atteso: 26 profili distinti ma 25 posizioni distinte nel piano; 35 campioni
# sulla posizione "solo Cervus elaphus" e 35 su "solo Capreolus capreolus"
# sommando le due aree

species_scores <- as.data.frame(
  wascores(site_scores[, c("PCoA1", "PCoA2")], matrice_prede))
species_scores$Taxon <- rownames(species_scores)

plot_pcoa <- ggplot() +
  stat_ellipse(data = site_scores,
               aes(x = PCoA1, y = PCoA2,
                   colour = Popolazione, fill = Popolazione),
               geom = "polygon", alpha = 0.12, level = 0.95, linewidth = 0.8) +
  geom_point(data = site_points,
             aes(x = x_plot, y = y, colour = Popolazione,
                 fill = Popolazione, size = n_campioni),
             alpha = 0.65, shape = 21, stroke = 0.4) +
  geom_point(data = species_scores, aes(x = PCoA1, y = PCoA2),
             colour = "#B22222", size = 3) +
  geom_text_repel(data = species_scores,
                  aes(x = PCoA1, y = PCoA2, label = Taxon),
                  fontface = "bold.italic", size = 4, seed = 123,
                  max.overlaps = 20) +
  scale_colour_manual(values = COL_AREA, labels = LAB_AREA) +
  scale_fill_manual(values = COL_AREA, labels = LAB_AREA) +
  scale_size_continuous(name = "Samples at the same position",
                        range = c(2.5, 14), breaks = c(1, 5, 10, 20, 35)) +
  labs(colour = NULL, fill = NULL,
       title = "Dietary niche space (Principal Coordinates Analysis)",
       subtitle = paste0("Bray-Curtis dissimilarity on Hellinger-transformed ",
                         "RRA proportions, 10 prey taxa | ",
                         var_asse1, "% + ", var_asse2, "% of positive eigenvalues"),
       x = paste0("Coordinate 1 (", var_asse1, "%)"),
       y = paste0("Coordinate 2 (", var_asse2, "%)")) +
  tema_tesi +
  # tema_tesi azzera i titoli di legenda: qui serve quello della dimensione
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = 9, colour = "grey25"))

mostra(plot_pcoa)
ggsave("Grafico_PCoA_final.png", plot_pcoa, width = 10, height = 8, dpi = 300)

# --- Figura 3.6: tipi di dieta ------------------------------------------------
tipi_dieta <- grezzi %>%
  filter(RRA_percentuale > 0) %>%
  group_by(Sample_ID, Popolazione) %>%
  summarise(n_taxa = n(),
            dominante = scientific_name[which.max(RRA_percentuale)],
            .groups = "drop") %>%
  mutate(tipo = case_when(
    n_taxa > 1 ~ "Mixed (>1 taxon)",
    dominante %in% c("Bos", "Capra", "Ovis", "Ovis aries") ~ "Domestic livestock only",
    TRUE ~ paste(dominante, "only")))

ordine_tipi <- c("Cervus elaphus only", "Capreolus capreolus only",
                 "Mixed (>1 taxon)", "Sus scrofa only",
                 "Caprinae only", "Domestic livestock only")

riepilogo_tipi <- expand_grid(Popolazione = c("Slovenia", "Croazia"),
                              tipo = ordine_tipi) %>%
  left_join(count(tipi_dieta, Popolazione, tipo, name = "n"),
            by = c("Popolazione", "tipo")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(Popolazione) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(tipo = factor(tipo, levels = rev(ordine_tipi)),
         Popolazione = factor(Popolazione, levels = c("Slovenia", "Croazia")))

print(riepilogo_tipi, n = Inf)
mostra_tabella(riepilogo_tipi, "Tipi di dieta")
# atteso: Cervus only 15 SLO / 20 CRO; Capreolus only 20 / 15; misti 14 / 6;
# Sus only 1 / 7; Caprinae only 2 / 1; domestici 0 / 2

etichette_tipi <- function(x) {
  parse(text = vapply(x, function(v) {
    if (v == "Mixed (>1 taxon)")        return("plain('Mixed (>1 taxon)')")
    if (v == "Caprinae only")           return("plain('Caprinae only')")
    if (v == "Domestic livestock only") return("plain('Domestic livestock only')")
    sp <- sub(" only$", "", v)
    sprintf("italic('%s')~plain('only')", sp)
  }, character(1)))
}

plot_tipi <- ggplot(riepilogo_tipi, aes(x = pct, y = tipo, fill = Popolazione)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = ifelse(n > 0, sprintf("%.1f%%", pct), "0")),
            position = position_dodge(width = 0.75), hjust = -0.15, size = 3.2) +
  scale_fill_manual(values = COL_AREA, labels = LAB_AREA) +
  scale_y_discrete(labels = etichette_tipi) +
  scale_x_continuous(limits = c(0, 48), expand = c(0, 0)) +
  labs(title = "Distribution of samples by diet type",
       subtitle = "83 of 103 samples contained a single prey taxon",
       x = "Percentage of samples in the geographic area (%)", y = NULL) +
  tema_tesi + theme(legend.position = "bottom")

mostra(plot_tipi)
ggsave("Figura_TipiDieta.png", plot_tipi, width = 9.5, height = 5.5, dpi = 300)


# ==============================================================================
# 10. MODELLI DIETA-AMBIENTE
# ==============================================================================

ambiente <- read_csv2(FILE_AMBIENTE, locale = locale(decimal_mark = ","),
                      show_col_types = FALSE)

df_env <- df_mv %>%
  inner_join(ambiente %>% dplyr::select(Sample, ELEV_mean, NDVI_mean, NDVI_sd),
             by = c("Sample_ID" = "Sample")) %>%
  rename(Geographic_Area = Popolazione)
stopifnot(nrow(df_env) == N_TOT)   # se si ferma qui, mancano i 3 campioni nel file ambientale

df_env <- df_env %>%
  mutate(ELEV_mean_z = as.numeric(scale(ELEV_mean)),
         NDVI_mean_z = as.numeric(scale(NDVI_mean)),
         NDVI_sd_z   = as.numeric(scale(NDVI_sd)),
         Geographic_Area = relevel(factor(Geographic_Area), ref = "Slovenia"),
         ricchezza = rowSums(across(all_of(ordine_taxa)) > 0),
         Dieta_mista = as.integer(ricchezza > 1))

cat("\nCampioni a preda singola / misti:",
    sum(df_env$Dieta_mista == 0), "/", sum(df_env$Dieta_mista == 1), "\n")
stopifnot(sum(df_env$Dieta_mista) == 20)


# --- 10.1 HABITAT HETEROGENEITY HYPOTHESIS ------------------------------------
modello_hhh <- glm(Dieta_mista ~ NDVI_sd_z + ELEV_mean_z + Geographic_Area,
                   data = df_env, family = binomial(link = "logit"))

cat("\n=== HHH ===\n"); print(summary(modello_hhh))
cat("\nOdds ratio e IC 95% da profilo di verosimiglianza (richiede MASS):\n")
print(exp(cbind(OR = coef(modello_hhh), confint(modello_hhh))))
# valori del 22/09 su 100 campioni, da rileggere: NDVI_sd_z OR = 2.860, IC [1.641, 5.448], p = 0.00049
# anteprima Python 23/09, 103 campioni e NDVI nuovi: NDVI_sd_z OR = 2.487, p = 0.0010 (da confermare in R)

# Test del rapporto di verosimiglianza sull'interazione (riportato in tesi)
modello_hhh_int <- update(modello_hhh, . ~ . + NDVI_sd_z:Geographic_Area)
cat("\nTest dell'interazione NDVI_sd x area:\n")
print(anova(modello_hhh, modello_hhh_int, test = "LRT"))
# valori del 22/09 su 100 campioni, da rileggere: chi2(1) = 0.184, p = 0.668
# anteprima Python 23/09, 103 campioni e NDVI nuovi: chi2(1) = 0.008, p = 0.930

# Curve predette DAL MODELLO RIPORTATO IN TESI, per area, a quota media.
# (La versione precedente usava geom_smooth(), che rifitta un modello
#  univariato su NDVI_sd e quindi disegnava una curva diversa da quella del
#  modello discusso nel testo.)
griglia_hhh <- expand_grid(
    NDVI_sd = seq(min(df_env$NDVI_sd), max(df_env$NDVI_sd), length.out = 200),
    Geographic_Area = factor(c("Slovenia", "Croazia"),
                             levels = levels(df_env$Geographic_Area))) %>%
  mutate(NDVI_sd_z  = (NDVI_sd - mean(df_env$NDVI_sd)) / sd(df_env$NDVI_sd),
         ELEV_mean_z = 0)

pred_hhh <- predict(modello_hhh, newdata = griglia_hhh,
                    type = "link", se.fit = TRUE)
griglia_hhh <- griglia_hhh %>%
  mutate(p  = plogis(pred_hhh$fit),
         lo = plogis(pred_hhh$fit - 1.96 * pred_hhh$se.fit),
         hi = plogis(pred_hhh$fit + 1.96 * pred_hhh$se.fit))

plot_hhh <- ggplot() +
  geom_ribbon(data = griglia_hhh,
              aes(x = NDVI_sd, ymin = lo, ymax = hi, fill = Geographic_Area),
              alpha = 0.18) +
  geom_line(data = griglia_hhh,
            aes(x = NDVI_sd, y = p, colour = Geographic_Area), linewidth = 1) +
  geom_point(data = df_env,
             aes(x = NDVI_sd, y = Dieta_mista, colour = Geographic_Area),
             position = position_jitter(width = 0, height = 0.025),
             size = 2.3, alpha = 0.7) +
  scale_colour_manual(values = COL_AREA, labels = LAB_AREA) +
  scale_fill_manual(values = COL_AREA, labels = LAB_AREA) +
  scale_y_continuous(breaks = c(0, 1),
                     labels = c("Single prey\n(1 taxon)", "Mixed\n(>1 taxon)")) +
  labs(title = "Mixed-prey occurrence and habitat heterogeneity",
       subtitle = paste("Model-adjusted probabilities at mean standardised",
                        "elevation; shaded bands are 95% confidence intervals"),
       x = "Habitat heterogeneity (NDVI standard deviation)", y = NULL) +
  tema_tesi + theme(legend.position = "top")

mostra(plot_hhh)
ggsave("HHH_Logistic_Plot_final.png", plot_hhh, width = 10, height = 6.5, dpi = 300)


# --- 10.2 MODELLI GERARCHICI DI PREDA DOMINANTE -------------------------------
df_sub <- df_env %>%
  filter(`Cervus elaphus` > 0 | `Capreolus capreolus` > 0 | `Sus scrofa` > 0) %>%
  rowwise() %>%
  mutate(dominante = c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa")[
      which.max(c_across(c(`Cervus elaphus`, `Capreolus capreolus`, `Sus scrofa`)))]) %>%
  ungroup()

stopifnot(nrow(df_sub) == 93)
cat("\nPreda dominante:\n"); print(table(df_sub$dominante))
# atteso: Capreolus 42, Cervus 39, Sus 12 (nessun pareggio)

df_sub$is_boar_dominant <- as.integer(df_sub$dominante == "Sus scrofa")
modello_livello1 <- glm(is_boar_dominant ~ ELEV_mean_z + NDVI_mean_z,
                        data = df_sub, family = binomial(link = "logit"))
cat("\n=== LIVELLO 1: cinghiale vs cervidi ===\n")
print(summary(modello_livello1))
print(exp(cbind(OR = coef(modello_livello1), confint(modello_livello1))))
# valori del 22/09 su 100 campioni, da rileggere: ELEV_mean_z OR = 0.381, IC [0.161, 0.739], p = 0.0109
# anteprima Python 23/09, 103 campioni e NDVI nuovi: ELEV_mean_z OR = 0.369, p = 0.0085

df_cervidi <- df_sub %>% filter(dominante != "Sus scrofa")
stopifnot(nrow(df_cervidi) == 81)
df_cervidi$is_cervo_dominant <- as.integer(df_cervidi$dominante == "Cervus elaphus")
modello_livello2 <- glm(is_cervo_dominant ~ ELEV_mean_z + NDVI_mean_z,
                        data = df_cervidi, family = binomial(link = "logit"))
cat("\n=== LIVELLO 2: cervo vs capriolo ===\n")
print(summary(modello_livello2))
print(exp(cbind(OR = coef(modello_livello2), confint(modello_livello2))))
# valori del 22/09 su 100 campioni, da rileggere: ELEV_mean_z OR = 0.809, IC [0.376, 1.716], p = 0.581
# anteprima Python 23/09, 103 campioni e NDVI nuovi: ELEV_mean_z OR = 0.579, p = 0.144

# ATTENZIONE: i modelli sono stati stimati su ELEV_mean_z calcolato sui 103
# campioni (df_env). La griglia di previsione deve usare LE STESSE costanti di
# centratura e scala, non la media e la deviazione standard del sottoinsieme,
# altrimenti la curva disegnata non e' quella del modello riportato in tesi.
ELEV_MU <- mean(df_env$ELEV_mean)
ELEV_SD <- sd(df_env$ELEV_mean)

curva_quota <- function(modello, dati) {
  griglia <- tibble(ELEV_mean = seq(min(dati$ELEV_mean), max(dati$ELEV_mean),
                                    length.out = 200)) %>%
    mutate(ELEV_mean_z = (ELEV_mean - ELEV_MU) / ELEV_SD,
           NDVI_mean_z = 0)
  pr <- predict(modello, newdata = griglia, type = "link", se.fit = TRUE)
  griglia %>% mutate(p  = plogis(pr$fit),
                     lo = plogis(pr$fit - 1.96 * pr$se.fit),
                     hi = plogis(pr$fit + 1.96 * pr$se.fit))
}

pannello_gerarchico <- function(curva, dati, risposta, titolo, sottotitolo, etichette) {
  ggplot() +
    geom_ribbon(data = curva, aes(x = ELEV_mean, ymin = lo, ymax = hi),
                fill = "grey70", alpha = 0.4) +
    geom_line(data = curva, aes(x = ELEV_mean, y = p), linewidth = 1) +
    geom_point(data = dati,
               aes(x = ELEV_mean, y = .data[[risposta]], colour = Geographic_Area),
               position = position_jitter(width = 0, height = 0.03),
               size = 2.2, alpha = 0.7) +
    scale_colour_manual(values = COL_AREA, labels = LAB_AREA) +
    scale_y_continuous(breaks = c(0, 1), labels = etichette) +
    labs(title = titolo, subtitle = sottotitolo, x = "Elevation (m)", y = NULL) +
    tema_tesi
}

# p della quota letti dai modelli, cosi' il sottotitolo non resta indietro
p_quota <- function(modello) {
  p <- summary(modello)$coefficients["ELEV_mean_z", "Pr(>|z|)"]
  paste0("Elevation: p = ", formatC(p, digits = 2, format = "g"))
}

fig_gerarchica <-
  pannello_gerarchico(curva_quota(modello_livello1, df_sub), df_sub,
                      "is_boar_dominant",
                      "A  Wild boar vs cervids", p_quota(modello_livello1),
                      c("Cervid-dominated", "Boar-dominated")) +
  pannello_gerarchico(curva_quota(modello_livello2, df_cervidi), df_cervidi,
                      "is_cervo_dominant",
                      "B  Red deer vs roe deer (within cervids)",
                      p_quota(modello_livello2),
                      c("Roe deer-dominated", "Red deer-dominated")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

mostra(fig_gerarchica)
ggsave("Grafico_Gerarchico_Cinghiale_Cervidi_final.png", fig_gerarchica,
       width = 12, height = 6, dpi = 300)


# --- 10.3 REDUNDANCY ANALYSIS -------------------------------------------------
# La RDA usa le 6 categorie aggregate (i domestici raggruppati) e i predittori
# sulla scala originale non standardizzata, come dichiarato in tesi.
df_rda <- df_env %>%
  mutate(Domestic = `Ovis aries` + Bos + Capra + Ovis)

categorie_rda <- c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa",
                   "Caprinae", "Rupicapra rupicapra", "Domestic")

species_matrix <- df_rda %>% dplyr::select(all_of(categorie_rda)) %>% as.data.frame()
rownames(species_matrix) <- df_rda$Sample_ID
env <- df_rda %>% dplyr::select(ELEV_mean, NDVI_mean, NDVI_sd) %>% as.data.frame()

species_hel <- decostand(species_matrix, method = "hellinger")
modello_rda <- rda(species_hel ~ ELEV_mean + NDVI_mean + NDVI_sd, data = env)

cat("\n=== RDA ===\n"); print(summary(modello_rda))
set.seed(123); print(anova(modello_rda, permutations = 999))
set.seed(123); print(anova(modello_rda, by = "margin", permutations = 999))
print(vif.cca(modello_rda))
# valori del 22/09 su 100 campioni, da rileggere: R2 non aggiustato = 8.65%, RDA1 = 4.09%, RDA2 = 3.04%,
# F(3,96) = 3.031, p = 0.002; VIF 1.014 / 3.371 / 3.368
# anteprima Python 23/09, 103 campioni e NDVI nuovi: R2 = 7.79%, RDA1 = 3.79%, RDA2 = 2.49%,
# F(3,99) = 2.790; VIF 1.024 / 3.839 / 3.793; margini: quota p ~0.02, NDVI_mean p ~0.05,
# NDVI_sd p ~0.05 (p per permutazione: il valore esatto lo da' R)

site_scores_rda <- as.data.frame(scores(modello_rda, display = "sites", scaling = 2))
site_scores_rda$Geographic_Area <- df_rda$Geographic_Area

species_scores_rda <- as.data.frame(scores(modello_rda, display = "species", scaling = 2))
species_scores_rda$Taxon <- rownames(species_scores_rda)

env_scores <- as.data.frame(scores(modello_rda, display = "bp", scaling = 2))
env_scores$Variable <- recode(rownames(env_scores),
                              "ELEV_mean" = "Elevation",
                              "NDVI_mean" = "NDVI mean",
                              "NDVI_sd"   = "NDVI sd")

var_explained <- round(100 * eigenvals(modello_rda) / sum(eigenvals(modello_rda)), 1)
arrow_mult <- 2

plot_rda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(data = site_scores_rda,
             aes(x = RDA1, y = RDA2, colour = Geographic_Area),
             position = position_jitter(width = 0.03, height = 0.03),
             size = 2.5, alpha = 0.35) +
  geom_segment(data = species_scores_rda,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
  geom_text_repel(data = species_scores_rda,
                  aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = Taxon),
                  fontface = "italic", size = 3.2, seed = 42, max.overlaps = 20) +
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = RDA1 * arrow_mult, yend = RDA2 * arrow_mult),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "darkred", linewidth = 0.7) +
  geom_text_repel(data = env_scores,
                  aes(x = RDA1 * arrow_mult * 1.15,
                      y = RDA2 * arrow_mult * 1.15, label = Variable),
                  colour = "darkred", fontface = "bold", size = 3.5,
                  seed = 42, max.overlaps = 20) +
  scale_colour_manual(values = COL_AREA, labels = LAB_AREA) +
  labs(title = "Redundancy analysis: wolf diet and environmental gradients",
       subtitle = paste(N_TOT, "faecal samples, Slovenian and Croatian geographic areas"),
       x = paste0("RDA1 (", var_explained[1], "%)"),
       y = paste0("RDA2 (", var_explained[2], "%)")) +
  tema_tesi + theme(legend.position = "top")

mostra(plot_rda)
ggsave("RDA_triplot_final.png", plot_rda, width = 10, height = 8, dpi = 200)


# ==============================================================================
# 11. ANALISI ESPLORATIVA PER REGIONE (Discussion, Sezione 4.2)
# Su indicazione della relatrice (22/09/2026) l'attribuzione dei campioni ai
# cluster genetici di Snjegota et al. (2021) e il confronto fra regioni sono
# trattati nella Discussion. Le regioni croate sono attribuite dalle coordinate
# dei campioni, seguendo i raggruppamenti delle mappe del report di laboratorio.
# Le liste sono scritte per esteso, cosi' l'attribuzione e' verificabile campione
# per campione.
# ==============================================================================

regioni_croazia <- list(
  "Zumberak"     = c("HRV005", "HRV00C", "HRV00E", "HRV00J", "HRV00K", "HRV012",
                     "HRV014", "HRV016", "HRV017", "HRV01K", "HRV02P", "HRV02U"),
  "Gorski kotar" = c("HRV003", "HRV004", "HRV011", "HRV015", "HRV018", "HRV01C",
                     "HRV01E", "HRV01F", "HRV01H", "HRV01J", "HRV01X", "HRV020",
                     "HRV021", "HRV022", "HRV025", "HRV026", "HRV028", "HRV02A",
                     "HRV02K", "HRV060", "HRV061", "HRV062"),
  "Lika N"       = c("HRV01U", "HRV02E", "HRV02H", "HRV032", "HRV033", "HRV05A",
                     "HRV05C", "HRV05E", "HRV05F", "HRV05H", "HRV05J", "HRV05U"),
  "Lika S"       = c("HRV007", "HRV008", "HRV02J"),   # Lika meridionale, zona di Gracac
  "Dalmatia"     = c("HRV05X"),                       # area del cluster 2
  "HRV027"       = c("HRV027")                        # non attribuibile
)

tab_regioni <- tibble(Sample_ID = unlist(regioni_croazia, use.names = FALSE),
                      Regione   = rep(names(regioni_croazia),
                                      lengths(regioni_croazia)))
stopifnot(nrow(tab_regioni) == 51, !anyDuplicated(tab_regioni$Sample_ID))

individui <- read_csv(FILE_INDIVIDUI, show_col_types = FALSE) %>%
  dplyr::select(Sample, Individual)

df_reg <- df_env %>%
  left_join(tab_regioni, by = "Sample_ID") %>%
  mutate(Regione = if_else(Geographic_Area == "Slovenia", "Slovenia", Regione)) %>%
  left_join(individui, by = c("Sample_ID" = "Sample")) %>%
  rowwise() %>%
  mutate(Preda_dominante = ordine_taxa[which.max(c_across(all_of(ordine_taxa)))]) %>%
  ungroup()
stopifnot(!any(is.na(df_reg$Regione)))

riepilogo_regioni <- df_reg %>%
  group_by(Regione) %>%
  summarise(n = n(),
            individui      = n_distinct(Individual, na.rm = TRUE),  # MSV164 e HRV018 senza genotipo
            quota_media    = round(mean(ELEV_mean)),
            NDVI_sd_medio  = round(mean(NDVI_sd), 3),
            campioni_misti = sum(ricchezza > 1),
            .groups = "drop")
cat("\n=== Riepilogo per regione ===\n"); print(riepilogo_regioni)
mostra_tabella(riepilogo_regioni, "Riepilogo per regione")
# atteso: Slovenia 52 campioni, 22 individui, 14 misti (quota e NDVI_sd da
#         rileggere: cambiano con i 3 campioni recuperati)
#         Gorski kotar 22, 13, 850 m, 0.050, 1 misto
#         Lika N 12, 8, 926 m, 0.056, 3 misti; Lika S 3, 3, 704 m, 0.141, 1 misto
#         Zumberak 12, 3, 740 m, 0.048, 0 misti

cat("\n=== Preda dominante per regione ===\n")
print(table(df_reg$Regione, df_reg$Preda_dominante))
# atteso: Gorski kotar 21 cervo su 22; Zumberak 11 capriolo su 12 (+1 Ovis);
#         Lika N 8 cinghiale su 12; Lika S 2 Caprinae + 1 Capra

# PERMANOVA esplorativa fra le tre regioni croate principali (n = 46)
sel_reg  <- df_reg$Regione %in% c("Zumberak", "Gorski kotar", "Lika N")
dati_reg <- df_reg[sel_reg, ]
mat_reg  <- dati_reg %>% dplyr::select(all_of(ordine_taxa)) %>% as.data.frame()
dist_reg <- vegdist(decostand(mat_reg, method = "hellinger"), method = "bray")

set.seed(123)
permanova_regioni <- adonis2(dist_reg ~ Regione, data = dati_reg,
                             permutations = 9999)
cat("\n=== PERMANOVA fra regioni croate (esplorativa) ===\n")
print(permanova_regioni)
# ottenuto (22/09/2026): R2 = 0.67182, F(2,43) = 44.014, p = 1e-04, cioe' il
# minimo ottenibile con 9999 permutazioni (1/10000).
# ATTENZIONE: i 12 campioni dello Zumberak vengono da 3 individui. L'R2 e'
# gonfiato dalla non indipendenza e va presentato come descrittivo.


# ==============================================================================
# 13. TRE AREE DI STUDIO PER LE MAPPE (Risultati) E FILE PER QGIS
# Le tre aree seguono i raggruppamenti della presentazione del laboratorio
# (De Barba et al., EVMC 2025): Slovenia; Northern Croatia = Zumberak + Gorski
# kotar (+ HRV027, Banovina, che Octenjak et al. 2020 mettono nella regione
# North); Southern Croatia = Lika + Dalmazia. Servono SOLO per le mappe e le
# tabelle descrittive: nessun test statistico usa questa divisione.
# ==============================================================================

leggi_punti <- function(f) {
  # Il file puo' essere stato risalvato da Excel: separatore ";" oppure ",",
  # decimali con il punto oppure con la virgola. Si legge tutto come testo e si
  # convertono le coordinate a mano, cosi' il punto decimale non viene mai
  # scambiato per un separatore delle migliaia (errore del 24/09: 46.14 -> 4614).
  prima_riga <- readLines(f, n = 1, encoding = "UTF-8")
  separatore <- if (grepl(";", prima_riga)) ";" else ","
  a_numero <- function(x) as.numeric(gsub(",", ".", x, fixed = TRUE))
  dati <- read_delim(f, delim = separatore,
                     col_types = cols(.default = col_character()),
                     locale = locale(encoding = "UTF-8")) %>%
    mutate(Latitude = a_numero(Latitude), Longitude = a_numero(Longitude))
  # controllo: tutti i punti devono cadere fra Slovenia e Croazia
  stopifnot(all(between(dati$Latitude, 42, 47.5)),
            all(between(dati$Longitude, 13, 18.5)))
  dati
}
punti <- leggi_punti(FILE_PUNTI)
stopifnot(nrow(punti) == N_TOT, all(punti$Sample %in% df_env$Sample_ID))

branchi <- read_csv(FILE_INDIVIDUI, show_col_types = FALSE) %>%
  dplyr::select(Sample, Pack)

# Categorie delle torte (uguali in tutte le mappe e nella Figura G)
mappa_categorie <- list(
  ROE      = "Capreolus capreolus",
  RED      = "Cervus elaphus",
  BOAR     = "Sus scrofa",
  CHAMOIS  = "Rupicapra rupicapra",
  CAPRINAE = "Caprinae",                              # non risolto, mai assegnato
  DOMESTIC = c("Ovis aries", "Bos", "Capra", "Ovis"),  # taxa assigned to domestic livestock
  LEPUS    = "Lepus"
)
# Palette verificata (tutte le coppie: CVD dE minimo 9.0, visione normale 15.7);
# non usa il blu e il rosso delle due aree geografiche.
COL_TAXA <- c(ROE = "#eda100", RED = "#8c510a", BOAR = "#4a3aa7",
              CHAMOIS = "#00a5a8", CAPRINAE = "#c51b7d", DOMESTIC = "#7570b3",
              LEPUS = "#4d9221")
LAB_TAXA <- parse(text = c("'Roe deer'", "'Red deer'", "'Wild boar'",
                           "'Northern chamois'", "'Caprinae (unresolved)'",
                           "'Domestic livestock'", "italic('Lepus')~plain('sp.')"))

ORDINE_AREE <- c("Slovenia", "Northern Croatia", "Southern Croatia")

aree_studio <- df_reg %>%
  mutate(Study_area = case_when(
           Geographic_Area == "Slovenia"                          ~ "Slovenia",
           Regione %in% c("Zumberak", "Gorski kotar", "HRV027")   ~ "Northern Croatia",
           TRUE                                                   ~ "Southern Croatia"),
         Study_area = factor(Study_area, levels = ORDINE_AREE))
for (k in names(mappa_categorie)) {
  aree_studio[[k]] <- rowSums(as.matrix(aree_studio[, mappa_categorie[[k]], drop = FALSE]))
}
stopifnot(all(abs(rowSums(aree_studio[, names(mappa_categorie)]) - 100) < 1e-6))
print(table(aree_studio$Study_area))
stopifnot(identical(as.vector(table(aree_studio$Study_area)), c(52L, 35L, 16L)))

# --- 13.1 CSV per QGIS: una riga per campione (torte delle mappe di dettaglio) --
torte_campioni <- aree_studio %>%
  left_join(branchi, by = c("Sample_ID" = "Sample")) %>%
  left_join(punti %>% dplyr::select(Sample, Latitude, Longitude),
            by = c("Sample_ID" = "Sample")) %>%
  transmute(Sample = Sample_ID,
            Area   = if_else(Geographic_Area == "Slovenia", "Slovenia", "Croatia"),
            Study  = as.character(Study_area),
            Region = if_else(Geographic_Area == "Slovenia",
                             str_remove(Pack, " pack$"), Regione),
            Latitude, Longitude,
            across(all_of(names(mappa_categorie)), ~ round(.x, 3)),
            NTAXA = ricchezza)
stopifnot(!any(is.na(torte_campioni$Latitude)))
write_csv(torte_campioni, "Torte_campioni_103.csv")

# --- 13.2 CSV per QGIS: una riga per area (torte della mappa generale) --------
torte_aree <- aree_studio %>%
  left_join(punti %>% dplyr::select(Sample, Latitude, Longitude),
            by = c("Sample_ID" = "Sample")) %>%
  group_by(Study = Study_area) %>%
  summarise(N         = n(),
            N_IND     = n_distinct(Individual, na.rm = TRUE),  # MSV164 e HRV018 senza genotipo
            Latitude  = round(mean(Latitude), 5),
            Longitude = round(mean(Longitude), 5),
            across(all_of(names(mappa_categorie)), ~ round(mean(.x), 2)),
            .groups = "drop")
print(torte_aree)
mostra_tabella(torte_aree, "Torte per area di studio")
write_csv(torte_aree, "Torte_aree_3.csv")
# atteso (anteprima Python):
#   Slovenia         52 22 | 48.48 36.29  3.90 2.39  8.87  0.08 0
#   Northern Croatia 35 17 | 35.43 58.85  2.86 0     0     2.86 0
#   Southern Croatia 16 12 | 18.75  0    53.81 0    12.84 13.73 0.87

# --- 13.3 Tabella FOO e RRA per area di studio (10 taxa) ----------------------
tab_3aree <- aree_studio %>%
  dplyr::select(Study_area, all_of(ordine_taxa)) %>%
  pivot_longer(-Study_area, names_to = "Taxon", values_to = "RRA") %>%
  group_by(Study_area, Taxon) %>%
  summarise(FOO = round(100 * mean(RRA > 0), 2),
            RRA = round(mean(RRA), 2), .groups = "drop") %>%
  pivot_wider(names_from = Study_area, values_from = c(FOO, RRA)) %>%
  mutate(Taxon = factor(Taxon, levels = ordine_taxa)) %>%
  arrange(Taxon)
cat("\n=== FOO e RRA per area di studio (descrittivo) ===\n"); print(tab_3aree, width = Inf)
mostra_tabella(tab_3aree, "FOO e RRA per area di studio")
write.csv2(tab_3aree, "Tabella_FOO_RRA_3aree.csv", row.names = FALSE)


# ==============================================================================
# 14. FIGURA G - COMPOSIZIONE DEI SINGOLI CAMPIONI (come nella presentazione)
# Una barra per campione, divisa per categoria secondo la RRA del campione.
# ==============================================================================

comp <- aree_studio %>%
  dplyr::select(Sample_ID, Study_area, all_of(names(mappa_categorie)))
mat_comp <- as.matrix(comp[, names(mappa_categorie)])
comp$dominante <- names(mappa_categorie)[max.col(mat_comp, ties.method = "first")]
comp$quota_dom <- apply(mat_comp, 1, max)
ordine_campioni <- comp %>%
  arrange(Study_area, factor(dominante, levels = names(COL_TAXA)), desc(quota_dom)) %>%
  pull(Sample_ID)

comp_long <- comp %>%
  pivot_longer(all_of(names(mappa_categorie)), names_to = "Categoria",
               values_to = "RRA") %>%
  filter(RRA > 0) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = ordine_campioni),
         Categoria = factor(Categoria, levels = rev(names(COL_TAXA))))

n_per_area <- table(comp$Study_area)
etichette_aree <- setNames(paste0(names(n_per_area), "\n(n = ", n_per_area, ")"),
                           names(n_per_area))

fig_composizione <- ggplot(comp_long, aes(x = Sample_ID, y = RRA, fill = Categoria)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  facet_grid(~ Study_area, scales = "free_x", space = "free_x",
             labeller = labeller(Study_area = etichette_aree)) +
  scale_fill_manual(values = COL_TAXA, breaks = names(COL_TAXA), labels = LAB_TAXA) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Prey composition of individual samples",
       subtitle = paste("Each bar is one faecal sample; segments are relative read",
                        "abundance (proportion of prey sequences, not biomass)"),
       x = "Samples, ordered by dominant prey", y = "Relative read abundance") +
  tema_tesi +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

mostra(fig_composizione)
ggsave("Figura_ComposizioneCampioni.png", fig_composizione,
       width = 13, height = 5.5, dpi = 300)


# ==============================================================================
# 15. TORTE DEL TAXONOMIC COVERAGE PER FAMIGLIA (dati grezzi, non filtrati)
# Si usa count_total di ogni MOTU (totale del laboratorio su tutti i campioni
# della corsa). Le colonne dei singoli campioni NON si usano: nel file sloveno
# l'intestazione delle ultime colonne e' spostata (MSV1MU compare come "0" e c'e'
# una colonna spuria "7524"), mentre count_total e' integro (totale 839,437).
# Nomi dopo la riassegnazione BLAST (scientific_name_edited, dove presente).
# Le sequenze riassegnate al predatore (Canis lupus) sono escluse.
# ==============================================================================

leggi_coverage <- function(f, salta) {
  read_delim(f, delim = ";", skip = salta, col_types = cols(.default = col_character()),
             locale = locale(encoding = "UTF-8")) %>%
    rename_with(~ str_remove(.x, "^﻿")) %>%
    mutate(nome = if ("scientific_name_edited" %in% names(.))
                    coalesce(na_if(str_trim(scientific_name_edited), ""),
                             str_trim(scientific_name))
                  else str_trim(scientific_name),
           reads = as.numeric(count_total)) %>%
    filter(!is.na(reads))
}
cov_slo <- leggi_coverage(FILE_COV_SLO, salta = 1) %>% mutate(Area = "Slovenia")
cov_cro <- leggi_coverage(FILE_COV_CRO, salta = 0) %>% mutate(Area = "Croatia")
stopifnot(sum(cov_slo$reads) == 839437, sum(cov_cro$reads) == 3980866)

famiglia <- function(n) case_when(
  n %in% c("Cervus elaphus", "Capreolus capreolus", "Cervinae", "Dama dama") ~ "Cervidae",
  n %in% c("Caprinae", "Capra", "Bos", "Bos taurus", "Ovis", "Ovis aries",
           "Rupicapra", "Rupicapra rupicapra", "Bovidae")                   ~ "Bovidae",
  n == "Sus scrofa"                                                           ~ "Suidae",
  n %in% c("Vulpes vulpes", "Canidae")                                        ~ "Canidae",
  n == "Lepus"                                                                ~ "Leporidae",
  n %in% c("Gallus gallus", "Phasianinae")                                    ~ "Phasianidae",
  n == "Canis lupus"                                                          ~ "PREDATORE",
  TRUE                                                                        ~ "Not assigned to family")

ORDINE_FAMIGLIE <- c("Cervidae", "Bovidae", "Suidae", "Canidae", "Leporidae",
                     "Phasianidae", "Not assigned to family")
# Colori scelti per le coppie di fette adiacenti nelle due torte (CVD dE >= 15):
# Suidae e Leporidae hanno lo stesso colore di Sus scrofa e Lepus nelle mappe.
COL_FAMIGLIE <- c(Cervidae = "#882255", Bovidae = "#999933", Suidae = "#4a3aa7",
                  Canidae = "#e87ba4", Leporidae = "#4d9221",
                  Phasianidae = "#d95f02", "Not assigned to family" = "#9a9a9a")

coverage <- bind_rows(cov_slo, cov_cro) %>%
  mutate(Famiglia = famiglia(nome))
cat("\nReads escluse perche' riassegnate al predatore:\n")
print(coverage %>% filter(Famiglia == "PREDATORE") %>% count(Area, wt = reads))
cat("Nomi rimasti sopra il livello di famiglia:\n")
print(coverage %>% filter(Famiglia == "Not assigned to family") %>%
        count(Area, nome, wt = reads))

coverage_fam <- coverage %>%
  filter(Famiglia != "PREDATORE") %>%
  group_by(Area, Famiglia) %>%
  summarise(reads = sum(reads), MOTU = n(), .groups = "drop") %>%
  group_by(Area) %>%
  mutate(pct = 100 * reads / sum(reads)) %>%
  ungroup() %>%
  mutate(Famiglia = factor(Famiglia, levels = ORDINE_FAMIGLIE),
         Area = factor(Area, levels = c("Slovenia", "Croatia"))) %>%
  arrange(Area, Famiglia)
cat("\n=== Taxonomic coverage per famiglia (reads, dati grezzi) ===\n")
print(coverage_fam %>% mutate(pct = round(pct, 2)), n = Inf)
mostra_tabella(coverage_fam, "Coverage per famiglia")
write.csv2(coverage_fam %>% mutate(pct = round(pct, 2)),
           "Tabella_coverage_famiglie.csv", row.names = FALSE)
# atteso: Slovenia Cervidae 85.38, Bovidae 8.40, Suidae 5.98, Canidae 0.22,
#         Phasianidae 0.02, non assegnati 0.01 (839,285 reads, 152 del predatore escluse)
#         Croazia Cervidae 62.38, Suidae 20.82, Canidae 12.41, Bovidae 4.30,
#         Leporidae 0.09 (3,980,866 reads)

totali_cov <- coverage_fam %>% group_by(Area) %>% summarise(reads = sum(reads), .groups = "drop")
etichette_cov <- setNames(paste0(totali_cov$Area, "\n",
                                 formatC(totali_cov$reads, format = "d", big.mark = ","),
                                 " reads"),
                          totali_cov$Area)

fig_coverage <- ggplot(coverage_fam, aes(x = 1, y = pct, fill = Famiglia)) +
  geom_col(width = 1, colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(pct >= 3, sprintf("%.1f%%", pct), "")),
            position = position_stack(vjust = 0.5), colour = "white",
            fontface = "bold", size = 3.6) +
  coord_polar(theta = "y", direction = -1) +
  facet_wrap(~ Area, labeller = labeller(Area = etichette_cov)) +
  scale_fill_manual(values = COL_FAMIGLIE, breaks = ORDINE_FAMIGLIE, drop = TRUE) +
  labs(title = "Taxonomic coverage of the unfiltered sequence data",
       subtitle = paste("Share of sequence reads by family, all sequenced samples;",
                        "slices below 3% are listed in the table")) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35"),
        strip.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(), legend.position = "bottom")

mostra(fig_coverage)
ggsave("Figura_TaxonomicCoverage.png", fig_coverage, width = 10, height = 5.5, dpi = 300)


# ==============================================================================
# 16. COVARIATE AMBIENTALI PER AREA GEOGRAFICA (Risultati, tabella descrittiva)
# Sostituisce la colonna della quota della tabella della Sezione 2.1: nei Metodi
# resta il disegno di campionamento (130 campioni), qui i dati usati (103).
# ==============================================================================

tab_covariate <- df_env %>%
  mutate(Area = if_else(Geographic_Area == "Slovenia", "Slovenia", "Croatia")) %>%
  group_by(Area) %>%
  summarise(n = n(),
            across(c(ELEV_mean, NDVI_mean, NDVI_sd),
                   list(min = min, mediana = median, max = max)),
            .groups = "drop")
cat("\n=== Covariate ambientali per area ===\n"); print(tab_covariate, width = Inf)
mostra_tabella(tab_covariate, "Covariate ambientali")
write.csv2(tab_covariate, "Tabella_covariate_ambientali.csv", row.names = FALSE)
# atteso: quota Slovenia 634-1,509 m, Croazia 103-1,257 m (media nel buffer di 2 km)


# ==============================================================================
# 17. FIGURA NDVI: TRE BUFFER DI ESEMPIO (Metodi 2.9)
# Richiede i pacchetti terra e sf e i tre GeoTIFF esportati con
# GEE_export_per_mappe.js. Se mancano, la sezione viene saltata.
# Scelta con regola fissa, non in base alla dieta: NDVI_sd minimo (HRV033),
# mediano (MSV04C) e massimo (MSV17C) fra i 103 campioni.
# ==============================================================================

FILE_NDVI <- c(HRV033 = "NDVI_HRV033_2022.tif",
               MSV04C = "NDVI_MSV04C_2020.tif",
               MSV17C = "NDVI_MSV17C_2020.tif")
RUOLO_NDVI <- c(HRV033 = "minimum", MSV04C = "median", MSV17C = "maximum")

if (all(file.exists(FILE_NDVI)) &&
    requireNamespace("terra", quietly = TRUE) &&
    requireNamespace("sf", quietly = TRUE)) {

  centri <- punti %>%
    filter(Sample %in% names(FILE_NDVI)) %>%
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    sf::st_transform(32633)
  xy_centri <- sf::st_coordinates(centri)
  rownames(xy_centri) <- centri$Sample

  pannelli <- purrr::map_dfr(names(FILE_NDVI), function(s) {
    r <- terra::rast(FILE_NDVI[[s]])
    names(r) <- "NDVI"
    r <- terra::classify(r, cbind(-Inf, -1.0001, NA))   # nodata (-9999) -> NA
    # controllo: media e dev. standard nel buffer devono tornare con GEE
    buf <- sf::st_buffer(centri[centri$Sample == s, ], 2000, nQuadSegs = 16)
    v <- terra::extract(r, terra::vect(buf))$NDVI
    rif <- df_env %>% filter(Sample_ID == s)
    cat(sprintf("%s: NDVI_sd dal GeoTIFF %.4f (GEE %.4f); NDVI_mean %.4f (GEE %.4f)\n",
                s, sd(v, na.rm = TRUE), rif$NDVI_sd, mean(v, na.rm = TRUE), rif$NDVI_mean))
    as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
      mutate(Sample = s,
             x_km = (x - xy_centri[s, "X"]) / 1000,
             y_km = (y - xy_centri[s, "Y"]) / 1000)
  })

  lim_basso <- floor(quantile(pannelli$NDVI, 0.02) * 20) / 20
  lim_alto  <- ceiling(max(pannelli$NDVI) * 20) / 20
  cat("Scala dei colori NDVI:", lim_basso, "-", lim_alto,
      "(valori sotto il limite inferiore mostrati col colore piu' chiaro)\n")

  etichette_ndvi <- sapply(names(FILE_NDVI), function(s) {
    sprintf("%s: NDVI sd = %.3f (%s)", s,
            df_env$NDVI_sd[df_env$Sample_ID == s], RUOLO_NDVI[[s]])
  })
  cerchio <- tibble(t = seq(0, 2 * pi, length.out = 361),
                    x_km = 2 * cos(t), y_km = 2 * sin(t))

  fig_ndvi <- ggplot(pannelli %>% mutate(Sample = factor(Sample, levels = names(FILE_NDVI))),
                     aes(x = x_km, y = y_km)) +
    geom_raster(aes(fill = NDVI)) +
    geom_path(data = cerchio, colour = "black", linewidth = 0.6) +
    geom_point(data = tibble(x_km = 0, y_km = 0), shape = 21, fill = "white",
               colour = "black", size = 2.2, stroke = 0.8) +
    facet_wrap(~ Sample, labeller = labeller(Sample = etichette_ndvi)) +
    scale_fill_gradientn(colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b"),
                         limits = c(lim_basso, lim_alto), oob = scales::squish,
                         name = "NDVI") +
    coord_equal(xlim = c(-2.4, 2.4), ylim = c(-2.4, 2.4), expand = FALSE) +
    labs(title = "Within-buffer variation in NDVI",
         subtitle = paste("Sentinel-2 median composite, 1 June-31 August of the",
                          "collection year; circle = 2-km buffer, point = sampling location"),
         x = "Distance from sampling point (km)", y = "Distance from sampling point (km)") +
    tema_tesi +
    theme(legend.title = element_text(), strip.text = element_text(face = "bold"),
          panel.grid = element_blank())

  mostra(fig_ndvi)
  ggsave("Figura_NDVI_buffer_esempio.png", fig_ndvi, width = 12, height = 4.8, dpi = 300)
} else {
  cat("\nSezione 17 saltata: servono i pacchetti terra e sf e i file",
      paste(FILE_NDVI, collapse = ", "), "nella cartella di lavoro.\n")
}


# ==============================================================================
# 18. AMBIENTE DI ESECUZIONE (per la sezione Software and Reproducibility)
# ==============================================================================
cat("\n\n=== sessionInfo() ===\n")
print(sessionInfo())
# Copia da qui i numeri di versione di R, tidyverse, vegan, MASS, ggplot2 e
# patchwork e inseriscili nella Sezione 2.12 della tesi.
