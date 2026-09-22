# ==============================================================================
# DIET COMPOSITION AND ENVIRONMENTAL CORRELATES OF THE GREY WOLF
# IN SLOVENIA AND CROATIA - SCRIPT UNICO DI ANALISI
#
# Alfio Tomarchio - Universita' di Bologna
# Ultima revisione: 22 settembre 2026
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
FILE_INDIVIDUI <- "individui_e_branchi_100_campioni.csv"  # solo per la sezione 11

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

# Controlli: 100 campioni, 49 + 51
stopifnot(n_distinct(grezzi$Sample_ID) == 100)
stopifnot(n_distinct(grezzi$Sample_ID[grezzi$Popolazione == "Slovenia"]) == 49)
stopifnot(n_distinct(grezzi$Sample_ID[grezzi$Popolazione == "Croazia"])  == 51)

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
write.csv2(foo_rra_finale, OUT_FOO_RRA, row.names = FALSE)

# Valori attesi in tesi (Tabella 3.2):
#   Slovenia: Capreolus 55.10 / 45.41 - Cervus 48.98 / 38.42 - Caprinae 10.20 / 9.41
#   Croazia : Cervus 41.18 / 40.39 - Capreolus 31.37 / 30.20 - Sus 19.61 / 18.84


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
             "FOO_percentuale" = "Frequency of Occurrence",
             "RRA_percentuale" = "Relative Read Abundance"))
  df_long$Metrica <- factor(df_long$Metrica,
    levels = c("Frequency of Occurrence", "Relative Read Abundance"))

  ggplot(df_long, aes(x = Valore, y = Preda, fill = Categoria)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = round(Valore, 1)), hjust = -0.2, size = 3.5) +
    facet_wrap(~ Metrica, scales = "free_x") +
    scale_fill_manual(values = COL_CATEGORIA) +
    scale_y_discrete(labels = etichette_taxa_plotmath) +
    scale_x_continuous(limits = c(0, max(df_long$Valore) * 1.15),
                       expand = c(0, 0)) +
    labs(title = paste("Diet composition:", LAB_AREA[[nazione]]),
         x = "Percentage (%)", y = NULL, fill = NULL) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 11, colour = "black"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                margin = margin(b = 15)),
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
ggsave("Grafico_2Pannelli_Slovenia_final.png", plot_slo,
       width = 12, height = 6, dpi = 300)
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

stopifnot(nrow(community_matrix) == 100)
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
# Atteso: 81 campioni su 100 con un solo taxon (36 Slovenia, 45 Croazia);
# mediana 8,647 reads in Slovenia (2,110-110,813) e 58,579 in Croazia
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
    subtitle = paste("Percentages of the 49 Slovenian and 51 Croatian samples",
                     "retained after filtering"),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) &
  theme(legend.position = "bottom")

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
print(levins_finale)      # atteso: Slovenia 0.193, Croazia 0.268

p_slo <- as.numeric(rra_matrice[rra_matrice$Popolazione == "Slovenia", -1]) / 100
p_cro <- as.numeric(rra_matrice[rra_matrice$Popolazione == "Croazia",  -1]) / 100
pianka_finale <- sum(p_slo * p_cro) / sqrt(sum(p_slo^2) * sum(p_cro^2))
cat("\nIndice di Pianka:", round(pianka_finale, 3), "\n")   # atteso 0.930


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
# atteso: R2 = 0.02349, F = 2.357, p = 0.0731

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
# atteso: F = 3.744, p = 0.063; distanze medie 0.486 Slovenia, 0.571 Croazia.
# Nota: vegan avvisa che alcune distanze al quadrato sono negative e le porta a
# zero. E' il motivo per cui questi valori non si riproducono con un calcolo
# fatto a mano che tratti diversamente gli autovalori negativi della matrice di
# Bray-Curtis. I numeri da usare in tesi sono quelli stampati qui.

png("Betadisper_Plot_final.png", width = 2000, height = 1600, res = 300)
plot(dispersion_mod, hull = FALSE, ellipse = TRUE,
     main = "Multivariate dispersion of diet composition",
     sub  = "Bray-Curtis dissimilarity on Hellinger-transformed RRA proportions",
     col = unname(COL_AREA),
     lwd = 2, seg.col = "grey80", seg.lwd = 0.5)
legend("topleft", legend = levels(gruppo),
       col = unname(COL_AREA), pch = 16, bty = "n", cex = 1.1)
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
stopifnot(abs(var_asse1 - 52.70) < 0.01, abs(var_asse2 - 26.72) < 0.01)

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
# atteso: 25 profili distinti ma 24 posizioni distinte nel piano; 35 campioni
# sulla posizione "solo Cervus elaphus" e 33 su "solo Capreolus capreolus"
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
# atteso: Cervus only 15 SLO / 20 CRO; Capreolus only 18 / 15; misti 13 / 6;
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
       subtitle = "81 of 100 samples contained a single prey taxon",
       x = "Percentage of samples in the geographic area (%)", y = NULL) +
  tema_tesi + theme(legend.position = "bottom")

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
stopifnot(nrow(df_env) == 100)

df_env <- df_env %>%
  mutate(ELEV_mean_z = as.numeric(scale(ELEV_mean)),
         NDVI_mean_z = as.numeric(scale(NDVI_mean)),
         NDVI_sd_z   = as.numeric(scale(NDVI_sd)),
         Geographic_Area = relevel(factor(Geographic_Area), ref = "Slovenia"),
         ricchezza = rowSums(across(all_of(ordine_taxa)) > 0),
         Dieta_mista = as.integer(ricchezza > 1))

cat("\nCampioni a preda singola / misti:",
    sum(df_env$Dieta_mista == 0), "/", sum(df_env$Dieta_mista == 1), "\n")
stopifnot(sum(df_env$Dieta_mista) == 19)


# --- 10.1 HABITAT HETEROGENEITY HYPOTHESIS ------------------------------------
modello_hhh <- glm(Dieta_mista ~ NDVI_sd_z + ELEV_mean_z + Geographic_Area,
                   data = df_env, family = binomial(link = "logit"))

cat("\n=== HHH ===\n"); print(summary(modello_hhh))
cat("\nOdds ratio e IC 95% da profilo di verosimiglianza (richiede MASS):\n")
print(exp(cbind(OR = coef(modello_hhh), confint(modello_hhh))))
# atteso: NDVI_sd_z OR = 2.860, IC [1.641, 5.448], p = 0.00049

# Test del rapporto di verosimiglianza sull'interazione (riportato in tesi)
modello_hhh_int <- update(modello_hhh, . ~ . + NDVI_sd_z:Geographic_Area)
cat("\nTest dell'interazione NDVI_sd x area:\n")
print(anova(modello_hhh, modello_hhh_int, test = "LRT"))
# atteso: chi2(1) = 0.184, p = 0.668

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

ggsave("HHH_Logistic_Plot_final.png", plot_hhh, width = 10, height = 6.5, dpi = 300)


# --- 10.2 MODELLI GERARCHICI DI PREDA DOMINANTE -------------------------------
df_sub <- df_env %>%
  filter(`Cervus elaphus` > 0 | `Capreolus capreolus` > 0 | `Sus scrofa` > 0) %>%
  rowwise() %>%
  mutate(dominante = c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa")[
      which.max(c_across(c(`Cervus elaphus`, `Capreolus capreolus`, `Sus scrofa`)))]) %>%
  ungroup()

stopifnot(nrow(df_sub) == 90)
cat("\nPreda dominante:\n"); print(table(df_sub$dominante))
# atteso: Capreolus 39, Cervus 39, Sus 12 (nessun pareggio)

df_sub$is_boar_dominant <- as.integer(df_sub$dominante == "Sus scrofa")
modello_livello1 <- glm(is_boar_dominant ~ ELEV_mean_z + NDVI_mean_z,
                        data = df_sub, family = binomial(link = "logit"))
cat("\n=== LIVELLO 1: cinghiale vs cervidi ===\n")
print(summary(modello_livello1))
print(exp(cbind(OR = coef(modello_livello1), confint(modello_livello1))))
# atteso: ELEV_mean_z OR = 0.381, IC [0.161, 0.739], p = 0.0109

df_cervidi <- df_sub %>% filter(dominante != "Sus scrofa")
stopifnot(nrow(df_cervidi) == 78)
df_cervidi$is_cervo_dominant <- as.integer(df_cervidi$dominante == "Cervus elaphus")
modello_livello2 <- glm(is_cervo_dominant ~ ELEV_mean_z + NDVI_mean_z,
                        data = df_cervidi, family = binomial(link = "logit"))
cat("\n=== LIVELLO 2: cervo vs capriolo ===\n")
print(summary(modello_livello2))
print(exp(cbind(OR = coef(modello_livello2), confint(modello_livello2))))
# atteso: ELEV_mean_z OR = 0.809, IC [0.376, 1.716], p = 0.581

# ATTENZIONE: i modelli sono stati stimati su ELEV_mean_z calcolato sui 100
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

fig_gerarchica <-
  pannello_gerarchico(curva_quota(modello_livello1, df_sub), df_sub,
                      "is_boar_dominant",
                      "A  Wild boar vs cervids", "Elevation: p = 0.011",
                      c("Cervid-dominated", "Boar-dominated")) +
  pannello_gerarchico(curva_quota(modello_livello2, df_cervidi), df_cervidi,
                      "is_cervo_dominant",
                      "B  Red deer vs roe deer (within cervids)",
                      "Elevation: p = 0.581",
                      c("Roe deer-dominated", "Red deer-dominated")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

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
# atteso: R2 non aggiustato = 8.65%, RDA1 = 4.09%, RDA2 = 3.04%,
# F(3,96) = 3.031, p = 0.002; VIF 1.014 / 3.371 / 3.368

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
       subtitle = "100 faecal samples, Slovenian and Croatian geographic areas",
       x = paste0("RDA1 (", var_explained[1], "%)"),
       y = paste0("RDA2 (", var_explained[2], "%)")) +
  tema_tesi + theme(legend.position = "top")

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
            individui      = n_distinct(Individual),
            quota_media    = round(mean(ELEV_mean)),
            NDVI_sd_medio  = round(mean(NDVI_sd), 3),
            campioni_misti = sum(ricchezza > 1),
            .groups = "drop")
cat("\n=== Riepilogo per regione ===\n"); print(riepilogo_regioni)
# atteso: Slovenia 49 campioni, 21 individui, 1047 m, 0.083, 13 misti
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
# atteso: R2 = 0.672, F = 44.0; p <= 0.0002 (il valore esatto dipende dalle
# permutazioni: riporta in tesi quello stampato qui).
# ATTENZIONE: i 12 campioni dello Zumberak vengono da 3 individui. L'R2 e'
# gonfiato dalla non indipendenza e va presentato come descrittivo.


# ==============================================================================
# 12. AMBIENTE DI ESECUZIONE (per la sezione Software and Reproducibility)
# ==============================================================================
cat("\n\n=== sessionInfo() ===\n")
print(sessionInfo())
# Copia da qui i numeri di versione di R, tidyverse, vegan, MASS, ggplot2 e
# patchwork e inseriscili nella Sezione 2.12 della tesi.
