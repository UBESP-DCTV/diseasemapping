################################################################################
############# SCRIPT A SCOPO ILLUSTRATIVO PER IL DISEASE MAPPING ###############
################################################################################
# Prima di eseguire questo script per la prima volta:
#   1. Aprire il progetto tramite il file .Rproj
#   2. Eseguire in console: renv::restore()
# Questo installerà automaticamente le versioni esatte di tutti i pacchetti
# necessari, registrate in renv.lock

################################################################################
# SETUP ------------------------------------------------------------------------

# Settare seed (per risultati riproducibili)
set.seed(123)
# Il seed garantisce la riproducibilità dei risultati MCMC (che sono stocastici)

# Librerie 
library(sf)
library(tidyverse)
library(nimble)
library(coda)
library(spdep)
library(patchwork)
library(ggmcmc)
library(RColorBrewer)
library(classInt)
library(INLA)
library(gt)
library(writexl)

# Tutti i pacchetti sopra sono necessari per: gestione dati spaziali (sf, spdep, INLA),
# manipolazione dati (dplyr, tidyverse), modelli bayesiani (nimble, coda, ggmcmc),
# visualizzazione (ggplot2, patchwork, RColorBrewer, gt) e classificazione dei break (classInt)

# Se non esiste già una cartella output all'interno del progetto, ne crea una
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/diagnostics")) dir.create("output/diagnostics") #per check di convergenza dei modelli


#funzione per generare le mappe
plot_disease_map <- function(data, fill_var, title, palette) {
  ggplot() + 
    geom_sf(data = data, aes(fill = .data[[fill_var]]), col = "black", show.legend = TRUE) + 
    ggtitle(title) + 
    theme_bw() +
    theme(
      plot.title = element_text(size = 12),
      legend.key.size = unit(0.5, "cm"),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.text = element_text(size = 10),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    scale_fill_manual(
      values = palette,
      name = NULL,
      drop = FALSE,
      guide = guide_legend(override.aes = list(colour = NA))
    ) + 
    coord_sf(datum = NULL)
}

# funzione per creare le etichette delle legende nelle mappe 
make_labels <- function(breaks, decimals = 2) {
  b <- round(breaks, decimals)
  n <- length(b) - 1
  labels <- character(n)
  labels[1] <- paste0("< ", b[2])
  if (n > 2) for (i in 2:(n - 1)) labels[i] <- paste0(b[i], " - ", b[i + 1])
  labels[n] <- paste0("> ", b[n])
  labels
}

# Shp file 
# shape file dei comuni della regione x (tutti gli shape files delle unità amministrative
# italiane possono essere scaricati qui https://www.istat.it/notizia/confini-delle-unita-amministrative-a-fini-statistici-al-1-gennaio-2018-2/)
comuni <- read_sf("data/regionex.shp") %>% rename("COD_COM"=PRO_COM) %>% mutate(COD_COM=as.numeric(COD_COM))
N <- nrow(comuni)
nomi_comuni <- comuni %>% sf::st_drop_geometry() %>% dplyr::select(COD_COM, COMUNE)
# NOTA: sempre tenere nella stessa cartella i file shp e shx per una corretta lettura
 
# Importare il dataset sui casi di tumore nei comuni della regione x
data <- readxl::read_excel("data/data.xlsx")
data %>% head(12)
# il dataset contiene:
# x righe = (numero comuni * numero classi età) 
# 4 colonne: id comune, classe di età, numero di casi, popolazione

# NOTA: # Le stesse analisi possono essere condotte separatamente per uomini e donne, 
# se i dati a disposizione lo permettono.
# es. data_M <- data %>% filter(sex=="M")
#     data_F <- data %>% filter(sex=="F")
# Usare i due dataframe generati per svolgere le analisi separatamente.

################################################################################
# ANALISI ----------------------------------------------------------------------
 
# STANDARDIZZAZIONE INDIRETTA E STANDARDISED MORTALITY RATIO (SMR) -------------
 
# Nota: le classi di età utilizzate sono quinquennali (0-4, 5-9, ..., 80-84, 85+)
# Step 1: estrarre tassi=decessi/popolazione collassando sulla classe d'età 
rates <- data %>% 
  group_by(age_grp) %>%
  summarise(
    total_deaths = sum(cases),
    total_pop = sum(pop),
    rate = total_deaths / total_pop
  )
# Calcolo dei tassi di riferimento specifici per età (standard "interno":
# calcolati sull'intera popolazione in studio, sommando su tutti i comuni)

# Step 2: Aggiungere la colonna tassi al dataset e calcolare attesi=popolazione*tasso
data <- data %>%
  left_join(rates %>% dplyr::select(age_grp, rate), by = "age_grp") %>%
  mutate(expected = pop * rate)
# Numero di casi attesi in ciascun comune/classe d'età, calcolato applicando 
# il tasso di riferimento (rate) alla popolazione locale (metodo della 
# standardizzazione indiretta)

# Step 3: Calcolo SMR=osservati/attesi collassando sul comune
smr <- data %>%
  group_by(COD_COM) %>%
  summarise(
    population=sum(pop),
    observed = sum(cases),
    expected = sum(expected),
    SMR = observed / expected,
    SMR_low  = qchisq(0.025, df = 2*observed) / (2*expected), #intervallo al 95%, modificabile
    SMR_high =qchisq(0.975, df = 2*(observed+1)) / (2*expected),
    p_value = pmin(1, 2 * pmin(ppois(observed, expected), 1 - ppois(observed - 1, expected)))
  ) %>% 
  ungroup() %>% 
  mutate(q_value = p.adjust(p_value, method = "fdr")) #q-value di Benjamini-Hochberg
# Collassando su tutte le classi d'età si ottiene, per ciascun comune,
# il totale osservato, il totale atteso e il loro rapporto (SMR)

# Allineamento esplicito di smr all'ordine di comuni 
smr <- comuni %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(COD_COM) %>% 
  left_join(smr, by = "COD_COM")

# Step 4: Tabella, Istogramma e Funnel plot degli SMR

# Genera tabella con: nome comune, casi osservati, attesi, SMR, intervallo di confidenza
table_SMR <- nomi_comuni %>%
  left_join(smr, by = "COD_COM") %>% 
  dplyr::select(COMUNE, observed, expected, SMR, SMR_low, SMR_high) %>% 
  mutate(expected=round(expected,2), SMR=round(SMR,2),
         SMR_low=round(SMR_low,2), SMR_high=round(SMR_high, 2)) %>% 
  gt()
print(table_SMR)
#Salva la tabella in formato Word nella cartella output
f_table_SMR <- "output/table_SMR.docx"
if (file.exists(f_table_SMR)) stopifnot(file.remove(f_table_SMR))  # se fallisce, il file è aperto in Word: chiudilo
gtsave(table_SMR, f_table_SMR)

# Genera istogramma degli SMR
hist_SMR <- smr %>% ggplot()+
  geom_histogram(mapping = aes(x=SMR), bins = 30, fill = "grey70", col = "black")+
  labs(y="")
print(hist_SMR)
#Salva l'istogramma nella cartella output (formato modificabile)
ggsave("output/hist_SMR.tiff", hist_SMR, width=4, height=3)

# Genera funnel plot con p-value e q-value
band <- local({
  #si può essere più tolleranti sulle bande del q-value, es. qui 0.20 
  p_crit <- max(c(0, sort(smr$p_value)[sort(smr$p_value) <= seq_along(smr$p_value)/nrow(smr) * 0.20])) 
  E <- seq(min(smr$expected), max(smr$expected), length.out = 200)
  data.frame(E,
             low_p = qchisq(0.025, 2*E)/(2*E),         high_p = qchisq(0.975, 2*(E+1))/(2*E),
             low_q = qchisq(p_crit/2, 2*E)/(2*E),      high_q = qchisq(1 - p_crit/2, 2*(E+1))/(2*E))
})

funnel_plot <- ggplot() +
  geom_line(data = band, aes(E, low_p,  linetype = "p-value"), col = "grey40") +
  geom_line(data = band, aes(E, high_p, linetype = "p-value"), col = "grey40") +
  geom_line(data = band, aes(E, low_q,  linetype = "q-value"), col = "black") +
  geom_line(data = band, aes(E, high_q, linetype = "q-value"), col = "black") +
  geom_hline(yintercept = 1, col = "red") +
  geom_point(data = smr, aes(expected, SMR))+
  coord_cartesian(ylim = c(0, quantile(smr$SMR, 0.99, na.rm = TRUE))) +
  scale_linetype_manual(values = c("p-value" = "dashed", "q-value" = "solid"), name = NULL) +
  labs(x = "Expected", y = "SMR") +
  theme_bw()
print(funnel_plot)
#Salva il funnel plot nella cartella output (formato modificabile)
ggsave("output/funnel_plot.tiff", funnel_plot, width=6, height=4)

# Step 5: Mappare gli SMR
com_dat <- comuni %>%
  left_join(smr, by = "COD_COM")

# break assoluti
abs.brks <- c(0, 0.5, 0.75, 1, 1.25, Inf)
abs.brks.labels <- c("< 0.50", "0.50 - 0.75", "0.75 - 1.00", "1.00 - 1.25", "> 1.25")
# SMR=1 indica rischio pari allo standard di riferimento; 
# SMR>1 indica un eccesso di rischio locale, SMR<1 un rischio inferiore
# break relativi
rel.brks <- classIntervals(smr$SMR, n = 5, style = "quantile", intervalClosure = "right")
# NOTA: I break assoluti definiscono intervalli fissi, consentendo di confrontare più mappe
# (ad esempio, mappe nel tempo o da diversi modelli) con la stessa scala. 
# I break relativi, invece, sono calcolati specificamente sul dataset corrente
# e sono utili per evidenziare meglio le zone con valori estremi, 
# ma non permettono un confronto diretto tra mappe diverse.

# genera palette in scala di grigi con 5 tonalità (costumizzabile)
palette <- brewer.pal(n = 5, name = "Greys")

com_dat <- com_dat |>
  dplyr::mutate(abs.brks = cut(SMR, abs.brks, include.lowest = TRUE, labels = abs.brks.labels),
                rel.brks = cut(SMR, rel.brks$brks, include.lowest = TRUE, labels = make_labels(rel.brks$brks)))

# Richiama la funzione creata a inizio script per generare tutte le mappe
#break assoluti
SMR_map_abs <- plot_disease_map(com_dat, "abs.brks", "SMR", palette)
print(SMR_map_abs)
#break relativi
SMR_map_rel <- plot_disease_map(com_dat, "rel.brks", "SMR", palette)
print(SMR_map_rel)
#Salva la mappa nella cartella output (formato modificabile)
ggsave("output/SMR_map_abs.tiff", SMR_map_abs, width=5, height=5, bg="white")
ggsave("output/SMR_map_rel.tiff", SMR_map_rel, width = 5, height = 5, bg="white")

# NOTA: quando si tratta di casi incidenti (come in questo caso), gli SMR
# diventano SIR (Standardised Incidence Ratio), e in caso di ricoveri ospedalieri SHR
# (standardised hospitalised ratio)

# MODELLI BAYESIANI GERARCHICI --------------------------------------------
# NOTA: Esistono diversi software e pacchetti per modellare dati spaziali con approccio bayesiano.
# I risultati possono variare leggermente a seconda degli algoritmi implementati. 
# In questo esempio verrà utilizzato il pacchetto Nimble.

# Definizione matrice delle adiacenze 
com_nb <- poly2nb(comuni)
adiacenze <- nb2WB(com_nb)
sumNumNeigh <- sum(adiacenze$num)

# Dati e costanti in formato per nimble
data_mod = list(O = smr$observed, E = smr$expected)

# scale (usato solo nel BYM2) normalizza la varianza del grafo di adiacenza
W.scale <- nb2mat(com_nb, zero.policy = TRUE, style = "B")           
W.scale <- -W.scale
diag(W.scale) <- abs(apply(W.scale, 1, sum))
Q = inla.scale.model(W.scale, constr=list(A=matrix(1, nrow=1, ncol=nrow(W.scale)), e=0))
scale = exp((1/nrow(W.scale))*sum(log(1/diag(Q))))

constants <- list(
  N = N, L = sumNumNeigh,
  adj = adiacenze[["adj"]], num = adiacenze[["num"]],
  weights = adiacenze[["weights"]],
  scale = scale
)

### MODELLO POISSON-GAMMA ###

# Step 1: Definizione del codice per il modello
PGcode <- nimbleCode(
  {
    for (i in 1:N) {
      O[i] ~ dpois(lambda[i]) #O[i]: numero di casi osservati nel comune i, modellato 
      #come Poisson poichè si tratta di conteggi di eventi rari
      lambda[i] <- theta[i]*E[i]  #lambda=RR*attesi
      #a seconda della parametrizzazione del modello si può scegliere se usare
      # i casi attesi (rischio relativo) o la popolazione (tasso di rischio) come offset
      
      theta[i] ~ dgamma(alpha, beta) # theta[i] è il rischio relativo (RR) del comune i: 
      #la Gamma è la prior coniugata della Poisson, il che rende il modello analiticamente 
      #trattabile ed è la ragione del nome "Poisson-Gamma"
      
      post[i] <- step(theta[i]-1) #probabilità a posteriori che il RR>1 
      #(ovvero che ci sia un eccesso di rischio)
    }
    # Iperparametri della Gamma stimati dai dati stessi (approccio empirico-bayesiano 
    #esteso, dato che qui sono trattati come nodi stocastici e non fissati a priori)
    alpha ~ dexp(0.1) #hyperprior
    beta ~ dexp(0.1) #hyperprior
  }
)


# Step 2: Valori iniziali per i nodi stocastici 
inits_PG <- list(alpha=1.60233, beta=3.59906, theta = rgamma(N, shape = 1.60233, scale = 3.59906))
# NOTA: qui sono stati utilizzati i valori standard del modello PG, ma possono 
# essere utilizzati altri valori

# Step 3: Iterazioni MCMC
ni <- 10000 #numero di iterazioni 
nt <- 10 #thinning interval (si considera una iterazione ogni 10) 
nb <- 5000 #burn-in (scartare prime 5000 iterazioni) 
nc <- 1 #numero di catene di initial values 
# (per un check più accurato della convergenza usare più catene)

# Step 4: Fit del modello
modelPG <- nimbleMCMC(
  code = PGcode, data = data_mod, constants = constants, inits = inits_PG,
  monitors = c("theta", "post"), # con monitors definiamo i parametri che vogliamo monitorare
  niter = ni, nburnin = nb, thin = nt, nchains = nc, 
  setSeed = 123, progressBar = FALSE,
  samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC =TRUE
)

# Il WAIC (Watanabe Akaike Information Criterion) permette di confrontare 
# modelli diversi: valori più bassi indicano un migliore compromesso 
# tra bontà di adattamento e complessità del modello
modelPG$WAIC$WAIC

# Dataframe con i parametri monitorati
results_PG <- data.frame(
  COD_COM = smr$COD_COM,
  rr = modelPG$summary[paste0("theta[", 1:N, "]"), "Mean"],
  #intervallo di credibilità al 95%
  cri_low  = modelPG$summary[paste0("theta[", 1:N, "]"), "95%CI_low"], 
  cri_high = modelPG$summary[paste0("theta[", 1:N, "]"), "95%CI_upp"],
  post = modelPG$summary[paste0("post[", 1:N, "]"), "Mean"] 
)
# NOTA: il modello restituisce, volendo, tutti i valori percentili

# Strumenti di diagnostica della convergenza del modello
# Traceplot
# Il traceplot verifica visivamente la convergenza delle catene MCMC 
# (assenza di trend, buon "mixing"); il plot di autocorrelazione verifica 
# quanto le iterazioni successive sono correlate tra loro 
# (un'alta autocorrelazione richiederebbe più thinning o più iterazioni)
sample_ids <- sample(1:N, min(12, N)) #selezioniamo 12 comuni casuali per cui mostrare i grafici
ggs(modelPG$samples) %>% filter(Parameter %in% paste0("theta[", sample_ids, "]"))  %>% 
  ggs_traceplot() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 4) +
  theme_bw() -> PG_traceplot
# Plot autocorrelazione
ggs(modelPG$samples) %>% filter(Parameter %in% paste0("theta[", sample_ids, "]"))  %>% 
  ggs_autocorrelation()+
  facet_wrap(~ Parameter, scales = "free_y", ncol = 4) +
  theme_bw() -> PG_autocorr

print(PG_traceplot)
print(PG_autocorr)

#Salva plot nella cartella diagnostics dentro la cartella output
ggsave("output/diagnostics/PG_traceplot.tiff", PG_traceplot, width=8, height=6)
ggsave("output/diagnostics/PG_autocorr.tiff", PG_autocorr, width=8, height=6)

# Step 5: Mappare i risultati
com_dat <- comuni %>%
  left_join(results_PG, by = "COD_COM")

# break assoluti
# per gli RR riutilizziamo gli stessi break assoluti degli SMR
abs.brks.post <- c(0, 0.7, 0.9, 0.95, 0.99, 1) #break per P(RR>1) 
abs.brks.post.labels <- c("< 0.70", "0.70 - 0.90", "0.90 - 0.95", "0.95 - 0.99", "> 0.99")
# break relativi
rel.brks <- classIntervals(results_PG$rr, n = 5, style = "quantile", intervalClosure = "right")

com_dat <- com_dat |>
  dplyr::mutate(abs.brks = cut(rr, abs.brks, include.lowest = TRUE, labels = abs.brks.labels),
                rel.brks = cut(rr, rel.brks$brks, include.lowest = TRUE, labels = make_labels(rel.brks$brks)),
                abs.brks.post = cut(post, abs.brks.post, include.lowest = TRUE, labels = abs.brks.post.labels))

# Mappa dei rischi relativi (o tassi di rischio, a seconda della parametrizzazione)
#break assoluti
PG_map_abs <- plot_disease_map(com_dat, "abs.brks", "RR from Poisson-Gamma Model", palette)
print(PG_map_abs) 
#break relativi
PG_map_rel <- plot_disease_map(com_dat, "rel.brks", "RR from Poisson-Gamma Model", palette)
print(PG_map_rel)

# Mappa delle probabilita a posteriori
PGpost_map <- plot_disease_map(com_dat, "abs.brks.post", "P(RR>1) from Poisson-Gamma Model", palette)
print(PGpost_map)

#Salva mappe dentro output (formato modificabile)
ggsave("output/PG_rrmap_abs.tiff", PG_map_abs, width=5, height=5, bg="white")
ggsave("output/PG_rrmap_rel.tiff", PG_map_rel, width = 5, height = 5, bg="white")
ggsave("output/PG_postmap.tiff", PGpost_map, width=5, height=5, bg="white")

### MODELLO BESAG-YORK-MOLLIE (BYM) ###

# Step 1: Definizione del codice per il modello
# Modello log-lineare: log(E[i]) è l'offset (casi attesi), alpha0 è il 
# rischio medio globale, b.ete cattura la variabilità casuale non spaziale 
# (eterogeneità), b.clu cattura la componente di rischio spazialmente 
# strutturata (clustering tra comuni vicini)
Bymcode <- nimbleCode(
  {
    for (i in 1:N) {
      O[i] ~ dpois(mu[i]) #likelihood (poisson)
      log(mu[i]) <- log(E[i]) + alpha0 + b.ete[i] + b.clu[i] #mu=offset+intercetta+eterogeneità+clustering
      
      rr[i] <- exp(alpha0+b.ete[i]+b.clu[i])
      
      postrr[i] <- step(rr[i]-1)
      
      # Effetto casuale non strutturato spazialmente (iid tra i comuni): 
      # rappresenta variabilità/rumore locale non spiegata dalla struttura spaziale
      b.ete[i] ~ dnorm(0,tau.ete)
    } # <-- chiusura del for, FUORI da qualsiasi commento
    
    # Effetto casuale spaziale CAR (Conditional Autoregressive): il rischio 
    # di ciascun comune è modellato come dipendente dai rischi dei comuni 
    # adiacenti (definiti dalla matrice di vicinanza 'adiacenze'), 
    # assumendo che comuni vicini tendano ad avere rischi simili
    b.clu[1:N] ~ dcar_normal(adj[1:L], weights[1:L], 
                             num[1:N], tau = tau.clu)
    
    # Prior debolmente informative sulle precisioni (inverso della varianza) 
    # degli effetti casuali; questi iperparametri (0.5, 0.0005) sono una 
    # scelta convenzionale in letteratura di disease mapping
    alpha0  ~ dnorm(0, 0.001)  
    tau.ete ~ dgamma(0.5,0.0005)
    sigma.ete <- sqrt(1/tau.ete)
    tau.clu ~ dgamma(0.5,0.0005)
    sigma.clu <- sqrt(1/tau.clu)
  } # <-- chiusura del blocco principale del modello
)

# Step 2: Valori iniziali per i nodi stocastici 
inits_BYM <- list(tau.ete = 1, tau.clu = 1, alpha0 = 0, b.ete = rep(0, N), b.clu = rep(0, N))
# i valori iniziali sono tutti costumizzabili, ma e importante che le precisioni siano >0

# Step 3: Iterazioni MCMC, possiamo tenere i valori assegnati in precedenza

# Step 4: Fit del modello
modelBYM <- nimbleMCMC(
  code = Bymcode, data = data_mod, constants = constants,
  inits = inits_BYM, monitors = c("rr", "postrr"),
  niter = ni, nburnin = nb, thin = nt, nchains = nc, WAIC=TRUE,
  setSeed = 123, progressBar = FALSE, samplesAsCodaMCMC = TRUE, summary = TRUE
)

# WAIC per confronto tra modelli
modelBYM$WAIC$WAIC

# Dataframe con risultati estratti 
results_BYM <- data.frame(
  COD_COM = smr$COD_COM,
  rr = modelBYM$summary[paste0("rr[", 1:N, "]"), "Mean"],
  #intervallo di credibilità al 95%
  cri_low  = modelBYM$summary[paste0("rr[", 1:N, "]"), "95%CI_low"],
  cri_high = modelBYM$summary[paste0("rr[", 1:N, "]"), "95%CI_upp"],
  post = modelBYM$summary[paste0("postrr[", 1:N, "]"), "Mean"]
)

# Strumenti di diagnostica del modello
# Traceplot
ggs(modelBYM$samples) %>% filter(Parameter %in% paste0("rr[", sample_ids, "]"))  %>% 
  ggs_traceplot() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 4) +
  theme_bw()-> BYM_traceplot
# Plot autocorrelazione
ggs(modelBYM$samples) %>% filter(Parameter %in% paste0("rr[", sample_ids, "]"))  %>% 
  ggs_autocorrelation()+
  facet_wrap(~ Parameter, scales = "free_y", ncol = 4) +
  theme_bw() -> BYM_autocorr

print(BYM_traceplot)
print(BYM_autocorr)

#Salva plot nella cartella diagnostics dentro la cartella output
ggsave("output/diagnostics/BYM_traceplot.tiff", BYM_traceplot, width=8, height=6)
ggsave("output/diagnostics/BYM_autocorr.tiff", BYM_autocorr, width=8, height=6)

# Step 5: Mappare i risultati
com_dat <- comuni %>%
  left_join(results_BYM, by = "COD_COM")

# non è necessario ridefinire i break assoluti
# break relativi
rel.brks <- classIntervals(results_BYM$rr, n = 5, style = "quantile", intervalClosure = "right")

com_dat <- com_dat |>
  dplyr::mutate(abs.brks = cut(rr, abs.brks, include.lowest = TRUE, labels = abs.brks.labels),
                rel.brks = cut(rr, rel.brks$brks, include.lowest = TRUE, labels = make_labels(rel.brks$brks)),
                abs.brks.post = cut(post, abs.brks.post, include.lowest = TRUE, labels = abs.brks.post.labels))

# Mappe dei rischi relativi (o tassi di rischio, a seconda della parametrizzazione)
#break assoluti
BYM_map_abs <- plot_disease_map(com_dat, "abs.brks", "RR from BYM Model", palette)
print(BYM_map_abs)
#break relativi
BYM_map_rel <- plot_disease_map(com_dat, "rel.brks", "RR from BYM Model", palette)
print(BYM_map_rel)

# Mappa delle probabilita a posteriori
BYMpost_map <- plot_disease_map(com_dat, "abs.brks.post", "P(RR>1) from BYM Model", palette)
print(BYMpost_map)

#Salva mappe nella cartella output (formato modificabile)
ggsave("output/BYM_rrmap_abs.tiff", BYM_map_abs, width=5, height=5, bg="white")
ggsave("output/BYM_rrmap_rel.tiff", BYM_map_rel, width=5, height=5, bg="white")
ggsave("output/BYM_postmap.tiff", BYMpost_map, width=5, height=5, bg="white")

### MODELLO BYM2 ###
# Riparametrizzazione BYM2: un unico effetto casuale 
# combinato che pesa la componente spaziale (b.clu) e quella non strutturata 
# (b.ete) tramite il parametro phi, con tau che ne controlla la scala globale. 
# Questo evita problemi di identificabilità presenti nel BYM classico, dove 
# b.ete e b.clu non sono separatamente identificabili 
# Step 1: Definizione dl codice per il modello

Bymcode2 <- nimbleCode(
  {
    for (i in 1:N) {
      O[i] ~ dpois(mu[i]) #likelihood
      
      log(mu[i]) <- log(E[i]) + alpha0 + beta[i] 
      
      beta[i] <- (1/sqrt(tau))*((b.clu[i]*sqrt(phi/scale)) + (b.ete[i]*sqrt(1-phi))) #combinazione pesata
      
      rr[i] <- exp(alpha0+beta[i])
      
      postrr[i] <-step(rr[i]-1)
      
      b.ete[i]~dnorm(0, tau.ete)
    }
    b.clu[1:N] ~ dcar_normal(adj[1:L], weights[1:L], 
                             num[1:N], tau = 1)
    
    #prior non informative
    alpha0  ~ dnorm(0, 0.001) 
    tau.ete ~ dgamma(1,0.01)
    phi ~ dbeta(1, 1) # phi è vincolato tra 0 e 1: rappresenta la proporzione di varianza 
    #attribuibile alla componente spaziale rispetto al totale; la prior Beta(1,1) 
    #è debolmente informativa e simmetrica
    tau ~ dgamma(1,0.01) # tau è la precisione complessiva della combinazione spaziale+non strutturata;
    sigma <- 1/tau
  }
)

# Step 2: Valori iniziali per i nodi stocastici 
inits_BYM2 <- list(tau=1, alpha0 = 0, tau.ete=1, b.ete = rep(0, N), b.clu = rep(0, N), phi=0.5)

# Step 3: Iterazioni MCMC, possiamo tenere i valori assegnati in precedenza

# Step 4: Fit del modello
modelBYM2 <- nimbleMCMC(
  code = Bymcode2, data = data_mod, constants = constants,
  inits = inits_BYM2, monitors = c("rr", "postrr", "phi"),
  niter = ni, nburnin = nb, thin = nt, nchains = nc, WAIC=TRUE,
  setSeed = 123, progressBar = FALSE, samplesAsCodaMCMC = TRUE, summary = TRUE
)

# WAIC per confronto tra modelli
modelBYM2$WAIC$WAIC

# Dataframe con risultati estratti 
results_BYM2 <- data.frame(
  COD_COM = smr$COD_COM,
  rr = modelBYM2$summary[paste0("rr[", 1:N, "]"), "Mean"],
  #intervallo di credibilità al 95%
  ci_low  = modelBYM2$summary[paste0("rr[", 1:N, "]"), "95%CI_low"],
  ci_high = modelBYM2$summary[paste0("rr[", 1:N, "]"), "95%CI_upp"],
  post = modelBYM2$summary[paste0("postrr[", 1:N, "]"), "Mean"]
)
# il parametro phi ci dice quanto l’effetto spaziale strutturato (clustering) 
# contribuisce al RR rispetto all’overdispersion non strutturata (eterogeneità)
# phi vicino a 1 --> clustering contribuisce maggiormente 
# phi vicino a 0 --> eterogeneità contribuisce maggiormente 
modelBYM2$summary["phi", "Mean"]

#Salvare i risultati di SMR e RR dai 3 modelli in un unico file excel
risultati_PG   <- nomi_comuni %>% left_join(results_PG,   by = "COD_COM")
risultati_BYM  <- nomi_comuni %>% left_join(results_BYM,  by = "COD_COM")
risultati_BYM2 <- nomi_comuni %>% left_join(results_BYM2, by = "COD_COM")

writexl::write_xlsx(
  list(
    SMR             = nomi_comuni %>% left_join(smr, by = "COD_COM"),
    Poisson_Gamma   = risultati_PG,
    BYM             = risultati_BYM,
    BYM2            = risultati_BYM2
  ),
  path = "output/risultati_modelli.xlsx"
)

# Strumenti di diagnostica del modello
# Traceplot
ggs(modelBYM2$samples) %>% filter(Parameter %in% paste0("rr[", sample_ids, "]"))  %>% 
  ggs_traceplot() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 4) +
  theme_bw() -> BYM2_traceplot
# Plot autocorrelazione
ggs(modelBYM2$samples) %>% filter(Parameter %in% paste0("rr[", sample_ids, "]"))  %>% 
  ggs_autocorrelation()+
  facet_wrap(~ Parameter, scales = "free_y", ncol = 4) +
  theme_bw() -> BYM2_autocorr

print(BYM2_traceplot)
print(BYM2_autocorr)

#Salva plot nella cartella diagnostics dentro la cartella output
ggsave("output/diagnostics/BYM2_traceplot.tiff", BYM2_traceplot, width=8, height=6)
ggsave("output/diagnostics/BYM2_autocorr.tiff", BYM2_autocorr, width=8, height=6)

# Step 5: Mappare i risultati
com_dat <- comuni %>%
  left_join(results_BYM2, by = "COD_COM")

# non è necessario ridefinire i break assoluti
# break relativi
rel.brks <- classIntervals(results_BYM2$rr, n = 5, style = "quantile", intervalClosure = "right")

com_dat <- com_dat |>
  dplyr::mutate(abs.brks = cut(rr, abs.brks, include.lowest = TRUE, labels = abs.brks.labels),
                rel.brks = cut(rr, rel.brks$brks, include.lowest = TRUE, labels = make_labels(rel.brks$brks)),
                abs.brks.post = cut(post, abs.brks.post, include.lowest = TRUE, labels = abs.brks.post.labels))

# Mappa dei rischi relativi (o tassi di rischio, a seconda della parametrizzazione)
#break assoluti
BYM2_map_abs <- plot_disease_map(com_dat, "abs.brks", "RR from BYM2 Model", palette)
print(BYM2_map_abs)
#break relativi
BYM2_map_rel <- plot_disease_map(com_dat, "rel.brks", "RR from BYM2 Model", palette)
print(BYM2_map_rel)

# Mappa delle probabilita a posteriori
BYM2post_map <- plot_disease_map(com_dat, "abs.brks.post", "P(RR>1) from BYM2 Model", palette)
print(BYM2post_map)

#Salva mappe nella cartella output (formato modificabile)
ggsave("output/BYM2_rrmap_abs.tiff", BYM2_map_abs, width=5, height=5, bg="white")
ggsave("output/BYM2_rrmap_rel.tiff", BYM2_map_rel, width=5, height=5, bg="white")
ggsave("output/BYM2_postmap.tiff", BYM2post_map, width=5, height=5, bg="white")


# Confronto visivo diretto tra le 4 mappe: ci si aspetta che i modelli 
# gerarchici (PG, BYM, BYM2) producano stime più "lisciate" (shrinkage) 
# rispetto al dato grezzo SMR, specialmente nei comuni con popolazione piccola
(SMR_map_abs|PG_map_abs)/(BYM_map_abs|BYM2_map_abs) -> mappe
print(mappe)
ggsave("output/mappe_confronto.tiff", mappe, width=8, height=8, bg="white")


