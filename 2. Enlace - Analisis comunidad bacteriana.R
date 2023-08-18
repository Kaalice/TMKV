############################################
## Enlace - Analisis comunidad bacteriana ##
############################################

##############
## PAQUETES ##
##############

library("ggplot2")
library("vegan")
library("devtools")
library("pairwiseAdonis")
library(phyloseq)
library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(microbiome)
library(genefilter)
library(SpiecEasi)
library(seqtime)
library(igraph)
library(qgraph)
library(ggnet)
library(RColorBrewer)
library(tidyverse)
library(grid)
library(gridExtra)
library(network)
library(sna)
library(ggplot2)
 
###########################
##### CATEGORIZACION ASV ##
###########################
########################################
## PHYLOSEQs Core, Pan, Peripheral #####
########################################

# Directorio de trabajo
setwd("~/1. TMKV/5. Phyloseq")

# Llamar objeto ps ej. ps.Liq.COY19.F2
ps.Liq.COY19<- readRDS("ps.Liq.COY19.F2.rds")

# Crear funcion para calculo de prevalencia
prev = apply(X = otu_table(ps.Liq.COY19), MARGIN=ifelse(taxa_are_rows(ps.Liq.COY19), yes=1, no=2),
             FUN=function(x){sum(x > 0)})

# Transformar la tax_table en un dataframe
tax_table_df=as.data.frame(ps.Liq.COY19@tax_table@.Data)

# Ordenar en un dataframe
df.prev = data.frame(Prevalencia=prev, AbundanciaTotal=taxa_sums(ps.Liq.COY19),
                     tax_table_df)

prevdf1 = subset(df.prev, Phylum %in% get_taxa_unique(ps.Liq.COY19, "Phylum"))

# Definir valores de corte
prevalence.core = 8
prevalence.pan= 4

# Executar filtro de prevalencia
keepTaxa.core = rownames(prevdf1)[(prevdf1$Prevalencia >= prevalence.core)]
keepTaxa.pan1 = rownames(prevdf1)[(prevdf1$Prevalencia >= prevalence.pan)]
keepTaxa.pan2 = rownames(prevdf1)[(prevdf1$Prevalencia < prevalence.core)]
keepTaxa.peripheral = rownames(prevdf1)[(prevdf1$Prevalencia < prevalence.pan)]

# Filtrar
ps.Liq.COY19.core = prune_taxa(keepTaxa.core, ps.Liq.COY19)
ps.Liq.COY19.pan = prune_taxa(keepTaxa.pan1, ps.Liq.COY19)
ps.Liq.COY19.pan2 = prune_taxa(keepTaxa.pan2, ps.Liq.COY19.pan)
ps.Liq.COY19.pan =ps.Liq.COY19.pan2
ps.Liq.COY19.peripheral = prune_taxa(keepTaxa.peripheral, ps.Liq.COY19)

# Guardar los archivos como rds
saveRDS(ps.Liq.COY19.core, "ps.Liq.COY19.core.rds")
saveRDS(ps.Liq.COY19.pan, "ps.Liq.COY19.pan.rds")
saveRDS(ps.Liq.COY19.peripheral, "ps.Liq.COY19.peripheral.rds")

##########################
## GRAFICOS PREVALENCIA ##
##########################

# Referencia: Sierra et al. 2020
# https://github.com/mariaasierra/Lichen_Microbiome/blob/master/Core-microbiome/Prevalence.R

# Directorio de trabajo 
setwd("~/1. TMKV/5. Phyloseq")

# Crear la función "lichen_melt"

lichen_melt = function(physeq,
                       includeSampleVars = character(),
                       omitZero = FALSE){
  
  ## Se indican los argumentos necesarios para que la función.
  ## Y despues del simbolo { comienza el cuerpo de la función, es decir, las operaciones que se realizaran. 
  
  
  # Se indica a otu_table como `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  
  ## funciones utilizadas: 
  ### otu_table: accede a la tabla de OTU de un archivo phyloseq  
  ### as: Coaccionar (forzar a ser) un objeto a una clase determinada. En este caso, matrix.
  
  ## transponer la tabla
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  
  
  ## como data.frame() pero mejorado
  otudt = data.table(otutab, keep.rownames = TRUE)
  ### el argumento keep.rownames en TRUE retendrá los nombres de fila de ese objeto en una columna denominada rn.
  
  ## cambio de nombres en data frame. "rn" nombre viejo, "TaxaID" nombre nuevo. 
  setnames(otudt, "rn", "TaxaID")
  
  ## función ":=" : Para agregar columna TaxaIDchar igual a taxaID forzada a ser  caracter.
  otudt[, TaxaIDchar := as.character(TaxaID)]
  
  ## se elimina la primera columna TaxaID. Entonces, queda la misma tabla pero con las secuencias al final. 
  otudt[, TaxaID := NULL]
  
  ## cambio de nombres en data frame, de "TaxaIDchar" por "TaxaID"
  setnames(otudt, "TaxaIDchar", "TaxaID")
  
  
  # fusión de tablas de conteo
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  
  ## funciones: 
  ### melt.data.table permite la fusión en múltiples columnas simultáneamente.
  ### Resulta una tabla con 3 columnas: TaxaID, SampleID y count. 
  ### Nueva tabla: La secuencia se encuentra en "x" muestra "y" veces. Ahora el "nombre de la secuencia" aparecerá repetido tanta cantidad de veces como en tantas muestras se encuentre. 
  
  
  ### Se calcula la abundancia relativa de cada secuencia en cada muestra. 
  mdt[, RelativeAbundance := count / sum(count), by = SampleID] 
  
  ### "Se agregará la nueva columna "RelativeAbundace" que será la división entre el conteo de esa secuencia y la suma del conteo total, por cada grupo de muestras
  
  # Si existe una tax_table, acceder a ella.
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    
    # hacer que la tabla de taxonomia sea matrix, y luego data frame. 
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    
    # Cambiar el nombre de la tabla de taxas de la variable rn por "TaxaID"
    setnames(taxdt, "rn", "TaxaID")
    
    # forzar que TaxaID sea de clase caracter, y se cambia a ultima columna de la tabla (antes primera)
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    
    # setkey ordena la columna en orden ascendente, entonces ordeno ambas tablas igual
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    
    mdt <- taxdt[mdt]
  }
  setkey(mdt, "TaxaID")
  
  #Entrega tabla mdt
  return(mdt)
}
# Termina de crearse la función.

# Ejemplo con Sue TAM19  ##

# Leer objeto phyloseq
ps<-readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F1.rds")

# Fusionar objetos phyloseq
physeq<- merge_phyloseq(ps)

# Aplicar la función creada al archivo phyloseq. Resulta una tabla 
mdt <- lichen_melt(physeq)

# Revisar el nombre para cada rango
rank_names(ps)

## Resulta una tabla que muestra: Cada secuencia, su taxonomia, a que muestra pertenece, la cantidad de lecturas y la abundancia relativa. 

# Construir tabla con la prevalencia y total de lecturas de cada secuencia.
lichen.prev = mdt[, list(Prevalence = ((mean(count > 0)*9)), #Mean of samples with counts greater than cero
                         TotalCounts = sum(count)),
                  by = TaxaID]
# Construir tabla con la secuencia y su filo correspondiente. 
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)])) #Color taxa by phylum

# Ordenar las tablas, y se unen. Queda la tabla lichen.prev con el filo para cada secuencia.  
setkey(lichen.prev, TaxaID)
setkey(addPhylum, TaxaID)
lichen.prev <- addPhylum[lichen.prev]

# Cambiar numero de filos segun set de datos
get_taxa_unique(ps, "Phylum")

# Sumar el total de conteos por filo, ajustar en la parte []$filo, el numero de filos deseado
showPhyla = lichen.prev[, sum(TotalCounts), by = Phylum][order(-V1)][1:26]$Filo
# Ordenar
setkey(lichen.prev, Phylum)

#view(lichen.prev)
#Eliminar NA
# lichen.prev <- lichen.prev[!is.na(lichen.prev$Phylum),]

# Crar dataframe colores
Filo <- c("Proteobacteria","Actinobacteriota","Acidobacteriota","Bacteroidota",
          "Planctomycetota","Verrucomicrobiota","Myxococcota","Armaimonadota",
          "Cyanobacteria","Chloroflexota","Gemmatimonadota","Patescibacteria",
          "Eremiobacterota","Bdellovibrionota","Nitrospirota","Desulfobacterota_B",
          "Firmicutes","Dormibacterota","Methylomirabilota","Elusimicrobiota",
          "Omnitrophota","Firmicutes_A","Firmicutes_C","Dependentiae",
          "Desulfobacterota","Eisenbacteria")
cod <- c("#69AFB2","#E4B664","#6184BE","#E0E168","#C95F6F","#64B0D1","#1E5664","#57599D","#A1BE57",
         "#7A609F","#BB68A1","#67AA8C","#7AAD59","#CD6926","#C96299","#C23234","#4D7AAF","#9B66A0",
         "#918854","#355B8A","#040006","#7F7F7F","#DAD6C1","#738E3C","#8F3734","#E9E337")

colours <- data.frame(Filo, cod)

### Grafico "termo"
p <- ggplot(lichen.prev, 
            mapping = aes(Prevalence, TotalCounts, color=Phylum)) + 
  geom_jitter(size = 5, alpha = 0.3) + 
  geom_vline(xintercept=1.35, color='gray87',linetype="dashed") +
  geom_vline(xintercept=2.25, color='gray87',linetype="dashed") +
  geom_vline(xintercept=3.15, color='red',linetype="dashed") +
  geom_vline(xintercept=4.05, color='gray87',linetype="dashed") +
  geom_vline(xintercept=4.95, color='gray87',linetype="dashed") +
  geom_vline(xintercept=5.85, color='gray87',linetype="dashed") +
  geom_vline(xintercept=6.75, color='red',linetype="dashed") +
  geom_vline(xintercept=7.65, color='gray87',linetype="dashed") +
  geom_vline(xintercept=8.55, color='gray87',linetype="dashed") +
  geom_hline(yintercept=9.55, color='blue',linetype="dashed") +
  annotate(geom='label', x=1.5, y=9999, label="Pan",angle=45, fontface=2) +
  annotate(geom='label', x=5, y=9999, label="Peripheral", fontface=2) +
  annotate(geom='label', x=8.5, y=9999, label="Core", fontface=2) +
  scale_y_log10() +  theme_minimal() + ylab("Total Counts") +
  theme(strip.text = element_text(size=20),
        axis.text.x = element_text(colour = "black", size=11 ), 
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title.x = element_text(size=13), 
        axis.title.y = element_text(size=13),
        legend.text  = element_text(size=12, colour = "gray40"),
        legend.title = element_text(size = 14)) + ggtitle("Prevalencia Sue TAMahique")+
  scale_x_continuous(breaks = seq(1,9,1)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_color_manual(values = colours$cod, limits=colours$Filo)

p

setwd("~/1. TMKV/2. Results/7. Delimitación por prevalencia (Sierra)")
ggsave("Prevalencia Sue TAM19.png", width = 20, height = 10, dpi = 300)


################################
## Medidas de diversidad alfa ##
################################

# En phyloseq simplemente llamamos la función plot_richness y podemos visualizar las medidas de diversidad.
setwd("~/1. TMKV/5. Phyloseq")

ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")

ps.Liq.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.core.rds")
ps.Sus.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.core.rds")
ps.Liq.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.core.rds")
ps.Sus.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.core.rds")
ps.Sue.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.core.rds")
ps.Sue.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.core.rds")

ps.Liq.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.pan.rds")
ps.Sus.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.pan.rds")
ps.Liq.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.pan.rds")
ps.Sus.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.pan.rds")
ps.Sue.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.pan.rds")
ps.Sue.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.pan.rds")

ps.Liq.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.peripheral.rds")
ps.Sus.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.peripheral.rds")
ps.Liq.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.peripheral.rds")
ps.Sus.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.peripheral.rds")
ps.Sue.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.peripheral.rds")
ps.Sue.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.peripheral.rds")
# Armar Phyloseqs


physeq.F2<- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2, ps.Sue.COY19.F2)

physeq.core <- merge_phyloseq(ps.Liq.TAM19.core,ps.Liq.COY19.core,ps.Sus.TAM19.core,ps.Sus.COY19.core,ps.Sue.TAM19.core, ps.Sue.COY19.core)
physeq.pan<- merge_phyloseq(ps.Liq.TAM19.pan,ps.Liq.COY19.pan,ps.Sus.TAM19.pan,ps.Sus.COY19.pan,ps.Sue.TAM19.pan, ps.Sue.COY19.pan)
physeq.peripheral <- merge_phyloseq(ps.Liq.TAM19.peripheral,ps.Liq.COY19.peripheral,ps.Sus.TAM19.peripheral,ps.Sus.COY19.peripheral,ps.Sue.TAM19.peripheral, ps.Sue.COY19.peripheral)

plot.rich.alfa.F2<-plot_richness(physeq.F2, color = "Tipo_de_muestra_sitio", x = "Tipo_de_muestra_sitio", 
                                 measures = c("Shannon")) + 
  geom_boxplot(aes(fill = Tipo_de_muestra_sitio), alpha=.7) +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))

setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
ggsave("Diversidad alfa F2 F2.png", width = 5, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")

plot.rich.alfa.core<-plot_richness(physeq.core, color = "Tipo_de_muestra_sitio", x = "Tipo_de_muestra_sitio", 
                                   measures = c("Shannon")) + 
  geom_boxplot(aes(fill = Tipo_de_muestra_sitio), alpha=.7) +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))

setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
ggsave("Diversidad alfa core.png", width = 5, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")

plot.rich.alfa.pan<-plot_richness(physeq.pan, color = "Tipo_de_muestra_sitio", x = "Tipo_de_muestra_sitio", 
                                  measures = c("Shannon")) + 
  geom_boxplot(aes(fill = Tipo_de_muestra_sitio), alpha=.7) +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))

setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
ggsave("Diversidad alfa pan.png", width = 5, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")

plot.rich.alfa.peripheral<-plot_richness(physeq.peripheral, color = "Tipo_de_muestra_sitio", x = "Tipo_de_muestra_sitio", 
                                         measures = c("Shannon")) + 
  geom_boxplot(aes(fill = Tipo_de_muestra_sitio), alpha=.7) +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))

setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
ggsave("Diversidad alfa peripheral.png", width = 5, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")

## Medida de diversidad alfa con diferencia  


install.packages('devtools')
devtools::install_github('twbattaglia/btools')
library(btools)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(ggpubr)

# En phyloseq simplemente llamamos la función plot_richness y podemos visualizar las medidas de diversidad.
setwd("~/1. TMKV/5. Phyloseq")

physeq.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.rds")

# Guardamos un dataframe con las medidas de diversidad alfa
alpha_pd <- estimate_pd(physeq.F2)
# Combinamos la metadata con alpha.diversity
data <- cbind(sample_data(physeq.F2), alpha_pd) 

# Y calculamos un ANOVA
psd5.anova <- aov(PD ~ Tipo_de_muestra_sitio, data) 
# install.packages("xtable")
library(xtable)
psd5.anova.table <- xtable(psd5.anova)
View(psd5.anova.table)

setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
write.table(psd5.anova.table, file = "Anova-Diversidadalfa.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# La función global nos da 26 medidas de diversidad que nos ayudan a entender la estructura de las comunidades microbianas. En general, estas medidas se dividen en riqueza, diversidad, dominancia, rareza, cobertura y uniformidad.

tab <- global(physeq.F2, index = "all")
head(tab)
setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
write.table(tab, file = "Indices de diversidad, rareza y dominancia.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# Generamos un objeto `phyloseq` sin taxa que sume 0 reads
physeq.F2.2 <- prune_taxa(taxa_sums(physeq.F2) > 0, psd5.anova)
# Calculamos los índices de diversidad
tab <- diversities(physeq.F2.2, index = "all")
# Y finalmente visualizamos la tabla de resultados
head(tab)

physeq.F2.meta <- meta(physeq.F2)
head(physeq.F2.meta)
physeq.F2.meta$Shannon <- tab$shannon 
physeq.F2.meta$InverseSimpson <- tab$inverse_simpson
# Obtenemos las variables desde nuestro objeto `phyloseq`
spps <- levels(physeq.F2.meta$Tipo_de_muestra_sitio)
# Creamos una lista de lo que queremos comparar
pares.spps <- combn(seq_along(spps), 2, simplify = FALSE, FUN = function(i)spps[i])
# Imprimimos en pantalla el resultado
print(pares.spps)

plot.rich.alfa.F2<-plot_richness(physeq.F2, color = "Tipo_de_muestra_sitio", x = "Tipo_de_muestra_sitio", 
                                 measures = c("Shannon")) + 
  geom_boxplot(aes(fill = Tipo_de_muestra_sitio), alpha=.7) +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) + stat_compare_means(comparisons = pares.spps, label = "p.signif",hide.ns = FALSE,tip.length = 0)


setwd("~/1. TMKV/2. Results/8. Diversidad alfa")
ggsave("diversidad alfa.png", width = 10, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")

##################################
## Diversidad beta PCoA Unifrac ##
##################################

ps.TMKV<- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.rds")

# Para set completo 
ps.TMKV.mds.unifrac <- ordinate(ps.TMKV, method = "PCoA", distance = "unifrac")
evals <- ps.TMKV.mds.unifrac$values$Eigenvalues
UNIFRAC.TMKV <- plot_ordination(ps.TMKV, ps.TMKV.mds.unifrac, color = "Tipo_de_muestra_sitio") +
  labs(col = "Tipo de muestra") +
  coord_fixed(sqrt(evals[2] / evals[1]))

# https://github.com/joey711/phyloseq/issues/323
UNIFRAC.TMKV <- UNIFRAC.TMKV + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))


setwd("~/1. TMKV/2. Results/9. Diversidad beta")
ggsave("PCoA_Unifrac_TMKV.png", width = 20, height = 10, dpi = 300)

# Permanova 
ps.TMKV<- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.rds") #F2
meta(ps.TMKV)

unifrac.dist <- UniFrac(ps.TMKV, weighted = TRUE, normalized = TRUE,  
                        parallel = FALSE, fast = TRUE)
permanova <- adonis(unifrac.dist ~ Tipo_de_muestra_sitio, data = meta(ps.TMKV))

pwadonis = pairwise.adonis(unifrac.dist, meta(ps.TMKV)$Tipo_de_muestra_sitio)
pwadonis.dat=as.data.frame(pwadonis)
View(pwadonis.dat)
write.csv2(pwadonis.dat,"~/1. TMKV/2. Results/9. Diversidad beta/pwadonis.csv", row.names = FALSE)

#########################
# Diversidad beta NMDS ##
#########################


## F2 ##

setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")

physeq.F2 <- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2,ps.Sue.COY19.F2)

physeq.mds.bray <- ordinate(physeq.F2, method = "MDS", distance = "bray")
evals <- physeq.mds.bray$values$Eigenvalues
pord.F2 <- plot_ordination(physeq.F2, physeq.mds.bray, color = "Tipo_de_muestra_sitio") +
  labs(col = "Tipo_de_muestra_sitio")  +
  coord_fixed(sqrt(evals[2] / evals[1]))+ 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))
pord.F2

setwd("~/1. TMKV/2. Results/9. Diversidad beta")
ggsave("NMDS F2.png", width = 20, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")



## Core ##


setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Liq.TAM19.core.rds")
ps.Sus.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sus.TAM19.core.rds")
ps.Liq.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Liq.COY19.core.rds")
ps.Sus.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sus.COY19.core.rds")
ps.Sue.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sue.TAM19.core.rds")
ps.Sue.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sue.COY19.core.rds")

## Reconstrucción de archivos phyloseq (sacar arbol) ##

## Liquen TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.TAM19.core)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.TAM19.core)

#Tercero, tabla de metadatos
meta.data= meta(ps.Liq.TAM19.core)
# La reconstruccion del objeto phyloseq:
ps.Liq.TAM19.core <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                              tax_table(Tabla.Taxonomia),
                              sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.TAM19.core, "ps.Liq.TAM19.core.rds")

## Sustrato TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sus.TAM19.core)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sus.TAM19.core)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sus.TAM19.core)
# La reconstruccion del objeto phyloseq:
ps.Sus.TAM19.core <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                              tax_table(Tabla.Taxonomia),
                              sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sus.TAM19.core, "ps.Sus.TAM19.core.rds")

## Suelo TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sue.TAM19.core)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sue.TAM19.core)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sue.TAM19.core)
# La reconstruccion del objeto phyloseq:
ps.Sue.TAM19.core <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                              tax_table(Tabla.Taxonomia),
                              sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sue.TAM19.core, "ps.Sue.TAM19.core.rds")

## Liquen COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.COY19.core)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.COY19.core)

#Tercero, tabla de metadatos
meta.data= meta(ps.Liq.COY19.core)
# La reconstruccion del objeto phyloseq:
ps.Liq.COY19.core <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                              tax_table(Tabla.Taxonomia),
                              sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.COY19.core, "ps.Liq.COY19.core.rds")

## Sustrato COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sus.COY19.core)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sus.COY19.core)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sus.COY19.core)
# La reconstruccion del objeto phyloseq:
ps.Sus.COY19.core <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                              tax_table(Tabla.Taxonomia),
                              sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sus.COY19.core, "ps.Sus.COY19.core.rds")

## Suelo COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sue.COY19.core)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sue.COY19.core)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sue.COY19.core)
# La reconstruccion del objeto phyloseq:
ps.Sue.COY19.core <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                              tax_table(Tabla.Taxonomia),
                              sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sue.COY19.core, "ps.Sue.COY19.core.rds")


## Diversidad beta core ##

physeq.core <- merge_phyloseq(ps.Liq.TAM19.core,ps.Liq.COY19.core,ps.Sus.TAM19.core,ps.Sus.COY19.core,ps.Sue.TAM19.core,ps.Sue.COY19.core)

physeq.mds.bray <- ordinate(physeq.core, method = "MDS", distance = "bray")
evals <- physeq.mds.bray$values$Eigenvalues
pord.core <- plot_ordination(physeq.core, physeq.mds.bray, color = "Tipo_de_muestra_sitio") +
  labs(col = "Tipo_de_muestra_sitio") +
  coord_fixed(sqrt(evals[2] / evals[1]))+ 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))


pord.core


setwd("~/1. TMKV/2. Results/9. Diversidad beta")
ggsave("NMDS core.png", width = 20, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")


## Pan ##


setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Liq.TAM19.pan.rds")
ps.Sus.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sus.TAM19.pan.rds")
ps.Liq.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Liq.COY19.pan.rds")
ps.Sus.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sus.COY19.pan.rds")
ps.Sue.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sue.TAM19.pan.rds")
ps.Sue.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sue.COY19.pan.rds")

## Reconstrucción de archivos phyloseq (sacar arbol) ##

## Liquen TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.TAM19.pan)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.TAM19.pan)

#Tercero, tabla de metadatos
meta.data= meta(ps.Liq.TAM19.pan)
# La reconstruccion del objeto phyloseq:
ps.Liq.TAM19.pan <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                             tax_table(Tabla.Taxonomia),
                             sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.TAM19.pan, "ps.Liq.TAM19.pan.rds")

## Sustrato TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sus.TAM19.pan)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sus.TAM19.pan)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sus.TAM19.pan)
# La reconstruccion del objeto phyloseq:
ps.Sus.TAM19.pan <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                             tax_table(Tabla.Taxonomia),
                             sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sus.TAM19.pan, "ps.Sus.TAM19.pan.rds")

## Suelo TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sue.TAM19.pan)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sue.TAM19.pan)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sue.TAM19.pan)
# La reconstruccion del objeto phyloseq:
ps.Sue.TAM19.pan <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                             tax_table(Tabla.Taxonomia),
                             sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sue.TAM19.pan, "ps.Sue.TAM19.pan.rds")

## Liquen COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.COY19.pan)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.COY19.pan)

#Tercero, tabla de metadatos
meta.data= meta(ps.Liq.COY19.pan)
# La reconstruccion del objeto phyloseq:
ps.Liq.COY19.pan <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                             tax_table(Tabla.Taxonomia),
                             sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.COY19.pan, "ps.Liq.COY19.pan.rds")

## Sustrato COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sus.COY19.pan)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sus.COY19.pan)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sus.COY19.pan)
# La reconstruccion del objeto phyloseq:
ps.Sus.COY19.pan <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                             tax_table(Tabla.Taxonomia),
                             sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sus.COY19.pan, "ps.Sus.COY19.pan.rds")

## Suelo COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sue.COY19.pan)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sue.COY19.pan)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sue.COY19.pan)
# La reconstruccion del objeto phyloseq:
ps.Sue.COY19.pan <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                             tax_table(Tabla.Taxonomia),
                             sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sue.COY19.pan, "ps.Sue.COY19.pan.rds")
## Diversidad beta ##

physeq.pan <- merge_phyloseq(ps.Liq.TAM19.pan,ps.Liq.COY19.pan,ps.Sus.TAM19.pan,ps.Sus.COY19.pan,ps.Sue.TAM19.pan,ps.Sue.COY19.pan)


physeq.mds.bray <- ordinate(physeq.pan, method = "MDS", distance = "bray")
evals <- physeq.mds.bray$values$Eigenvalues
pord.pan <- plot_ordination(physeq.pan, physeq.mds.bray, color = "Tipo_de_muestra_sitio") +
  labs(col = "Tipo_de_muestra_sitio") +
  coord_fixed(sqrt(evals[2] / evals[1]))+ 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))

pord.pan
setwd("~/1. TMKV/2. Results/9. Diversidad beta")
ggsave("NMDS pan.png", width = 20, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")


## peripheral ##


setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Liq.TAM19.peripheral.rds")
ps.Sus.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sus.TAM19.peripheral.rds")
ps.Liq.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Liq.COY19.peripheral.rds")
ps.Sus.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sus.COY19.peripheral.rds")
ps.Sue.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sue.TAM19.peripheral.rds")
ps.Sue.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.Sue.COY19.peripheral.rds")

## Reconstrucción de archivos phyloseq (sacar arbol) ##

# Liquen TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.TAM19.peripheral)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.TAM19.peripheral)

#Tercero, tabla de metadatos
meta.data= meta(ps.Liq.TAM19.peripheral)
# La reconstruccion del objeto phyloseq:
ps.Liq.TAM19.peripheral <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                                    tax_table(Tabla.Taxonomia),
                                    sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.TAM19.peripheral, "ps.Liq.TAM19.peripheral.rds")

## Sustrato TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sus.TAM19.peripheral)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sus.TAM19.peripheral)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sus.TAM19.peripheral)
# La reconstruccion del objeto phyloseq:
ps.Sus.TAM19.peripheral <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                                    tax_table(Tabla.Taxonomia),
                                    sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sus.TAM19.peripheral, "ps.Sus.TAM19.peripheral.rds")

## Suelo TAM19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sue.TAM19.peripheral)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sue.TAM19.peripheral)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sue.TAM19.peripheral)
# La reconstruccion del objeto phyloseq:
ps.Sue.TAM19.peripheral <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                                    tax_table(Tabla.Taxonomia),
                                    sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sue.TAM19.peripheral, "ps.Sue.TAM19.peripheral.rds")

## Liquen COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.COY19.peripheral)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.COY19.peripheral)

#Tercero, tabla de metadatos
meta.data= meta(ps.Liq.COY19.peripheral)
# La reconstruccion del objeto phyloseq:
ps.Liq.COY19.peripheral <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                                    tax_table(Tabla.Taxonomia),
                                    sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.COY19.peripheral, "ps.Liq.COY19.peripheral.rds")

## Sustrato COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sus.COY19.peripheral)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sus.COY19.peripheral)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sus.COY19.peripheral)
# La reconstruccion del objeto phyloseq:
ps.Sus.COY19.peripheral <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                                    tax_table(Tabla.Taxonomia),
                                    sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sus.COY19.peripheral, "ps.Sus.COY19.peripheral.rds")

## Suelo COY19 ##

# Primero, la tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Sue.COY19.peripheral)

# Segundo, la tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Sue.COY19.peripheral)

#Tercero, tabla de metadatos
meta.data= meta(ps.Sue.COY19.peripheral)
# La reconstruccion del objeto phyloseq:
ps.Sue.COY19.peripheral <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                                    tax_table(Tabla.Taxonomia),
                                    sample_data(meta.data))

# Guardamos los archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Sue.COY19.peripheral, "ps.Sue.COY19.peripheral.rds")
## Diversidad beta ##

physeq.peripheral <- merge_phyloseq(ps.Liq.TAM19.peripheral,ps.Liq.COY19.peripheral,ps.Sus.TAM19.peripheral,ps.Sus.COY19.peripheral,ps.Sue.TAM19.peripheral,ps.Sue.COY19.peripheral)


physeq.mds.bray <- ordinate(physeq.peripheral, method = "MDS", distance = "bray")
evals <- physeq.mds.bray$values$Eigenvalues
pord.peripheral <- plot_ordination(physeq.peripheral, physeq.mds.bray, color = "Tipo_de_muestra_sitio") +
  labs(col = "Tipo_de_muestra_sitio") +
  coord_fixed(sqrt(evals[2] / evals[1]))+ 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43")) +
  scale_fill_manual(values=c("#6ECB63","#FF9A4A","#1F965E","#FF7044","#1E7179","#FE2D43"))

pord.peripheral
setwd("~/1. TMKV/2. Results/9. Diversidad beta")
ggsave("NMDS peripheral.png", width = 20, height = 10, dpi = 300)
setwd("~/1. TMKV/5. Phyloseq")


#########################
## Abundancia relativa ##
#########################

###################################
## Tabla con abundancia relativa ##
###################################

# Objetos ps
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.COY19<- readRDS("ps.Liq.COY19.F2.rds")
ps.Sus.COY19<- readRDS("ps.Sus.COY19.F2.rds")
ps.Sue.COY19<- readRDS("ps.Sue.COY19.F2.rds")
ps.Liq.TAM19<- readRDS("ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19<- readRDS("ps.Sus.TAM19.F2.rds")
ps.Sue.TAM19<- readRDS("ps.Sue.TAM19.F2.rds")

# Merge ps
ps.completo <- merge_phyloseq(ps.Liq.COY19,ps.Sus.COY19,ps.Sue.COY19,ps.Liq.TAM19,ps.Sus.TAM19,ps.Sue.TAM19)
ps.completo
meta(ps.completo)


## Filos ##


# Agrupamos a nivel de filo
ps.completo.PsNA <- subset_taxa(ps.completo, Phylum!="NA")
ps.completo.PsNA.agg <- aggregate_taxa(ps.completo.PsNA, "Phylum")

# Creamos un objeto con las abundancias relativas de reads por cada filo y por cada muestra
Ab.completo.PsNA<-abundances(ps.completo.PsNA.agg, "compositional")
Ab.completo.PsNA <- as.data.frame(Ab.completo.PsNA)
View(Ab.completo.PsNA)

# Guardamos abundancias relativas de filo por muestra
setwd("~/1. TMKV/2. Results/11. Abundancia diferencial")
write.csv(Ab.completo.PsNA, "Ab.rel.Phylum.csv")


## Proteobacteria ##


# Todos los Generos
ps.completo.Proteobacteria <- subset_taxa(ps.completo, Phylum=="Proteobacteria")
ps.completo.Proteobacteria.sNA <- subset_taxa(ps.completo.Proteobacteria, Genus!="NA")
ps.completo.Proteobacteria.sNA.agg <- aggregate_taxa(ps.completo.Proteobacteria.sNA, "Genus")

Ab.completo.Proteobacteria.sNA <-abundances(ps.completo.Proteobacteria.sNA.agg, "compositional")
Ab.completo.Proteobacteria.sNA <- as.data.frame(Ab.completo.Proteobacteria.sNA)

write.csv(Ab.completo.Proteobacteria.sNA, "Ab.rel.Proteobacteria.csv")


## Actinobacteriota ##


ps.completo.Actinobacteriota <- subset_taxa(ps.completo, Phylum=="Actinobacteriota")
ps.completo.Actinobacteriota.sNA <- subset_taxa(ps.completo.Actinobacteriota, Genus!="NA")
ps.completo.Actinobacteriota.sNA.agg <- aggregate_taxa(ps.completo.Actinobacteriota.sNA, "Genus")

Ab.completo.Actinobacteriota.sNA <-abundances(ps.completo.Actinobacteriota.sNA.agg, "compositional")
Ab.completo.Actinobacteriota.sNA <- as.data.frame(Ab.completo.Actinobacteriota.sNA)

write.csv(Ab.completo.Actinobacteriota.sNA, "Ab.rel.Actinobacteriota.csv")


## Acidobacteriota ##


ps.completo.Acidobacteriota <- subset_taxa(ps.completo, Phylum=="Acidobacteriota")
ps.completo.Acidobacteriota.sNA <- subset_taxa(ps.completo.Acidobacteriota, Genus!="NA")
ps.completo.Acidobacteriota.sNA.agg <- aggregate_taxa(ps.completo.Acidobacteriota.sNA, "Genus")

Ab.completo.Acidobacteriota.sNA <-abundances(ps.completo.Acidobacteriota.sNA.agg, "compositional")
Ab.completo.Acidobacteriota.sNA <- as.data.frame(Ab.completo.Acidobacteriota.sNA)

write.csv(Ab.completo.Acidobacteriota.sNA, "Ab.rel.Acidobacteriota.csv")


## Bacteroidota ##


ps.completo.Bacteroidota <- subset_taxa(ps.completo, Phylum=="Bacteroidota")
ps.completo.Bacteroidota.sNA <- subset_taxa(ps.completo.Bacteroidota, Genus!="NA")
ps.completo.Bacteroidota.sNA.agg <- aggregate_taxa(ps.completo.Bacteroidota.sNA, "Genus")

Ab.completo.Bacteroidota.sNA <-abundances(ps.completo.Bacteroidota.sNA.agg, "compositional")
Ab.completo.Bacteroidota.sNA <- as.data.frame(Ab.completo.Bacteroidota.sNA)

write.csv(Ab.completo.Bacteroidota.sNA, "Ab.rel.Bacteroidota.csv")


## Planctomycetota ##

ps.completo.Planctomycetota <- subset_taxa(ps.completo, Phylum=="Planctomycetota")
ps.completo.Planctomycetota.sNA <- subset_taxa(ps.completo.Planctomycetota, Genus!="NA")
ps.completo.Planctomycetota.sNA.agg <- aggregate_taxa(ps.completo.Planctomycetota.sNA, "Genus")

Ab.completo.Planctomycetota.sNA <-abundances(ps.completo.Planctomycetota.sNA.agg, "compositional")
Ab.completo.Planctomycetota.sNA <- as.data.frame(Ab.completo.Planctomycetota.sNA)

write.csv(Ab.completo.Planctomycetota.sNA, "Ab.rel.Planctomycetota.csv")


## Armatimonadota ##

ps.completo.Armatimonadota <- subset_taxa(ps.completo, Phylum=="Armatimonadota")
ps.completo.Armatimonadota.sNA <- subset_taxa(ps.completo.Armatimonadota, Genus!="NA")
ps.completo.Armatimonadota.sNA.agg <- aggregate_taxa(ps.completo.Armatimonadota.sNA, "Genus")

Ab.completo.Armatimonadota.sNA <-abundances(ps.completo.Armatimonadota.sNA.agg, "compositional")
Ab.completo.Armatimonadota.sNA <- as.data.frame(Ab.completo.Armatimonadota.sNA)

write.csv(Ab.completo.Armatimonadota.sNA, "Ab.rel.Armatimonadota.csv")



##################################
# Grafico de barras apiladas AR ##
##################################

setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.peripheral.rds")
ps.Sus.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.peripheral.rds")
ps.Liq.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.peripheral.rds")
ps.Sus.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.peripheral.rds")
ps.Sue.TAM19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.peripheral.rds")
ps.Sue.COY19.peripheral <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.peripheral.rds")

physeq.peripheral <- merge_phyloseq(ps.Liq.TAM19.peripheral,ps.Liq.COY19.peripheral,ps.Sus.TAM19.peripheral,ps.Sus.COY19.peripheral,ps.Sue.TAM19.peripheral,ps.Sue.COY19.peripheral)
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.pan.rds")
ps.Sus.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.pan.rds")
ps.Liq.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.pan.rds")
ps.Sus.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.pan.rds")
ps.Sue.TAM19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.pan.rds")
ps.Sue.COY19.pan <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.pan.rds")

physeq.pan <- merge_phyloseq(ps.Liq.TAM19.pan,ps.Liq.COY19.pan,ps.Sus.TAM19.pan,ps.Sus.COY19.pan,ps.Sue.TAM19.pan,ps.Sue.COY19.pan)
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.core.rds")
ps.Sus.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.core.rds")
ps.Liq.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.core.rds")
ps.Sus.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.core.rds")
ps.Sue.TAM19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.core.rds")
ps.Sue.COY19.core <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.core.rds")

physeq.core <- merge_phyloseq(ps.Liq.TAM19.core,ps.Liq.COY19.core,ps.Sus.TAM19.core,ps.Sus.COY19.core,ps.Sue.TAM19.core,ps.Sue.COY19.core)
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")

physeq.F2 <- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2,ps.Sue.COY19.F2)


rank_names(physeq.peripheral)
rank_names(physeq.pan)
rank_names(physeq.core)
rank_names(physeq.F2)

## Comunidad completa 
#podemos ver cuantos filos hay en total
length(get_taxa_unique(physeq.F2, "Phylum"))
#23 filos (cony); a mi me dan 29 filos (nayla); 27 filos (KV)

#vemos los filos que hay y cuantos ASVs tiene cada uno
tabla_filos <- table(tax_table(physeq.F2)[, "Phylum"], exclude = NULL)
tabla_filos <- as.data.frame(tabla_filos)

#Ahora creamos un archvo phyloseq pero a nivel de filo
#aglomeramos nivel de filo y relativizamos abundancia
filo.ps.F2.para.grafico = tax_glom(physeq.F2, "Phylum", NArm = FALSE)#aglomeramos a nivel de filo e incluimos NA para ver su abundancia
#Ahora una vez aglomerados lo relativisamos
rel.filo.F2.para.grafico = transform_sample_counts(filo.ps.F2.para.grafico, function(x) x/sum(x)) #abundancia relativa

# Tabla de codigo de colores
# Los nombres de los grupos que quieres colorear (TAL CUAL como los arroja R (fijarse en las mayusculas))
Filo <- c("Proteobacteria","Actinobacteriota","Acidobacteriota","Bacteroidota",
          "Planctomycetota","Verrucomicrobiota","Myxococcota","Armaimonadota",
          "Cyanobacteria","Chloroflexota","Gemmatimonadota","Patescibacteria",
          "Eremiobacterota","Bdellovibrionota","Nitrospirota","Desulfobacterota_B",
          "Firmicutes","Dormibacterota","Methylomirabilota","Elusimicrobiota",
          "Omnitrophota","Firmicutes_A","Firmicutes_C","Dependentiae",
          "Desulfobacterota","Eisenbacteria")
# Los codigos de los colores asociados
cod <- c("#69AFB2","#E4B664","#6184BE","#E0E168","#C95F6F","#64B0D1","#1E5664","#57599D","#A1BE57",
         "#7A609F","#BB68A1","#67AA8C","#7AAD59","#CD6926","#C96299","#C23234","#4D7AAF","#9B66A0",
         "#918854","#355B8A","#040006","#7F7F7F","#DAD6C1","#738E3C","#8F3734","#E9E337")

# El data frame con ambos vectores asociados
color <- data.frame(Filo, cod)

meta <- meta(physeq.F2)
es <- row.names(meta)

# Graficamos
Filo.F2 = plot_bar(rel.filo.F2.para.grafico, x= "Réplica", fill="Phylum") + guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~Tipo_de_muestra_sitio, scales = "free_x", ncol = 5) + 
  scale_fill_manual(values = color$cod, limits=color$Filo) + geom_bar(stat="identity", color="black", size=0) +
  labs( ylab="Frecuencia relativa", xlab("")) + #definimos nombre de ejes
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = ",")) + #comas para decimales
  theme(text = element_text(size = 10)) + theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(strip.text = element_text(size = 8 , color="black")) +
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 6))+ 
  labs(fill = "Filo")


setwd("~/1. TMKV/2. Results/10. Abundancia relativa")
ggsave("ARFilos_F2.png", width = 12.5, height = 7.3, dpi = 300)

## Core microbioma

#podemos ver cuantos filos hay en total
length(get_taxa_unique(physeq.core, "Phylum"))
#23 filos (cony); a mi me dan 29 filos (nayla); 27 filos (KV)

#vemos los filos que hay y cuantos ASVs tiene cada uno
tabla_filos <- table(tax_table(physeq.core)[, "Phylum"], exclude = NULL)
tabla_filos <- as.data.frame(tabla_filos)

#Ahora creamos un archvo phyloseq pero a nivel de filo
#aglomeramos nivel de filo y relativizamos abundancia
filo.ps.core.para.grafico = tax_glom(physeq.core, "Phylum", NArm = FALSE)#aglomeramos a nivel de filo e incluimos NA para ver su abundancia
#Ahora una vez aglomerados lo relativisamos
rel.filo.core.para.grafico = transform_sample_counts(filo.ps.core.para.grafico, function(x) x/sum(x)) #abundancia relativa

# Tabl de codigo de colores
Filo <- c("Proteobacteria","Actinobacteriota","Acidobacteriota","Bacteroidota",
          "Planctomycetota","Verrucomicrobiota","Myxococcota","Armaimonadota",
          "Cyanobacteria","Chloroflexota","Gemmatimonadota","Patescibacteria",
          "Eremiobacterota","Bdellovibrionota","Nitrospirota","Desulfobacterota_B",
          "Firmicutes","Dormibacterota","Methylomirabilota","Elusimicrobiota",
          "Omnitrophota","Firmicutes_A","Firmicutes_C","Dependentiae",
          "Desulfobacterota","Eisenbacteria")
cod <- c("#69AFB2","#E4B664","#6184BE","#E0E168","#C95F6F","#64B0D1","#1E5664","#57599D","#A1BE57",
         "#7A609F","#BB68A1","#67AA8C","#7AAD59","#CD6926","#C96299","#C23234","#4D7AAF","#9B66A0",
         "#918854","#355B8A","#040006","#7F7F7F","#DAD6C1","#738E3C","#8F3734","#E9E337")

color <- data.frame(Filo, cod)

meta <- meta(physeq.core)
es <- row.names(meta)

# Graficamos
Filo.core = plot_bar(rel.filo.core.para.grafico, x= "Réplica", fill="Phylum") + guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~Tipo_de_muestra_sitio, scales = "free_x", ncol = 5) + 
  scale_fill_manual(values = color$cod, limits=color$Filo) + geom_bar(stat="identity", color="black", size=0) +
  labs( ylab="Frecuencia relativa", xlab("")) + #definimos nombre de ejes
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = ",")) + #comas para decimales
  theme(text = element_text(size = 10)) + theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(strip.text = element_text(size = 8 , color="black")) +
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 6))+ 
  labs(fill = "Filo")


setwd("~/1. TMKV/2. Results/10. Abundancia relativa")
ggsave("ARFilos_core.png", width = 12.5, height = 7.3, dpi = 300)

# pan

#podemos ver cuantos filos hay en total
length(get_taxa_unique(physeq.pan, "Phylum"))
#23 filos (cony); a mi me dan 29 filos (nayla); 27 filos (KV)

#vemos los filos que hay y cuantos ASVs tiene cada uno
tabla_filos <- table(tax_table(physeq.pan)[, "Phylum"], exclude = NULL)
tabla_filos <- as.data.frame(tabla_filos)

#Ahora creamos un archvo phyloseq pero a nivel de filo
#aglomeramos nivel de filo y relativizamos abundancia
filo.ps.pan.para.grafico = tax_glom(physeq.pan, "Phylum", NArm = FALSE)#aglomeramos a nivel de filo e incluimos NA para ver su abundancia
#Ahora una vez aglomerados lo relativisamos
rel.filo.pan.para.grafico = transform_sample_counts(filo.ps.pan.para.grafico, function(x) x/sum(x)) #abundancia relativa

# Tabl de codigo de colores
Filo <- c("Proteobacteria","Actinobacteriota","Acidobacteriota","Bacteroidota",
          "Planctomycetota","Verrucomicrobiota","Myxococcota","Armaimonadota",
          "Cyanobacteria","Chloroflexota","Gemmatimonadota","Patescibacteria",
          "Eremiobacterota","Bdellovibrionota","Nitrospirota","Desulfobacterota_B",
          "Firmicutes","Dormibacterota","Methylomirabilota","Elusimicrobiota",
          "Omnitrophota","Firmicutes_A","Firmicutes_C","Dependentiae",
          "Desulfobacterota","Eisenbacteria")
cod <- c("#69AFB2","#E4B664","#6184BE","#E0E168","#C95F6F","#64B0D1","#1E5664","#57599D","#A1BE57",
         "#7A609F","#BB68A1","#67AA8C","#7AAD59","#CD6926","#C96299","#C23234","#4D7AAF","#9B66A0",
         "#918854","#355B8A","#040006","#7F7F7F","#DAD6C1","#738E3C","#8F3734","#E9E337")

color <- data.frame(Filo, cod)

meta <- meta(physeq.pan)
es <- row.names(meta)

# Graficamos
Filo.pan = plot_bar(rel.filo.pan.para.grafico, x= "Réplica", fill="Phylum") + guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~Tipo_de_muestra_sitio, scales = "free_x", ncol = 5) + 
  scale_fill_manual(values = color$cod, limits=color$Filo) + geom_bar(stat="identity", color="black", size=0) +
  labs( ylab="Frecuencia relativa", xlab("")) + #definimos nombre de ejes
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = ",")) + #comas para decimales
  theme(text = element_text(size = 10)) + theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(strip.text = element_text(size = 8 , color="black")) +
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 6))+ 
  labs(fill = "Filo")


setwd("~/1. TMKV/2. Results/10. Abundancia relativa")
ggsave("ARFilos_pan.png", width = 12.5, height = 7.3, dpi = 300)
# peripheral

#podemos ver cuantos filos hay en total
length(get_taxa_unique(physeq.core, "Phylum"))
#23 filos (cony); a mi me dan 29 filos (nayla); 27 filos (KV)

#vemos los filos que hay y cuantos ASVs tiene cada uno
tabla_filos <- table(tax_table(physeq.core)[, "Phylum"], exclude = NULL)
tabla_filos <- as.data.frame(tabla_filos)

#Ahora creamos un archvo phyloseq pero a nivel de filo
#aglomeramos nivel de filo y relativizamos abundancia
filo.ps.core.para.grafico = tax_glom(physeq.core, "Phylum", NArm = FALSE)#aglomeramos a nivel de filo e incluimos NA para ver su abundancia
#Ahora una vez aglomerados lo relativisamos
rel.filo.core.para.grafico = transform_sample_counts(filo.ps.core.para.grafico, function(x) x/sum(x)) #abundancia relativa

# Tabl de codigo de colores
Filo <- c("Proteobacteria","Actinobacteriota","Acidobacteriota","Bacteroidota",
          "Planctomycetota","Verrucomicrobiota","Myxococcota","Armaimonadota",
          "Cyanobacteria","Chloroflexota","Gemmatimonadota","Patescibacteria",
          "Eremiobacterota","Bdellovibrionota","Nitrospirota","Desulfobacterota_B",
          "Firmicutes","Dormibacterota","Methylomirabilota","Elusimicrobiota",
          "Omnitrophota","Firmicutes_A","Firmicutes_C","Dependentiae",
          "Desulfobacterota","Eisenbacteria")
cod <- c("#69AFB2","#E4B664","#6184BE","#E0E168","#C95F6F","#64B0D1","#1E5664","#57599D","#A1BE57",
         "#7A609F","#BB68A1","#67AA8C","#7AAD59","#CD6926","#C96299","#C23234","#4D7AAF","#9B66A0",
         "#918854","#355B8A","#040006","#7F7F7F","#DAD6C1","#738E3C","#8F3734","#E9E337")

color <- data.frame(Filo, cod)

meta <- meta(physeq.core)
es <- row.names(meta)

# Graficamos
Filo.core = plot_bar(rel.filo.core.para.grafico, x= "Réplica", fill="Phylum") + guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~Tipo_de_muestra_sitio, scales = "free_x", ncol = 5) + 
  scale_fill_manual(values = color$cod, limits=color$Filo) + geom_bar(stat="identity", color="black", size=0) +
  labs( ylab="Frecuencia relativa", xlab("")) + #definimos nombre de ejes
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = ",")) + #comas para decimales
  theme(text = element_text(size = 10)) + theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(strip.text = element_text(size = 8 , color="black")) +
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 6))+ 
  labs(fill = "Filo")


setwd("~/1. TMKV/2. Results/10. Abundancia relativa")
ggsave("ARFilos_core.png", width = 12.5, height = 7.3, dpi = 300)

#########################################
## Heat Map abundancia relativa phylum ##
#########################################

setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")

physeq.F2 <- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2,ps.Sue.COY19.F2)

# Revisión
sample_variables(physeq.F2)
rank_names(physeq.F2)
sample_names(physeq.F2)
meta(physeq.F2)

# Heat Map Filos 

physeq.F2.sNA <- subset_taxa(physeq.F2, Phylum!="NA")
physeq.F2.sNA.Phylum = tax_glom(physeq.F2.sNA, "Phylum", NArm = FALSE)

unique
taxa_names(physeq.F2.sNA.Phylum)

p<-plot_heatmap(physeq.F2.sNA.Phylum, "NMDS", "bray", "Tipo_de_muestra_sitio", "Phylum", 
                low="#FFFFCC", high="#000033", na.value="white", 
                sample.order = Sample.names_F2)

setwd("~/1. TMKV/2. Results/10. Abundancia relativa")
ggsave("HM-Filos-sNA-AbRel.png",plot=p, width = 15, height = 10)


######################
## Diagrama de Venn ##
######################

### Archivo Phyloseq 

# Cargar el archivo phyloseq
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")

ps.F2 <- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2,ps.Sue.COY19.F2)

# Más adelante se presentan problemas para graficar la red debido a los NA, entonces se modifica la tabla de taxonomia:
tax_table(ps.F2)[is.na(tax_table(ps.F2))] <- "NoId"
tax_table(ps.Liq.TAM19.F2)[is.na(tax_table(ps.Liq.TAM19.F2))] <- "NoId"
tax_table(ps.Sus.TAM19.F2)[is.na(tax_table(ps.Sus.TAM19.F2))] <- "NoId"
tax_table(ps.Liq.COY19.F2)[is.na(tax_table(ps.Liq.COY19.F2))] <- "NoId"
tax_table(ps.Sus.COY19.F2)[is.na(tax_table(ps.Sus.COY19.F2))] <- "NoId"
tax_table(ps.Sue.TAM19.F2)[is.na(tax_table(ps.Sue.TAM19.F2))] <- "NoId"
tax_table(ps.Sue.COY19.F2)[is.na(tax_table(ps.Sue.COY19.F2))] <- "NoId"

### Perfil taxonomico del conjunto de datos

# Perfil taxonomico: tabla de datos disponibles (util para consultar durante el análisis).
# Extraer y tabular los datos
taxotutable.F2 <- phyloseq::psmelt(ps.F2)
taxotutable.Liq.TAM19.F2 <- phyloseq::psmelt(ps.Liq.TAM19.F2)
taxotutable.Sus.TAM19.F2 <- phyloseq::psmelt(ps.Sus.TAM19.F2)
taxotutable.Liq.COY19.F2 <- phyloseq::psmelt(ps.Liq.COY19.F2)
taxotutable.Sus.COY19.F2 <- phyloseq::psmelt(ps.Sus.COY19.F2)
taxotutable.Sue.TAM19.F2 <- phyloseq::psmelt(ps.Sue.TAM19.F2)
taxotutable.Sue.COY19.F2 <- phyloseq::psmelt(ps.Sue.COY19.F2)

# Exportar/guardar tabla
setwd("~/1. TMKV/2. Results/12. Diagramas de Venn")
write.table(taxotutable.F2, file = "SummarizedTaxAssigments_ps.F2.tsv", sep = "\t", quote = FALSE)
write.table(taxotutable.Liq.TAM19.F2, file = "SummarizedTaxAssigments_ps.Liq.TAM19.F2.tsv", sep = "\t", quote = FALSE)
write.table(taxotutable.Sus.TAM19.F2, file = "SummarizedTaxAssigments_ps.Sus.TAM19.F2.tsv", sep = "\t", quote = FALSE)
write.table(taxotutable.Liq.COY19.F2, file = "SummarizedTaxAssigments_ps.Liq.COY19.F2.tsv", sep = "\t", quote = FALSE)
write.table(taxotutable.Sus.COY19.F2, file = "SummarizedTaxAssigments_ps.Sus.COY19.F2.tsv", sep = "\t", quote = FALSE)
write.table(taxotutable.Sue.TAM19.F2, file = "SummarizedTaxAssigments_ps.Sue.TAM19.F2.tsv", sep = "\t", quote = FALSE)
write.table(taxotutable.Sue.COY19.F2, file = "SummarizedTaxAssigments_ps.Sue.COY19.F2.tsv", sep = "\t", quote = FALSE)


#### El encabezado debe ser desplazado una columna a la derecha en el excel. 

TaxOTUTable.F2<-as.data.frame(taxotutable.F2)
TaxOTUTable.Liq.TAM19.F2<-as.data.frame(taxotutable.Liq.TAM19.F2)
TaxOTUTable.Sus.TAM19.F2<-as.data.frame(taxotutable.Sus.TAM19.F2)
TaxOTUTable.Liq.COY19.F2<-as.data.frame(taxotutable.Liq.COY19.F2)
TaxOTUTable.Sus.COY19.F2<-as.data.frame(taxotutable.Sus.COY19.F2)
TaxOTUTable.Sue.TAM19.F2<-as.data.frame(taxotutable.Sue.TAM19.F2)
TaxOTUTable.Sue.COY19.F2<-as.data.frame(taxotutable.Sue.COY19.F2)

colnames(TaxOTUTable)


## Diagrama de Venn para OTUs (ASV) ##

OTU.list <- list(LiqTAM=TaxOTUTable.Liq.TAM19.F2$OTU,
                 SusTAM=TaxOTUTable.Sus.TAM19.F2$OTU,
                 SueTAM=TaxOTUTable.Sue.TAM19.F2$OTU,
                 LiqCOY=TaxOTUTable.Liq.COY19.F2$OTU,
                 SusCOY=TaxOTUTable.Sus.COY19.F2$OTU,
                 SueCOY=TaxOTUTable.Sue.COY19.F2$OTU)

ggVennDiagram(OTU.list,label_alpha=0) + 
  scale_fill_gradient(low="blue",high = "red")

venn <- Venn(OTU.list)
data <- process_data(venn)
data

# https://cran.r-project.org/web/packages/ggVennDiagram/readme/README.html
#https://r-charts.com/part-whole/ggvenndiagram/

Colors<- c("#F8766D","#OOBA38","#B79F00","#00BFC4","#619CFF","#F564E3")

VennDiagram_OTU <- ggVennDiagram(OTU.list,label_alpha=0, label = "count",label_color = "white", color = 1, lwd = 0.7 ) + 
  scale_fill_gradient(low="blue",high = "red")

setwd("~/1. TMKV/2. Results/12. Diagramas de Venn")
ggsave("VennDiagram_OTU_F2.png", width = 20, height = 10, dpi = 300)


## Diagrama de Venn para Phylum ##


Phylum.list <- list(LiqTAM=TaxOTUTable.Liq.TAM19.F2$Phylum,
                    SusTAM=TaxOTUTable.Sus.TAM19.F2$Phylum,
                    SueTAM=TaxOTUTable.Sue.TAM19.F2$Phylum,
                    LiqCOY=TaxOTUTable.Liq.COY19.F2$Phylum,
                    SusCOY=TaxOTUTable.Sus.COY19.F2$Phylum,
                    SueCOY=TaxOTUTable.Sue.COY19.F2$Phylum)
# Saber elementos unicos en una lista
unique(unlist(Phylum.list))

# https://cran.r-project.org/web/packages/ggVennDiagram/readme/README.html
#https://r-charts.com/part-whole/ggvenndiagram/

VennDiagram_Phylum <- ggVennDiagram(Phylum.list,label_alpha=0, label = "count",label_color = "white", color = 1, lwd = 0.7 ) + 
  scale_fill_gradient(low="blue",high = "red")

setwd("~/1. TMKV/2. Results/12. Diagramas de Venn")
ggsave("VennDiagram_Phylum_F2.png", width = 20, height = 10, dpi = 300)


## Diagrama de Venn para Genus ##


Genus.list <- list(LiqTAM=TaxOTUTable.Liq.TAM19.F2$Genus,
                   SusTAM=TaxOTUTable.Sus.TAM19.F2$Genus,
                   SueTAM=TaxOTUTable.Sue.TAM19.F2$Genus,
                   LiqCOY=TaxOTUTable.Liq.COY19.F2$Genus,
                   SusCOY=TaxOTUTable.Sus.COY19.F2$Genus,
                   SueCOY=TaxOTUTable.Sue.COY19.F2$Genus)
# Saber elementos unicos en una lista
unique(unlist(Genus.list))

# https://cran.r-project.org/web/packages/ggVennDiagram/readme/README.html
#https://r-charts.com/part-whole/ggvenndiagram/

VennDiagram_Genus <- ggVennDiagram(Genus.list,label_alpha=0, label = "count",label_color = "white", color = 1, lwd = 0.7 ) + 
  scale_fill_gradient(low="blue",high = "red")

setwd("~/1. TMKV/2. Results/12. Diagramas de Venn")
ggsave("VennDiagram_Genus_F2.png", width = 20, height = 10, dpi = 300)






