
#########################################
## Enlace - Analisis de secuencias 16S ##
#########################################

##############
## PAQUETES ##
##############

library("dada2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
library(gridExtra)
library(ggplot2)
library(DECIPHER)
library(grid)
library(Biostrings)
library(phyloseq)
library(vegan)
library(readr)
library(devtools)
library(iNEXT)
library(ggplot2)
library(tibble)
library(readxl)
library(ggrepel)
library(car)
library(ggpubr)
library(rstatix)
library(phangorn)
library(tidyverse)
library(plyr)
library(gridExtra)
library(kableExtra)
library(xtable)
library(ggpubr)
library(BiocManager)
library(DESeq2)
library(microbiome)
library(philr)
library(btools)
library(fantaxtic)
library(ampvis2)
library(tsnemicrobiota)


######################################
##### ANALISIS SECUENCIAS BRUTAS #####
######################################

# Análisis Rawdata TMKV

# Fijar directorio de trabajo
setwd("~/1. TMKV/2. Results/0. Rawdata")

# Leer archivo rawdata
Dataframe_rawdata <- read.csv("~/1. TMKV/2. Results/0. Rawdata/Dataframe_rawdata.csv", sep=";")

# Revisar
View(Dataframe_rawdata)

# Cambiar nombre de columna
colnames(Dataframe_rawdata)
Dataframe_rawdata <- tibble::column_to_rownames(Dataframe_rawdata,var="Tesis")
View(Dataframe_rawdata)

###################### 
## 1. BASES TOTALES ##
######################

# Seleccionar datos de interés en el dataframe
dat.Bases.totales <- Dataframe_rawdata %>% tibble::rownames_to_column(var="outlier") %>% group_by("Tipo.de.muestra") %>% 
  mutate(is_outlier=ifelse(is_outlier(Bases.totales), Bases.totales, as.numeric(NA)))

# COMPARAR
# ANOVA de una vía (independiente)
one.way <- aov(Bases.totales ~ Tipo.de.muestra, data = dat.Bases.totales)
summary(one.way)

# 1. Homoestacidad 
LT=leveneTest(Bases.totales ~ Tipo.de.muestra, data = dat.Bases.totales)
setwd("~/1. TMKV/2. Results/0. Rawdata")
write.table(LT, file = "LT.Bases.totales.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
# 2. Normalidad

aov_residuals <- residuals(object = one.way )
ST = shapiro.test(x = aov_residuals )
setwd("~/1. TMKV/2. Results/0. Rawdata")
ST.list=sapply(ST, unlist)
write.table(ST.list, file = "ST.Bases.totales.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# 3. Comparación por pares
stat.test=tukey_hsd(one.way)
# 4. Grafico 
bxp.Bases.totales <- ggboxplot(dat.Bases.totales, x = "Tipo.de.muestra", y = "Bases.totales", fill = "Tipo.de.muestra", 
                               palette = c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))+
  theme_grey(base_size = 10) +
  theme(legend.position = "none") +
  xlab("Tipo de muestra") + ylab("Bases totales") +
  stat_pvalue_manual(
    stat.test, 
    y.position = 35, step.increase = 0.1, size =4,
    label = "p.adj.signif", hide.ns = TRUE,color = "black", linetype = 9,tip.length = 0)
bxp.Bases.totales

######################### 
## 2. LECTURAS TOTALES ##
#########################

# Seleccionar datos de interés en el dataframe
dat.Lecturas.totales <- Dataframe_rawdata %>% tibble::rownames_to_column(var="outlier") %>% group_by("Tipo.de.muestra") %>% 
  mutate(is_outlier=ifelse(is_outlier(Lecturas.totales), Lecturas.totales, as.numeric(NA)))

# COMPARAR
# ANOVA de una vía (independiente)
one.way <- aov(Lecturas.totales ~ Tipo.de.muestra, data = dat.Lecturas.totales)
summary(one.way)

# 1. Homoestacidad
LT=leveneTest(Lecturas.totales ~ Tipo.de.muestra, data = dat.Lecturas.totales)
setwd("~/1. TMKV/2. Results/0. Rawdata")
write.table(LT, file = "LT.Lecturas.totales.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
# 2. Normalidad
aov_residuals <- residuals(object = one.way )
ST = shapiro.test(x = aov_residuals )
setwd("~/1. TMKV/2. Results/0. Rawdata")
ST.list=sapply(ST, unlist)
write.table(ST.list, file = "ST.Lecturas.totales.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# 3. Comparación por pares
stat.test=tukey_hsd(one.way)
# 4. Grafico 
bxp.Lecturas.totales <- ggboxplot(dat.Lecturas.totales, x = "Tipo.de.muestra", y = "Lecturas.totales", fill = "Tipo.de.muestra", 
                                  palette = c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))+
  theme_grey(base_size = 10) +
  theme(legend.position = "none") +
  xlab("Tipo de muestra") + ylab("Lecturas totales") +
  stat_pvalue_manual(
    stat.test, 
    y.position = 35, step.increase = 0.1, size =4,
    label = "p.adj.signif", hide.ns = TRUE,color = "black", linetype = 9,tip.length = 0)
bxp.Lecturas.totales

######################
## 3. PORCENTAJE GC ##
######################

# Seleccionar datos de interés en el dataframe
dat.Porcentaje.GC <- Dataframe_rawdata %>% tibble::rownames_to_column(var="outlier") %>% group_by("Tipo.de.muestra") %>% 
  mutate(is_outlier=ifelse(is_outlier(Porcentaje.GC), Porcentaje.GC, as.numeric(NA)))

# COMPARAR
# ANOVA de una vía (independiente)
one.way <- aov(Porcentaje.GC ~ Tipo.de.muestra, data = dat.Porcentaje.GC)
summary(one.way)

# 1. Homoestacidad
LT=leveneTest(Porcentaje.GC ~ Tipo.de.muestra, data = dat.Porcentaje.GC)
setwd("~/1. TMKV/2. Results/0. Rawdata")
write.table(LT, file = "LT.Porcentaje.GC.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
# 2. Normalidad
aov_residuals <- residuals(object = one.way )
ST = shapiro.test(x = aov_residuals )
setwd("~/1. TMKV/2. Results/0. Rawdata")
ST.list=sapply(ST, unlist)
write.table(ST.list, file = "ST.Porcentaje.GC.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# 3. Comparación por pares
stat.test=tukey_hsd(one.way)
# 4. Grafico 
bxp.Porcentaje.GC <- ggboxplot(dat.Porcentaje.GC, x = "Tipo.de.muestra", y = "Porcentaje.GC", fill = "Tipo.de.muestra", 
                               palette = c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))+
  theme_grey(base_size = 10) +
  theme(legend.position = "none") +
  xlab("Tipo de muestra") + ylab("Porcentaje GC") +
  stat_pvalue_manual(
    stat.test, 
    y.position = 35, step.increase = 0.1, size =4,
    label = "p.adj.signif", hide.ns = TRUE,color = "black", linetype = 9,tip.length = 0)
bxp.Porcentaje.GC

#######################
## 4. PORCENTAJE Q30 ##
#######################

# Seleccionar datos de interés en el dataframe
dat.Porcentaje.Q30 <- Dataframe_rawdata %>% tibble::rownames_to_column(var="outlier") %>% group_by("Tipo.de.muestra") %>% 
  mutate(is_outlier=ifelse(is_outlier(Porcentaje.Q30), Porcentaje.Q30, as.numeric(NA)))

# COMPARAR
# ANOVA de una vía (independiente)
one.way <- aov(Porcentaje.Q30 ~ Tipo.de.muestra, data = dat.Porcentaje.Q30)
summary(one.way)

# 1. Homoestacidad
LT=leveneTest(Porcentaje.Q30 ~ Tipo.de.muestra, data = dat.Porcentaje.Q30)
setwd("~/1. TMKV/2. Results/0. Rawdata")
write.table(LT, file = "LT.Porcentaje.Q30.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
# 2. Normalidad
aov_residuals <- residuals(object = one.way )
ST = shapiro.test(x = aov_residuals )
setwd("~/1. TMKV/2. Results/0. Rawdata")
ST.list=sapply(ST, unlist)
write.table(ST.list, file = "ST.Porcentaje.Q30.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# 3. Comparación por pares
stat.test=tukey_hsd(one.way)
# 4. Grafico 
bxp.Porcentaje.Q30 <- ggboxplot(dat.Porcentaje.Q30, x = "Tipo.de.muestra", y = "Porcentaje.Q30", fill = "Tipo.de.muestra", 
                                palette = c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))+
  theme_grey(base_size = 10) +
  theme(legend.position = "none") +
  xlab("Tipo de muestra") + ylab("Porcentaje Q30") +
  stat_pvalue_manual(
    stat.test, 
    y.position = 35, step.increase = 0.1, size =4,
    label = "p.adj.signif", hide.ns = TRUE,color = "black", linetype = 9,tip.length = 0)
bxp.Porcentaje.Q30

#######################
## 5. PORCENTAJE Q20 ##
#######################

# Seleccionar datos de interés en el dataframe
dat.Porcentaje.Q20 <- Dataframe_rawdata %>% tibble::rownames_to_column(var="outlier") %>% group_by("Tipo.de.muestra") %>% 
  mutate(is_outlier=ifelse(is_outlier(Porcentaje.Q20), Porcentaje.Q20, as.numeric(NA)))

# COMPARAR
# ANOVA de una vía (independiente)
one.way <- aov(Porcentaje.Q20 ~ Tipo.de.muestra, data = dat.Porcentaje.Q20)
summary(one.way)

# 1. Homoestacidad
LT=leveneTest(Porcentaje.Q20 ~ Tipo.de.muestra, data = dat.Porcentaje.Q20)
setwd("~/1. TMKV/2. Results/0. Rawdata")
write.table(LT, file = "LT.Porcentaje.Q20.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
# 2. Normalidad
aov_residuals <- residuals(object = one.way )
ST = shapiro.test(x = aov_residuals )
setwd("~/1. TMKV/2. Results/0. Rawdata")
ST.list=sapply(ST, unlist)
write.table(ST.list, file = "ST.Porcentaje.Q20.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# 3. Comparación por pares
stat.test=tukey_hsd(one.way)

# 4. Grafico 
bxp.Porcentaje.Q20 <- ggboxplot(dat.Porcentaje.Q20, x = "Tipo.de.muestra", y = "Porcentaje.Q20", fill = "Tipo.de.muestra", 
                                palette = c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))+
  theme_grey(base_size = 10) +
  theme(legend.position = "none") +
  xlab("Tipo de muestra") + ylab("Porcentaje Q20") +
  stat_pvalue_manual(
    stat.test, 
    y.position = 35, step.increase = 0.1, size =4,
    label = "p.adj.signif", hide.ns = TRUE,color = "black", linetype = 9,tip.length = 0)
bxp.Porcentaje.Q20

###############################
## GRAFICO SECUENCIAS BRUTAS ##
###############################

ggarrange(bxp.Bases.totales, bxp.Lecturas.totales,bxp.Porcentaje.GC,bxp.Porcentaje.Q20,bxp.Porcentaje.Q30, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3) %>%
  ggexport(filename = "Boxplots.rawdata.pdf")

############################################
##### ANALISIS DE SECUENCIAS CON DADA2 #####
############################################
####################################
## REVISION CALIDAD DE SECUENCIAS ##
####################################

# Se deben revisar previamente los largos de secuencia, tanto para forward 
# como de reverse. Las secuencias pueden revisarse en el programa fastaqc.                                                             

# Graficos en R
# Definir path con la dirección donde se encuentran las secuencias a analizar. 

setwd("/home/karla/TMKV/Samples")
path <- "/home/karla/TMKV/Samples" 
list.files (path) # Corroborar los archivos
length(list.files (path)) #Corroborar la cantidad de archivos

# Crear lista con archivos de secuencia, segun un patrón:

Sue <- sort(list.files(path, pattern = "Sue-", full.names=TRUE))
fnFs <- sort(list.files(path, pattern="_F", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R", full.names = TRUE))

# Grafico de visión general:

## Gráficos de calidad para Liquenes Forward:
plotQualityProfile(fnFs[c(1:20)]) +
  scale_x_continuous(name="pb",breaks=seq(0,300,100)) +
  scale_y_continuous(name="Calidad phred", breaks=seq(0, 40, 10)) 

## Gráficos de calidad para Sustrato Forward:
plotQualityProfile(fnFs[c(31:50)]) +
  scale_x_continuous(name="pb",breaks=seq(0,300,100)) +
  scale_y_continuous(name="Calidad phred", breaks=seq(0, 40, 10)) 

## Gráficos de calidad para Liquenes Reverse:
plotQualityProfile(fnRs[c(1:20)]) +
  scale_x_continuous(name="pb",breaks=seq(0,300,100)) +
  scale_y_continuous(name="Calidad phred", breaks=seq(0, 40, 10)) 

## Gráficos de calidad para Sustrato Reverse:
plotQualityProfile(fnRs[c(31:50)]) +
  scale_x_continuous(name="pb",breaks=seq(0,300,100)) +
  scale_y_continuous(name="Calidad phred", breaks=seq(0, 40, 10)) 

## Gráficos de calidad para Suelos Forward y reverse:
plotQualityProfile(Sue[c(1,3,5,7,9,11,13,15,17,19,2,4,6,8,10,12,14,16,18,20)]) +
  scale_x_continuous(name="pb",breaks=seq(0,300,100)) +
  scale_y_continuous(name="Calidad phred", breaks=seq(0, 40, 10)) 

# Gráfico con zoom para definir el corte:

## Gráficos de calidad para Liquenes Forward:

plotQualityProfile(fnFs[c(1:20)])+ 
  scale_x_continuous(name="pb",breaks=seq(285,301,2), limits=c(285,301))+
  scale_y_continuous(name="Calidad phred", breaks=seq(24, 40, 2), limits=c(24,40))+
  geom_hline(yintercept=28, linetype="dotted", color = "blue", size=0.5)+
  geom_hline(yintercept=30, linetype="dotted", color = "blue", size=0.5)+
  geom_vline(xintercept=287, linetype="dotted", color = "red", size=0.5)+
  theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.5))

## Gráficos de calidad para Sustratos Forward:

plotQualityProfile(fnFs[c(31:50)])+ 
  scale_x_continuous(name="pb",breaks=seq(285,301,2), limits=c(285,301))+
  scale_y_continuous(name="Calidad phred", breaks=seq(24, 40, 2), limits=c(24,40))+
  geom_hline(yintercept=28, linetype="dotted", color = "blue", size=0.5)+
  geom_hline(yintercept=30, linetype="dotted", color = "blue", size=0.5)+
  geom_vline(xintercept=287, linetype="dotted", color = "red", size=0.5)+
  theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.5))

## Gráficos de calidad para Liquenes Reverse:

plotQualityProfile(fnRs[c(1:20)])+ 
  scale_x_continuous(name="pb",breaks=seq(195,210,2), limits=c(195,210))+
  scale_y_continuous(name="Calidad phred", breaks=seq(24, 40, 2), limits=c(24,40))+
  geom_hline(yintercept=28, linetype="dotted", color = "blue", size=0.5)+
  geom_hline(yintercept=30, linetype="dotted", color = "blue", size=0.5)+
  geom_vline(xintercept=201, linetype="dotted", color = "red", size=0.5)+
  theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.5))

## Gráficos de calidad para Sustratos Reverse:

plotQualityProfile(fnRs[c(31:50)])+ 
  scale_x_continuous(name="pb",breaks=seq(195,210,2), limits=c(195,210))+
  scale_y_continuous(name="Calidad phred", breaks=seq(24, 40, 2), limits=c(24,40))+
  geom_hline(yintercept=28, linetype="dotted", color = "blue", size=0.5)+
  geom_hline(yintercept=30, linetype="dotted", color = "blue", size=0.5)+
  geom_vline(xintercept=201, linetype="dotted", color = "red", size=0.5)+
  theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.5))

## Gráficos de calidad para Suelos Forward y reverse:

plotQualityProfile(Sue[c(1,3,5,7,9,11,13,15,17,19)]) +
  scale_x_continuous(name="pb",breaks=seq(285,301,2), limits=c(285,301))+
  scale_y_continuous(name="Calidad phred", breaks=seq(24, 40, 2), limits=c(24,40))+
  geom_hline(yintercept=28, linetype="dotted", color = "blue", size=0.5)+
  geom_hline(yintercept=30, linetype="dotted", color = "blue", size=0.5)+
  geom_vline(xintercept=287, linetype="dotted", color = "red", size=0.5)+
  theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.5)) 


plotQualityProfile(Sue[c(2,4,6,8,10,12,14,16,18,20)])+ 
  scale_x_continuous(name="pb",breaks=seq(195,210,2), limits=c(195,210))+
  scale_y_continuous(name="Calidad phred", breaks=seq(24, 40, 2), limits=c(24,40))+
  geom_hline(yintercept=28, linetype="dotted", color = "blue", size=0.5)+
  geom_hline(yintercept=30, linetype="dotted", color = "blue", size=0.5)+
  geom_vline(xintercept=201, linetype="dotted", color = "red", size=0.5)+
  theme(axis.text.x=element_text(size=10, angle = 90, vjust=0.5))


###########
## DADA2 ##
###########

# Ejemplo con Liq_COY19 - La tuberia se ejecuto para cada set de muestras (LiqCOY, LiqTAM, SusCOY, SusTAM, SueCOY, SueTAM)

# Definir el path con la direccion donde se encuentran las secuencias.
path.Liq_COY19 <-"~/1. TMKV/1. Samples/Liq_COY19"

# Ver los archivos que hay en "path"
list.files (path.Liq_COY19)

# Listar los archivos forward y reverse
fnFs.Liq_COY19 <- sort(list.files(path.Liq_COY19, pattern="_F.fastq.gz")) 
fnRs.Liq_COY19 <- sort(list.files(path.Liq_COY19, pattern="_R.fastq.gz")) 

# Revisar que contengan las 10 secuencias en cada objeto (o cantidad de secuencias a trabajar)
length(fnFs.Liq_COY19)
length(fnRs.Liq_COY19)

# Extraer los nombres de los archivos
sample.names.Liq_COY19 <- paste0((sapply(strsplit(basename(fnFs.Liq_COY19), "_"), `[`, 1))) 

# Revisar los nombres obtenidos
sample.names.Liq_COY19

# Especificar el path completo para evitar errores de ambiguedad. 
fnFs.Liq_COY19 <- file.path(path.Liq_COY19, fnFs.Liq_COY19)
fnRs.Liq_COY19 <- file.path(path.Liq_COY19, fnRs.Liq_COY19)

# Crear directorio donde dejar las secuencias editadas
filtRs.Liq_COY19 <- file.path(path.Liq_COY19, "filtered_287-201", paste0(sample.names.Liq_COY19, "_R_filt.fastq.gz"))  
filtFs.Liq_COY19 <- file.path(path.Liq_COY19, "filtered_287-201", paste0(sample.names.Liq_COY19, "_F_filt.fastq.gz"))

# Revisar que el nombre de los archivos este correcto
filtFs.Liq_COY19
filtRs.Liq_COY19

# Revisar el numero de archivos que se a generado
length(filtFs.Liq_COY19)
length(filtRs.Liq_COY19)

# Filtrar las secuencias
out <- filterAndTrim(fnFs.Liq_COY19, filtFs.Liq_COY19, fnRs.Liq_COY19, filtRs.Liq_COY19, truncLen=c(287,201), trimLeft=c(19,25), 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Guardar el espacio de trabajo cada tanto.
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")   

# Revisar
head(out) #head muestra las primeras 6 columnas

# Desreplicar
derepFs.Liq_COY19 <- derepFastq(filtFs.Liq_COY19, verbose=TRUE) 
derepRs.Liq_COY19 <- derepFastq(filtRs.Liq_COY19, verbose=TRUE)

# Guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")    

# Nombrar los dereplicados segun los nombres de muestras 
names(derepFs.Liq_COY19) <- sample.names.Liq_COY19 
names(derepRs.Liq_COY19) <- sample.names.Liq_COY19

# Correr el metodo "learnErrors".
errF.Liq_COY19 <- learnErrors(filtFs.Liq_COY19, multithread=10) 
errR.Liq_COY19 <- learnErrors(filtRs.Liq_COY19, multithread=10)

# Graficar los errores para cada par
pErrF.Liq_COY19 <-plotErrors(errF.Liq_COY19, nominalQ=TRUE) + 
  ggtitle("Secuencias forward") + theme(plot.title = element_text(hjust = 0.5, size=12))
pErrR.Liq_COY19 <-plotErrors(errR.Liq_COY19, nominalQ=TRUE) + 
  ggtitle("Secuencias reverse") + theme(plot.title = element_text(hjust = 0.5, size=12))

f1.Liq_COY19<-grid.arrange(pErrF.Liq_COY19, pErrR.Liq_COY19, ncol=2,
                           top = textGrob("Método learnErrors para liquenes Rufescens Coyhaique", 
                                          gp=gpar(fontsize=17,font=3))) 
g.Liq_COY19 <- arrangeGrob(f1.Liq_COY19)

# Guardar grafico
setwd("~/1. TMKV/2. Results/2. pErr")
ggsave("pErr_Liq-COY19.png",plot=g.Liq_COY19, width = 20, height = 14)
setwd("~/1. TMKV")

# Calcular la inferencia de las Amplicon Sequence Variants (ASVs).
dadaFs.Liq_COY19 <- dada(derepFs.Liq_COY19, err=errF.Liq_COY19, pool=TRUE, multithread=10)
dadaRs.Liq_COY19 <- dada(derepRs.Liq_COY19, err=errR.Liq_COY19, pool=TRUE, multithread=10)

# Guardar el espacio de trabajo cada tanto.
save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/DADA16S_TMKV_LiqCOY19.RData")   

# Revisar
dadaFs.Liq_COY19
dadaRs.Liq_COY19

# Fusionar las lecturas Forward y Reverse para obtener las secuencias sin ruido completas. 
mergers.Liq_COY19 <- mergePairs(dadaFs.Liq_COY19, derepFs.Liq_COY19, dadaRs.Liq_COY19, derepRs.Liq_COY19, minOverlap = 20, verbose=TRUE)

# Guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")    

# Revisar
head(mergers.Liq_COY19[[1]])

# Generar una tabla de variantes de secuencia de amplicones (ASV)
seqtab.Liq_COY19 <- makeSequenceTable(mergers.Liq_COY19) 

# Revisar dimensiones 
dim(seqtab.Liq_COY19) 

# Generar dataframe para graficar
dist.Liq_COY19 <- table(nchar(getSequences(seqtab.Liq_COY19)))
dist.Liq_COY19 <-as.data.frame(dist.Liq_COY19)
View(dist.Liq_COY19)

dist.Liq_COY19$Var1 <- as.numeric(as.character(dist.Liq_COY19$Var1))
dist.Liq_COY19$Freq <- as.numeric(as.character(dist.Liq_COY19$Freq))

write.csv(dist.Liq_COY19,"~/1. TMKV/2. Results/3. Dist. largo seq/dist.Liq_COY19.csv", row.names = FALSE)

# Graficar
p.Liq_COY19.dls <- ggplot(data=dist.Liq_COY19, aes(x=Var1, y=Freq)) + geom_histogram(stat='identity') + 
  scale_x_continuous(breaks = seq(355, 385, 5),limits=c(355,385), name = "Largo secuencia (pb)") + scale_y_continuous(name= "Frecuencia")

# Guardar grafico
setwd("~/1. TMKV/2. Results/3. Dist. largo seq")
ggsave("Largo_secuencia_Liq_COY19.png", width = 7, height = 5, device= "png")
setwd("~/1. TMKV")

# Identificar y eliminar secuencias quimericas
seqtab.nochim.Liq_COY19 <- removeBimeraDenovo(seqtab.Liq_COY19, method="consensus", multithread=20, verbose=TRUE)

# Guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")

# Revisar
dim(seqtab.nochim.Liq_COY19)

# Porcentaje de quimeras
PorcentajeQuimeras.Liq_COY19 <- as.data.frame(sum(seqtab.nochim.Liq_COY19)/sum(seqtab.Liq_COY19))*100
rownames(PorcentajeQuimeras.Liq_COY19) <- "% Quimeras"

PorcentajeQuimeras.Liq_COY19 <- rename(PorcentajeQuimeras.Liq_COY19, "Liquenes Rufescens Coyhaique"="sum(seqtab.nochim.Liq_COY19)/sum(seqtab.Liq_COY19)")


# Verificación final de nuestro progreso, observemos la cantidad de lecturas que lograron en cada paso del proceso:

getN <- function(x) sum(getUniques(x))
track.Liq_COY19 <- cbind(out, sapply(dadaFs.Liq_COY19, getN), sapply(dadaRs.Liq_COY19, getN), sapply(mergers.Liq_COY19, getN), rowSums(seqtab.nochim.Liq_COY19))
colnames(track.Liq_COY19) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.Liq_COY19) <- sample.names.Liq_COY19
head(track.Liq_COY19)

# Se recomienda guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")

# Exportar la tabla
setwd("~/1. TMKV/2. Results/4. Resumen Merge+Quimeras")
write.csv(track.Liq_COY19, "Resumen-Secuencias_Liq-COY19.csv")
setwd("~/1. TMKV")

# Guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")

#########################################################
## ASIGNACIÓN TAXONOMICA CON DECIPHER, MEDIANTE IdTaxa ##
#########################################################

# Se descargo desde la pagina de decipher el trainingset GTDB_r95
# Referencia: http://www2.decipher.codes/Downloads.html

# Cargamos el ambiente de trabajo de GTDB_r95 
load("~/1. TMKV/3. GTDB/GTDB_r95-mod_August2020.RData")

# Crear un DNAStringSet a partir de los ASV
dna.Liq_COY19 <- DNAStringSet(getSequences(seqtab.nochim.Liq_COY19)) 

#  Clasificar secuencias de acuerdo al trainingset
ids.Liq_COY19 <- IdTaxa(dna.Liq_COY19, trainingSet, strand="top", threshold = 60, processors=NULL, verbose=FALSE)

# Niveles de identificación
ranks.Liq_COY19 <- c("domain", "phylum", "class", "order", "family", "genus", "species") 

# Convertir el objeto de salida de la clase "Taxa" en una matriz analoga a la salida de assignTaxonomy
taxid.Liq_COY19 <- t(sapply(ids.Liq_COY19, function(x) {
  m.Liq_COY19 <- match(ranks.Liq_COY19, x$rank)
  taxa.Liq_COY19 <- x$taxon[m.Liq_COY19]
  taxa.Liq_COY19[startsWith(taxa.Liq_COY19, "unclassified_")] <- NA
  taxa.Liq_COY19
}))
colnames(taxid.Liq_COY19) <- ranks.Liq_COY19; rownames(taxid.Liq_COY19) <- getSequences(seqtab.nochim.Liq_COY19)

taxa.Liq_COY19 <- taxid.Liq_COY19

# Guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")

# Revisar
taxa.print.Liq_COY19 <- taxa.Liq_COY19 #removemos nombres de filas solo para mostrar
rownames(taxa.print.Liq_COY19) <- NULL
head(taxa.print.Liq_COY19)

######################
## AGREGAR METADATA ##
######################

#definir tema
theme_set(theme_light())

# Añadir la metadata
metadata.Liq_COY19 <- read.table("~/1. TMKV/4. Metadata/Metadata_Liq_COY19.csv", sep = ";")

# Agregar nombres a las columnas
colnames(metadata.Liq_COY19) <- c("Muestra", "Tipo_de_Muestra", "Tipo_de_muestra_sitio", "Origen", "Sitio","Nombre_Tesis","Especie","Cianobionte","Mes_Recoleccion","Año_Recoleccion","Latitud","Longitud","Altura", "Ambiente","Distribución_global","Nombre_Terreno","Nombre_GBIF","Numero_orden_secuenciación","Réplica","Número")
metadata.Liq_COY19 <-metadata.Liq_COY19[1:10,]
rownames(metadata.Liq_COY19) <- sample.names.Liq_COY19

# Agregar nombre a los niveles taxonomicos a la tabla de la taxonomia
colnames(taxa.Liq_COY19) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")

# Importar a phyloseq como tal
ps.Liq.COY19 <- phyloseq(otu_table(seqtab.nochim.Liq_COY19, taxa_are_rows=FALSE), tax_table(taxa.Liq_COY19), sample_data(metadata.Liq_COY19))

# Revisar
ps.Liq.COY19

# Exportar archivo phyloseq

saveRDS(ps.Liq.COY19, "~/1. TMKV/5. Phyloseq/ps.Liq.COY19.rds")

# Guardar el espacio de trabajo cada tanto. 
save.image("~/1. TMKV/0. Scripts&Environment/DADA16S_TMKV_LiqCOY19.RData")

########################
## ARBOL FILOGENETICO ##
########################

# Llamar phyloseqs
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")

# Unir Phyloseq
physeq.F2 <- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2,ps.Sue.COY19.F2)
physeq.F2

# Generar data.frame desde OTU table
OTU1 = as(otu_table(physeq.F2), "matrix")

if(taxa_are_rows(physeq.F2)){OTU1 <- t(OTU1)}

OTUdf = as.data.frame(OTU1)

# Revisar data.frame 
str(OTUdf)

# Extraer secuencias desde el nombre de las columnas
seq<-c(colnames(OTUdf))
str(seq)
class(seq)
length(seq)

# Sumar las columnas (abundancia de ASV)
csum <- data.frame(colSums(OTUdf))
str(csum)
class(csum)
length(csum)
head(csum)

# Extraer abundancias 
abundance<- c(csum[,1])
str(abundance)
class(abundance)
length(abundance)

df.s_a<-(as.data.frame(abundance,seq))
class(df.s_a)
str(df.s_a)
head(df.s_a)
df.s_a[1,]
df.s_a[,1]

rownames_to_column(df.s_a, var = "sequence")   -> df.s_a
hdf<-head(df.s_a)
view(hdf)

seqs<-getSequences(df.s_a)

# Propagar los nombres de las secuencias al árbol
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")

# Guardar ambiente de trabajo cada tanto
save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Árblfilogenetico_TMKV.RData") 

dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# Guardar ambiente de trabajo cada tanto
save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Árblfilogenetico_TMKV.RData") 

# Tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(physeq.F2)

# Tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(physeq.F2)
colnames(Tabla.Taxonomia) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") 

# Tabla de metadatos
meta.data<-meta(physeq.F2)

# Reconstruccion del objeto phyloseq:
physeq.F2 <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                      tax_table(Tabla.Taxonomia),
                      sample_data(meta.data),
                      phy_tree(fitGTR$tree))
physeq.F2
rank_names(physeq.F2)
view(meta(physeq.F2))
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(physeq.F2, "ps.TMKV.rds")

# Guardar ambiente de trabajo 
save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Árblfilogenetico_TMKV.RData") 


#####################################
##### OTROS ANALISIS Y GRAFICOS #####
#####################################
############################
## FILTRADO DE SECUENCIAS ##
############################

## Ejemplo con Liq TAM19 ##

setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19 <- readRDS("ps.Liq.TAM19.rds")

# 1. Filtro por dominio: Eliminar lo que NO es bacteria 

# Revisar si se detecto algo que NO sea bacterias:
get_taxa_unique(ps.Liq.TAM19, "Domain")
# Filtrar
ps.Liq.TAM19.F1 <- subset_taxa(ps.Liq.TAM19, Domain =="Bacteria")
# Guardar phyloseq filtrado
saveRDS(ps.Liq.TAM19.F1, "ps.Liq.TAM19.F1.rds")

# 2. Filtrado por abundancia

# Grafico filtro por abundacia 

setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.COY19.F1 <- readRDS("ps.Liq.COY19.F1.rds")

# Primero, determinar el mejor corte por abundancia

# Sumar las taxas (abundancia de cada seq)
Suma_taxas<-data.frame(taxa_sums(ps.Liq.COY19.F1))
# Revisar
Head_ST<-head(Suma_taxas)

# Convertir las secuencias en una variable
Suma_taxas_var <- rownames_to_column(Suma_taxas, var = "sequence")

# Reordenar la tabla de manera ascendente 
ST_asc<-arrange(Suma_taxas_var, taxa_sums.ps.Liq.COY19.F1.)

# Revisar
Head_ST_asc<-head(ST_asc)

# Ordenar la tabla para graficar
str(ST_asc)

# Agregar columna con secuencia 
ST<-mutate(ST_asc, order=1:2614)   ## Cambiar segun la cantidad de taxas del archivo phyloseq

# Graficamos 

# Grafico completo
q <-ggplot(ST, aes(x=order, y=taxa_sums.ps.Liq.COY19.F1.)) + geom_point()
pS <- q + scale_y_continuous(name="Reads") + scale_x_continuous(name="Secuencia") + geom_hline(yintercept=10, linetype="dotted", color = "red", size=1)
pS

# Zoom final 
ST_zoom2 <-ST[1:1900, ]
Zoom_pS <-ggplot(ST_zoom2, aes(x=order, y=taxa_sums.ps.Liq.COY19.F1.)) +
  geom_point() + 
  scale_y_continuous(limits = c(0,50), breaks=seq(0,50,10))+
  scale_x_continuous(limits=c(0,1900), breaks=seq(0,1900,100))  + 
  geom_hline(yintercept=10, linetype="dotted", color = "red", size=1) + 
  scale_y_continuous(name="Reads") + 
  scale_x_continuous(name="Secuencia")

Zoom_pS

# Juntar los graficos: Original + zoom final
pAbF <- grid.arrange(pS, Zoom_pS, ncol=2, top = "Filtro abundancia liquen Coyhaique")
pAbF

# Guardar
setwd("~/1. TMKV/2. Results/5. Filtro abundancia")
ggsave("Filtro abundancia Liq COY19.png", plot=pAbF,width = 12, height = 6, device= "png",path = "~/1. TMKV/2. Results/5. Filtro abundancia")


# Corte de abundancia
ps.Liq.TAM19.F2<-prune_taxa(taxa_sums(ps.Liq.TAM19.F1)>10,ps.Liq.TAM19.F1)

# Guardar Phyloseq filtrado
saveRDS(ps.Liq.TAM19.F2, "ps.Liq.TAM19.F2.rds")


#########################################
## RECONSTRUCCION DE ARCHIVOS PHYLOSEQ ##
#########################################

## Ejemplo Liquen TAM19 
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19 <- readRDS("ps.Liq.TAM19.rds")

# Tabla de abundancias de cada ASV:
Tabla.Conteos.ASV <- otu_table(ps.Liq.TAM19)

# Tabla de taxonomia asignada a cada ASV:
Tabla.Taxonomia <- tax_table(ps.Liq.TAM19)
colnames(Tabla.Taxonomia) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") 

# Tabla de metadatos
setwd("~/1. TMKV/4. Metadata")
samdf<- read.csv("Metadata_Liq_TAM19.csv", header=TRUE,sep = ";", row.names = 1)
sample_names <- c("Liq-TAM19-003","Liq-TAM19-009","Liq-TAM19-010","Liq-TAM19-022","Liq-TAM19-023",
                  "Liq-TAM19-026","Liq-TAM19-028","Liq-TAM19-033","Liq-TAM19-036","Liq-TAM19-037")
rownames(Tabla.Conteos.ASV)<- sample_names
rownames(samdf)<- sample_names
colnames(samdf) <- c("Tipo_de_Muestra", "Tipo_de_muestra_sitio", "Origen", "Sitio","Nombre_Tesis","Especie","Cianobionte","Mes_Recoleccion","Año_Recoleccion","Latitud","Longitud","Altura", "Ambiente","Distribución_global","Nombre_Terreno","Nombre_GBIF","Numero_orden_secuenciación","Réplica","Número")

# Reconstruccion del objeto phyloseq:
ps.Liq.TAM19 <- phyloseq(otu_table(Tabla.Conteos.ASV, taxa_are_rows=FALSE),
                         tax_table(Tabla.Taxonomia),
                         sample_data(samdf))
rank_names(ps.Liq.TAM19)
view(meta(ps.Liq.TAM19))

# Guardar archivos como rds
setwd("~/1. TMKV/5. Phyloseq")
saveRDS(ps.Liq.TAM19, "ps.Liq.TAM19.rds")

# Guardar ambiente de trabajo 
save.image("~/1. TMKV/0. Scripts&Environment/Reconstruccionphyloseq-ps.Liq.TAM19.RData")


###########################################
## GRAFICO FRECUENCIA LARGO DE SECUENCIA ##
###########################################

# Llamar dataframe para graficos de frecuencia de largo de secuencia
dist_Liq_COY19 <- read_csv("~/1. TMKV/2. Results/3. Dist. largo seq/dist.Liq_COY19.csv")
#View(dist_Liq_COY19)

dist_Sus_COY19 <- read_csv("~/1. TMKV/2. Results/3. Dist. largo seq/dist.Sus_COY19.csv")
#View(dist_Sus_COY19)

dist_Sue_COY19 <- read_csv("~/1. TMKV/2. Results/3. Dist. largo seq/dist.Sue_COY19.csv")
#View(dist_Sue_COY19)

dist_Liq_TAM19 <- read_csv("~/1. TMKV/2. Results/3. Dist. largo seq/dist.Liq_TAM19.csv")
#View(dist_Liq_TAM19)

dist_Sus_TAM19 <- read_csv("~/1. TMKV/2. Results/3. Dist. largo seq/dist.Sus_TAM19.csv")
#View(dist_Sus_TAM19)

dist_Sue_TAM19 <- read_csv("~/1. TMKV/2. Results/3. Dist. largo seq/dist.Sue_TAM19.csv")
#View(dist_Sue_TAM19)

# Construir graficos 
# Se puede reducir 360-380
p.Liq_TAM19.dls <- ggplot(data=dist_Liq_COY19, aes(x=Var1, y=Freq)) + geom_histogram(stat='identity') + 
  scale_x_continuous(breaks = seq(355, 385, 5),limits=c(355,385), name = "Largo secuencia (pb)") + scale_y_continuous(name= "Frecuencia")


################################
## GRAFICO RESUMEN SECUENCIAS ##
################################

setwd("~/1. TMKV/2. Results/4. Resumen Merge+Quimeras")

Resumen_Secuencias<- read_excel("Resumen-Secuencias_DADA2_para grafico.xlsx")

p <- ggplot(Resumen_Secuencias, aes(x=reorder(Analysis,-Sequences), y=Sequences, fill=Type)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#FF9A4A","#FF7044","#FE2D43","#6ECB63","#1F965E","#1E7179")) +
  scale_fill_manual(values=c("#FF9A4A","#FF7044","#FE2D43","#6ECB63","#1F965E","#1E7179")) +
  xlab("Filtros") + ylab("Secuencias") 

p 

#########################
## GRAFICO RARREFACION ##
#########################

setwd("~/1. TMKV/5. Phyloseq")

# cargar los archivos phyloseq
ps.TMKV<-readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.rds")

# Generar una copia del archivo 
lich<-ps.TMKV

# Pegar "ASV" antes del numero de cada secuencia
taxa_names(lich) <- paste0("ASV", seq(ntaxa(lich)))
# Generar una matriz 
OTU1 = as(otu_table(lich), "matrix")
if(taxa_are_rows(lich)){OTU1 <- t(OTU1)}
OTU2 = as.data.frame(OTU1)
# Generar un data frame con la informacion de origen para cada muestra
origin<-as.data.frame(lich@sam_data[["Tipo_de_muestra_sitio"]])
colnames(origin)<-"Origin"
# Juntar dataframes
OTU<-cbind(origin,OTU2)

# Generar dataframe individuales para cada origen

OTULiqCOY<-OTU[OTU$Origin=="Talo P. rufescens COY",]
OTUSusCOY<-OTU[OTU$Origin=="Sustrato P. rufescens COY",]
OTULiqTAM<-OTU[OTU$Origin=="Talo P. rufescens TAM",]
OTUSusTAM<-OTU[OTU$Origin=="Sustrato P. rufescens TAM",]
OTUSueCOY<-OTU[OTU$Origin=="Suelo COY",]
OTUSueTAM<-OTU[OTU$Origin=="Suelo TAM",]

# Generar los nombres con cursiva

LiqCOY <- expression(paste("Talo", italic("P. rufescens"), COY))
SusCOY <- expression(paste("Sustrato ", italic("P. rufescens"), COY))
SueCOY <- expression(paste("Suelo COY"))
LiqTAM <- expression(paste("Talo", italic("P. rufescens"), TAM))
SusTAM <- expression(paste("Sustrato ", italic("P. rufescens"), TAM))
SueTAM <- expression(paste("Suelo TAM"))

# Generar el promedio de cada muestra
OTULiqCOY_mean<-as.data.frame(colMeans(OTULiqCOY[2:ncol(OTULiqCOY)]))
colnames(OTULiqCOY_mean)<-"Liquen P. rufescens COY"

OTUSusCOY_mean<-as.data.frame(colMeans(OTUSusCOY[2:ncol(OTUSusCOY)]))
colnames(OTUSusCOY_mean)<-"Sustrato P. rufescens COY"

OTUSueCOY_mean<-as.data.frame(colMeans(OTUSueCOY[2:ncol(OTUSueCOY)]))
colnames(OTUSueCOY_mean)<-"Suelo COY"

OTULiqTAM_mean<-as.data.frame(colMeans(OTULiqTAM[2:ncol(OTULiqTAM)]))
colnames(OTULiqTAM_mean)<-"Liquen P. rufescens TAM"

OTUSusTAM_mean<-as.data.frame(colMeans(OTUSusTAM[2:ncol(OTUSusTAM)]))
colnames(OTUSusTAM_mean)<-"Sustrato P. rufescens TAM"

OTUSueTAM_mean<-as.data.frame(colMeans(OTUSueTAM[2:ncol(OTUSueTAM)]))
colnames(OTUSueTAM_mean)<-"Suelo TAM"

# Juntar la informacion de los promedios
OTU_ave<-cbind(OTULiqCOY_mean, OTUSusCOY_mean, OTUSueCOY_mean,
               OTULiqTAM_mean, OTUSusTAM_mean, OTUSueTAM_mean)

# Transformar de un argumento numero a un vector numerico
OTU_ave_c<-ceiling(OTU_ave)

acc.todo<-iNEXT(OTU_ave_c, q = 0, datatype = "abundance", nboot=200, se=TRUE, conf=0.95)
acc.todo$DataInfo

# Graficar
Rarefaccion <- ggiNEXT(acc.todo, facet.var="none", color.var="site", grey=FALSE, se=TRUE) +
  scale_color_manual(values=c("#FF9A4A","#6ECB63","#FE2D43","#1E7179","#FF7044","#1F965E")) +
  scale_fill_manual(values=c("#FF9A4A","#6ECB63","#FE2D43","#1E7179","#FF7044","#1F965E")) +
  theme_bw(base_size = 18) + 
  theme(legend.position="left") +
  xlab("Nº de secuencias") + ylab("ASVs")

Rarefaccion

# Guardar
dev.off()
tiff("Rarefaccion.tiff", units="in", width=15, height=12, res=300)
Rarefaccion
dev.off()

ggsave("Rarefaccion.png", width = 12.5, height = 7.3, dpi = 300)

View(acc.todo)
View(acc.todo[["AsyEst"]])
View(acc.todo[["iNextEst"]])

# Guardar
setwd("~/1. TMKV/0. Scripts&Environment/2. Enviroment")
save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Rarefaccion.RData")

