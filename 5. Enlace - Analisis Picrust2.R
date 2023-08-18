###################################
## 5. Enlace - Analisis Picrust2 ##
###################################

##############
## Paquetes ##
##############

library(car)
library(tidyverse)
library(fs)
library(ade4)
library("factoextra")
library(ggrepel)
library(ggplot2)
library(vegan)
library(microbiome)
library(pairwiseAdonis)
library(factoextra)

###########
## Datos ##
###########

#Importar data de contribucion de cada ASV
KO_contrib<-read.delim("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv.gz")

# Revisar archivo
#View(KO_contrib)
colnames(KO_contrib)
levels(as.factor(KO_contrib$function.))
levels(as.factor(KO_contrib$sample))

# Agregar columas
KO_contrib = mutate(KO_contrib,Sampletype_site = case_when(
  startsWith(sample, "Liq-TAM") ~ "LiqTAM",
  startsWith(sample, "Sus-TAM") ~ "SusTAM",
  startsWith(sample, "Sue-TAM") ~ "SueTAM",
  startsWith(sample, "Liq-COY") ~ "LiqCOY",
  startsWith(sample, "Sus-COY") ~ "SusCOY",
  startsWith(sample, "Sue-COY") ~ "SueCOY"))

KO_contrib = mutate(KO_contrib,Enzyme = case_when(
  startsWith(function., "K00117") ~ "gcd/mGDH",
  startsWith(function., "K00937") ~ "ppk1",
  startsWith(function., "K01077") ~ "phoA/phoB",
  startsWith(function., "K01078") ~ "olpA",
  startsWith(function., "K01093") ~ "appA",
  startsWith(function., "K01113") ~ "phoD",
  startsWith(function., "K01114") ~ "plc",
  startsWith(function., "K01126") ~ "ugpQ/glpQ",
  startsWith(function., "K01507") ~ "ppa",
  startsWith(function., "K01524") ~ "ppx",
  startsWith(function., "K02039") ~ "phoU",
  startsWith(function., "K02040") ~ "pstS",
  startsWith(function., "K02043") ~ "phnF",
  startsWith(function., "K02044") ~ "phnD",
  startsWith(function., "K02445") ~ "GlpT",
  startsWith(function., "K03306") ~ "pit",
  startsWith(function., "K03430") ~ "phnW",
  startsWith(function., "K05306") ~ "phnX",
  startsWith(function., "K05774") ~ "phnN",
  startsWith(function., "K05813") ~ "ugpB",
  startsWith(function., "K06163") ~ "phnJ",
  startsWith(function., "K06167") ~ "phnP",
  startsWith(function., "K06193") ~ "phnA",
  startsWith(function., "K07048") ~ "opd",
  startsWith(function., "K07093") ~ "phoX",
  startsWith(function., "K07636") ~ "phoR",
  startsWith(function., "K07657") ~ "phoB",
  startsWith(function., "K07658") ~ "PhoP/PhoB1",
  startsWith(function., "K09474") ~ "phoN",
  startsWith(function., "K11081") ~ "phnS",
  startsWith(function., "K03788") ~ "aphA",
  startsWith(function., "K09994") ~ "PhnO",
  startsWith(function., "K11751") ~ "UshA",
  startsWith(function., "K19669") ~ "pphA/Pal",
))

# Definir niveles
KO_contrib$sample <- factor(KO_contrib$sample,levels = c("Liq-COY19-002","Liq-COY19-008","Liq-COY19-016","Liq-COY19-019","Liq-COY19-022",
                                                         "Liq-COY19-023","Liq-COY19-027","Liq-COY19-028","Liq-COY19-030","Liq-COY19-031",
                                                         "Liq-TAM19-003","Liq-TAM19-009","Liq-TAM19-010","Liq-TAM19-022","Liq-TAM19-023",
                                                         "Liq-TAM19-026","Liq-TAM19-028","Liq-TAM19-033","Liq-TAM19-036","Liq-TAM19-037",
                                                         "Sus-COY19-002","Sus-COY19-008","Sus-COY19-016","Sus-COY19-019","Sus-COY19-022",
                                                         "Sus-COY19-023","Sus-COY19-027","Sus-COY19-028","Sus-COY19-030","Sus-COY19-031",
                                                         "Sus-TAM19-003","Sus-TAM19-009","Sus-TAM19-010","Sus-TAM19-022","Sus-TAM19-023",
                                                         "Sus-TAM19-026","Sus-TAM19-028","Sus-TAM19-033","Sus-TAM19-036","Sus-TAM19-037",
                                                         "Sue-COY19-01","Sue-COY19-02","Sue-COY19-03","Sue-COY19-04",
                                                         "Sue-COY19-05","Sue-COY19-06","Sue-COY19-07","Sue-COY19-08",
                                                         "Sue-TAM19-001","Sue-TAM19-002","Sue-TAM19-003","Sue-TAM19-004","Sue-TAM19-005",
                                                         "Sue-TAM19-007","Sue-TAM19-008","Sue-TAM19-009","Sue-TAM19-011","Sue-TAM19-014"))


KO_contrib$Enzyme <- factor(KO_contrib$Enzyme,levels = c("ppa","gcd/mGDH","appA","phoN","aphA","olpA",
                                                         "phoX","phoA/phoB","phoD","ushA","phnW","phnX",
                                                         "phnA","pphA/Pal","phnF","phnJ","phnN","phnO","phnP",
                                                         "ugpQ/glpQ","plc","opd","ppx","ppk1",
                                                         "pit","phoB","phoR","phoU","PhoP/PhoB1",
                                                         "pstS","phnD","phnS","ugpB","GlpT"))

# Se revisaron 34 enzimas, 4 fueron descartadas porque su detección fue insignificante 

# Revisamos la tabla
#View(KO_contrib)

KO_contrib=subset(KO_contrib, Enzyme != "ushA")
KO_contrib=subset(KO_contrib, Enzyme != "aphA")
KO_contrib=subset(KO_contrib, Enzyme != "pphA/Pal")
KO_contrib=subset(KO_contrib, Enzyme != "phnO")


# Importar data de taxonomia
Tax_contrib<-read.delim("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline/taxonomia.tsv")
#View(Tax_contrib)
# Cambiar el nombre de columnas para hacer merge 
colnames(Tax_contrib)=c("taxon","Domain","Phylum","Class","Order","Family","Genus","Species")

# Unimos
merge_contrib=merge(KO_contrib,Tax_contrib)
#View(merge_contrib)

# Creamos dataframe para colores en graficos 
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
#View(color)

##########################################
## Graficos enzimas por tipo de muestra ##
##########################################

##############
## Heat map ##
##############

#Reestructurar dataframe

Heatmap_AB <- merge_contrib %>%
  group_by(Enzyme, sample) %>%
  summarise(
    Ab = sum(taxon_rel_abun)
  )
# Revisar
#View(Heatmap_AB)

# Definir niveles
Heatmap_AB$sample <- factor(Heatmap_AB$sample,levels = c("Liq-COY19-002","Liq-COY19-008","Liq-COY19-016","Liq-COY19-019","Liq-COY19-022",
                                                         "Liq-COY19-023","Liq-COY19-027","Liq-COY19-028","Liq-COY19-030","Liq-COY19-031",
                                                         "Liq-TAM19-003","Liq-TAM19-009","Liq-TAM19-010","Liq-TAM19-022","Liq-TAM19-023",
                                                         "Liq-TAM19-026","Liq-TAM19-028","Liq-TAM19-033","Liq-TAM19-036","Liq-TAM19-037",
                                                         "Sus-COY19-002","Sus-COY19-008","Sus-COY19-016","Sus-COY19-019","Sus-COY19-022",
                                                         "Sus-COY19-023","Sus-COY19-027","Sus-COY19-028","Sus-COY19-030","Sus-COY19-031",
                                                         "Sus-TAM19-003","Sus-TAM19-009","Sus-TAM19-010","Sus-TAM19-022","Sus-TAM19-023",
                                                         "Sus-TAM19-026","Sus-TAM19-028","Sus-TAM19-033","Sus-TAM19-036","Sus-TAM19-037",
                                                         "Sue-COY19-01","Sue-COY19-02","Sue-COY19-03","Sue-COY19-04",
                                                         "Sue-COY19-05","Sue-COY19-06","Sue-COY19-07","Sue-COY19-08",
                                                         "Sue-TAM19-001","Sue-TAM19-002","Sue-TAM19-003","Sue-TAM19-004","Sue-TAM19-005",
                                                         "Sue-TAM19-007","Sue-TAM19-008","Sue-TAM19-009","Sue-TAM19-011","Sue-TAM19-014"))


Heatmap_AB$Enzyme <- factor(Heatmap_AB$Enzyme,levels = c("ppa","gcd/mGDH","appA","phoN","olpA",
                                                         "phoX","phoA/phoB","phoD","phnW","phnX",
                                                         "phnA","phnF","phnJ","phnN","phnP","ugpQ/glpQ",
                                                         "plc","opd","ppx","ppk1","pit","phoB","phoR","phoU",
                                                         "PhoP/PhoB1","pstS","phnD","phnS","ugpB","GlpT"))

# Graficar
Heatmap_ES = ggplot(Heatmap_AB, aes(x = sample, y = Enzyme , fill=Ab)) + 
  geom_tile() + ylim(rev(levels(Heatmap_AB$Enzyme))) + scale_fill_gradient(low = "white", high = "#ad0303")+
  theme(axis.text.x = element_text(angle = 90))

Heatmap_ES

# Guardar

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
#png("Heatmap_ES.#png", width = 1000 ,height = 800)
Heatmap_ES
dev.off()

#####################################################
## Diferencia significativa - Analisis estadistico ##
#####################################################

# Adaptación dataframe

# dataframe a adaptar
#View(Heatmap_AB)

Dif_AB= spread(data = Heatmap_AB, key = sample, value = Ab)
# Revisar
#View(Dif_AB)

Dif_AB=as.data.frame(Dif_AB[c(1:30),])


Dif_AB = column_to_rownames(Dif_AB,"Enzyme")
Dif_AB=as.data.frame(t(Dif_AB))

#revisar
head(Dif_AB)
nrow(Dif_AB)

Dif_AB=mutate(Dif_AB,Tipo=c(rep("LiqCOY",10),rep("LiqTAM",10),
                            rep("SusCOY",10),rep("SusTAM",10),
                            rep("SueCOY",8),rep("SueTAM",10)))
# Revisar
head(Dif_AB)
tail(Dif_AB)
View(Dif_AB)

# Edicion de nombres de columna (problemas con "/")
Dif_AB_f=Dif_AB
colnames(Dif_AB)

colnames(Dif_AB_f) = c("ppa","gcdmGDH","appA","phoN","olpA","phoX",
                       "phoAphoB","phoD","phnW","phnX","phnA","phnF",
                       "phnJ","phnN","phnP","ugpQglpQ","plc","opd",
                       "ppx","ppk1","pit","phoB","phoR","phoU",
                       "PhoPPhoB1","pstS","phnD","phnS","ugpB","GlpT","Tipo")

#revisar
#View(Dif_AB_f)
as.factor(Dif_AB_f$Tipo)

##########################################
## Revisar parametrico o no parametrico ##
##########################################

# Normalidad 


for (i in seq(1:30)) {
  ST <- shapiro.test(Dif_AB_f[,i])
  print(ST)
}

# Siendo la hipótesis nula que la población está distribuida normalmente, si el p-valor es menor a alfa (nivel de significancia) entonces la hipótesis nula es rechazada (se concluye que los datos no vienen de una distribución normal).
# Resultados: no normales

# Levene test

for (i in seq(1:30)) {
  a <- colnames(Dif_AB_f)
  LT <- leveneTest(data=Dif_AB_f, as.formula(paste(a[i], a[31], sep="~")))
  print(LT)
}

# si 0.05 > Pr, NO homogeneos (varianza dif entre al menos 2 grupos). Resultado: solo 3 homogeneas. 

# Datos no parametricos, entonces se usa el test  Kruskal-Wallis 

#################
## Comparación ##
#################

# Kruskal-Wallis 

for (i in seq(1:30)) {
  a <- colnames(Dif_AB_f)
  KT <- kruskal.test(data=Dif_AB_f, as.formula(paste(a[i], a[31], sep="~")))
  KT=print(KT)
}

# As the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the treatment groups.

# Comparación por pares: Wilcox test

for (i in seq(1:30)) {
  a <- colnames(Dif_AB_f)
  WT <- pairwise.wilcox.test(Dif_AB_f[,i], Dif_AB_f$Tipo,
                             p.adjust.method = "BH")$p.value
  
  df <- data.frame(expand.grid(dimnames(WT)),array(WT))
  na.omit(df)
  
  colnames(df)=c("V1","V2",paste(a[i], "WT", sep="-"))
  df= mutate(df,comparation=paste(df$V1,df$V2, sep = " vs "))
  print(df)
  setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/WT")
  write.table(df, file = paste(a[i], "WT.txt", sep="-"), sep = "\t",
              row.names = TRUE, col.names = NA)
}

#########################################
## Manipulacion de datos para graficar ##
#########################################

files <- dir_ls("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/WT/", regexp = ".txt")  
files
WT.all = map(files, read.table) %>%
  reduce(inner_join, by = "V5") 
#View(WT.all)

colnames(WT.all)=WT.all[1,]

WT.df=WT.all[,c("comparation","ppa-WT","gcdmGDH-WT","appA-WT","phoN-WT","olpA-WT","phoX-WT",
                "phoAphoB-WT","phoD-WT","phnW-WT","phnX-WT","phnA-WT","phnF-WT",
                "phnJ-WT","phnN-WT","phnP-WT","ugpQglpQ-WT","plc-WT","opd-WT",
                "ppx-WT","ppk1-WT","pit-WT","phoB-WT","phoR-WT","phoU-WT",
                "PhoPPhoB1-WT","pstS-WT","phnD-WT","phnS-WT","ugpB-WT","GlpT-WT")]

colnames(WT.df)=WT.df[1,]
#WT.df=WT.df[c(2:6,8:11,14:16,20:21,26),]
WT.df=WT.df[c(2,26,14,5,3,11,9,15,21),]
#View(WT.df)

long.WT.df=gather(WT.df, key="Enzima",value="Pvalue",2:31)
#View(long.WT.df)

long.WT.df$Pvalue=as.numeric(long.WT.df$Pvalue)
class(long.WT.df$Pvalue)

long.WT.PV = mutate(long.WT.df, Significance =
                      case_when(Pvalue > 0.05 ~ "ns",
                                Pvalue <= 0.05 & Pvalue > 0.005 ~ "*",
                                Pvalue <= 0.005 & Pvalue > 0.0005~ "**",
                                Pvalue <= 0.0005 & Pvalue > 0.00005~ "***",
                                Pvalue <= 0.00005 ~ "****"))
#View(long.WT.PV)

long.WT.PV$Enzima=as.factor(long.WT.PV$Enzima)

long.WT.PV$Enzima<- factor(long.WT.PV$Enzima,levels = c("ppa-WT","gcdmGDH-WT","appA-WT","phoN-WT","olpA-WT",
                                                        "phoX-WT","phoAphoB-WT","phoD-WT","phnW-WT","phnX-WT",
                                                        "phnA-WT","phnF-WT","phnJ-WT","phnN-WT","phnP-WT",
                                                        "ugpQglpQ-WT","plc-WT","opd-WT","ppx-WT","ppk1-WT",
                                                        "pit-WT","phoB-WT","phoR-WT","phoU-WT","PhoPPhoB1-WT",
                                                        "pstS-WT","phnD-WT","phnS-WT","ugpB-WT","GlpT-WT"))

long.WT.PV$comparation=as.factor(long.WT.PV$comparation)
levels(long.WT.PV$comparation)

long.WT.PV$comparation<- factor(long.WT.PV$comparation,levels = c("LiqTAM vs LiqCOY",
                                                                  "SusTAM vs SusCOY",
                                                                  "SueTAM vs SueCOY",
                                                                  "SusCOY vs LiqCOY",
                                                                  "SueCOY vs LiqCOY",
                                                                  "SusTAM vs LiqTAM",
                                                                  "SueTAM vs LiqTAM",
                                                                  "SusCOY vs SueCOY",
                                                                  "SusTAM vs SueTAM"))
##########################################
## Graficar Heatmap dif. significativas ##
##########################################

Heatmap_WT= ggplot(long.WT.PV, aes(x = comparation, y = Enzima , fill=Significance)) + 
  geom_tile() + ylim(rev(levels(long.WT.PV$Enzima)))  + 
  scale_fill_manual(values=c("#ffda00","#ffcd01","#ffba01","#ffa701","white"))+
  theme(axis.text.x = element_text(angle = 90))

# Guardar

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
#png("Heatmap_WT.#png", width = 400 ,height = 800)
Heatmap_WT
dev.off()

#########
## PCA ## 
#########
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

View(Dif_AB)
as.factor(Dif_AB$Tipo)

Dif_AB$Tipo<- factor(Dif_AB$Tipo,levels = c("LiqCOY","LiqTAM","SusCOY","SusTAM","SueCOY","SueTAM"))



pca <- dudi.pca(Dif_AB[,1:30], scannf = F, nf = 5)

# Análisis de contribución de variables

# Revisar
class(pca)
print(pca)

# info variables 

# Contribución en cada dimención
pca_var=get_pca_var(pca)
pca_var_con=pca_var$contrib
pca_var_con_D1D2= pca_var_con[,1:2]

colnames(pca_var_con_D1D2)=c("D1","D2")
pca_var_con_D1D2=as.data.frame(pca_var_con_D1D2)
pca_var_con_D1order <- pca_var_con_D1D2[order(-pca_var_con_D1D2$D1),]
pca_var_con_D2order <- pca_var_con_D1D2[order(-pca_var_con_D1D2$D2),]

# Contribución de cada variable 
sumarize_pca=facto_summarize(pca,"var")
sumarize_pca <- sumarize_pca[order(-sumarize_pca$contrib),]
View(sumarize_pca)



setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
write.table(pca_var_con_D1order, file = "pca_gen_var_con_D1order.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(pca_var_con_D2order, file = "pca_gen_var_con_D2order.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(sumarize_pca, file = "sumarize_pca_gen.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
# Grafico
PCA <- fviz_pca_biplot(pca, 
                       geom.ind = "point",
                       fill.ind = Dif_AB$Tipo,
                       col.ind = Dif_AB$Tipo,
                       pointshape = 21, pointsize = 2,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       repel = TRUE,
                       geom = c("text", "arrow"),
                       col.var = "Grey20") + 
  scale_color_manual(values=c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179")) + 
  scale_fill_manual(values= c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179")) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

# Guardar

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
#png("PCA_ES.#png", width = 1500 ,height = 850)
PCA
dev.off()

##############
# Permanova ##
##############

ps.TMKV<- readRDS("~/1. TMKV/5. Phyloseq/ps.TMKV.rds") #F2
meta(ps.TMKV)
str(Dif_AB)
dist_EP <- vegdist((Dif_AB[,1:30]), method="bray")

permanova_EP <- adonis2((Dif_AB[,1:30]) ~ Tipo, data = (Dif_AB), permutations = 999, method="bray")

pwadonis_EP = pairwise.adonis(dist_EP, meta(ps.TMKV)$Tipo_de_muestra_sitio)
pwadonis.dat_EP=as.data.frame(pwadonis_EP)
#View(pwadonis.dat_EP)

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
write.table(pwadonis.dat_EP, file = "PERMANOVA-EnzimasP.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#########################################
## Graficos enzimas por tipo taxonomia ##
#########################################

# Creamos dataframe
bar_AB <- merge_contrib %>%
  group_by(Sampletype_site, Phylum, Enzyme) %>%
  summarise(
    Ab = sum(taxon_rel_abun)
  )

#View(bar_AB)

bar_AR = mutate(bar_AB, AR= case_when(Sampletype_site == "LiqCOY"  ~ Ab/10,
                                      Sampletype_site == "LiqTAM"  ~ Ab/10,
                                      Sampletype_site == "SusCOY"  ~ Ab/10,
                                      Sampletype_site == "SusTAM"  ~ Ab/10,
                                      Sampletype_site == "SueCOY"  ~ Ab/8,
                                      Sampletype_site == "SueTAM"  ~ Ab/10))

levels(as.factor(bar_AR$Sampletype_site))


bar_AR$Sampletype_site <- factor(bar_AR$Sampletype_site,levels = c("LiqCOY","LiqTAM","SusCOY","SusTAM","SueCOY","SueTAM"))



#View(bar_AR)

# Grafico filos
#png("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/taxcontribution.#png", width = 2000, height = 1000)
ggplot(bar_AR, aes(x = Sampletype_site, y = AR, fill = Phylum))  + 
  geom_bar(stat = 'identity') + facet_grid(~ Enzyme) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = color$cod, limits=color$Filo)
dev.off()

# Creamos dataframe
bar_AB_G <- merge_contrib %>%
  group_by(Sampletype_site, Genus, Enzyme) %>%
  summarise(
    Ab = sum(taxon_rel_abun)
  )

levels(as.factor(bar_AB_G$Genus))

bar_AR_G = mutate(bar_AB_G, AR= case_when(Sampletype_site == "LiqCOY"  ~ Ab/10,
                                          Sampletype_site == "LiqTAM"  ~ Ab/10,
                                          Sampletype_site == "SusCOY"  ~ Ab/10,
                                          Sampletype_site == "SusTAM"  ~ Ab/10,
                                          Sampletype_site == "SueCOY"  ~ Ab/8,
                                          Sampletype_site == "SueTAM"  ~ Ab/10))


bar_AR_G$Sampletype_site <- factor(bar_AR_G$Sampletype_site,levels = c("LiqCOY","LiqTAM","SusCOY","SusTAM","SueCOY","SueTAM"))



# Grafico para generos: para mostrar NA

#png("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/taxcontributionGenus.#png", width = 2000, height = 1000)
ggplot(bar_AR_G, aes(x = Sampletype_site, y = AR, fill = Genus)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Enzyme) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90))
dev.off()

#######################
## PCA Para Gen-Filo ##
#######################

# Adaptacion dataframe

#View(merge_contrib)

# Diseño dataframe
colnames(merge_contrib)
merge_contrib_PCA=mutate(merge_contrib,EP=paste(Enzyme, Phylum, sep=" - "))
#View(merge_contrib_PCA)

merge_contrib_PCA =merge_contrib_PCA[,c("EP","sample","taxon_rel_abun")]

merge_contrib_PCA_comp <- merge_contrib_PCA %>%
  group_by(EP, sample) %>%
  summarise(
    Ab = sum(taxon_rel_abun)
  )
#View(merge_contrib_PCA_comp)

contrib_PCA_comp = spread(merge_contrib_PCA_comp,key = sample,value = Ab)
#View(contrib_PCA_comp)


contrib_PCA_comp[is.na(contrib_PCA_comp)] <- 0
contrib_PCA_comp_df=as.data.frame(contrib_PCA_comp)

contrib_PCA_comp_T=t(contrib_PCA_comp_df)
#View(contrib_PCA_comp_T)

colnames(contrib_PCA_comp_T)=contrib_PCA_comp_T[1,]

contrib_PCA_comp_T=contrib_PCA_comp_T[2:59,]

contrib_PCA_comp_T=as.data.frame(contrib_PCA_comp_T)

contrib_PCA_comp_T=mutate(contrib_PCA_comp_T,Tipo=c(rep("LiqCOY",10),rep("LiqTAM",10),
                                                    rep("SusCOY",10),rep("SusTAM",10),
                                                    rep("SueCOY",8),rep("SueTAM",10)))

# PCA
# contrib_PCA_comp_T con tipo
# contrib_PCA_comp_df para PCA

as.factor(contrib_PCA_comp_T$Tipo)

contrib_PCA_comp_T$Tipo<- factor(contrib_PCA_comp_T$Tipo,levels = c("LiqCOY","LiqTAM","SusCOY","SusTAM","SueCOY","SueTAM"))

#View(contrib_PCA_comp_df)
rownames(contrib_PCA_comp_df)=c(contrib_PCA_comp_df$EP)

contrib_PCA_comp_df=contrib_PCA_comp_df[,2:59]

# PCA con texto variables principales

options(ggrepel.max.overlaps = Inf)

pca2 <- dudi.pca(t(contrib_PCA_comp_df), scannf = F, nf = 5)


# info variables 

pca_var2=get_pca_var(pca2)
pca_var_con2=pca_var2$contrib
pca_var_con_D1D2_2= pca_var_con2[,1:2]

colnames(pca_var_con_D1D2_2)=c("D1","D2")
pca_var_con_D1D2_2=as.data.frame(pca_var_con_D1D2_2)
pca_var_con_D1order_2 <- pca_var_con_D1D2_2[order(-pca_var_con_D1D2_2$D1),]
pca_var_con_D2order_2 <- pca_var_con_D1D2_2[order(-pca_var_con_D1D2_2$D2),]

sumarize_pca2=facto_summarize(pca2,"var")
sumarize_pca2 <- sumarize_pca2[order(-sumarize_pca2$contrib),]
sumarize_pca2_60=sumarize_pca2[1:60,]
View(sumarize_pca2_60)


setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
write.table(pca_var_con_D1order_2, file = "pca_genfilo_var_con_D1order.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(pca_var_con_D2order_2, file = "pca_genfilo_var_con_D2order.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

write.table(sumarize_pca2_60, file = "sumarize_pca2_60_genfilo.txt", sep = "\t",
            row.names = TRUE, col.names = NA)



PCA <- fviz_pca_biplot(pca2,geom.var = c("point","text"),
                       select.var = list(contrib = 60),
                       geom.ind = "point",
                       fill.ind = contrib_PCA_comp_T$Tipo,
                       col.ind = contrib_PCA_comp_T$Tipo,
                       pointshape = 21, pointsize = 2,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       repel = TRUE,
                       col.var = "Grey20")+ 
  scale_color_manual(values=c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179")) + 
  scale_fill_manual(values= c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))+
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

# Guardar

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
#png("PCA_EF_60.#png", width = 1500 ,height = 850)
PCA
dev.off()

# PCA con texto "todas"

options(ggrepel.max.overlaps = 10)
pca2 <- dudi.pca(t(contrib_PCA_comp_df), scannf = F, nf = 5)


PCA <- fviz_pca_biplot(pca2,geom.var = c("point","text"),
                       geom.ind = "point",
                       fill.ind = contrib_PCA_comp_T$Tipo,
                       col.ind = contrib_PCA_comp_T$Tipo,
                       pointshape = 21, pointsize = 2,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       repel = TRUE,
                       col.var = "Grey20")+ 
  scale_color_manual(values=c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179")) + 
  scale_fill_manual(values= c("#FF9A4A","#6ECB63","#FF7044","#1F965E","#FE2D43","#1E7179"))

# Guardar

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
##png("PCA_EF_todo.##png", width = 1500 ,height = 850)
PCA
dev.off()


# Permanova PCA filos

#View(contrib_PCA_comp_T)
class(contrib_PCA_comp_T)
class(contrib_PCA_comp_T$`appA - Acidobacteriota`)

contrib_PCA_comp_Perm=contrib_PCA_comp_T[,1:415]
contrib_PCA_comp_Perm_num <- as.data.frame(sapply(contrib_PCA_comp_Perm, as.numeric))

contrib_PCA_comp_Perm_num_t=as.data.frame(t(contrib_PCA_comp_Perm_num))



str(contrib_PCA_comp_Perm_num_t)


permanova_EPF <- adonis2(contrib_PCA_comp_Perm_num_t ~ Tipo, data = (contrib_PCA_comp_T), permutations = 999, method="bray")
dist_EPF <- vegdist((contrib_PCA_comp_Perm_num_t), method="bray")


pwadonis_EPF = pairwise.adonis(dist_EPF, meta(ps.TMKV)$Tipo_de_muestra_sitio)
pwadonis.dat_EPF=as.data.frame(pwadonis_EPF)
View(pwadonis.dat_EPF)

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
write.table(pwadonis.dat_EPF, file = "PERMANOVA-EnzimasP-Filo.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


###############################################
## Heat map para Gen-Filo-Generos(Conocidos) ##
###############################################

# Adaptacion dataframe

#View(merge_contrib)

# Diseño dataframe

#Revisar colnames
colnames(merge_contrib)
#Crear nueva columna con 
merge_contrib_PG=mutate(merge_contrib,EP=paste(Phylum, Genus, Enzyme, sep=" - "))
#View(merge_contrib_PG)

merge_contrib_PG$Enzyme=as.character(merge_contrib_PG$Enzyme)
# Problemas con /, entonces creamos columna nueva
merge_contrib_PG = mutate(merge_contrib_PG,Enzimas = case_when(
  startsWith(Enzyme, "gcd/mGDH") ~ "gcdmGDH",
  startsWith(Enzyme, "ppk1") ~ "ppk1",
  startsWith(Enzyme, "phoA/phoB") ~ "phoAphoB",
  startsWith(Enzyme, "olpA") ~ "olpA",
  startsWith(Enzyme, "appA") ~ "appA",
  startsWith(Enzyme, "phoD") ~ "phoD",
  startsWith(Enzyme, "plc") ~ "plc",
  startsWith(Enzyme, "ugpQ/glpQ") ~ "ugpQglpQ",
  startsWith(Enzyme, "ppa") ~ "ppa",
  startsWith(Enzyme, "ppx") ~ "ppx",
  startsWith(Enzyme, "phoU") ~ "phoU",
  startsWith(Enzyme, "pstS") ~ "pstS",
  startsWith(Enzyme, "phnF") ~ "phnF",
  startsWith(Enzyme, "phnD") ~ "phnD",
  startsWith(Enzyme, "GlpT") ~ "GlpT",
  startsWith(Enzyme, "pit") ~ "pit",
  startsWith(Enzyme, "phnW") ~ "phnW",
  startsWith(Enzyme, "phnX") ~ "phnX",
  startsWith(Enzyme, "phnN") ~ "phnN",
  startsWith(Enzyme, "ugpB") ~ "ugpB",
  startsWith(Enzyme, "phnJ") ~ "phnJ",
  startsWith(Enzyme, "phnP") ~ "phnP",
  startsWith(Enzyme, "phnA") ~ "phnA",
  startsWith(Enzyme, "opd") ~ "opd",
  startsWith(Enzyme, "phoX") ~ "phoX",
  startsWith(Enzyme, "phoR") ~ "phoR",
  startsWith(Enzyme, "phoB") ~ "phoB",
  startsWith(Enzyme, "PhoP/PhoB1") ~ "PhoPPhoB1",
  startsWith(Enzyme, "phoN") ~ "phoN",
  startsWith(Enzyme, "phnS") ~ "phnS",
))
#View(merge_contrib_PG)

# Sacar NA 
merge_contrib_GsNA=subset(merge_contrib_PG, Genus != "NA")

Enzymes <- c("ppa","gcdmGDH","appA","phoN","olpA",
             "phoX","phoAphoB","phoD","phnW","phnX",
             "phnA","phnF","phnJ","phnN","phnP",
             "ugpQglpQ","plc","opd","ppx","ppk1",
             "pit","phoB","phoR","phoU","PhoPPhoB1",
             "pstS","phnD","phnS","ugpB","GlpT")


# Loop para obtener los 10 generos más abundantes por enzima segun el promedio en todas las muestras. 

for (i in Enzymes) {
  
  # Seleccionar enzima
  merge_contrib_GsNA_ppa=subset(merge_contrib_GsNA, Enzimas == i)
  #View(merge_contrib_GsNA_ppa)
  
  
  #Seleccionar las 3 columnas de interés
  merge_contrib_GsNA_ppa_3c=merge_contrib_GsNA_ppa[,c("EP","sample","taxon_rel_abun")]
  # Renombrar las columnas
  colnames(merge_contrib_GsNA_ppa_3c)=c("EP","sample","Ab")
  #View(merge_contrib_GsNA_ppa_3c)
  
  # Agrupar por ID y muestra, y sumar las abundancias relativas
  merge_contrib_GsNA_ppa_3c_group <- merge_contrib_GsNA_ppa_3c %>%
    group_by(EP, sample) %>%
    summarise(
      Ab = sum(Ab)
    )
  
  #View(merge_contrib_GsNA_ppa_3c_group)
  
  
  # Cambiar estructura de los datos
  Contrib_ppa_order_sp = spread(merge_contrib_GsNA_ppa_3c_group,key = sample,value = Ab)
  #View(Contrib_ppa_order_sp)
  
  # Todo NA = 0 
  Contrib_ppa_order_sp[is.na(Contrib_ppa_order_sp)] <- 0
  # nombrar filas
  Contrib_ppa_order_sp=column_to_rownames(Contrib_ppa_order_sp, "EP")
  #Calcular el promedio
  Contrib_ppa_order_sp$mean <- rowMeans(Contrib_ppa_order_sp)
  #colnames(Contrib_ppa_order_sp)
  #View(Contrib_ppa_order_sp)
  # ordenamos por promedio
  Contrib_ppa_order_sp_meanorder=Contrib_ppa_order_sp[order(Contrib_ppa_order_sp$mean, decreasing = TRUE),]
  #View(Contrib_ppa_order_sp_meanorder)
  #Contrib_ppa_order_sp_meanorder$mean
  
  # Seleccionamos 10 más abundantes
  Contrib_ppa_sp_meanorder=Contrib_ppa_order_sp_meanorder[c(1:2),]
  
  
  #View(Contrib_ppa_sp_meanorder)
  
  setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/10 Genus + Ab x Gen")
  write.csv(Contrib_ppa_sp_meanorder, file = paste(i, "10G.csv", sep="_"),
            row.names = TRUE, col.names = TRUE,sep = "_")
  
}

# Reestructurar GlpT_10G porque 2 muestras desaparecen dps del spread (Liq-COY19-022 y Liq-COY19-028)

#GlpT_10G=read.csv("GlpT_10G.csv")
#View(GlpT_10G)

#GlpT_10G = mutate(GlpT_10G,"Liq-COY19-022"=rep(0,2))
#GlpT_10G = mutate(GlpT_10G,"Liq-COY19-028"=rep(0,2))
#GlpT_10G=GlpT_10G[,c(1,2,3,4,5,59,6,7,60,8:58)]
#colnames(GlpT_10G)
#GlpT_10G[is.na(GlpT_10G)] <- 0

#GlpT_10G=column_to_rownames(GlpT_10G,var="X")

#setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/10 Genus + Ab x Gen")
#write.csv(GlpT_10G, file = paste("GlpT", "10G.csv", sep="_"),
 #         row.names = TRUE, col.names = TRUE,sep = "_")

# Unir data.frames

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO/10 Genus + Ab x Gen")

files  <- list.files(pattern = '\\.csv')
tables <- lapply(files, read.csv, header = TRUE)
#dfs <- lapply(tables, as.data.frame, header = TRUE)
combined.df <- do.call(rbind , tables)
View(combined.df)


combined.df_g=gather(combined.df, key="Sample",value="Ab",2:60)
View(combined.df_g)
x=c(combined.df$X)

combined.df_g$X=as.factor(combined.df_g$X)
levels(combined.df_g$X)
combined.df_g$X <- factor(combined.df_g$X,levels = c ("Actinobacteriota - Actinoplanes - ppa","Acidobacteriota - Gp6-AA56 - ppa",
                                                      "Acidobacteriota - Gp6-AA56 - gcd/mGDH","Acidobacteriota - Terriglobus - gcd/mGDH",
                                                      "Acidobacteriota - Bog-209 - appA","Proteobacteria - Reyranella - appA",  
                                                      "Proteobacteria - Phenylobacterium - phoN","Proteobacteria - Methylobacterium - phoN",
                                                      "Actinobacteriota - Mycobacterium - olpA","Proteobacteria - Reyranella - olpA",
                                                      "Actinobacteriota - Actinoplanes - phoX","Proteobacteria - Reyranella - phoX",
                                                      "Proteobacteria - LB1R16 - phoA/phoB","Actinobacteriota - Solirubrobacter - phoA/phoB",
                                                      "Actinobacteriota - Solirubrobacter - phoD","Acidobacteriota - Gp6-AA56 - phoD",
                                                      "Proteobacteria - Caballeronia - phnW","Proteobacteria - Reyranella - phnW",
                                                      "Proteobacteria - Caulobacter - phnX","Proteobacteria - PMMR1 - phnX",
                                                      "Proteobacteria - Reyranella - phnA","Proteobacteria - LB1R16 - phnA" ,
                                                      "Proteobacteria - Z2-YC6860 - phnF","Proteobacteria - Mesorhizobium - phnF",
                                                      "Proteobacteria - Reyranella - phnJ","Proteobacteria - Mesorhizobium - phnJ",
                                                      "Proteobacteria - Mesorhizobium - phnN","Proteobacteria - Reyranella - phnN",
                                                      "Acidobacteriota - Terriglobus - phnP","Proteobacteria - Reyranella - phnP" ,
                                                      "Acidobacteriota - Gp6-AA56 - ugpQ/glpQ","Proteobacteria - Cronobacter - ugpQ/glpQ",
                                                      "Actinobacteriota - Actinoplanes - plc","Actinobacteriota - Solirubrobacter - plc",
                                                      "Actinobacteriota - Pseudonocardia - opd","Proteobacteria - Reyranella - opd",
                                                      "Actinobacteriota - Solirubrobacter - ppx","Acidobacteriota - Gp6-AA56 - ppx", 
                                                      "Actinobacteriota - Solirubrobacter - ppk1","Actinobacteriota - Actinoplanes - ppk1",
                                                      "Acidobacteriota - Gp6-AA56 - pit","Proteobacteria - Cronobacter - pit",
                                                      "Actinobacteriota - Solirubrobacter - phoB","Acidobacteriota - Gp6-AA56 - phoB", 
                                                      "Acidobacteriota - Gp6-AA56 - phoR","Proteobacteria - Cronobacter - phoR",
                                                      "Proteobacteria - Cronobacter - phoU","Acidobacteriota - Gp6-AA56 - phoU",
                                                      "Acidobacteriota - Gp6-AA56 - PhoP/PhoB1","Proteobacteria - Cronobacter - PhoP/PhoB1",
                                                      "Proteobacteria - Cronobacter - pstS","Acidobacteriota - Gp6-AA56 - pstS",
                                                      "Proteobacteria - Reyranella - phnD","Proteobacteria - Phenylobacterium - phnD",
                                                      "Actinobacteriota - Streptomyces - phnS","Proteobacteria - Caballeronia - phnS",
                                                      "Proteobacteria - Reyranella - ugpB","Proteobacteria - Mesorhizobium - ugpB",
                                                      "Myxococcota - Labilithrix - GlpT","Proteobacteria - Bradyrhizobium - GlpT"))                 




# Graficar
Heatmap_10G= ggplot(combined.df_g, aes(x = Sample, y = X , fill=Ab)) + 
  geom_tile()+ scale_fill_gradient(low = "white", high = "black") + ylim(rev(levels(combined.df_g$X)))

Heatmap_10G


setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/Graficos Picrust2 TMKV TODO")
##png("Heatmap_10G.#png", width = 1000 ,height = 800)
Heatmap_10G
dev.off()

# Se recomienda guardar el espacio de trabajo cada tanto. ##### Cambiar segun espacio de trabajo ##### 
##save.image("~/1. TMKV/0. Scripts&Environment/Analisis Picrust2_TMKV.RData")    ############################################


############################
## Radarchart taxa claves ##
############################

# Data predicciones
#View(KO_contrib)

# Data predicciones solo LiqCOY o solo LiqTAM
KO_contrib_LiqCOY=subset(KO_contrib, Sampletype_site == "LiqCOY")

KO_contrib_LiqCOY_sum <- KO_contrib_LiqCOY %>%
  group_by(taxon, Enzyme) %>%
  summarise(
    GFC = sum(genome_function_count)
  )

KO_contrib_LiqCOY_sumsd= spread(data = KO_contrib_LiqCOY_sum, key = Enzyme, value = GFC)

KO_contrib_LiqTAM=subset(KO_contrib, Sampletype_site == "LiqTAM")
KO_contrib_LiqTAM_sum <- KO_contrib_LiqTAM %>%
  group_by(taxon, Enzyme) %>%
  summarise(
    GFC = sum(genome_function_count)
  )

KO_contrib_LiqTAM_sumsd= spread(data = KO_contrib_LiqTAM_sum, key = Enzyme, value = GFC)



# tablas taxones claves
setwd("~/1. TMKV/2. Results/14. Red")
Keystonetaxa_LC <- read_csv("keystone_Liq.COY.csv")
str(Keystonetaxa_LC)
View(Keystonetaxa_LC)
Keystonetaxa_LT <- read_csv("keystone_Liq.TAM.csv")
str(Keystonetaxa_LT)

# seleccionar solo los 30 taxa claves de COY
Keystonetaxa_LC=Keystonetaxa_LC[order(Keystonetaxa_LC$bn_sort, decreasing = TRUE),]
Keystonetaxa_LC=Keystonetaxa_LC[c(1:30),]

colnames(Keystonetaxa_LC)= c("taxon","deg_sort","bn_sort","sitio") 
colnames(Keystonetaxa_LT)= c("taxon","deg_sort","bn_sort","sitio") 

# Merge de dataframes

KTGC_COY=merge(Keystonetaxa_LC,KO_contrib_LiqCOY_sumsd)
KTGC_COY=KTGC_COY[,c(1,4:33)]
KTGC_COY[is.na(KTGC_COY)] <- 0
View(KTGC_COY)
KTGC_COY=column_to_rownames(KTGC_COY,"taxon")

KTGC_TAM=merge(Keystonetaxa_LT,KO_contrib_LiqTAM_sumsd)
KTGC_TAM=KTGC_TAM[,c(1,4:33)]
KTGC_TAM[is.na(KTGC_TAM)] <- 0
View(KTGC_TAM)
KTGC_TAM=column_to_rownames(KTGC_TAM,"taxon")


setwd("~/1. TMKV/2. Results/14. Red")
write.table(KTGC_COY, file = "KTGC_COY.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

write.table(KTGC_TAM, file = "KTGC_TAM.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# radarchart

library("fmsb")

KTGC_TAM_copy=KTGC_TAM
rownames(KTGC_TAM_copy)=c("ASV1","ASV2","ASV3","ASV4","ASV5","ASV6","ASV7","ASV8","ASV9","ASV10",
                          "ASV11","ASV12","ASV13","ASV14","ASV15","ASV16","ASV17","ASV18","ASV19","ASV20",
                          "ASV21","ASV22","ASV23","ASV24","ASV25","ASV26","ASV27","ASV28","ASV29","ASV30")
View(KTGC_TAM_copy)

radarchart(KTGC_TAM_copy,cglcol = "grey", axislabcol = "grey",caxislabels = seq(min(KTGC_TAM_copy),max(KTGC_TAM_copy),1), 
           seg=length(seq(min(KTGC_TAM_copy),max(KTGC_TAM_copy),1))-1, pty=32)





create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title
  )
}


create_beautiful_radarchart(
  data = KTGC_COY, caxislabels = c(0, 2,4,6,8,10),
  color = c("#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9",
            "#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9",
            "#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9",
            "#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9",
            "#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9",
            "#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9")
)

