#################################################
## Red comunidad microbiana y Keystone species ##
#################################################

###############
### Paquetes ##
###############

library(phyloseq)
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
library(VennDiagram)

######################
## Archivo Phyloseq ## 
######################

# Cargar el archivo phyloseq
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.COY<- readRDS("ps.Liq.COY19.F2.rds")

# Revisión
ps.Liq.COY
sample_variables(ps.Liq.COY)
rank_names(ps.Liq.COY)

# Usaremos solo CORE así que se filtrará. 

# Creamos una funcion
prev = apply(X = otu_table(ps.Liq.COY), MARGIN=ifelse(taxa_are_rows(ps.Liq.COY), yes=1, no=2),
             FUN=function(x){sum(x > 0)})

# Transformar la tax_table en un dataframe para que no salga el error
tax_table_df=as.data.frame(ps.Liq.COY@tax_table@.Data)

# y ordenamos en un dataframe
df.prev = data.frame(Prevalencia=prev, AbundanciaTotal=taxa_sums(ps.Liq.COY),
                     tax_table_df)

prevdf1 = subset(df.prev, Phylum %in% get_taxa_unique(ps.Liq.COY, "Phylum"))

# Definimos el umbral de prevalencia a un 5%
#Se puede colocar el porcentaje deseado 0.9 * nsamples(ps) o el numero
prevalence.core = 8

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa.core = rownames(prevdf1)[(prevdf1$Prevalencia >= prevalence.core)]

#Filtramos ahora creando un archivo phyloseq que tendra el filtro de abundancia
ps.Liq.COY.core = prune_taxa(keepTaxa.core, ps.Liq.COY)

#Guardamos los archivos como rds
saveRDS(ps.Liq.COY.core, "ps.Liq.COY.core.rds")

# Más adelante se presentan problemas para graficar la red debido a los NA, entonces se modifica la tabla de taxonomia:
tax_table(ps.Liq.COY.core)[is.na(tax_table(ps.Liq.COY.core))] <- "NoId"

#############################################
## Perfil taxonomico del conjunto de datos ##
#############################################

# Perfil taxonomico: tabla de datos disponibles (util para consultar durante el análisis).
# Extraer y tabular los datos
taxotutable <- phyloseq::psmelt(ps.Liq.COY.core)

# Exportar/guardar tabla
setwd("~/1. TMKV/2. Results/14. Red")
write.table(taxotutable, file = "SummarizedTaxAssigments_ps.Liq.COY.tsv", sep = "\t", quote = FALSE)
#### El encabezado debe ser desplazado una columna a la derecha en el excel. 

#####################
### INFERIR LA RED ##
#####################

# SpiecEasi: inferencia de redes
se_mb <- spiec.easi(ps.Liq.COY.core, method='mb', 
                    lambda.min.ratio=1e-3, nlambda=100, 
                    pulsar.params=list(rep.num=100, ncores=6))

# Puedes exportar/guardar el objeto 'se_mb'
saveRDS(se_mb, "se_mb_ps.Liq.COY.RDS")

save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Keystonetaxa_LiqCOY_TMKV.RData")

# Leer red

setwd("~/1. TMKV/2. Results/14. Red")
se_mb_ps.Liq.COY <- readRDS(file = "se_mb_ps.Liq.COY.RDS")

# Construir red a partir de la 'sparse adjacency matrix' o 'refit matrix'
se_net <- adj2igraph(getRefit(se_mb_ps.Liq.COY), 
                     rmEmptyNodes = TRUE, diag = FALSE, 
                     vertex.attr = list(name = taxa_names(ps.Liq.COY.core))) # Usamos el ID de las taxas para nombrar los vértices o nodos de la red


# Creamos dataframe con codigo de colores
colour <- c("NoId"="#fff6ed","Proteobacteria"="#69AFB2","Actinobacteriota"="#E4B664","Acidobacteriota"="#6184BE","Bacteroidota"="#E0E168",
            "Planctomycetota"="#C95F6F","Verrucomicrobiota"="#64B0D1","Myxococcota"="#1E5664","Armatimonadota"="#57599D",
            "Cyanobacteria"="#A1BE57","Chloroflexota"="#7A609F","Gemmatimonadota"="#BB68A1","Patescibacteria"="#67AA8C",
            "Eremiobacterota"="#7AAD59","Bdellovibrionota"="#CD6926","Nitrospirota"="#C96299","Desulfobacterota_B"="#C23234",
            "Firmicutes"="#4D7AAF","Dormibacterota"="#9B66A0","Methylomirabilota"="#918854","Elusimicrobiota"="#355B8A",
            "Omnitrophota"="#040006","Firmicutes_A"="#7F7F7F","Firmicutes_C"="#DAD6C1","Dependentiae"="#738E3C",
            "Desulfobacterota"="#8F3734","Eisenbacteria"="#E9E337")

# Guardar Red 1

setwd("~/1. TMKV/2. Results/14. Red")
png("Red1_Liq.COY.png",width = 900, height = 500)
plot_network(se_net, ps.Liq.COY.core, type = "taxa",  color = "Phylum", label = NULL, alpha = .8) + 
  scale_color_manual(values = colour)
dev.off()

# Extraer matriz de adyacencia
net_class <- as_adjacency_matrix(se_net, type = "both")
# Generar objeto de clase 'network'
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(ps.Liq.COY.core), 
                     matrix.type = "adjacency", directed = F)
# Visualizamos la red usando ggnet2
ggnet2(net_class)

###############################
## Características de la red ##
###############################

# Copiar el objeto 'se_net' a 'net'
net <- se_net
# Características básicas de una red: vértices o nodos (nodes or vertex), bordes (edges), nombre de los nodos, y número total de nodos y bordes.
nodes <- V(net) # Nodos
edges <- E(net) # Bordes
node.names <- V(net)$name # Nombre de nodos
num.nodes <- vcount(net) # Número total de nodos
num.edges <- ecount(net) # Número total de bordes
# V() es por vértices o nodos
# E() es por edges

# Extraer los coeficientes de regresión del objeto 'se_mb'
beCOYat <- as.matrix(symBeta(getOptBeta(se_mb_ps.Liq.COY)))

# Calcular el número de edges positivos y negativos en la red
positive <- length(beCOYat[beCOYat>0])/2 
negative <- length(beCOYat[beCOYat<0])/2 
total <- length(beCOYat[beCOYat!=0])/2

# Asignarle color a los edges positivos y negativos y a visualizar la red.
tax_ids <- taxa_names(ps.Liq.COY.core)
edges <- E(net) # edges
edge_colors <- c()
for(e_index in 1:length(edges)){
  adj_nodes <- ends(net,edges[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- beCOYat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"forestgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red") # negative
  }
}
E(net)$color <- edge_colors
# Generamos el objeto de clase network y graficamos usando ggnet2.
# Extraer matriz de adyacencia
net_class <- as_adjacency_matrix(net, type = "both")
# Generar objeto de clase 'network'
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(ps.Liq.COY.core), 
                     matrix.type = "adjacency", directed = F)

# Graficar red
ggnet2(net_class, edge.color = E(net)$color)

# Extraer tabla de taxonomía
tax_tbl <- as.data.frame(ps.Liq.COY.core@tax_table@.Data)
# Reemplazar el nombre de cada nodo por su Phylum
V(net)$name <- as.character(getTaxonomy(V(net)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Guardar la lista de Phylum por nodo en 'nodenames'
nodenames <- V(net)$name

# Extraer matriz de adyacencia
net_class <- as_adjacency_matrix(net, type = "both")
# Generar objeto de clase 'network'
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(ps.Liq.COY.core), 
                     matrix.type = "adjacency", directed = F)
# Graficar red
ggnet2(net_class, color = nodenames, edge.color = E(net)$color)

# Definir paleta de colores, transparencias, etc. 
ggnet2(net_class, color = nodenames, edge.color = E(net)$color, palette = colour, alpha = 1, node.size = 4, edge.alpha = 0.4)

## Asignar coordenadas al diseño de la red
net$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
net$layout <- layout.fruchterman.reingold(net)

# Graficar red
ggnet2(net_class, mode = net$layout, color = nodenames, edge.color = E(net)$color, palette = colour, alpha = 1, node.size = 3, edge.alpha = 0.2)

# Guardar Red 2 manual

############################
### Estadística de la red ##
############################

# Calcular node degree, node centrality, y transitivity de la red.

# Copiar el objeto 'se_net' a 'net'
net <- se_net
# Asignar coordenadas al diseño de la red
net$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
net$layout <- layout.fruchterman.reingold(net)

#############
### Degree ##
#############

# Calcular node degree
deg <- igraph::degree(net, mode = "all")
# Calcular degree distribution
deg.dist <- degree_distribution(net, mode = "all", cumulative = F)

Degree.Liq.COY=plot(deg.dist, xlab = "Nodes degree", ylab = "Probability")
Degree.Liq.COY =print(lines(deg.dist))

######################
### Node centrality ##
######################

# Closeness
clos <- igraph::closeness(net, mode = "all")

# Betweenness
betw <- igraph::betweenness(net, v = V(net))

centrality.Liq.COY=centralityPlot(net, include = c("Betweenness", "Closeness", "Degree")) + 
  theme(axis.text.y = element_blank())

png("Centrality.Liq.COY.png",width = 500, height = 900)
centrality.Liq.COY
dev.off()

###################
### Transitivity ##
###################

clustering_coeff_global <- transitivity(net, type = "global")
clustering_coeff_local <- transitivity(net, type = "local")

net.knn <- knn(net, vids = V(net))
head(net.knn$knn)

Transitividad.Liq.COY=as.data.frame(net.knn$knn)
View(Transitividad.Liq.COY)
Transitividad.Liq.COY=rownames_to_column(Transitividad.Liq.COY)

write_csv(Transitividad.Liq.COY,"Transitividad.Liq.COY.csv")

###########################
### Estructura de la red ##
###########################

# Copiar el objeto 'se_net' a 'net'
net <- se_net
# Asignar coordenadas al diseño de la red
net$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
net$layout <- layout.fruchterman.reingold(net)
# Detección de módulos
wt <- walktrap.community(net)
membership(wt) %>% head()

igraph::plot_dendrogram(wt)

##################
### Modularidad ##
##################

# Calcular modularidad
modularity(net, membership(wt))

plot(wt, net)

net_class <- as_adjacency_matrix(net, type = "both")
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(ps.Liq.COY.core), 
                     matrix.type = "adjacency", directed = F)
ggnet2(net_class, mode = net$layout, color = wt$membership)

ggnet2(net_class, mode = net$layout, color = wt$membership) + 
  stat_ellipse(aes(color = factor(wt$membership)), type = "norm")

nodenames <- as.character(getTaxonomy(V(net)$name, tax_tbl, level = "phylum", useRownames = TRUE))

ggnet2(net_class, mode = net$layout, color = nodenames, edge.color = edge_colors, palette = colour, alpha = 1, node.size = 4, edge.alpha = 0.4) + 
  stat_ellipse(aes(group = factor(wt$membership)), type = "norm")


######################################
## Graficos finales de redes a usar ##
######################################

ggnet2(net_class, mode = net$layout, color = nodenames, edge.color = edge_colors, palette = colour, alpha = 1, node.size = 4, edge.alpha = 0.4) + theme(legend.position="none")

ggnet2(net_class, mode = net$layout, color = wt$membership, edge.color = edge_colors, alpha = 1, node.size = 4, edge.alpha = 0.4)

# Guardado manual 

#################################
# Búsqueda de keystone species ##
#################################


# Calcular degree
deg <- igraph::degree(net)
# Ordenar nodos según degree de mayor a menor
deg_sort <- sort(deg, decreasing = TRUE)
head(deg_sort) # Mira las primeras 6 filas del objeto 'deg_sort'
deg_df <- as.data.frame(deg_sort)
# Graficar histograma

ggplot(deg_df, aes(x = deg_sort)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60)) + 
  geom_vline(aes(xintercept=mean(deg_sort)),
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Degree", y = "Node count")

## Node centrality
# Calcular betweenness
bn <- igraph::betweenness(net)
# Ordenar nodos según betweenness de mayor a menor
bn_sort <- sort(bn, decreasing = TRUE)
head(bn_sort) # Mira las primeras 6 filas del objeto 'bn_sort'
bn_df <- as.data.frame(bn_sort)

# Graficar histograma
ggplot(bn_df, aes(x = bn_sort)) + 
  geom_histogram(binwidth = 30) +  
  geom_vline(aes(xintercept=mean(bn_sort)), 
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Node count")
pdeg <- ggplot(deg_df, aes(x = deg_sort)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60)) + 
  geom_vline(aes(xintercept=mean(deg_sort)),
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Degree", y = "Node count", title = "Network Nodes Degree")+  
  geom_vline(aes(xintercept=19), 
             color="red", linetype="dashed", size=1)

pbn <- ggplot(bn_df, aes(x = bn_sort)) + 
  geom_histogram() +  
  geom_vline(aes(xintercept=mean(bn_sort)), 
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Node count", title = "Network Nodes Centrality: Betweenness")+  
  geom_vline(aes(xintercept=1000), 
             color="red", linetype="dashed", size=1)

grid.arrange(pdeg, pbn, ncol = 2)

# guardado manual Degree_Betweenness_Keystonetaa_Liq.COY 1200x665


head(deg_df)

deg_df$TaxID <- row.names(deg_df) 
head(deg_df)

deg_high_df <- dplyr::filter(deg_df, deg_sort > 19)
head(deg_high_df)

head(bn_df)

bn_df$TaxID <- row.names(bn_df)
bn_high_df <- dplyr::filter(bn_df, bn_sort > 1000)
head(bn_high_df)

keystone <- merge(deg_high_df, bn_high_df, all.x = FALSE)
head(keystone)


keystone <- keystone[order(keystone$deg_sort, decreasing = TRUE),]
row.names(keystone) <- keystone$TaxID

View(keystone)
write_csv(keystone, "keystone_Liq.COY.csv")
save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Keystonetaxa_LiqCOY_TMKV.RData")

############################
## Graficos transitividad ##
############################

# Distribucion
setwd("~/1. TMKV/2. Results/14. Red")
T.Liq.TAM=read_csv("Transitividad.Liq.TAM.csv")

colnames(T.Liq.TAM)=c("OTU","Transitividad")

TTAM=ggplot(T.Liq.TAM, aes(x = Transitividad)) + 
  geom_histogram(binwidth = 1) + 
  theme_minimal()

T.Liq.COY=read_csv("Transitividad.Liq.COY.csv")
colnames(T.Liq.COY)=c("OTU","Transitividad")

TCOY=ggplot(T.Liq.COY, aes(x = Transitividad)) + 
  geom_histogram(binwidth = 1) + 
  theme_minimal()

grid.arrange(TTAM, TCOY, ncol = 2)

##########################
## Graficos taxa claves ##
##########################

setwd("~/1. TMKV/2. Results/14. Red")
library(readxl)
library(ggrepel)
library(readr)

# Importar data de taxonomia
Tax_contrib<-read.delim("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV Core/Picrust2 TMKV/taxonomia.tsv")
# Cambiar el nombre de columnas para hacer merge 
colnames(Tax_contrib)=c("taxon","Domain","Phylum","Class","Order","Family","Genus","Species")

# importar datos taxaclaves

Keystonetaxa_LC <- read_csv("keystone_Liq.COY.csv")
str(Keystonetaxa_LC)
Keystonetaxa_LT <- read_csv("keystone_Liq.TAM.csv")
str(Keystonetaxa_LT)

# seleccionar solo los 30 taxa claves de COY
Keystonetaxa_LC=Keystonetaxa_LC[order(Keystonetaxa_LC$bn_sort, decreasing = TRUE),]
Keystonetaxa_LC=Keystonetaxa_LC[c(1:30),]

# Agregamos variable por sitio y unimos tablas

Keystonetaxa_LC=mutate(Keystonetaxa_LC,sitio=rep("COY",30))
Keystonetaxa_LT=mutate(Keystonetaxa_LT,sitio=rep("TAM",30))
Keystonetaxa= rbind(Keystonetaxa_LC,Keystonetaxa_LT)

View(Keystonetaxa)

# Unimos taxa claves con taxonomia

colnames(Tax_contrib)
colnames(Keystonetaxa)= c("taxon","deg_sort","bn_sort","sitio") 
Keystonetaxa_ID=merge(Keystonetaxa,Tax_contrib)
View(Keystonetaxa_ID)

# Usamos ggplot2

Filo <- c("Proteobacteria","Actinobacteriota","Acidobacteriota","Bacteroidota",
          "Planctomycetota","Cyanobacteria","NA")
cod <- c("#69AFB2","#E4B664","#6184BE","#E0E168","#C95F6F","#A1BE57","gray")

color <- data.frame(Filo, cod)



Keystonetaxa_plot = ggplot(Keystonetaxa_ID, aes(x = bn_sort, y = deg_sort, label = Genus, color=Phylum, shape=sitio)) + 
  scale_y_continuous(limits = c(15, 30), breaks = c(15,20,25,30)) + 
  scale_x_continuous(limits = c(700, 2700), breaks = c(700,1200,1700,2200,2700)) + 
  geom_jitter(alpha = 0.5, size = 4)+ geom_text_repel(size = 6, max.overlaps = 20)+
  theme(axis.text=element_text(size=22)) +
  theme_minimal() +
  labs(x = "Betweenness", y = "Degree", 
       title = "Nodes With Highest Degree And Betweenness: Keystone Nodes") +
  scale_color_manual(values = color$cod, limits=color$Filo)

Keystonetaxa_plot

setwd("~/1. TMKV/2. Results/14. Red")
ggsave("Keystonetaxa_plot.png", plot=Keystonetaxa_plot,width = 15, height = 13, device= "png",path = "~/1. TMKV/2. Results/14. Red")


save.image("~/1. TMKV/0. Scripts&Environment/2. Enviroment/Keystonetaxa_LT.RData")

# diagrama de venn

v2 <- venn.diagram(list(COY=Keystonetaxa_LC$TaxID, TAM=Keystonetaxa_LT$TaxID),
                   fill = c("#ff9a4a", "#6ecb63"),
                   alpha = c(0.5, 0.5),
                   filename=NULL)

grid.newpage()
grid.draw(v2)

# guardado manual 