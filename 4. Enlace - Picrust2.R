#######################
## Enlace - Picrust2 ##
#######################

# Paginas utiles 
# https://www.nature.com/articles/s41587-020-0548-6
# https://github-wiki-see.page/m/picrust/picrust2/wiki/PICRUSt2-Tutorial-%28v2.3.0-beta%29

#################################################
## Obtencion de datos para ingresar en picrust ##
#################################################

# Paquetes
library(Biostrings)
library(dada2)
library(readr)
library(readxl)
library(phyloseq)

# Llamar phyloseq 
setwd("~/1. TMKV/5. Phyloseq")
ps.Liq.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.TAM19.F2.rds")
ps.Sus.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.TAM19.F2.rds")
ps.Liq.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Liq.COY19.F2.rds")
ps.Sus.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sus.COY19.F2.rds")
ps.Sue.COY19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.COY19.F2.rds")
ps.Sue.TAM19.F2 <- readRDS("~/1. TMKV/5. Phyloseq/ps.Sue.TAM19.F2.rds")


ps.picrust_F2<- merge_phyloseq(ps.Liq.TAM19.F2,ps.Liq.COY19.F2,ps.Sus.TAM19.F2,ps.Sus.COY19.F2,ps.Sue.TAM19.F2, ps.Sue.COY19.F2)


# Descargar tabla de taxonomia para cambiar los nombres a los ASV
Tabla.Taxonomia.ps <-tax_table(ps.picrust_F2)
setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo")
write.table(Tabla.Taxonomia.ps, "taxonomia.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Creamos un nuevo ps para generar el objeto refseq
dna <- Biostrings::DNAStringSet(taxa_names(ps.picrust_F2))
# Nota: Revisar que genera taxa_names(ps.picrust_F2), tiene que generar esto:
#  [2] "TACGAAGGGGGCTAGCGTTGCTCGGAATGACTGGGCGTAAAGGGCGCGTAGGCGGATCTGTCAGTCAGACGTGAAATTCCTGGGCTTAACCTGGGGGCTGCGTTTGAGACGGTGGGTCTAGAGTTTGGAAGAGGGTCGTGGAATTCCCAGTGTAGAGGTGAAATTCGTAGATATTGGGAAGAACACCGGTGGCGAAGGCGGCGACCTGGTCCTGGACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGTGCGCTGGATGTTGGGGAGCCTAGCTTTTCAGTGTCGTAGCTAACGCGGTAAGCGCACCGCCTGGGGAGTACGGCCGCAAGG"            
#  [3] "TACGAAGGGGGCTAGCGTTGCTCGGAATGACTGGGCGTAAAGGGCGCGTAGGCGGACATGTCAGTCAGACGTGAAATTCCTGGGCTTAACCTGGGGGCTGCGTTTGAGACGGTGTGTCTAGAGTTTGGAAGAGGGTCGTGGAATTCCCAGTGTAGAGGTGAAATTCGTAGATATTGGGAAGAACACCGGTGGCGAAGGCGGCGACCTGGTCCTGGACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGTGCGCTGGATGTTGGGGAGCCTAGCTTTTCAGTGTCGTAGCTAACGCGGTAAGCGCACCGCCTGGGGAGTACGGCCGCAAGG"            
#  [4] "TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGTGCGTAGGCGGATAAGTAAGTCAGTGGTGAAATCTCCAGGCTTAACCTGGAAACTGCCATTGATACTATTTATCTTGAATTACGTGGAGGTGAGCGGAATATGTCATGTAGCGGTGAAATGCTTAGATATGACATAGAACACCAATTGCGAAGGCAGCTCACTACACGTTGATTGACGCTGAGGCACGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGGATACTCGACATACGCGATATACAGTGTGTGTCTGAGCGAAAGCATTAAGTATCCCACCTGGGAAGTACGACCGCAAGG" 

names(dna) <- taxa_names(ps.picrust_F2)
ps2 <- merge_phyloseq(ps.picrust_F2, dna)
ps2

saveRDS(ps2, "~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/PhyloseqPicrust2")

# Tu nuevo phyloseq tiene que tener el objeto DNAStringSet
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4723 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 4723 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 4723 reference sequences ]

# Generar objeto con tabla de taxonomia
Tabla.Taxonomia.ps <-tax_table(ps2)
head(Tabla.Taxonomia.ps)

#ahora hacemos el archivo fasta con las secuencias de los ASVs, ocupando el refseq del ps2 generado anteriormente
refseq16S_DNAStringSet <- phyloseq::refseq(ps2, errorIfNULL=TRUE) 
# Guardamos el objeto que generamos como archivo fasta 
setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo")
writeXStringSet(refseq16S_DNAStringSet, "./refseq16S.fna", append=FALSE, compress=FALSE, compression_level=NA, format="fasta") 
# Se genera archivo fasta con cada secuencia del ps2

# ahora hacemos el archivo con la abundancia de cada ASV en las muestras, ocupando la otu_table del ps2

ASV16S_matrix <- t(as(otu_table(ps2),"matrix")) 
setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo")
write.table(ASV16S_matrix, "./seqtab.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# Va a generar un archivo de texto que puede ser visto en excel (delimitado por tab), y muestra la abundancia de cada ASV en cada muestra. 
#más info en https://forum.qiime2.org/t/importing-dada2-and-phyloseq-objects-to-qiime-2/4683

####################
## EN EL TERMINAL ##
####################

# Entrar a servidor
# Moverse a carpeta de trabajo donde estan los archivos generados usando los comandos:
#  ls: listar las carpetas del directorio donde estoy.
#  cd: ir al directorio

# Activar conda con los siguientes comandos (SOLO HACER UNA VEZ y en caso de que el comando "conda" no exista: "conda: command not found)
# source /opt/miniconda3/bin/activate
# conda init

# Activar picrust con el comando:
# conda activate picrust2
# https://programmerclick.com/article/60401346449/

#En el terminal con PICRUST2 cargado correr el siguiente comando:


#picrust2_pipeline.py -s refseq16S.fna -i seqtab.tsv -o picrust2_out_pipeline -p 15
# Puede tardar varios minutos

# https://github-wiki-see.page/m/picrust/picrust2/wiki/Full-pipeline-script
# https://qiime2.org/
# En Picrust, la precisión de cualquier tipo de muestra dependerá en gran medida de la 
# disponibilidad de genomas de referencia apropiados. Se puede evaluar parcialmente
# este problema calculando los valores por ASV y el índice de taxón secuenciado más
# cercano ponderado por muestra (NSTI), lo que le dará una idea aproximada de qué tan
# bien representados están sus ASV en la base de datos de referencia.

# PARA TODO
# 7 of 11306 ASVs were above the max NSTI cut-off of 2.0 and were removed.

# Genera tablas con todos los KO

#Pre procesamiento para volver a ingresar a picrust2

#Importar data (cuantas enzimas por microorganismo)

########################
## Filtrado de tablas ##
########################
#Pre procesamiento para volver a ingresar a picrust2

#Importar data (cuantas enzimas por microorganismo)

setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline")
EC_predicted_by_ASV <- read.delim("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline/EC_predicted.tsv.gz", row.names=1)
KO_predicted_by_ASV <- read.delim("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline/KO_predicted.tsv.gz", row.names=1)

#View(EC_predicted_by_ASV)
#View(KO_predicted_by_ASV)

# Aquí se seguirá trabajando solo con KO (si se usan ambos, se debe hacer para ambos archivos)
# Revisar si alguna fila es 0

conteo_0=0
conteo_total=0
for (i in colnames(KO_predicted_by_ASV)) {
  suma=sum(KO_predicted_by_ASV[,i])
  if (suma==0){
    conteo_0=conteo_0+1
  }
  conteo_total=conteo_total+1 }

#Reducir numeros a 1
reduccion<-function(x){
  for (i in 1:length(x)){
    if (x[i]>1) x[i]=1
  }
  x
}

KO_predicted_by_ASV_red<-apply(KO_predicted_by_ASV, 2, reduccion)
EC_predicted_by_ASV_red<-apply(EC_predicted_by_ASV, 2, reduccion)
#View(KO_predicted_by_ASV_red)

#Importar lista de enzimas Formato del excel segun "Fosforo.xlsx"
KOmetadat<- read_excel("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline/Fosforo.xlsx", sheet = "KO")

#KO_sNa
KO0<-KOmetadat$KO
KO1<-KO0[!(is.na(KO0))]
KO_sNa<-gsub(":",".",KO1)
KO_sNa
#¿Eliminan los NA? Probablemente sea para conjunto de datos más grandes donde se generan columnas con NA. Preguntar a Claudio

#Filtrar listas y juntar como vector
EC_in_data<-c()
KO_in_data<-c()
Code_in_data<-c()

# Este es un vector vacio donde quedan almacenadas los KO 
KO_sNa_not_in_data<-c()

counter=0
# Primero, crea un variable, en este caso "ec", que sería cada variable dentro del vector KO_sNA (que son la lista de nuestros KO)
for (ec in KO_sNa) {
  counter=counter+1
  if (ec %in% colnames(EC_predicted_by_ASV_red)) {
    EC_in_data<-c(EC_in_data,ec)
    Code_in_data<-c(Code_in_data,ec)
  }
  else {
    ko<-KO_sNa[counter]
    if (ko %in% colnames(KO_predicted_by_ASV_red)){
      KO_in_data<-c(KO_in_data,ko)
      Code_in_data<-c(Code_in_data,ko)
    }
    else{
      if (ec=="-") {
        KO_sNa_not_in_data<-c(KO_sNa_not_in_data,ko)
      }
      else {
        KO_sNa_not_in_data<-c(KO_sNa_not_in_data,ec) 
      }
    }
  }
}

#Todos los EC y KO en la data
EC_in_data<-c()
EC_not_in_data0<-c()
for (ec in KO_sNa) {
  if (ec %in% colnames(EC_predicted_by_ASV_red)) {
    EC_in_data<-c(EC_in_data,ec)
  }
  else {
    
  }
}
KO_in_data<-c()
for (ko in KO_sNa) {
  if (ko %in% colnames(KO_predicted_by_ASV_red)){
    KO_in_data<-c(KO_in_data,ko)
  }
  else {
    
  }
}
full_in_data<-c(EC_in_data,KO_in_data)



# Filtrar
data_filt<-data.frame(row.names = rownames(EC_predicted_by_ASV_red))
for (i in full_in_data) {
  if (i %in% colnames(EC_predicted_by_ASV_red)) {
    data_filt<-cbind.data.frame(data_filt,EC_predicted_by_ASV_red[,i])
  }
  if (i %in% colnames(KO_predicted_by_ASV_red)) {
    data_filt<-cbind.data.frame(data_filt,KO_predicted_by_ASV_red[,i])
  }
}
colnames(data_filt)<-c(full_in_data)
head(data_filt)

#Exportar en .tsv

EC_KO_predicted_filt0<-data.frame(sequence=rownames(data_filt))
KO_predicted_filt<-cbind.data.frame(EC_KO_predicted_filt0,data_filt)

colnames(KO_predicted_filt)

#colnames(KO_predicted_filt)<- Gen
setwd("~/1. TMKV/2. Results/13. Picrust/Picrust TMKV todo/picrust2_out_pipeline")
write.table(KO_predicted_filt,sep = "\t",file = "KO_predicted_filt.tsv",row.names = FALSE)

# En este caso se estudiará "que bacterias (o proporción de estas) cumple con un rol solubilizador de P" (No el "potencial solubilizador"), por eso se lleva a codigo binario la tabla (Presencia/ausencia)

############################
## Segunda parte picrust2 ##
############################

# Volver abrir picrust en la terminal_ conda activate picrust 
# Guardar todos los archivos en la carpeta "picrust2_out_pipeline" desde R
# Hacer correr lo siguiente desde el terminal, para poder crear nuestro archivo final a usar en PICRUST2
# metagenome_pipeline.py -i seqtab.tsv -m marker_predicted_and_nsti.tsv.gz -f KO_predicted_filt.tsv
# --strat_out -o KO_metagenome_out

#volver a la carpeta q se trabajara