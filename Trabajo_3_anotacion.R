###################################################################
# Nombre del Script           : Trabajo_3_Anotacion
# Descripción                 : Anotación del nombre de los genes de las sondas downreguladas de manera estadísticamente
#                               significativa y selección de genes compartidos por ambas líneas celulares.
# Autora                      : Lucía Sánchez García    
# Fecha última modificación   : 30 de junio de 2020
###################################################################

## 0. Establecimiento del directorio de trabajo
setwd("~/Documentos/GSE18198_RAW")

## Instalacion de la base de datos necesaria, en caso de que no se encuentre 
## ya instalada
## BiocManager::install("hgu95av2.db")


## 1. Cargar libreria de anotacion
# Se carga esta librera ya que es en la que se encuentra el nombre de los genes 
# de las sondas utilizadas, ya que en el articulo pone que se utilizo Affymetrix 
# Human Genome U133 Plus 2.0 Array annotation data (chip hgu133plus2)
## BiocManager::install("hgu1 33plus2.db") --> Instalacion de la base de datos necesaria
library("hgu133plus2.db")
 

## 2. Cargar los ficheros RData con los DEGs (FDR<0.05) previamente obtenidos.
# load("./MyResults_HPB_ALL_MOD.RData")
# load("./MyResults_KOPT_K1_MOD.RData")
toptableIQR_HPB_ALL <- readRDS(file = "toptableIQR_HPB_ALL.rds")
toptableIQR_KOPT_K1 <- readRDS(file = "toptableIQR_KOPT_K1.rds")

## --------------------------------------------------------------------------------- ##
## 3.-----------------------SELECCION DE GENES DOWNREGULADOS------------------------ ##
## ------------------------------------ HPB_ALL ------------------------------------ ##
## 3.1. Creacion de un subset de genes downregulados (logFC < 0) con p-valor <= 0.05
# Asi se seleccionan los genes significativos
ID.fdr.005.table_HPB_ALL_down <- subset(toptableIQR_HPB_ALL, toptableIQR_HPB_ALL$adj.P.Val<=0.05
                                   & toptableIQR_HPB_ALL$logFC < 0)
summary(ID.fdr.005.table_HPB_ALL_down)
dim(ID.fdr.005.table_HPB_ALL_down) # Muestra 682 sondas downreguladas significativamente

## 3.2. Adicion del simbolo a las sondas downreguladas del subset
#      Pasar a caracteres los nombres de las sondas
probenames.fdr.005_HPB_ALL_down <- as.character(rownames(ID.fdr.005.table_HPB_ALL_down))
#      Obtencion de los simbolos asociados a las sondas
list.GeneSymbol.fdr.005_HPB_ALL_down <- mget(probenames.fdr.005_HPB_ALL_down, hgu133plus2SYMBOL)
GeneSymbol_HPB_ALL_down <- as.character(list.GeneSymbol.fdr.005_HPB_ALL_down)
#      Adicion de la nueva columna obtenida
toptable.annotated_HPB_ALL_down <- cbind(ID.fdr.005.table_HPB_ALL_down,GeneSymbol_HPB_ALL_down)
#      Visualizacion de los primeros resultados
head(toptable.annotated_HPB_ALL_down) 

## ------------------------------------ KOPT_K1 ------------------------------------ ##

## 3.3 Creacion de un subset de genes downregulados (logFC < 0) con p-valor <= 0.05
ID.fdr.005.table_KOPT_K1_down <- subset(toptableIQR_KOPT_K1, toptableIQR_KOPT_K1$adj.P.Val<=0.05
                                   & toptableIQR_KOPT_K1$logFC < 0)
summary(ID.fdr.005.table_KOPT_K1_down)
dim(ID.fdr.005.table_KOPT_K1_down) # Muestra 5006 sondas downreguladas significativamente

## 3.4 Adicion del simbolo a las sondas downreguladas del subset
probenames.fdr.005_KOPT_K1_down <- as.character(rownames(ID.fdr.005.table_KOPT_K1_down))
list.GeneSymbol.fdr.005_KOPT_K1_down <- mget(probenames.fdr.005_KOPT_K1_down, hgu133plus2SYMBOL)
GeneSymbol_KOPT_K1_down <- as.character(list.GeneSymbol.fdr.005_KOPT_K1_down)
toptable.annotated_KOPT_K1_down <- cbind(ID.fdr.005.table_KOPT_K1_down,GeneSymbol_KOPT_K1_down)
head(toptable.annotated_KOPT_K1_down) 


## 3.5. Seleccion de los genes significativos en ambas lineas celulares
iteration <- 0
lista_genes_down <- c()
# Creacion de un loop de busqueda de genes significativos para HPB_ALL en la tabla de
# genes significativos en KOPT_K1
for (i in toptable.annotated_HPB_ALL_down$GeneSymbol_HPB_ALL_down) {
  if (i != "NA" && i %in% toptable.annotated_KOPT_K1_down$GeneSymbol_KOPT_K1_down) {
    # Si los genes no se encuentran previamente anotados en la lista "lista_genes",
    # seran introducidos -> evitamos la introduccion de genes repetidos, ya que sondas
    # diferentes pueden representar el mismo gen
    if (! i %in% lista_genes_down) {
      iteration <- iteration + 1; print(i);
      lista_genes_down[[iteration]] <- i} } }
print(iteration) # 171 genes downregulados significativamente en conjunto
sort(lista_genes_down) # Ordenación  alfabética


## 3.6. Creacion de un fichero .txt con los genes downregulados (171 genes)
write.table(lista_genes_down, "lista_genes_down.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


## ---------------------------- En caso de que se necesite ------------------------- ##
## 4.------------------------SELECCION DE GENES UPREGULADOS------------------------- ##
## ------------------------------------ HPB_ALL ------------------------------------ ##
## 4.1. Creacion de un subset de genes upregulados (logFC > 0) con p-valor <= 0.05
# Asi se seleccionan los genes Significativos
ID.fdr.005.table_HPB_ALL_up <- subset(toptableIQR_HPB_ALL, toptableIQR_HPB_ALL$adj.P.Val<=0.05
                                        & toptableIQR_HPB_ALL$logFC > 0)
summary(ID.fdr.005.table_HPB_ALL_up)
dim(ID.fdr.005.table_HPB_ALL_up) # Muestra 1423 sondas upreguladas significativamente

## 4.2. Adicion del simbolo a las sondas upreguladas del subset
#      Pasar a caracteres los nombres de las sondas
probenames.fdr.005_HPB_ALL_up <- as.character(rownames(ID.fdr.005.table_HPB_ALL_up))
#      Obtencion de los s?mbolos asociados a las sondas
list.GeneSymbol.fdr.005_HPB_ALL_up <- mget(probenames.fdr.005_HPB_ALL_up, hgu133plus2SYMBOL)
GeneSymbol_HPB_ALL_up <- as.character(list.GeneSymbol.fdr.005_HPB_ALL_up)
#      Adicion de la nueva columna obtenida
toptable.annotated_HPB_ALL_up <- cbind(ID.fdr.005.table_HPB_ALL_up,GeneSymbol_HPB_ALL_up)
#      Visualizacion de los primeros resultados
head(toptable.annotated_HPB_ALL_up) 

## ------------------------------------ KOPT_K1 ------------------------------------ ##

## 4.3 Creacion de un subset de genes upregulados (logFC > 0) con p-valor <= 0.05
ID.fdr.005.table_KOPT_K1_up <- subset(toptableIQR_KOPT_K1, toptableIQR_KOPT_K1$adj.P.Val<=0.05
                                        & toptableIQR_KOPT_K1$logFC > 0)
summary(ID.fdr.005.table_KOPT_K1_up)
dim(ID.fdr.005.table_KOPT_K1_up) # Muestra 2759 sondas upreguladas significativamente

## 4.4 Adicion del simbolo a las sondas downreguladas del subset
probenames.fdr.005_KOPT_K1_up <- as.character(rownames(ID.fdr.005.table_KOPT_K1_up))
list.GeneSymbol.fdr.005_KOPT_K1_up <- mget(probenames.fdr.005_KOPT_K1_up, hgu133plus2SYMBOL)
GeneSymbol_KOPT_K1_up <- as.character(list.GeneSymbol.fdr.005_KOPT_K1_up)
toptable.annotated_KOPT_K1_up <- cbind(ID.fdr.005.table_KOPT_K1_up,GeneSymbol_KOPT_K1_up)
head(toptable.annotated_KOPT_K1_up) 


## 4.5. Seleccion de los genes significativos en ambas lineas celulares
iteration <- 0
lista_genes_up <- c()
# Creacion de un loop de busqueda de genes significativos para HPB_ALL en la tabla de
# genes significativos en KOPT_K1
for (i in toptable.annotated_HPB_ALL_up$GeneSymbol_HPB_ALL_up) {
  if (i != "NA" && i %in% toptable.annotated_KOPT_K1_up$GeneSymbol_KOPT_K1_up) {
    # Si los genes no se encuentran previamente anotados en la lista "lista_genes_up",
    # seran introducidos -> evitamos la introducci?n de genes repetidos, ya que sondas
    # diferentes pueden representar el mismo gen
    if (! i %in% lista_genes_up) {
      iteration <- iteration + 1; print(i);
      lista_genes_up[[iteration]] <- i} } }
print(iteration) # 75 genes upregulados significativamente en conjunto
sort(lista_genes_up) # Ordenación  alfabética

## 4.6. Creacion de un fichero .txt con los genes upregulados (75 genes)
write.table(lista_genes_up, "lista_genes_up.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

