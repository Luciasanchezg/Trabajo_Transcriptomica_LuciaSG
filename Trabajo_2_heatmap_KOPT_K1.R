###################################################################
# Nombre del Script           : Trabajo_2_Heatmap_KOPT_K1
# Descripción                 : Realización de heatmap o mapa de calor para los genes downregulados con un logFC < -1 para la
#                               línea celular KOPT_k1
# Autora                      : Lucía Sánchez García    
# Fecha última modificación   : 30 de mayo de 2020
###################################################################

## 0. Establecimiento del directorio de trabajo
setwd("~/Documentos/GSE18198_RAW")

## Cargar las librerías necesarias
library(pheatmap)
library("genefilter")
library("affy")
library("hgu133plus2.db")
library(data.table)

## ------------------------------ REALIZACIÓN DE HEATMAP ------------------------------ ##
## ------------------------------------------------------------------------------------ ##
## ------------------------------------- KOPT_K1 -------------------------------------- ##
## ------------------------------------------------------------------------------------ ##
## 1. Cargar elementos
# Se carga el objeto eset
eset_KOPT_K1 <- readRDS(file = "eset_KOPT_K1.rds")

# Se realiza el filtrado por el rango intercuartílico
esetIQR_KOPT_K1 <- varFilter(eset_KOPT_K1, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
expression_total_KOPT_K1 <- exprs(esetIQR_KOPT_K1)

# Se cargan la tabla de KOPT_K1 con el pvalor y el logFC
toptableIQR_KOPT_K1 <- readRDS(file = "toptableIQR_KOPT_K1.rds")



## 2. Filtrado de las sondas de interés de toptableIQR_KOPT_K1
subgrupo_KOPT_K1_pvalor <- subset(toptableIQR_KOPT_K1, toptableIQR_KOPT_K1$adj.P.Val <= 0.05 
            & (toptableIQR_KOPT_K1$logFC < -1 ))

subgrupo_KOPT_K1_pvalor
dim(subgrupo_KOPT_K1_pvalor) # En total, 118 sondas downreguladas en KOPT_k1 (con LogFC < -1)



## 3. Adicion del simbolo a las sondas 
#      Pasar a caracteres los nombres de las sondas
probenames.KOPT_K1 <- as.character(rownames(subgrupo_KOPT_K1_pvalor))
#      Obtencion de los simbolos asociados a las sondas
list.GeneSymbol.KOPT_K1 <- mget(probenames.KOPT_K1, hgu133plus2SYMBOL)
GeneSymbol_KOPT_K1 <- as.character(list.GeneSymbol.KOPT_K1)
#      Adicion de la nueva columna obtenida
toptable.annotated_KOPT_K1 <- cbind(subgrupo_KOPT_K1_pvalor,GeneSymbol_KOPT_K1)
#      Visualizacion de los primeros resultados
head(toptable.annotated_KOPT_K1)



## 4. Eliminación de sondas que no tengan un gen asociado (en GeneSymbol pone NA)
nueva_toptable_KOPT_K1 <- droplevels(toptable.annotated_KOPT_K1[-which(
  toptable.annotated_KOPT_K1$GeneSymbol_KOPT_K1 == "NA"), ] )
dim(nueva_toptable_KOPT_K1)



## 5. Si los genes se repiten, se seleccionan aquellos con un pvalor ajustado mejor
filas_KOPT_K1_menor_pvalor <- lapply(split(nueva_toptable_KOPT_K1,nueva_toptable_KOPT_K1$GeneSymbol_KOPT_K1),
                function(x) {return(x[which.min(x$adj.P.Val),])})
df_DEG_menor_pvalor_KOPT_K1 <- rbindlist(filas_KOPT_K1_menor_pvalor)
df_DEG_menor_pvalor_KOPT_K1[,Sondas := unlist(lapply(filas_KOPT_K1_menor_pvalor , rownames))]
dim(df_DEG_menor_pvalor_KOPT_K1)



## 6. Obtención de las intensidades de las sondas filtradas en el paso anterior
filas_interes <- df_DEG_menor_pvalor_KOPT_K1$Sondas
df_intensidades_KOPT_K1 <- data.frame()
for (i in filas_interes){
  fila <- (which(rownames(expression_total_KOPT_K1) == i))
  filas_valores <- exprs(esetIQR_KOPT_K1[fila,])
  df_intensidades_KOPT_K1 <- rbind(df_intensidades_KOPT_K1,filas_valores)
  dim(df_intensidades_KOPT_K1)
}

# Cambio del nombre de columnas (sondas -> genes)
row.names(df_intensidades_KOPT_K1) <- df_DEG_menor_pvalor_KOPT_K1$GeneSymbol_KOPT_K1
dim(df_intensidades_KOPT_K1)



## 7. Realización de heatmap
pheatmap(df_intensidades_KOPT_K1)



## 8. Guardar el heatmap en una imagen
png("Heatmap_KOPT_K1.png", width=1000, height=1000)
pheatmap(df_intensidades_KOPT_K1)
dev.off()

