## 0. Establecimiento del directorio de trabajo
setwd("~/Documentos/Master_Lucia/Segundo_cuatri/Transcriptomica/Transcriptomica/Trabajo_Transcriptomica/GSE18198_RAW")


## Cargar las librerías necesarias
library(pheatmap)
library("genefilter")
library("affy")
library("hgu133plus2.db")
library(data.table)
  
## ------------------------------ REALIZACIÓN DE HEATMAP ------------------------------ ##
## ------------------------------------------------------------------------------------ ##
## ------------------------------------- HPB_ALL -------------------------------------- ##
## ------------------------------------------------------------------------------------ ##

## 1. Cargar elementos
# Se carga el objeto eset
eset_HPB_ALL <- readRDS(file = "eset_HPB_ALL.rds")

# Se realiza el filtrado por el rango intercuartílico
esetIQR_HPB_ALL <- varFilter(eset_HPB_ALL, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
expression_total_HPB_ALL <- exprs(esetIQR_HPB_ALL)

# Se cargan la tabla de HPB_ALL con el pvalor y el logFC
toptableIQR_HPB_ALL <- readRDS(file = "toptableIQR_HPB_ALL.rds")

## 2. Filtrado de las sondas de interés de toptableIQR_HPB_ALL
subgrupo_HPB_ALL_pvalor <- subset(toptableIQR_HPB_ALL, toptableIQR_HPB_ALL$adj.P.Val <= 0.05 
            & (toptableIQR_HPB_ALL$logFC < -1))

subgrupo_HPB_ALL_pvalor
dim(subgrupo_HPB_ALL_pvalor) # En total, 17 sondas downreguladas para HPB_ALL (con LogFC < -1)


## 3. Adicion del simbolo a las sondas 
#      Pasar a caracteres los nombres de las sondas
probenames.HPB_ALL <- as.character(rownames(subgrupo_HPB_ALL_pvalor))
#      Obtencion de los simbolos asociados a las sondas
list.GeneSymbol.HPB_ALL <- mget(probenames.HPB_ALL, hgu133plus2SYMBOL)
GeneSymbol_HPB_ALL <- as.character(list.GeneSymbol.HPB_ALL)
#      Adicion de la nueva columna obtenida
toptable.annotated_HPB_ALL <- cbind(subgrupo_HPB_ALL_pvalor,GeneSymbol_HPB_ALL)
#      Visualizacion de los primeros resultados
head(toptable.annotated_HPB_ALL)


## 4. Eliminación de sondas que no tengan un gen asociado (en GeneSymbol pone NA)
# nueva_toptable_HPB_ALL <- droplevels(toptable.annotated_HPB_ALL[-which(
#  toptable.annotated_HPB_ALL$GeneSymbol_HPB_ALL == "NA"), ] )
# dim(nueva_toptable_HPB_ALL) 
nueva_toptable_HPB_ALL <- toptable.annotated_HPB_ALL 

## 5. Si los genes se repiten, se seleccionan aquellos con un pvalor ajustado mejor
filas_HPB_ALL_menor_pvalor <- lapply(split(nueva_toptable_HPB_ALL,nueva_toptable_HPB_ALL$GeneSymbol_HPB_ALL),
                function(x) {return(x[which.min(x$adj.P.Val),])})
df_DEG_menor_pvalor <- rbindlist(filas_HPB_ALL_menor_pvalor)
df_DEG_menor_pvalor[,Sondas := unlist(lapply(filas_HPB_ALL_menor_pvalor , rownames))]
dim(df_DEG_menor_pvalor)


## 6. Obtención de las intensidades de las sondas filtradas en el paso anterior
filas_interes <- df_DEG_menor_pvalor$Sondas
df_intensidades_HPB_ALL <- data.frame()
for (i in filas_interes){
  fila <- (which(rownames(expression_total_HPB_ALL) == i))
  filas_valores <- exprs(esetIQR_HPB_ALL[fila,])
  df_intensidades_HPB_ALL <- rbind(df_intensidades_HPB_ALL,filas_valores)
  dim(df_intensidades_HPB_ALL)
}
df_intensidades_HPB_ALL

# Cambio del nombre de columnas (sondas -> genes)
row.names(df_intensidades_HPB_ALL) <- df_DEG_menor_pvalor$GeneSymbol_HPB_ALL
head(df_intensidades_HPB_ALL)
dim(df_intensidades_HPB_ALL)

## 7. Realización de heatmap
pheatmap(df_intensidades_HPB_ALL)


## 8. Guardar el heatmap en una imagen
png("Heatmap_HPB_ALL.png", width=500, height=500)
pheatmap(df_intensidades_HPB_ALL)
dev.off()
