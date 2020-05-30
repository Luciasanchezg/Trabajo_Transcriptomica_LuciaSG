###################################################################
# Nombre del Script           : Trabajo_5_preparacion_GSEA
# Descripción                 : Generación de un archivo csv con las intensidades de las muestras y controles para las líneas
#                               celulares HPB_ALL y KOPT_K1
# Autora                      : Lucía Sánchez García    
# Fecha última modificación   : 30 de mayo de 2020
###################################################################

## 0. Establecimiento del directorio de trabajo
setwd("~/Documentos/GSE18198_RAW")


## --------------------------------------------------------------------------------- ##
## ------------------------------------ HPB_ALL ------------------------------------ ##
## --------------------------------------------------------------------------------- ##
## Importe de los datos de expresión y generación de dataframe para HPB_ALL
eset_HPB_ALL <- readRDS(file = "eset_HPB_ALL.rds")
df_exp_HPB_ALL <- as.data.frame(exprs(eset_HPB_ALL))
dim(df_exp_HPB_ALL) # El dataframe contiene todas las sondas (54675 en total)

# Escritura del archivo con las intensidades normalizadas para HPB_ALL
write.csv(df_exp_HPB_ALL,file = "expression_HPB_ALL.csv", quote = FALSE)



## --------------------------------------------------------------------------------- ##
## ------------------------------------ KOPT_K1 ------------------------------------ ##
## --------------------------------------------------------------------------------- ##
## Importe de los datos de expresión y generación de dataframe para KOPT_K1 
eset_KOPT_K1 <- readRDS(file = "eset_KOPT_K1.rds")
df_exp_KOPT_K1 <- as.data.frame(exprs(eset_KOPT_K1))


# Escritura del archivo con las intensidades normalizadas para HPB_ALL
write.csv(df_exp_KOPT_K1,file = "expression_KOPT_K1.csv", quote = FALSE)
  
