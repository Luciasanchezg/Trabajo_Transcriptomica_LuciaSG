###################################################################
# Nombre del Script           : Trabajo_4_PCA
# Descripción                 : Realización de análisis de las componentes principales y visualización de las dos y tres
#                               componentes principales, para ambas líneas celulares.
# Autora                      : Lucía Sánchez García    
# Fecha última modificación   : 30 de mayo de 2020
###################################################################

## 0. Establecimiento del directorio de trabajo
setwd("~/Documentos/GSE18198_RAW")

## Cargar las librerías necesarias
library("genefilter")
library("affy")
library("limma")
library("scatterplot3d")

## ------------------------------------------------------------------------------------ ##
## ------------------------------------- HPB_ALL -------------------------------------- ##
## ------------------------------------------------------------------------------------ ##
## 1. Procesado de datos
# Lectura del fichero
eset_HPB_ALL <- readRDS(file = "eset_HPB_ALL.rds")
expression_total_HPB_ALL <- exprs(eset_HPB_ALL)

# Realización del análisis de componentes principales
expression_HPB_ALL.PC = prcomp(t(expression_total_HPB_ALL),scale.=TRUE)
PCAcomp_HPB_ALL <- expression_HPB_ALL.PC$x
head(PCAcomp_HPB_ALL)

# Elección de colores
colores <- c("cornflowerblue", "chocolate1")
colores.box <- rep(colores, each = 3)



## 2. PCA con las dos primeras componentes
# Localización de las muestras en una espacio 2D
plot(PCAcomp_HPB_ALL[,1],PCAcomp_HPB_ALL[,2],ylab="PC2",xlab="PC1",main=c("PCA 2D HPB_ALL"),
     pch=16,col=colores.box)

# División del espacio
abline(h=0,v=mean(PCAcomp_HPB_ALL[,1],h=0,col="grey"))

# Adición de etiquetas a las muestras
text(PCAcomp_HPB_ALL[,1],PCAcomp_HPB_ALL[,2],label=rownames(PCAcomp_HPB_ALL),font=0,
     col=colores.box,
     cex=0.75,pos=ifelse(PCAcomp_HPB_ALL[,1]<mean(PCAcomp_HPB_ALL[,1]),yes=4,no=2))



## 3. PCA con las tres primeras componentes
# Localización de las muestras en un espacio 3D
plot3d<-scatterplot3d(PCAcomp_HPB_ALL[,1],PCAcomp_HPB_ALL[,2],PCAcomp_HPB_ALL[,3],
                      bg="black", ylab="PC2", xlab="PC1",zlab="PC3",main=c("PCA 3D HPB_ALL"),
                      pch=20,type="h",angle=60, color = colores.box)

# Adición de las etiquetas a las muestras
text(plot3d$xyz.convert(PCAcomp_HPB_ALL[,1],PCAcomp_HPB_ALL[,2],PCAcomp_HPB_ALL[,3]),
     label=rownames(PCAcomp_HPB_ALL),font=0,col=colores.box,cex=0.75,
     pos=ifelse(PCAcomp_HPB_ALL[,1] < mean(PCAcomp_HPB_ALL[,1]), yes=4, no=2))



## ------------------------------------------------------------------------------------ ##
## ------------------------------------- KOPT_K1 -------------------------------------- ##
## ------------------------------------------------------------------------------------ ##
## 1. Procesado de datos
# Lectura del fichero
eset_KOPT_K1 <- readRDS(file = "eset_KOPT_K1.rds")
expression_total_KOPT_K1 <- exprs(eset_KOPT_K1)

# Realización de análisis de componentes principales
expression_KOPT_K1.PC = prcomp(t(expression_total_KOPT_K1),scale.=TRUE)
PCAcomp_KOPT_K1 <- expression_KOPT_K1.PC$x
head(PCAcomp_KOPT_K1)

# Elección de colores
colores <- c("coral", "darkolivegreen")
colores.box <- rep(colores, each = 3)



## 2. PCA con las dos primeras componentes
# Localización de las muestras en una espacio 2D
plot(PCAcomp_KOPT_K1[,1],PCAcomp_KOPT_K1[,2],ylab="PC2",xlab="PC1",main=c("PCA 2D KOPT_K1"),
     pch=16,col=colores.box)

# División del espacio
abline(h=0,v=mean(PCAcomp_KOPT_K1[,1],h=0,col="grey"))

# Adición de etiquetas a las muestras
text(PCAcomp_KOPT_K1[,1],PCAcomp_KOPT_K1[,2],label=rownames(PCAcomp_KOPT_K1),font=0,
     col=colores.box,
     cex=0.75,pos=ifelse(PCAcomp_KOPT_K1[,1]<mean(PCAcomp_KOPT_K1[,1]),yes=4,no=2))



## 3. PCA con las tres primeras componentes
# Localización de las muestras en un espacop 3D
plot3d<-scatterplot3d(PCAcomp_KOPT_K1[,1],PCAcomp_KOPT_K1[,2],PCAcomp_KOPT_K1[,3],
                      bg="black", ylab="PC2", xlab="PC1",zlab="PC3",main=c("PCA 3D KOPT_K1"),
                      pch=20,type="h",angle=60, color = colores.box)

# Adición de etiquetas a las muestras
text(plot3d$xyz.convert(PCAcomp_KOPT_K1[,1],PCAcomp_KOPT_K1[,2],PCAcomp_KOPT_K1[,3]),
     label=rownames(PCAcomp_KOPT_K1),font=0,col=colores.box,cex=0.75,
     pos=ifelse(PCAcomp_KOPT_K1[,1] < mean(PCAcomp_KOPT_K1[,1]), yes=4, no=2))

