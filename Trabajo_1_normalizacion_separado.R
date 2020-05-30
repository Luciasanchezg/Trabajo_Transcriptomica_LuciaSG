###################################################################
# Nombre del Script           : Trabajo_1_Normalización_separado
# Descripción                 : Normalización de los datos de microarrays de tecnología Affymetrix y realización de análisis
#                               de expresión diferencial con dos líneas celulares de leucemia linfoblástica aguda de linfocitos
#                               T, por separado (HPB_ALL y KOPT_K1). Visualización de datos en crudo y tras la normalización.
# Autora                      : Lucía Sánchez García    
# Fecha última modificación   : 30 de junio de 2020
###################################################################

## 0. Establecimiento del directorio de trabajo
setwd("~/Documentos/GSE18198_RAW")


## 1. Cargar librerias
library("affy")
library("limma")
library("genefilter")
library("ggplot2")
library("RColorBrewer") # Librería para elección de paleta de color



## 2. Importar los ficheros presentes en targets.txt
targets <- readTargets("targets.txt", row.names="FileName")
# Tratamiento por separado de las dos líneas celulares (KOPT-K1 y HPB-ALL). Para ello, se 
# realizan dos subsets diferentes, en funcion del tipo celular presente en Classes
targets_HPB_ALL <- subset(targets,grepl("^(HPB-ALL)",targets$Classes))
targets_KOPT_K1 <- subset(targets,grepl("^(KOPT-K1)",targets$Classes))



## 3. Importe de intensidades --> procedentes de los arrays de Affymetrix (.CEL), en bruto 
# para cada uno de los subsets
data_HPB_ALL <- ReadAffy(filenames=targets_HPB_ALL$FileName) 
data_KOPT_K1 <- ReadAffy(filenames=targets_KOPT_K1$FileName)

## Ambos subsets muestran:
##   Numero de muestras -> 6 samples
##   Numero de genes -> 54675 genes



## -------------------------------------------------------------------------------- ##
## --------------------------------- PARA HPB_ALL --------------------------------- ##
## -------------------------------------------------------------------------------- ##
## 4.1. Visualización de datos en crudo (HISTOGRAMA)
hist(data_HPB_ALL[,1:6], lwd = 2, which = 'pm', col = brewer.pal(n = 12, name = "Paired") , 
     ylab = 'Densidad', xlab ='Log2 intensidades', main ='Histograma de datos en crudo (HPB_ALL)')

## 4.2. Normalizacion con RMA
# Generación del objeto eset (clase ExprSet), 
# eset_HPB_ALL <- expresso(data_HPB_ALL,
#                        bg.correct = TRUE, # Corrige el fondo
#                        bgcorrect.method="rma",
#                        normalize = TRUE, 
#                        normalize.method="quantiles", # Normaliza por metodo de cuartiles
#                        pmcorrect.method="pmonly", # Metodo de aplicacion de correccion PM y MM
#                        summary.method="medianpolish", # Resume los resultados en un solo 
                                                        # valor, la mediana
#                        verbose = TRUE,
#                          )

# Se guardan los datos de la normalizacion en un fichero .rds
# saveRDS(eset_HPB_ALL, file = "eset_HPB_ALL.rds")



## 5. Generación de boxplots antes y después de la normalización
# Definición de un vector "colores", utilizado en la generación del boxplot para
# pintar las muestras en función de la clase (DMSO o SAHM1)
colores <- c("aquamarine1", "bisque1")
colores.box <- rep(colores, each = 3)

# Boxplot de los datos en crudo
boxplot(data_HPB_ALL,
        main = "Boxplot antes de normalización (HPB_ALL)",
        ylab = "Valores de expresión",
        xaxt = "n", # Con esto se indica la ausencia de etiquetas en el eje OX
        cex.axis = 0.5, # Tamaño de los valores de los ejes
        cex.lab = 1, # Tama?o de las etiquetas de los ejes
        col = colores.box)

# Establecimiento de los nombres del eje OX, rotados en un angulo de 45º
text(seq_along(data_HPB_ALL), par("usr")[3] - 0.5, labels = targets_HPB_ALL$FileName, srt = 45, adj = 1, xpd = TRUE, cex = 0.70)		

# Boxplot de los datos normalizados
exprseset_HPB_ALL <- as.data.frame(exprs(eset_HPB_ALL))
boxplot(data.frame(exprseset_HPB_ALL), # Boxplot para representar la normalizacion de los datos
        main="Boxplot tras la normalización (HPB_ALL)",
        ylab = "Valores de expresión",
        xaxt = "n",
        cex.axis = 0.5,
        cex.lab = 1,
        col = colores.box)

# Establecimiento de los nombres del eje OX, rotados en un ángulo de 45º
text(seq_along(data_HPB_ALL), par("usr")[3] - 0.5, labels = targets_HPB_ALL$FileName, srt = 45, adj = 1, xpd = TRUE, cex = 0.70)



## 6.Filtrado de datos con IQR
# Selección de genes con rango intercuartílico (IQR) mayor de 0.5 --> mayor varianza
# Resultado: descarte de patrones planos
esetIQR_HPB_ALL <- varFilter(eset_HPB_ALL, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

# eset_HPB_ALL    --> assayData: 54675 features, 6 samples
# esetIQR_HPB_ALL --> assayData: 27337 features, 6 samples
# Se ha reducido el número total de genes (de 56675 -> 27337)



## ----------------- Análisis de expresión diferencial para HPB_ALL ----------------- ##
## 7. Diseño de la matriz
design<-cbind(DMSO=c(1,1,1,0,0,0), SAHM1=c(0,0,0,1,1,1))
rownames(design) <- targets_HPB_ALL$FileName # Cambia los nombres de las filas a lo de
                                             # los ficheros (intuitivo)


## 8. Matriz de contraste
# makeContrasts prepara contrastes estadísticos, comparando SAHM1 frente a DMSO
cont.matrix_HPB_ALL <- makeContrasts(SAHM1vsDMSO = SAHM1-DMSO, levels = design) 



## 9. Obtención de genes diferencialmente expresados (DEGs)
# Ajuste de modelo lineal y realización de Bayes 
fit_HPB_ALL <- lmFit(esetIQR_HPB_ALL,design)  # Obtención de los DEGs con IQR 
fit2_HPB_ALL <- contrasts.fit(fit_HPB_ALL, cont.matrix_HPB_ALL) # Generación del contraste
fit2_HPB_ALL <- eBayes(fit2_HPB_ALL) # Realización el método bayesiano

# Tabla de los resultados de DEGs
toptableIQR_HPB_ALL <- topTable(fit2_HPB_ALL, number=dim(exprs(esetIQR_HPB_ALL))[1], adjust.method="BH", sort.by="p")
head(toptableIQR_HPB_ALL)



## 10. Guardar resultados
# save(toptableIQR_HPB_ALL,file="MyResults_HPB_ALL.RData")
saveRDS(toptableIQR_HPB_ALL, file = "toptableIQR_HPB_ALL.rds")



## 11. Creación de volcanoplot para la visualización de los datos
threshold_OE <- toptableIQR_HPB_ALL$adj.P.Val < 0.05
length(which(threshold_OE))
toptableIQR_HPB_ALL$p_valor_005 <- threshold_OE 

ggplot(toptableIQR_HPB_ALL) +
        geom_point(aes(x=logFC, y=-log10(adj.P.Val),colour=p_valor_005), size = 0.2) +
        geom_vline(xintercept=1, linetype = "longdash", color = "dimgrey") +
        geom_vline(xintercept=-1, linetype = "longdash", color = "dimgrey") +
        geom_vline(xintercept=0.5, linetype = "longdash", color = "darkgrey") +
        geom_vline(xintercept=-0.5, linetype = "longdash", color = "darkgrey") +
        geom_hline(yintercept=-log10(0.05), linetype = "longdash", color = "grey") +
        xlim(-3.3,3.3) +
        ggtitle("Volcanoplot HPB_ALL (Control Vs Tratamiento)") +
        xlab("Log2 Fold Change") + 
        ylab("-Log10 p-valor ajustado") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1), hjust = 0.5),
              axis.title = element_text(size = rel(0.9)))  



## -------------------------------------------------------------------------------- ##
## --------------------------------- PARA KOPT_K1 --------------------------------- ##
## -------------------------------------------------------------------------------- ##
## 4.1. Visualización de datos en crudo (HISTOGRAMA)
hist(data_KOPT_K1[,1:6], lwd = 2, which = 'pm', col = brewer.pal(n = 12, name = "Paired") , 
     ylab = 'Densidad', xlab ='Log2 intensidades', main ='Histograma de datos en crudo (KOPT_K1)')


## 4.2. Normalización con RMA
# eset_KOPT_K1 <- expresso(data_KOPT_K1,
#                         bg.correct = TRUE,
#                         bgcorrect.method="rma",
#                         normalize = TRUE, 
#                         normalize.method="quantiles",
#                         pmcorrect.method="pmonly",
#                         summary.method="medianpolish",
#                         verbose = TRUE,
#                          ) 

# Se guardan los datos de la normalizacion en un fichero .rds
# saveRDS(eset_KOPT_K1, file = "eset_KOPT_K1.rds")



## 5. Generación de Boxplots antes y después de la normalización
colores <- c("coral", "darkolivegreen1")
colores.box <- rep(colores, each = 3)

# Boxplot de los datos en crudo
boxplot(data_KOPT_K1,
        main = "Boxplot antes de la normalización (KOPT_K1)",
        ylab = "Valores de expresión",
        xaxt = "n", 
        cex.axis = 0.5, 
        cex.lab = 1,
        col = colores.box)

# Establecimiento de los nombres del eje OX, rotados en un ángulo de 45º
text(seq_along(data_KOPT_K1), par("usr")[3] - 0.5, labels = targets_KOPT_K1$FileName, srt = 45, adj = 1, xpd = TRUE, cex = 0.70)		

# Boxplot de los datos normalizados
exprseset_KOPT_K1 <- as.data.frame(exprs(eset_KOPT_K1))
boxplot(data.frame(exprseset_KOPT_K1),
        main="Boxplot tras la normalización (KOPT_K1)",
        ylab = "Valores de expresión",
        xaxt = "n",
        cex.axis = 0.5,
        cex.lab = 1,
        col = colores.box)

# Establecimiento de los nombres del eje OX, rotados en un ángulo de 45º
text(seq_along(data_KOPT_K1), par("usr")[3] - 0.5, labels = targets_KOPT_K1$FileName, srt = 45, adj = 1, xpd = TRUE, cex = 0.70)



## 6.Filtrado de datos con IQR
esetIQR_KOPT_K1 <- varFilter(eset_KOPT_K1, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

# eset_HPB_ALL    --> assayData: 54675 features, 6 samples
# esetIQR_HPB_ALL --> assayData: 27337 features, 6 samples
# Mismos resultados que con HPB_ALL



## ----------------- Análisis de expresión diferencial para KOPT_K1 ----------------- ##
## 7. Diseño de la matriz
design<-cbind(DMSO=c(1,1,1,0,0,0), SAHM1=c(0,0,0,1,1,1))
rownames(design) <- targets_KOPT_K1$FileName # Cambia los nombres de las filas a lo de los 
                                             # ficheros (intuitivo)


## 8. Matriz de constrate
# makeContrasts prepara contrastes estadísticos, comparando SAHM1 frente a DMSO
cont.matrix_KOPT_K1 <- makeContrasts(SAHM1vsDMSO = SAHM1-DMSO, levels = design) 



## 9. Obtención de genes diferencialmente expresados (DEGs)
# Ajuste de modelo lineal y realización de Bayes 
fit_KOPT_K1 <- lmFit(esetIQR_KOPT_K1,design)  # Obtención de los DEGs con IQR 
fit2_KOPT_K1 <- contrasts.fit(fit_KOPT_K1, cont.matrix_KOPT_K1) # Generación del contraste
fit2_KOPT_K1 <- eBayes(fit2_KOPT_K1) # Realización el método bayesiano

#Tabla de genes DEGs
toptableIQR_KOPT_K1 <- topTable(fit2_KOPT_K1, number=dim(exprs(esetIQR_KOPT_K1))[1], adjust.method="BH", sort.by="p")
head(toptableIQR_KOPT_K1)



## 10. Guardar resultados
# save(toptableIQR_KOPT_K1,file="MyResults_KOPT_K1.RData")
saveRDS(toptableIQR_KOPT_K1, file = "toptableIQR_KOPT_K1.rds")



## 11. Creación de volcanoplot para la visualización de los datos
threshold_OE <- toptableIQR_KOPT_K1$adj.P.Val < 0.05
length(which(threshold_OE))
toptableIQR_KOPT_K1$p_valor_005 <- threshold_OE 

ggplot(toptableIQR_KOPT_K1) +
        geom_point(aes(x=logFC, y=-log10(adj.P.Val),colour=p_valor_005), size = 0.2) +
        geom_vline(xintercept=1, linetype = "longdash", color = "dimgrey") +
        geom_vline(xintercept=-1, linetype = "longdash", color = "dimgrey") +
        geom_vline(xintercept=0.5, linetype = "longdash", color = "darkgrey") +
        geom_vline(xintercept=-0.5, linetype = "longdash", color = "darkgrey") +
        geom_hline(yintercept=-log10(0.05), linetype = "longdash", color = "grey") +
        xlim(-3,3) +
        ggtitle("Volcanoplot KOPT_K1 (Control Vs Tratamiento)") +
        xlab("Log2 Fold Change") + 
        ylab("-Log10 p-valor ajustado") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1), hjust = 0.5),
              axis.title = element_text(size = rel(0.9)))  
