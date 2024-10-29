BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)


# El dataset elegido es el del artículo "Metabotypes of response to bariatric surgery independent of the magnitude of weight loss". En el que se miden los metabolitos en sangre de pacientes con obesidad morbida sometidos a cirugia bariatrica. 
# Con estos podríamos saber que cambios metabólicos experimentan pacientes con este tipo de cirugia.

datos <- read.csv("metaboData/Datasets/2018-MetabotypingPaper/DataValues_S013.csv")
datos_variables <- read.csv("metaboData/Datasets/2018-MetabotypingPaper/DataInfo_S013.csv")

colnames(datos)
colnames(datos_variables)

#Explorando el archivo DataValues, veo que contiene tanto los metadatos de los pacientes como las mediciones, por otro lado DataInfo contiene información sobre las variables del estudio. Por ello, es necesario primero separar los metadatos de las mediciones del archivo "valores"

#Primero separamos los metadatos de los pacientes y los transformamos en el tipo "DataFrame", que es el reconocido por SummarizedExperiment
metadatos <- datos[,c("SUBJECTS", "SURGERY", "AGE", "GENDER", "Group")]
metadatos_df <- as(metadatos, "DataFrame")

datos <- subset(datos, select = -c(SUBJECTS, SURGERY, AGE, GENDER, Group,X))

#Las mediciones tienen bastante valores NA, para que no den problema después podemos realizar una imputación por kNN
library("VIM")

datos_imputados <- kNN(datos, k = 7, imp_var = FALSE)


#Además, las mediciones están separadas en las tomas en varios momentos, T0, T2, T4 y T5, en función del mes en el que se tomo la muestra, antes de la cirugia, al mes, a los 3 o a los 6. Por ello, aprovechando que Summarized experiment puede contener más de una matriz de mediciones, creo que es mejor separarlas.
#Dentro del bucle selecciono las columnas para cada T, para cada matriz, le quito del nombre el tiempo (ya se sabe por la matriz en la que está contenida) y la transformo en la matriz transpuesta (Se requieren las mediciones en las filas y las observaciones en las columnas).


for (t in c("T0","T2","T4","T5")){
  mediciones_t <- datos_imputados[,grepl(t,colnames(datos))]
  print(dim(mediciones_t))
  colnames(mediciones_t) <- gsub(t,"",colnames(mediciones_t))
  assign(paste0("mediciones_", t), t(as.matrix(mediciones_t)))
}

#El T2 tiene una fila de más, y todas deberían tener las mismas, la buscamos y la eliminamos.

print(setdiff(rownames(mediciones_T2),rownames(mediciones_T0)))

mediciones_T2 <- mediciones_T2[rownames(mediciones_T2)!="lysoPC.a.C14.0_",]

#Seleccionaremos uno de los tiempos, T0, por ejemplo y luego le quitaremos el nombre, para quedarnos con las mismas variables que en los archivos de mediciones. Después lo transformamos en un archivo DataFrame también.

datos_variables <- datos_variables[(grepl("T0",datos_variables$VarName)),]
datos_variables$VarName <- gsub("T0","",datos_variables$VarName)
datos_variables_df <- as(datos_variables,"DataFrame")


#Ahora podemos crear el contenedor SummarizedExperiment con los datos de las mediciones, los metadatos y la información sobre las variables.

se <- SummarizedExperiment(
  assays = list(Tprev = mediciones_T0,
                T1 = mediciones_T2,
                T3 = mediciones_T4,
                T6 = mediciones_T5),
  colData = metadatos_df,
  rowData = datos_variables_df)

#Ahora que tenemos el contenedor, podemos pasar a explorarlo.

#Podemos acceder al metadata del contenedor con colData()

head(colData(se))
dim(colData(se))

#La información de las variables está almacenada dentro de rowData()

head(rowData(se))

#La ventaja de SummarizedExperiment vs Expression sets es que SummarizedExperiment puede almacenar más de una matriz de mediciones, podemos ver todas con assays().

assays(se)

#Podriamos acceder a los datos de un ensayo en concreto

head(assay(se,"Tprev"),2)



#Los valores son recogidos en distintos tiempos, antes de la cirugia y a los meses de esta, por lo que puede ser más interesante ver la evolución. Podemos añadir una nueva matriz (en este caso la evolución entre antes de la cirugia y el ultimo mes del ensayo) al contenedor así:

assays(se)[["Tprev_T6"]] <- assay(se,"T6")-assay(se,"Tprev")


#Al ser tantisimas variables, podríamos hacer un análisis de componentes principales para ver la fuente de variabilidad de los datos:

pca_T6_Tprev <- prcomp(t(assay(se,"Tprev_T6")),scale. = TRUE)

summary(pca_T6_Tprev)

#Podemos observar que el primer y segundo componente explican el 34,4 y 12,4% de la variabilidad respectivamente, necesitando 20 componentes para explicar el 90% de la variabilidad de la muestra.

plot(pca_T6_Tprev$x[, 1:2], col = factor(colData(se)$SURGERY), pch = 19, 
     main = "PCA: Primeros dos componentes T6-T0", xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(factor(colData(se)$SURGERY)), 
       col = 1:length(levels(factor(colData(se)$SURGERY))), pch = 19)

#La mayor parte de la variabilidad está explicada por el primer componente principal, en el eje del segundo componente podemos observar dos outliers.

# Otra forma de hacer un análisis exploratorio de los datos es mediante un clustering de las observaciones.

dist_T6_Tprev <- dist(t(assay(se,"Tprev_T6")))
clust_T6_Tprev <- hclust(dist_T6_Tprev)
plot(clust_T6_Tprev,labels = factor(colData(se)$SURGERY), main = "Dendograma T6-T0")


#Podemos hacer estas exploraciones para las distintas mediciones. Ocultaré el código en este caso para no hacer el informe muy largo.

#OCULTAR CÓDIGO Y ESTA FRASE

#T3-T0 
assays(se)[["Tprev_T3"]] <- assay(se,"T3")-assay(se,"Tprev")
data_Tprev_T3 <- assay(se,"Tprev_T3")
data_Tprev_T3 <- data_Tprev_T3[, apply(data_Tprev_T3, 2, var) != 0]

pca_T3_Tprev <- prcomp(data_Tprev_T3,scale. = TRUE)
summary(pca_T3_Tprev)
plot(pca_T3_Tprev$x[, 1:2], col = factor(colData(se)$SURGERY), pch = 19, 
     main = "PCA: Primeros dos componentes T3 - T0", xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(factor(colData(se)$SURGERY)), 
       col = 1:length(levels(factor(colData(se)$SURGERY))), pch = 19)
dist_T3_Tprev <- dist(t(assay(se,"Tprev_T3")))
clust_T3_Tprev <- hclust(dist_T3_Tprev)
plot(clust_T3_Tprev,labels = factor(colData(se)$SURGERY), main = "Dendograma T3-T0")


#T1-T0 
assays(se)[["Tprev_T1"]] <- assay(se,"T1")-assay(se,"Tprev")
data_Tprev_T1 <- assay(se,"Tprev_T1")
data_Tprev_T1 <- data_Tprev_T1[, apply(data_Tprev_T1, 2, var) != 0]

pca_T1_Tprev <- prcomp(data_Tprev_T1,scale. = TRUE)
summary(pca_T1_Tprev)
plot(pca_T1_Tprev$x[, 1:2], col = factor(colData(se)$SURGERY), pch = 19, 
     main = "PCA: Primeros dos componentes T1 - T0", xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(factor(colData(se)$SURGERY)), 
       col = 1:length(levels(factor(colData(se)$SURGERY))), pch = 19)
dist_T1_Tprev <- dist(t(assay(se,"Tprev_T1")))
clust_T1_Tprev <- hclust(dist_T1_Tprev)
plot(clust_T1_Tprev,labels = factor(colData(se)$SURGERY), main = "Dendograma T1-T0")

#Parece que la cirugia no explica el mayor porcentaje de la variabilidad de los datos en las distintas mediciones. 


#Para guardar el objeto contenedor con los datos y metadatos en formato binario usamos la función save

save(se,file = "SummarizedExperiment_PSU.rda")

#Para crear el archivo markdown con los metadatos primero extraemos los metadatos de las muestras y las variables:

metadatos_muestras <- colData(se)

metadatos_variables <- rowData(se)


# Ahora creamos el archivo markdown

metadatos <- file("metadatos_dataset.md")

# Escribir la cabecera del archivo Markdown
writeLines(c("# Metadatos del Dataset", 
             "## Metadatos de las Muestras", 
             "```", 
             capture.output(print(metadatos_muestras)),  # Imprimir los metadatos de las muestras
             "```", 
             "## Metadatos de las Variables", 
             "```", 
             capture.output(print(metadatos_variables)),  # Imprimir los metadatos de las variables
             "```"), metadatos)

# Cerrar el archivo
close(metadatos)

#Por último, creamos el archivo de texto con los datos
write.table(datos_imputados, file = "Texto_datos.txt", sep = "\t", quote = FALSE)
