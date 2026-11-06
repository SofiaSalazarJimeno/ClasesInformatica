##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################

getwd()

### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)



##########################################################################################
#####################                                                        ##################
#####################       TRABAJO ORDINARIA BIOINFORMATICA R               ##################
#####################                                                        ##################
##########################################################################################

### 1. CARGA Y EXPLORACION DEL ARCHIVO

# Comprobamos el directorio de trabajo.
getwd()

# Cargamos el archivo de entrada.
datos <- read.delim("archivo-entrada.txt",
                    header = TRUE,
                    sep = "\t",
                    dec = ".")

# Examinamos los primeros y ultimos datos.
head(datos)
tail(datos)

# Comprobamos las dimensiones.
dim(datos)

# Comprobamos los nombres de las columnas.
names(datos)


### 2. CALCULO DE MEDIA Y DESVIACION ESTANDAR

# Seleccionamos las tres replicas del Control.
control <- datos[, c("C2.TPM", "C1.TPM", "C3.TPM")]

# Seleccionamos las dos replicas TPM de la Condicion.
condicion <- datos[, c("X1W5.TPM", "X1W6.TPM")]

# Calculamos las medias de cada grupo.
media_control <- apply(control, 1, mean, na.rm = TRUE)
media_condicion <- apply(condicion, 1, mean, na.rm = TRUE)

# Calculamos las desviaciones estandar.
DE_control <- apply(control, 1, sd, na.rm = TRUE)
DE_condicion <- apply(condicion, 1, sd, na.rm = TRUE)

# Añadimos las cuatro columnas a la tabla original.
datos$Media_Control <- media_control
datos$DE_Control <- DE_control
datos$Media_Condicion <- media_condicion
datos$DE_Condicion <- DE_condicion

# Comprobamos los resultados.
head(datos)


### 3. HISTOGRAMAS DE LAS MEDIAS

# Definimos los colores.
color_control <- "skyblue"
color_condicion <- "salmon"

# Mostramos los dos histogramas en una misma imagen.
par(mfrow = c(1, 2))

hist(datos$Media_Control,
     main = "Media del Control",
     xlab = "Expresion genica",
     ylab = "Frecuencia",
     col = color_control,
     border = "black")

hist(datos$Media_Condicion,
     main = "Media de la Condicion",
     xlab = "Expresion genica",
     ylab = "Frecuencia",
     col = color_condicion,
     border = "black")

par(mfrow = c(1, 1))


### 4. BOXPLOT DE CONTROL Y CONDICION

boxplot(list(Control = datos$Media_Control,
             Condicion = datos$Media_Condicion),
        main = "Comparacion de expresion genica",
        ylab = "Expresion genica",
        col = c(color_control, color_condicion),
        border = "black",
        names = c("Control", "Condicion"))

legend("topleft",
       legend = c("Control", "Condicion"),
       fill = c(color_control, color_condicion),
       title = "Grupos")


### 5. EVALUACION DE LA NORMALIDAD Y COMPARACION

# Test de Shapiro-Wilk para evaluar la normalidad.
normalidad_control <- shapiro.test(datos$Media_Control)
normalidad_condicion <- shapiro.test(datos$Media_Condicion)

# Mostramos los resultados.
normalidad_control
normalidad_condicion

# Si p < 0.05, rechazamos la hipotesis de normalidad.
# Si p >= 0.05, no rechazamos la hipotesis de normalidad.

# Comparamos cada gen mediante el test de Wilcoxon-Mann-Whitney,
# ya que las distribuciones de las medias no siguen una distribucion normal.

# Creamos un vector para guardar los p-valores.
p_valores <- numeric(nrow(datos))

# Realizamos la comparacion para cada gen.
for (i in 1:nrow(datos)) {

  x <- as.numeric(control[i, ])
  y <- as.numeric(condicion[i, ])

  p_valores[i] <- wilcox.test(x, y,
                              paired = FALSE,
                              exact = FALSE)$p.value
}

# Añadimos los p-valores a la tabla original.
datos$p_value <- p_valores

# Comprobamos los resultados.
head(datos)


### 6. IDENTIFICACION DE GENES DIFERENCIALMENTE EXPRESADOS

# Seleccionamos los genes con p-valor inferior a 0.05.
genes_diferenciales <- datos[datos$p_value < 0.05, ]

# Mostramos los genes diferencialmente expresados.
genes_diferenciales

# Calculamos el numero total de genes diferenciales.
numero_genes <- nrow(genes_diferenciales)

print(paste("Numero de genes diferenciales:", numero_genes))

# Identificamos los genes sobreexpresados.
sobreexpresados <- genes_diferenciales[
  genes_diferenciales$Media_Condicion >
    genes_diferenciales$Media_Control, ]

# Identificamos los genes reprimidos.
reprimidos <- genes_diferenciales[
  genes_diferenciales$Media_Condicion <
    genes_diferenciales$Media_Control, ]

# Mostramos los resultados.
print(paste("Genes sobreexpresados:", nrow(sobreexpresados)))
print(paste("Genes reprimidos:", nrow(reprimidos)))


### 7. MAPA DE CALOR (HEATMAP)

# Instalamos pheatmap si no esta instalado.
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# Comprobamos si existen genes diferencialmente expresados.
if (nrow(genes_diferenciales) > 0) {

  # Seleccionamos las columnas de expresion.
  matriz <- as.matrix(genes_diferenciales[
    , c("C2.TPM", "C1.TPM", "C3.TPM",
        "X1W5.TPM", "X1W6.TPM")])

  # Asignamos los nombres de los genes a las filas.
  rownames(matriz) <- genes_diferenciales$Gene.Name

  # Transformamos los valores mediante log2.
  matriz <- log2(matriz + 1)

  # Representamos el mapa de calor.
  pheatmap(matriz,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = "Genes diferencialmente expresados",
           color = colorRampPalette(c("blue", "white", "red"))(100))

} else {

  print("No se han encontrado genes diferencialmente expresados.")

}


