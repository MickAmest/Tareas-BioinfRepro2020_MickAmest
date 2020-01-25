### tarea 3.5
## 1. Generar un gráfico de puntos del data.frame iris, coloreando los 
## puntos por color
library(ggplot2)
ggplot(data=iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) + 
  geom_point()

## 2. Generar datos dat y generar un barplot, coloreado por especie

# Generar datos
set.seed(10)
dat <- data.frame(species=factor(rep(c("inventus", "otrus"), each=300)), 
                  size=c(rnorm(300, mean = 162), rnorm(300, mean=165)))

# Generar histograma
ggplot(data=dat, aes(x=size, fill= species)) + # En x va el tamaño para ver la frecuencia
  geom_histogram(binwidth=0.5, alpha=0.5, position = "identity") + # Para hacerlo translúcido por el sobrelape
  scale_fill_manual(values = (c("green3", "darkblue"))) + # cambiar manualmente los colores
  labs(fill= "Especie", x="Tamaño", y="Individuos") # cambiar el nombre de los ejes y la etiqueta de datos
    



