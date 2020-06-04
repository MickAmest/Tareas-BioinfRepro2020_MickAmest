### Tarea D. Script to plot iris data on separete boxes.
### Miguel Amaro, Mar 2020
library(ggplot2)
## See data
head(iris)
## Plot
ggplot(data=iris, aes(x= Sepal.Length, y=Sepal.Width)) +
  geom_point(aes(color=Species)) +
  facet_grid(cols = vars(Species)) + ## Vertical panels, by species.
  labs(x="Largo del sépalo", y="Ancho del sépalo")

