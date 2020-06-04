library(tidyr)
data <- iris
## Agreggate data and calculate mean
data_mean <- aggregate(data[, 1:4], list(data$Species), mean)
clean_mean <- data_mean %>% pivot_longer(cols = c("Sepal.Length", "Sepal.Width",
                                                         "Petal.Length","Petal.Width"),
                                                names_to = "Character")
## Agreggate data and calculate standard deviation
data_sd <- aggregate(data[, 1:4], list(data$Species), sd)
clean_sd <- data_sd %>% pivot_longer(cols = c("Sepal.Length", "Sepal.Width",
                                                  "Petal.Length","Petal.Width"),
                                         names_to = "Character")
## Bind into one dataframe
statistics <- cbind(clean_mean,clean_sd[,3])
colnames(statistics) <- c("Species","Character","Mean","SD")

## Plot mean with error bar
library(ggplot2)
p <-ggplot(statistics, aes(x=Character, y=Mean, fill=Species)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))
p     

## Plot dispersion between petal length and width
petal <- ggplot(data, aes(x=Petal.Length, y= Petal.Width, color=Species)) + 
  geom_point()
petal
## Plot dipsersion between sepal length and width
sepal <- ggplot(data, aes(x=Sepal.Length, y= Sepal.Width, color=Species)) + 
  geom_point()
sepal
