### Análisis de datos de Illumina MiSeq desde amptk

library(phyloseq)
library(vegan)
library(ggplot2)

# Cargar los datos
suelo <- import_biom("taxonomy.biom")
head(tax_table(suelo))
colnames(tax_table(suelo)) <- c("Domain", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species")
# Todas las muestras deben tener el mismo número de reads
no.reads = sort(sample_sums(suelo))
no.reads

# Diversidad alfa
p <-  plot_bar(suelo, "Host", fill="Phylum") + geom_bar(aes(color="Phylum"), stat = "identity",
                                                       position = "stack")
p + facet_wrap("Treatment")

# Indices de diversidad
diversidad <- estimate_richness(suelo, measures = c("Observed","Fisher"))
