#### TAREA 6.1
#### Marzo, 2020

## Contar con comandos de Unix el número de reads en cada archivo fastq
grep -o @M S19_R1_filter.fastq | wc
grep -o @M S19_R2_filter.fastq | wc
## Visualizar las 40 primeras líneas del fastq que corresponden a las 10
## primeras muestras
head -40 S19_R1_filter.fastq

## Generar reporte fastqc para las muestras crudas
fastqc S19_R1.fastq -o .
fastqc S19_R2.fastq -o .

## Generar reportes fastqc para las muestras filtradas
fastqc S19_R1_filter.fastq -o .
fastqc S19_R2_filter.fastq -o .
