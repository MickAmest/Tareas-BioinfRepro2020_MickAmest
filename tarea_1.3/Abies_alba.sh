## Descargar cinco secuencias cortas desde la base de NCBI y buscar y contar
## el número de veces que se encuentra la secuencia TGCA en cada una de ellas

# Crear directorio
mkdir Abies_alba
# Descargar las secuencias
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&id=378878883,378878882,378878881,378878880,378878879" > Abies_alba/secuencias_Abies_alba.txt
# Recortar los nombres de las secuencias
sed -i "s/ /_/g" Abies_alba/secuencias_Abies_alba.txt
cut -d_ -f 1 Abies_alba/secuencias_Abies_alba.txt > Abies_alba/secuencias_Abies_alba_cut.txt
# Crear archivos separados para cada secuencia
cd Abies_alba
grep ">" secuencias_Abies_alba_cut.txt > id.txt
for i in $(cat id.txt); do grep -A8 "$i" secuencias_Abies_alba_cut.txt > $i.fasta; done
# Contar en cada archivo fasta el número de veces que aparece la secuencia indicada
grep -cH TGCA *.fasta
