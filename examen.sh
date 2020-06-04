### Script para crear mil archivos llamados elefante1...elefante1000, cada uno
### de los cuales contiene el texto indicado.
## Miguel Amaro, Feb 2020

for i in {0..1000}; do echo "Hola mundo, este es el elefante numero $i" > elefante${i}.txt; done
