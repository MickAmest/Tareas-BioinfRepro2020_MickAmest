# Phylogeny of conifer CDPK genes
##### J. Miguel Amaro Estrada

### Retrieving gene sequences

The gene sequences of the CDPK family in conifers were obtained by
performing BLAST searches on the assembled transcriptomes of severeal
conifer species. The BLAST searches were performed locally to optimize
the time needed for the search. The sequences used as query, were the 14
sequences annotated as CDPK genes in the *Abies sachalinensis* transcriptome from
[TodoFirGene](http://plantomics.mind.meiji.ac.jp/todomatsu/help.html).

For this, the **blast+** utility has to be installed on the
machine. This can be done in Linux systems by running
`sudo apt-get install ncbi-blast+`. Alternatively, the executables can
be downloaded from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
to install blast+ on your preferred system.
It is also recommended that, if you have enough space on your machine, you
download the transcriptome assemblies and create local blast databases to
perform the search. Since they will be small databases, it won't take many
compututational resources to perform the search and you can even loop through
severeal transcriptomes on each search.

### Integrity and identification of the sequences
All the resulting sequences from the BLAST search have to be analyzed to find their longest
ORF. For this, the **ORFfinder** utility was used, which can be downloaded from the
[NCBI ORFfinder website](https://www.ncbi.nlm.nih.gov/orffinder/). Note that this
executable **only runs on Linux systems**. When running ORFfinder, make sure to use the appropiate relative path to the directory where you keep the executable. Alternatively, you can add this location to your `PATH` so it can be run from anywhere on your machine. Any sequence which longer ORF was shorter
than 1000 bp was eliminated for further analyses.

For the verification of the presence of the kinase and EF-hand domains on the sequences, InterPro Scan v.5 (Jones et al., 2014) was used, which can be downloaded from the
[InterPro website](https://www.ebi.ac.uk/interpro/download/). Note that Interpro Scan
**only runs on Linux systems**. Just as with ORFfinder, make sure to use the appropiate paths to the InterPro Scan directory when running this. Before running InterPro Scan, it is necessary to
translate the nucleotide sequences to aminoacid. This is done using the ape R package
v. 5.3 (Paradis & Schliep, 2018). All sequences that do not have any of the two
domains, were eliminated for further analyses.

Finally, all the retained sequences have to be blasted against two protein databases for
identification, this case: Araport11, downloaded from [The Arabidopsis Information Resource](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FProteins%2FAraport11_protein_lists)
and Swiss-Prot, downloaded from [UniProt](https://www.uniprot.org/downloads).
All the sequences have to be annotated according to their query sequence and their best hit.
of every identification blast. In this case, the `blastx` comand has to be used intead of `blastn`.

### Alingments and phylogenies

After the final sequences for each species were selected, a custom bash script can be written that first, asssmebles one-species matrixes with all sequences of every species, and then, assembles all individual matrixes into one big multispecies matrix. This multispecies matrix is then alingned by codons, using the DECIPHER package from Bioconductor (Wright, 2016).The resulting alingment in fasta format needs to be transformed to nexus format, which in this case was done using Mesquite 3.6 (Maddison & Maddison, 2018).

For the construction of a bayesian phylogeny,  MrBayes 3.2.7 was used and run on [The CIPRES Science Gateway V. 3.3](http://www.phylo.org/). The command line for MrBayes to run the selected model and parameters needs to be included at the bottom of the Nexus file, so it can be directly run in CIPRES without manually adjusting the parameters on their website. This same Nexus file can be run from the terminal in any cluster or personal computer with a pre-installed version of MrBayes by running the following:

```
### Assuming this is run from the bin directory:
mb
execute ../data/path_to_the/nexus_alignment.nex
```
For the plotting and visualization the MrBayes output tree, the `ggtree` package from Bioconductor (Yu et al., 2017), whichs has the advantage that it can read the outputs of several software.
