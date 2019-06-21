# R course 11

Het standaard protocol voor C. elegans is een voedingsbodem van E.coli. Tijdens dit project is onderzocht welke gen expressie veranderingen er te weeg worden gebracht als de voedingsbodem wordt veranderd naar B. subtillis.

### Data
De fastq bestanden zijn verkregen uit de NCBI database. De bestanden zijn aan de hand van een FastQC analyse getrimd en gefilterd op Q score.

### Mappen
Met als input de getrimde fastq bestanden is er gemapt met behulp van Bowtie2. Als eerst is er een reference gebouwd. 
Het referentie genoom is van de directory NAS gehaald.

```
bowtie2-build c_elegans.PRJEB28388.WS271.genomic.fa reference_genome
```

Als output zijn er een aantal BT2-bestanden uit voort gekomen.
Reference_genome is het referentie gen dat is gebouwd in het voorgaande commando.

```
bowtie2 -q -x reference_genome -U trimmed_SRR5832185.fastq -S SRR5832185_result.sam
```

Als output wordt er een sam bestand verkregen. Dit stam bestand wordt gesorteerd.

```
samtools sort -n SRR5832185_result.sam -o SRR5832185_result-sorted.sam
``` 

Deze gesorteerde bestanden worden gebruikt voor het maken van de countbestanden door middel van R samen een countmatrix vormen.

```
htseq-count -r name ./SRR5832185_result-sorted.sam ./c_elegans.PRJEB28388.WS271.annotations.gtf  > counts_SRR5832185.txt
```

### RNA-seq analyse
Met behulp van R zijn de significante genen bepaald. 

Het script vormt de losse count matrices om tot een volledige countmatrix met alle 6 samples en hun genen. Hierna is een low count filtering uitgevoerd waarbij alle genen met een count totaal lager dan 50 worden verwijderd uit de dataset. Normalisatie is uitgevoerd met behulp van de TMM methode, welke speciaal is ontwikkeld voor RNA-seq data. Dispersie is berekend en geplot in een MDS en BCV plot. Tevens is er een hierarchische clustering uitgevoerd om de afstanden tussen samples weer te geven.

Genen met een p-value lager dan 0.05 zijn tijdens dit project beschouwt als significant. Dit houd in dat het expressie niveau significant is veranderd tussen de condities. 
 
Genen samen met hun p-value en false discovery rate zijn weggeschreven naar een .txt bestand voor een latere GSEA analyse. 

