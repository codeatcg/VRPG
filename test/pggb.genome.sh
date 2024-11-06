mkdir pggb_genome
zcat  GCF_000146045.2_R64_genomic.fna.gz > pggb_genome/SGDref.fa
fastix -p "SGDref#1#" pggb_genome/SGDref.fa >> pggb_genome/all.fastix.fa
for i in DBVPG6765 Y12 SK1 DBVPG6044
do
    zcat  $i.genome.fa.gz > pggb_genome/$i.genome.fa
    fastix -p "${i}#1#" pggb_genome/$i.genome.fa >> pggb_genome/all.fastix.fa
done

samtools faidx pggb_genome/all.fastix.fa

