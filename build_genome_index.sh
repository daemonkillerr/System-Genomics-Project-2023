cd $SCRATCH/rawdata23/genome

cd genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gzip -d mm39.fa.gz

if [ ! -d "star-index" ]; then
    mkdir star-index
fi

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir star-index --genomeFastaFiles mm39.fa

