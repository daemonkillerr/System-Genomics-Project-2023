cd $SCRATCH/rawdata23
if [ ! -e genome ]; then
    mkdir genome
fi
wget -O genome/gencode.vM33.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz
gzip -d genome/gencode.vM33.annotation.gtf.gz