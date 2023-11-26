cd $SCRATCH/rawdata23
mkdir genome/rsem_vM33_gencode41
rsem-prepare-reference --gtf genome/gencode.vM33.annotation.gtf \
    genome/mm39.fa \
    genome/rsem_vM33_gencode41/rsem_vM33_gencode41