cd $SCRATCH/rawdata23

if [ ! -e mapping_transcriptome ]; then
    mkdir mapping_transcriptome
fi

if [ ! -e rsem ]; then
    mkdir rsem
fi

for id in `cat SRR_Acc_List.txt`; do
    echo "start to process sample $id"

    if [ ! -e mapping_transcriptome/$id ]; then
        echo " mapping started"
        mkdir mapping_transcriptome/$id
        STAR --genomeDir genome/star-index \
             --runThreadN 10 \
             --readFilesIn trimmed_dedup/${id}*.fastq.gz \
             --readFilesCommand zcat \
             --sjdbGTFfile genome/gencode.vM33.annotation.gtf \
             --quantMode TranscriptomeSAM \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix mapping_transcriptome/$id/
        echo " mapping done"
    fi

    if [ ! -e rsem/$id ]; then
        echo "rsem started"
        mkdir rsem/$id
        rsem-calculate-expression --alignments \
                                    -p 10 \
                                    mapping_transcriptome/$id/Aligned.toTranscriptome.out.bam \
                                    genome/rsem_vM33_gencode41/rsem_vM33_gencode41 \
                                    rsem/$id/$id
        echo " rsem done"
    fi
done