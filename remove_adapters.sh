cd {$SCRATCH}/rawdata23
mkdir trimmed_dedup
for file in *.fastq.gz; do
    fastp --dedup --adapter_sequence AAGCAGTGGTATCAACGCAGAGTAC -l 25 -i $file -o trimmed_dedup/$file
done