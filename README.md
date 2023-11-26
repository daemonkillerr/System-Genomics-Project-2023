# System-Genomics-Project-2023

## Preproccessing
For preprocessing we ran bash scripts in the following order:
- remove_adapters.sh, to remove adaptor sequences in the reads ->
- build_genome_index.sh, to download the reference genome as a fasta format and index it ->
- get_genome.sh, to download the genome annotations ->
- rsem_prepare_reference.sh, to create an RSEM index for the genome+annotation ->
- counting.sh, to compute transcription counts for differential expression analysis 
