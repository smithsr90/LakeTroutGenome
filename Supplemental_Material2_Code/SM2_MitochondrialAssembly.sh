# Assemble the lake trout mitome using Unicycler
conda create -n unicycler python=3
source activate unicycler
conda install -c bioconda unicycler

# Grab the arctic char mitome sequence
cd /nfs/westslope/1/storage/seth_smith/PacBio_LKTDoubledHaploid_GenomeAssembly/LT/mitome
samtools faidx /scratch/seth/Genomes/Arctic_Char/salpinus_v2_ref_norm.fa NC_000861.1 > achar_mitome.fa

# Map the arctic char mitome sequence to lake trout scaffolds
minimap2 -a -x map-pb ../reads/m54083_181211_211206.fakequal.fastq achar_mitome.fa > m54083_181211_211206_mitome_alignment.sam &

# Grab a subset of PacBio reads that map the the arctic char genome
minimap2 -a -x map-pb achar_mitome.fa ../reads/m54083_181211_211206.fakequal.fastq | grep -v "@" | grep "NC_000861.1" | awk '{print $1}' > pacbio_reads.list &
minimap2 -a -x sr achar_mitome.fa ../../../GLFC_LKT_WGSSeq_071119/SLW_52_F_1.fq.gz ../../../GLFC_LKT_WGSSeq_071119/SLW_52_F_2.fq.gz | grep -v "@" | grep "NC_000861.1" | awk '{print $1}' > illumina_reads.list &

# Grab mitome associated reads using sequence IDs from mapping
seqtk subseq ../reads/m54083_181211_211206.fakequal.fastq <(uniq pacbio_reads.list) > long_reads.fq &
seqtk subseq ../../../GLFC_LKT_WGSSeq_071119/SLW_52_F_1.fq.gz <(uniq illumina_reads.list) > illumina_reads.1.fq &
seqtk subseq ../../../GLFC_LKT_WGSSeq_071119/SLW_52_F_2.fq.gz <(uniq illumina_reads.list) > illumina_reads.2.fq &

# Check for proper pairing - this shouldn't return anything
diff <(cat illumina_reads.1.fq | grep "@" | awk '{print $1}') <(cat illumina_reads.2.fq | grep "@" | awk '{print $1}')

# Assemble using unicycler
unicycler -1 illumina_reads.1.fq -2 illumina_reads.2.fq -l long_reads.fq -o ./ --min_fasta_length 15000 --keep 0 -t 1 &