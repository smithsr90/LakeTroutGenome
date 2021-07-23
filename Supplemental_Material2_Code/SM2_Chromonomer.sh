### - Pipeline for integrating scaffolding information from the Lake Trout linkage map.

## First round of chromonomer - using existing gaps
# Map the linkage map to HiC scaffolds
minimap2 -t 8 -ax sr ./scaffolds_HiC.fasta ../map.fasta | samtools view -h -q60 > map.best.sam

# Run the first round of chromonomer
chromonomer --map map_info.txt --alns map.best.sam --fasta scaffolds_HiC.fasta --agp alt.agp --out_path ./chromonomer_std &

## Second round of chromonomer - gaps open at regions with depth that is greater than 2 SD from mean
#Map the complete PacBio dataset to initial chromonomer assembly
minimap2 -ax map-pb -t 20 chromonomer_std.fasta allreads.fa.gz > pacbio_mapped.sam &
samtools view -b -S -F2308 pacbio_mapped.sam > pacbio_mapped_filtered.bam &
samtools sort pacbio_mapped_filtered.bam -o pacbio_mapped_filtered_sorted.bam &

#Calculate per-base depth
samtools depth -aa pacbio_mapped_filtered_sorted.bam >  chromonomer_std.depth &

# Second round of chromonomer with --rescaffold option and gap opening at putative errors
minimap2 -t 8 -ax sr ./chromonomer_std.fa ../map.fasta | samtools view -h -q60 > chromonomer_std_lm.sam
chromonomer --map ../map_info.txt --alns ./chromonomer_std_lm.sam --agp ./chromonomer_std.agp --depth chromonomer_std.depth --depth_stdevs 2 --depth_win_size 1000 --rescaffold --fasta ./chromonomer_std.fa --out_path ./chromonomer_fin &


