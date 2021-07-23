# Polishing with polca
polca_mod.sh -a genome.fa -r 'SLW_52_F_paired.1.fq SLW_52_F_paired.2.fq' -t 36 &
mv genome.fa.PolcaCorrected.fa genome_p1.fa 

polca_mod.sh -a genome_p1.fa -r 'SLW_52_F_paired.1.fq SLW_52_F_paired.2.fq' -t 36 &
mv genome_p1.fa.PolcaCorrected.fa genome_p2.fa 

polca_mod.sh -a genome_p2.fa -r 'SLW_52_F_paired.1.fq SLW_52_F_paired.2.fq' -t 36 &
mv genome_p2.fa.PolcaCorrected.fa genome_p3.fa 