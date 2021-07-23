#!/bin/bash

time java -Xms100g -Xmx210g -jar 01_scripts/pilon-1.23.jar --genome 03_genome/genome.fasta --frags 02_raw/GLW_31_F_USPD16101428_H32YYCCX2_L4_trimmed_1.sorted.bam --frags 02_raw/SIW_25_M_USPD16101429_H32YYCCX2_L3_trimmed_1.sorted.bam --frags 02_raw/SLW_11_M_USPD16101426_H32YYCCX2_L3_trimmed_1.sorted.bam --frags 02_raw/SLW_52_F_USPD16101427_H32YYCCX2_L4_trimmed_1.sorted.bam --outdir 04_after_pilon/ --changes --tracks --diploid --fix indels,gaps,local --mindepth 5
