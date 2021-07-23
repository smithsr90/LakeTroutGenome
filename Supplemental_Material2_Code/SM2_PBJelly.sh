################################################################################
# - Commands used for genome polishing with PBJelly
################################################################################
# Map reads to gaps, get flanking reads, extract flanking reads, split into evenly sized fasta files for more efficient parallelization
minimap2 -ax map-pb -t 32 ./falcon_pilon_salsa2i3_chromostd_chromo2sd1000.fa allreads.fa.gz > scaf.sam
getFlanks.sh ./falcon_pilon_salsa2i3_chromostd_chromo2sd1000.fa &
samtools view -h scaf.sam -L flanks.bed -F2304 | grep -v "@" | awk '{print $1}' > informative_primary_hits.txt &
cat informative_primary_hits.txt | uniq | sort | uniq > informative_primary_reads.txt
seqkit grep -n -f informative_primary_reads.txt allreads.fa.gz  > informative_primary.fa &
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%42000==0){file=sprintf("split%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < informative_primary.fa &

# Run PBJelly
conda create -n pbsuite -c bioconda python=2
git clone https://github.com/esrice/PBJelly.git
# Modify path and conda environment in setup.py file
byobu
source activate pbsuite
source /home/seth_smith/software/PBJelly/setup.sh

# Make file called Protocol2.xml - File looks like this:
<jellyProtocol>                                    
    <reference>/nfs/westslope/1/storage/seth_smith/PacBio_LKTDoubledHaploid_GenomeAssembly/gapfilling/before/genome.fasta</reference>                                                                         
    <outputDir>/nfs/westslope/1/storage/seth_smith/PacBio_LKTDoubledHaploid_GenomeAssembly/gapfilling/jelly/</outputDir>                                                                                      
    <cluster>                                      
        <command notes="For single node, multi-core machines" >${CMD} ${JOBNAME} 2> ${STDERR} 1> ${STDOUT} &amp;</command>                                                                                    
        <nJobs>37</nJobs>                          
    </cluster>                                     
    <blasr>--minMatch 11 --minPctIdentity 75 --bestn 1 --nCandidates 10 --maxScore -500 --nproc 1 --fastSDP --noSplitSubreads</blasr>                                                                         
    <input baseDir="/nfs/westslope/1/storage/seth_smith/PacBio_LKTDoubledHaploid_GenomeAssembly/gapfilling/input/fastq/">                                                                                     
        <job>split0.fastq</job>                    
        <job>split1008000.fastq</job>              
        <job>split1050000.fastq</job>              
        <job>split1092000.fastq</job>              
        <job>split1134000.fastq</job>              
        <job>split1176000.fastq</job>              
        <job>split1218000.fastq</job>              
        <job>split1260000.fastq</job>              
        <job>split126000.fastq</job>               
        <job>split1302000.fastq</job>              
        <job>split1344000.fastq</job>              
        <job>split1386000.fastq</job>              
        <job>split1428000.fastq</job>              
        <job>split1470000.fastq</job>              
        <job>split168000.fastq</job>               
        <job>split210000.fastq</job>               
        <job>split252000.fastq</job>               
        <job>split294000.fastq</job>               
        <job>split336000.fastq</job>               
        <job>split378000.fastq</job>               
        <job>split420000.fastq</job>               
        <job>split42000.fastq</job>                
        <job>split462000.fastq</job>               
        <job>split504000.fastq</job>               
        <job>split546000.fastq</job>               
        <job>split588000.fastq</job>               
        <job>split630000.fastq</job>               
        <job>split672000.fastq</job>               
        <job>split714000.fastq</job>               
        <job>split756000.fastq</job>               
        <job>split798000.fastq</job>               
        <job>split840000.fastq</job>               
        <job>split84000.fastq</job>                
        <job>split882000.fastq</job>               
        <job>split924000.fastq</job>               
        <job>split966000.fastq</job>               
    </input>                                       
</jellyProtocol>              


# The main pipeline
Jelly.py setup Protocol1.xml 
Jelly.py mapping Protocol1.xml
Jelly.py support Protocol1.xml
Jelly.py extraction Protocol1.xml
Jelly.py assembly Protocol1.xml -x "--maxWiggle=100000"
Jelly.py output Protocol1.xml






