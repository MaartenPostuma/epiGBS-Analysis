#Bash code for repeatmasker epiGBS denovo mode (Do not know how this works for ref branch)

conda create -f epiGBSanalysis.yaml
conda activate epiGBSanalysis.yaml

wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz
gunzip Dfam.h5.gz

whereis RepeatMasker #Check where Repeatmasker is

#Add the first part of the path here (Hopefully you do not need sudo rights to adjust things in your conda evv)
mv Dfam.h5 /conda/env/epiGBSanalysis/share/RepeatMasker/Libraries


RepeatMasker -pa 8 -species Embryophyta -xsmall output/output_denovo/consensus_cluster.renamed.fa




#Get the blast database! Unless you already have it then do not do this.  
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -O /scratch/tmp//nr.gz 
diamond makedb -p 8 --in /scratch/tmp//nr.gz -d /scratch/tmp//nrdb.dmnd
rm /scratch/tmp//nr.gz

# search for matches to the proteins
diamond blastx -d /scratch/tmp//nrdb.dmnd -q output/output_denovo/consensus_cluster.renamed.fa -p 8 -f 5 -o /scratch/tmp//diamond.xml.gz --tmpdir /scratch/tmp/ --compress 1 -e 0.00001 -k 10

# replace some weird characters
zcat /scratch/tmp//diamond.xml.gz | sed "s/ &gt;/|/g" | sed "s/ & / and /g" | sed "s/&//g" | pigz -p 4 > GitIgnore_annotation/diamond_reformed.xml.gz

# parse blast and get the best hits
python scripts/parseBlastXMLbestHit.py GitIgnore_annotation/diamond_reformed.xml.gz
