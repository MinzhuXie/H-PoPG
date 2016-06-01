# H-PoPG
H-PoP and H-PoPG:  Heuristic Partitioning Algorithms for Single Individual Haplotyping of Polyploids

Usage: 
       java -jar H-PoPG.jar -p k-ploidy w weight (default 0.9) -f input_read_file_name -vcf vcf_file_name  -o output_phased_haplotypes_file_name

If -vcf option and the vcf file name are provided, H-PoPG is called, otherwise H-PoP is called.

For example: 
(1)  java -jar H-PoPG.jar -p 3 -f frags.txt -o phased_haplotypes.txt
   The command will call the algorithm H-PoP to reconstruct the haplotypes for a triploid (k=3). 
   The input SNP matrix is provided by the file frags.txt.  
   The reconstructed three haplotypes are stored in the output file phased_haplotypes.txt.

(2) java -jar H-PoPG.jar -p 3 -f frags.txt -vcf example.vcf -o phased_haplotypes.txt
   The command will call the algorithm H-PoPG to reconstruct the haplotypes for a triploid (k=3)
   The genotype informaion is provided by the input file example.vcf
   
Please see the files frags.txt, example.vcf and phased_haplotypes.txt
for the formats of the input read file, genotype file and the output file. 
