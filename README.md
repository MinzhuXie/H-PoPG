# H-PoPG
H-PoP and H-PoPG:  Heuristic Partitioning Algorithms for Single Individual Haplotyping of Polyploids

```
Usage: 
       java -jar H-PoPG.jar -p k-ploidy w weight (default 0.9) -f input_SNP_matrix_file_name -vcf vcf_file_name  -o output_phased_haplotypes_file_name
```

If -vcf option and the vcf file name are provided, H-PoPG is called, otherwise H-PoP is called.

For example: 

```
(1)  
java -jar H-PoPG.jar -p 3 -f frags.txt -o phased_haplotypes.txt

   The command will call the algorithm H-PoP to reconstruct the haplotypes for a triploid (k=3). 
   The input SNP matrix is provided by the file frags.txt.  
   The reconstructed three haplotypes are stored in the output file phased_haplotypes.txt.
```

```
(2) 
java -jar H-PoPG.jar -p 3 -f frags.txt -vcf example.vcf -o phased_haplotypes.txt

   The command will call the algorithm H-PoPG to reconstruct the haplotypes for a triploid (k=3)
   The genotype informaion is provided by the input file example.vcf
```

Input SNP matrix files uses a space(or tab) delimited file format with aligned SNP reads separated by new lines. 
The columns are specified as follows:

```
column 1: the number of read blocks (consecutive 0 1 sequence separated by a gap (i.e. a sequence of -))
column 2: read name
column 3: the start position of the first block
column 4: the 0 1 sequence of the first block
```

If there are more than 1 blocks, additional column pairs like columns 3 and 4 are needed for additional blocks. 

the last column: 
the sequence encoding the quality values for the sequence of all blocks, and must contain the same number of symbols as the total number of 0 and 1 of the read.

Example of input_SNP_matrix_file: 

```
1 read_1 23 0101 G,&G
1 read_2 24 0011 HGGG
4 read_3 25 0 26 10 29 0 31 1 GGGGG
2 read_4 26 00011 40 11100 GGGGGGGGGG
```

If a VCF file is provided, the start position of a read block _MUST_ match the index of the corrsponding SNP in the VCF file.

Please see the example files frags.txt, example.vcf and phased_haplotypes.txt
for the formats of the input SNP matrix file , genotype file and the output file . 

There is a discussion about the problems of using H-PoPG in stackoverflow (Thank Mel Carbajo for starting the discussion), and if you have any problems about using H-PoPG, please join the discussion. 
http://stackoverflow.com/questions/39859409/h-popg-haplotyper-nullpointerexception-error-at-algorithms-hbop2builder


# Version 0.2.0 updated on Oct. 18, 2019

H-PoPG version: 0.2.0. Author: Xie Minzhu. Email: xieminzhu@hotmail.com 

(1) Without -a option, the polyplotyping algorithm called is H-PoP or H-PoPG according to whether a vcf is provided.

  Usage: 
  ```
    java -jar H-PoPGv0.2.0.jar [-p ploidy (default is 3) -m maxsize (default is ploidy*ploidy*10) -w weight (default is 0.9)] -v[cf] vcffilename -b[am] bamfilenme [-f[rag] alignedFragsFileName].. -o outputPhasedFilename
      
    java -jar H-PoPGv0.2.0.jar [-p ploidy (default is 3) -m maxsize (default is ploidy*ploidy*10) -w weight (default is 0.9)] -f[rag] alignedFragsFileName [-f[rag] alignedFragsFileName].. [-v[cf] vcffilename] -o outputPhasedFilename
```

(2) With the option "-a HapTree", HapTree is called. (HapTree_v0 should be placed at the same directory as this java package.)

  Usage: 
```  
    java -jar H-PoPGv0.2.0.jar -a HapTree [-p ploidy (default is 3) -v[cf] vcffilename -b[am] bamfilenme [-f[rag] alignedFragsFileName].. -o outputPhasedFilename
      
    java -jar H-PoPGv0.2.0.jar -a HapTree [-p ploidy (default is 3)] -v[cf] vcffilename -f[rag] alignedFragsFileName [-f[rag] alignedFragsFileName].. -o outputPhasedFilename
```

(3) Only transfer a BAM File to a Frag file (the read covering less than 2 SNPs are discarded): 

```
    java -jar H-PoPGv0.2.0.jar -v[cf] vcffilename -b[am] bamfilenme -o outputFragFilename -b2f
```

(4) Output a File with detailed information about position in genome, SNP id, SNP alleles of the phased haplotypes, use the option -detailedoutput or -d as follows:

```
    java -jar H-PoPGv0.2.0.jar [-p ploidy (default is 3) -m maxsize (default is ploidy*ploidy*10) -w weight (default is 0.9)] -v[cf] vcffilename -b[am] bamfilenme [-f[rag] alignedFragsFileName].. -d[etailedoutput] outputdetailedPhasedFilename
```

(5) Print this help information: 
```
    java -jar H-PoPGv0.2.0.jar -h
```

 Notes: 
 
    (1) If a BAM File has been provide, a VCF file is needed. The input files can include multiple Frag Files, one Bam file and one VCF file.
    
    (2) For case (4), the vcf file must be provided.
