package io;

import java.io.*;
import java.util.*;

import main.PolyPlotyping;


import dataStructure.ComplexVariant;
import dataStructure.InDel;
import dataStructure.SNP;
import dataStructure.Variant;

//obtain the genotypes of the individual
public class VCFProcessor
{
	 PolyPlotyping pBop;
	 private int ploidy_number;
	 private int number_of_people;
	 public Map<String, List<List<byte[]>>> ref_to_population_genotypes;
	 public List<String> vcf_reference_names;
	 public Map <String, List<Variant>> ref_to_variant_list;
	 public Map <String, Map<String, Variant>> ref_to_id_to_SNP;	
	 public VCFProcessor(int ploidy_n, int n_of_people, PolyPlotyping alg)
	 {    
	     ref_to_population_genotypes = new HashMap<String, List<List<byte[]>>>(); 
	     vcf_reference_names = new ArrayList<String>();        
	     ref_to_variant_list = new HashMap<String, List<Variant>>();     
	     ref_to_id_to_SNP = new HashMap<String, Map<String, Variant>>();	        
	     ploidy_number = ploidy_n;
	     number_of_people = n_of_people;	
	     pBop = alg;
	 }	

	 public void loadVCF(String filename)
	     throws IOException
	 {
		 
	     BufferedReader vreader = new BufferedReader(new FileReader(new File(filename)));
	     String line = null;
	     int iter = 0;
	     while((line = vreader.readLine()) != null)
	     {
	    	 if(line.startsWith("#"))
	    	 {	    		
	    		 continue;
	    	 }	
	         if(!line.trim().equals(""))
	         {
	             int index_of_GT = -1;
	             int index_of_FT = -1;
	             String parts[] = line.split("\t");            
	             if(parts.length != 9 + number_of_people)
	             {                 
	                 System.out.println("Non-header line of the vcf file does not have the correct number of columns:\n");
	                 System.out.println("Expected 9 + number_of_people columns."); 
	                 System.exit(-1);
	             }
	             String person_parts[] = parts[8].split(":");
	             for(int i = 0; i < person_parts.length; i++)
	             {
	                 if(person_parts[i].equals("GT"))
	                     index_of_GT = i;
	                 if(person_parts[i].equals("FT"))
	                     index_of_FT = i;
	             }
	
	             if(index_of_GT == -1)
	             {
	                 System.out.println("Line in VCF file did not have a GT field. Genotype field required for haplotype assembly.\n");
	                 continue;
	             }          
	             
	             
	             boolean unknown = false;   
	             boolean nopass = true;
	             
	            
	             for(int i = 0; i < number_of_people; i++)
	             {            	 
	                 String[] tmpts = parts[9+i].split(":");                 
	                 if(index_of_FT!=-1)
	     //Xie added ||tmpts[index_of_FT].equals(".") to deal with the missing Filter value
	     //2017.08.12
	                	 if(tmpts[index_of_FT].equals("PASS")||tmpts[index_of_FT].equals("."))
	                		 nopass = false;
	                 
	                 String[] gt = tmpts[index_of_GT].split("/|\\|");
	                 for(int k= 0; k<gt.length; k++)
	                     if(gt[k].charAt(0) == '.')
	                         unknown = true;
	             }
	             
	 //Xie added ||parts[6].equals(".") to deal with the missing Filter value
	             //2017.08.12
	             if(parts[6].equals("PASS")||parts[6].equals("."))
	            	 nopass = false;
	            
	 //            if(unknown || nopass)
	 //            {
	 //           	 System.out.println("Skipping variant because an allele was missing. We don't assume we can call genotypes better than the genotype caller.\n");
	 //            } 
	 //            else
	             {
	                 String vcf_snp_chr = new String(parts[0]);
	                 if(!vcf_reference_names.contains(vcf_snp_chr))
	                 {
	                     ref_to_variant_list.put(vcf_snp_chr, new ArrayList<Variant>());
	                     ref_to_id_to_SNP.put(vcf_snp_chr, new HashMap<String, Variant>());	                   
	                     vcf_reference_names.add(vcf_snp_chr);
	                     ref_to_population_genotypes.put(vcf_snp_chr, new ArrayList<List<byte[]>>());
	                     for(int i = 0; i < number_of_people; i++)
	                         ref_to_population_genotypes.get(vcf_snp_chr).add(new ArrayList<byte[]>());
	
	                 }
	           
	                  for(int i = 0; i < number_of_people; i++)
	                  {
	                        byte genotype[] = getGenotype(parts[9 + i].split(":")[index_of_GT]);
	                        if(unknown || nopass) genotype = null;
	                        ref_to_population_genotypes.get(vcf_snp_chr).get(i).add(genotype);
	                  }                        
	                                           
	                  addVariant(parts, iter, vcf_snp_chr);
	                  iter++; 
	             }
	         }
	     }
	     vreader.close();
	     
	     Set<String> toRemove = new HashSet<String>();
	     for(String ref : ref_to_variant_list.keySet())
	     {
	         int position = 1;
	         for(Variant var : ref_to_variant_list.get(ref) )
	         {                        
	             var.reduced_genome_start_position = position;
	             position++;
	         }
	         if(ref_to_variant_list.get(ref).size() == 0)
	             toRemove.add(ref);
	     }
	
	     String ref;
	     for(Iterator<String> iterator1 = toRemove.iterator(); iterator1.hasNext(); ref_to_id_to_SNP.remove(ref))
	     {
	         ref = iterator1.next();
	         ref_to_variant_list.remove(ref);
	         vcf_reference_names.remove(ref);
	     }
	     
	 }

	 private void addVariant(String parts[], int iter, String ref)
	 {
	     String id = (new StringBuilder("VAR_POS_")).append(parts[1]).toString();
	     if(!parts[2].equals("."))
	         id = new String(parts[2]);
	     if(ref_to_variant_list.get(ref).size() > 0 
	    		 && 
	    		 ref.equals((ref_to_variant_list.get(ref).get(ref_to_variant_list.get(ref).size() - 1)).reference) 
	    		 && 
	    		 Integer.parseInt(parts[1]) < (ref_to_variant_list.get(ref).get(ref_to_variant_list.get(ref).size() - 1)).genome_start_position + (ref_to_variant_list.get(ref).get(ref_to_variant_list.get(ref).size() - 1)).getEffectiveLength())
	         System.out.println("Variant "+ (ref_to_variant_list.get(ref).get(ref_to_variant_list.get(ref).size() - 1)).toString()+" overlaps with variant: "+ parts[0]+"\t"+ parts[1]+"\t"+parts[2]+"\t"+parts[3]+"\t"+parts[4]+ "\t" + parts[5]+"\t"+parts[6]+"\n");
	     String alt_alleles[] = parts[4].split(",");
	     int alt_length = 0;
	    // String as[];
	     int j = alt_alleles.length;
	     for(int i = 0; i < j; i++)
	     {
	         String allele_string = alt_alleles[i];
	         alt_length += allele_string.length();
	     }
	
	     if(parts[3].length() == 1 && alt_length == alt_alleles.length)
	     {
	         String alleles[] = new String[1 + alt_alleles.length];
	         alleles[0] = new String(parts[3].substring(0, 1));
	         int index = 1;
	         String as1[];
	         int i1 = (as1 = alt_alleles).length;
	         for(int k = 0; k < i1; k++)
	         {
	             String str = as1[k];
	             alleles[index++] = new String(str.substring(0, 1));
	         }
	
	         SNP snp_var = new SNP(iter, Integer.parseInt(parts[1]), id, new String(parts[0]), alleles);
	         ref_to_variant_list.get(ref).add(snp_var);
	         ref_to_id_to_SNP.get(ref).put(snp_var.identifier, snp_var);
	     } 
	     else if(1 + alt_alleles.length == 2 && parts[3].length() != alt_alleles[0].length())
	     {
	         String alleles[] = new String[1 + alt_alleles.length];
	         alleles[0] = new String(parts[3]);
	         int index = 1;
	         String as2[];
	         int j1 = (as2 = alt_alleles).length;
	         for(int l = 0; l < j1; l++)
	         {
	             String str = as2[l];
	             alleles[index++] = new String(str);
	         }
	
	         InDel indel = new InDel(iter, Integer.parseInt(parts[1]), id, new String(parts[0]), alleles);
	         ref_to_variant_list.get(ref).add(indel);
	         ref_to_id_to_SNP.get(ref).put(indel.identifier, indel);
	     } else
	     {
	         String all_alleles[] = new String[alt_alleles.length + 1];
	         all_alleles[0] = new String(parts[3]);
	         System.arraycopy(alt_alleles, 0, all_alleles, 1, alt_alleles.length);         
	         ComplexVariant cv = new ComplexVariant(iter, Integer.parseInt(parts[1]), id, new String(parts[0]), all_alleles);
	         ref_to_variant_list.get(ref).add(cv);
	         ref_to_id_to_SNP.get(ref).put(cv.identifier, cv);
	     }
	 }
	 
	 private byte[] getGenotype(String combo)
	 {
	     byte alleles_byte[] = new byte[ploidy_number];
	     Set<Byte> allele_set = new HashSet<Byte>();
	     String alleles[] = combo.split("/|\\|");     
	     for(int i = 0; i < alleles.length; i++)
	     {
	         Byte allele_byte = Byte.valueOf(Byte.parseByte(alleles[i]));
	         allele_set.add(allele_byte);
	         alleles_byte[i] = allele_byte.byteValue();
	     }
	
	     if(allele_set.size() >= 2)
	         return alleles_byte;
	     else
	         return null;
	 }

}
