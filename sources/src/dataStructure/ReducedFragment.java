package dataStructure;

import java.util.ArrayList;
import java.util.List;

public class ReducedFragment
{
	 List<List<Byte>> alleles;
	 List<List<Integer>> positions;
	 List<String> references;
	 List<List<String>> qualities;
	 private int num_alleles;	 

	 public ReducedFragment()
	 {
	     alleles = null;
	     positions = null;
	     references = null;
	     qualities = null;
	     num_alleles = 0;
	     alleles = new ArrayList<List<Byte>> (2);
	     positions = new ArrayList<List<Integer>>(2);
	     references = new ArrayList<String>(2);
	     qualities = new ArrayList<List<String>>(2);
	 }
	
	 public ReducedFragment(int size_of_collections)
	 {
	     alleles = null;
	     positions = null;
	     references = null;
	     qualities = null;
	     num_alleles = 0;
	     alleles = new ArrayList<List<Byte>>(size_of_collections);
	     positions = new ArrayList<List<Integer>>(size_of_collections);
	     references = new ArrayList<String>(size_of_collections);
	     qualities = new ArrayList<List<String>>(size_of_collections);
	 }
	
	 public void addFragment(List<Byte> alleles_in, List<Integer> pos_in, List<String> qualities_in, String reference)
	 {
	     int added = 0;
	     alleles.add(new ArrayList<Byte>(2));	     
	     qualities.add(new ArrayList<String>(2));
	     positions.add(new ArrayList<Integer>(2));
	     for(int i = 0; i < pos_in.size(); i++)
	     {
	         boolean found = false;
	         for(int j = 0; j < positions.size(); j++)
	             if((positions.get(j)).contains(pos_in.get(i)))
	                 found = true;
	
	         if(!found)
	         {
	             added++;
	             (alleles.get(references.size())).add(alleles_in.get(i));
	             (qualities.get(references.size())).add(qualities_in.get(i));
	             (positions.get(references.size())).add(pos_in.get(i));
	         }
	     }
	
	     if(added == 0)
	     {
	         alleles.remove(references.size());
	         qualities.remove(references.size());
	         positions.remove(references.size());
	     } else
	     {
	         references.add(reference);
	         num_alleles += added;
	     }
	 }
	
	 public void addFragment(ReducedFragment rf)
	 {
	     if(!rf.references.isEmpty())
	         addFragment(rf.alleles.get(0), rf.positions.get(0), rf.qualities.get(0), rf.references.get(0));
	 }
	
	 public String toString()
	 {
	     if(alleles.size() == 0)
	         return new String();
	     StringBuilder toString = new StringBuilder();
	     for(int i = 0; i < references.size(); i++)
	     {
	         toString.append((String)references.get(i));
	         for(int j = 0; j < (alleles.get(i)).size(); j++)
	         {
	             toString.append(" ");
	             toString.append((positions.get(i)).get(j));
	             toString.append(" ");
	             toString.append((alleles.get(i)).get(j));
	             toString.append(" ");
	             toString.append(qualities.get(i).get(j));
	         }
	
	         toString.append("\t");
	     }
	
	     toString.replace(toString.length() - 1, toString.length(), System.getProperty("line.separator"));
	     return toString.toString();
	 }
	
	 public int getNumberOfReads()
	 {
	     return references.size();
	 }
	
	 public static ReducedFragment parseToString(String line)
	 {
	     String parts[] = line.split("\t");
	     int cnt = Integer.parseInt(parts[0]);
	     ReducedFragment rf = new ReducedFragment(cnt);
	     for(int i = 2; i < parts.length; i++)
	     {
	         String frag_parts[] = parts[i].split(" ");
	         rf.references.add(frag_parts[0]);
	         rf.alleles.add(new ArrayList<Byte>());
	         rf.qualities.add(new ArrayList<String>());
	         rf.positions.add(new ArrayList<Integer>());
	         for(int j = 1; j < frag_parts.length; j += 3)
	         {
	             (rf.positions.get(i - 2)).add(Integer.valueOf(Integer.parseInt(frag_parts[j])));
	             (rf.alleles.get(i - 2)).add(Byte.valueOf(Byte.parseByte(frag_parts[j + 1])));
	             (rf.qualities.get(i - 2)).add(frag_parts[j + 2]);
	             rf.num_alleles++;
	         }
	
	     }
	
	     return rf;
	 }
	
	 public List<List<Byte>> getAlleles()
	 {
	     return alleles;
	 }
	
	 public List<List<Integer>> getPositions()
	 {
	     return positions;
	 }
	
	 public List<String> getReferences()
	 {
	     return references;
	 }
	
	 public List<List<String>> getQualities()
	 {
	     return qualities;
	 }
	
	 public int getNum_alleles()
	 {
	     return num_alleles;
	 }
	
	 public void addFragment(List<List<Byte>> alleles2, List<List<Integer>> variant_pointers, List<List<String>> qualities2, List<String> refs)
	 {
	     for(int i = 0; i < refs.size(); i++)
	         addFragment(alleles2.get(i), variant_pointers.get(i), qualities2.get(i), refs.get(i));
	
	 }	
	 
	 public String HapCutFormatFragment(){
		 
		 
		 return null;
	 }
	 

	 
	 public String merge(String rid)
	 {
		 List<Integer> snp_pointers_1, snp_pointers_2;
		 
		 List<Byte> alleles_1, alleles_2 = null;	 
		
		 
		 if(alleles.size() == 0)
	         return new String();
		 snp_pointers_1 = positions.get(0);	
		 alleles_1 = alleles.get(0);
		 snp_pointers_2 = null;
		 if(alleles.size()==2)
		 {
			 snp_pointers_2 = positions.get(1);	
			 alleles_2 = alleles.get(1);
			 if(!references.get(0).equals(references.get(1)))
			 {
				 System.err.println("The references of two reads should be same: "+ references.get(0)+ " vs. " + references.get(1));
				 System.exit(-1);
				 
			 }
		 }
		 else if(alleles.size()>2)
		 {
			 System.err.println("The number " + alleles.size() + "  of read of a fragment should not greater than 2!");
			 System.exit(-1);
		 }
		 
		 	int len_2;
			if(snp_pointers_2 == null)
			{
				len_2 = 0;			
			}
			else
				len_2 = snp_pointers_2.size();
			
			List<Integer> snp_p_t = new ArrayList<Integer>();		
			List<Byte> snps = new ArrayList<Byte>();
			StringBuilder strQuality = new StringBuilder();
			int i, j;
			int len_1 = snp_pointers_1.size();
			for(i=0,j=0; i<len_1 && j<len_2; ){
				if(snp_pointers_1.get(i)<snp_pointers_2.get(j)){
					snp_p_t.add(snp_pointers_1.get(i));
					snps.add(alleles_1.get(i));
					strQuality.append(qualities.get(0).get(i));
					i++;
				}
				else if(snp_pointers_1.get(i)>snp_pointers_2.get(j)){
					snp_p_t.add(snp_pointers_2.get(j));
					snps.add(alleles_2.get(j));	
					strQuality.append(qualities.get(1).get(j));
					j++;
				}
				else{
					if(snp_pointers_1.get(i)==snp_pointers_2.get(j)){					
						snps.add(alleles_1.get(i));	
						snp_p_t.add(snp_pointers_1.get(i));	
						strQuality.append(qualities.get(0).get(i));
					}			
							
					i++; 
					j++;			
				}
			}
			while(i<len_1){
				snp_p_t.add(snp_pointers_1.get(i));
				snps.add(alleles_1.get(i));
				strQuality.append(qualities.get(0).get(i));
				i++;			
			}
			while(j<len_2){
				snp_p_t.add(snp_pointers_2.get(j));
				snps.add(alleles_2.get(j));	
				strQuality.append(qualities.get(1).get(j));
				j++;			
			}
			StringBuilder f_l = new StringBuilder();
			int n_b = 0;
			int n_snps = 0;
			len_1 = snp_p_t.size();	
			i = 0;
			int exp_p = -2;  
			StringBuilder qt = new StringBuilder();
			while(i<len_1)
			{
				if(snps.get(i)==(byte)'-')
				{
					i++;
					continue;
				}
				n_snps++;
				if(exp_p != snp_p_t.get(i)) //a new block
				{
					f_l.append("\t");
					f_l.append(snp_p_t.get(i)+1);
					f_l.append("\t");
					f_l.append(Integer.valueOf(snps.get(i)));	
					n_b++;
					exp_p = snp_p_t.get(i)+1;
				}
				else{
					f_l.append(Integer.valueOf(snps.get(i)));	
					exp_p++;
				}
				qt.append(strQuality.charAt(i));
				i++;
			}
			if(n_snps<2){
				return null;
			}
			else
			{
				String readid = rid;
				if(rid.startsWith(">"))
					readid = rid.substring(1);
				f_l.append("\t");			
				f_l.append(qt); 
				f_l.append("\n");
				return(""+n_b+"\t"+ readid +f_l.toString());
			}		
		}
	 
}


