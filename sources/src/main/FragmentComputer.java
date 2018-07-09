package main;

import java.io.*;
import java.util.*;

import shares.Tools;

import dataStructure.ReducedFragment;
import dataStructure.SAMFileDescriptor;

import net.sf.samtools.*;

public class FragmentComputer
{
	 private BufferedWriter frag_writer;
	 private int number_of_alleles_computed;
	 private PolyPlotyping ha;
	 private String suffix_to_remove;
	 public int processed_sam_files;

	 public FragmentComputer(String frag_writer_filename, String suffix_to_remove, PolyPlotyping ha)
     throws IOException
	 {
	     frag_writer = null;	    
	     number_of_alleles_computed = 0; 	    
	     this.suffix_to_remove = null;
	     processed_sam_files = 0;
	     frag_writer = new BufferedWriter(new FileWriter(frag_writer_filename));
	     this.suffix_to_remove = suffix_to_remove;
	     this.ha = ha;	    
	 }

	
	 protected int go()
	     throws IOException
	 {
	     Random r = new Random();
	     SAMFileReader sreader = null;	   
	     for(SAMFileDescriptor samDescriptor : (List<SAMFileDescriptor>)ha.sam_file_descriptors)
	     {	      
	         String sam_file = samDescriptor.getFilename();
	         processed_sam_files++;
	         sreader = new SAMFileReader(new File(sam_file));
	         sreader.setValidationStringency(net.sf.samtools.SAMFileReader.ValidationStringency.SILENT);
	         Map<String, ReducedFragment> name_to_fragments = new HashMap<String, ReducedFragment>();
	         SAMRecord samLine = null;
	         SAMRecordIterator samIterator = sreader.iterator();
	         String readName = null;
	         String readRef = null;
	         Set<String> start_or_end_seen = new HashSet<String>();
	         int fragments = 0;
	         while(samIterator.hasNext()) 
	         {
	             samLine = (SAMRecord)samIterator.next();
	             readName = Tools.getReadName(samLine, suffix_to_remove);
	             readRef = new String(samLine.getReferenceName());
	             if(Tools.is_primary(samLine))
	             {
	                 if(!ha.sam_processor.ref_to_variant_ptrs.containsKey(readRef))
	                 {
	                     ha.sam_processor.ref_to_variant_ptrs.put(readRef, Integer.valueOf(0));
	                     ha.sam_processor.ref_to_before_variant_ptrs.put(readRef, Integer.valueOf(0));
	                 }
	                 List<SAMRecord> samList = new ArrayList<SAMRecord>(1);
	                 samList.add(samLine);
	                 //obtain a ruduced fragment containing only SNPs
	                 ReducedFragment frag = ha.sam_processor.processRead(samList, 0);
	                 if(frag.getNum_alleles() > 0)
	                     if(!name_to_fragments.containsKey(readName))
	                         name_to_fragments.put(readName, frag);
	                     else
	                         ((ReducedFragment)name_to_fragments.get(readName)).addFragment(frag);
	                 if(samLine.getReadPairedFlag() && 
	                 		(samLine.getFirstOfPairFlag() || samLine.getSecondOfPairFlag()) && 
	                 		!start_or_end_seen.add(readName))
	                 {
	                     if(name_to_fragments.containsKey(readName))
	                     {
	                         if(((ReducedFragment)name_to_fragments.get(readName)).getNum_alleles() >= 2 && r.nextFloat() < samDescriptor.getSample_prob())
	                         {                   	 
	                        	 
	                        	 String lineHapCutFormat = name_to_fragments.get(readName).merge(readName);
	                          //   frag_writer.write(name_to_fragments.get(readName).getNumberOfReads()+"\t"+readName+"\t"+name_to_fragments.get(readName).toString());
	                        	 frag_writer.write(lineHapCutFormat);	                             
	                             
	                             number_of_alleles_computed += ((ReducedFragment)name_to_fragments.get(readName)).getNum_alleles();
	                         }
	                         name_to_fragments.remove(readName);
	                     }
	                     start_or_end_seen.remove(readName);
	                 }
	                 if(!samLine.getReadPairedFlag() && name_to_fragments.containsKey(readName))
	                 {
	                     if(((ReducedFragment)name_to_fragments.get(readName)).getNum_alleles() >= 2 && r.nextFloat() < samDescriptor.getSample_prob())
	                     {
	                    	 String lineHapCutFormat = name_to_fragments.get(readName).merge(readName);
	                       	 frag_writer.write(lineHapCutFormat);	
	                       //  frag_writer.write((new StringBuilder(String.valueOf(((ReducedFragment)name_to_fragments.get(readName)).getNumberOfReads()))).append("\t").append(readName).append("\t").append(((ReducedFragment)name_to_fragments.get(readName)).toString()).toString());
	                         number_of_alleles_computed += ((ReducedFragment)name_to_fragments.get(readName)).getNum_alleles();
	                     }
	                     name_to_fragments.remove(readName);
	                 }
	                 if(++fragments % 1000000 == 0)
	                     System.out.println((new StringBuilder("Processed ")).append(fragments).append(" fragments.").toString());
	             }
	         }
	         for(Iterator<String> iterator1 = name_to_fragments.keySet().iterator(); iterator1.hasNext();)
	         {
	             String key = (String)iterator1.next();
	             ReducedFragment frag = (ReducedFragment)name_to_fragments.get(key);
	             if(frag.getNum_alleles() >= 2 && r.nextFloat() < samDescriptor.getSample_prob())
	             {
	            	 String lineHapCutFormat = frag.merge(key);
                   	 frag_writer.write(lineHapCutFormat);	
	                // frag_writer.write((new StringBuilder(String.valueOf(((ReducedFragment)name_to_fragments.get(key)).getNumberOfReads()))).append("\t").append(key).append("\t").append(((ReducedFragment)name_to_fragments.get(key)).toString()).toString());
	                 number_of_alleles_computed += ((ReducedFragment)name_to_fragments.get(key)).getNum_alleles();
	             }
	         }
	
	         System.out.println("Read BAM File completed: "+sam_file);
	         ha.sam_processor.loaded_sam_files++;
	     }
	     if(sreader!=null)
	    	 sreader.close();
	
	     frag_writer.flush();
	     frag_writer.close();
	     System.out.println("BAM file --> Frag file completed. "+ "Total effective calles: "+ number_of_alleles_computed);
	     return number_of_alleles_computed;
	 }


}

