package io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import main.PolyPlotyping;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import dataStructure.ReducedFragment;
import dataStructure.Variant;

public class BAMfileReader
{
  private static int ERRORMASK = 1796;
  private static int PRIMARYMASK = 256;
  public VCFProcessor vcf_processor = null;
 // public DownsampleSam downsampler = null;
  
 // private static Logger logger = LogManager.getLogger(HaplotypeAssembler.class.getName());
  
  public Map<String, Integer> ref_to_variant_ptrs = new HashMap<String, Integer>();
  public Map<String, Integer> ref_to_before_variant_ptrs = new HashMap<String, Integer>();
  
  public int base_mismatch_counter = 0;
  public int base_match_counter = 0;
  public int base_match_good_read_counter = 0;
  public int number_reads_with_bad_flags = 0;
  public int number_reads_with_no_reference_alignment = 0;
  public int paired_reads_mapping_to_same_chromosome = 0;
  public int paired_reads_mapping_to_diff_chromosome = 0;
  public int skipped_non_primary_alignments = 0;
  public int readsPassingFilters = 0;
  public int readsPassingFiltersUsedInAssembly = 0;
  public int readsPassingFiltersWrongReference = 0;
  public int loaded_sam_files = 0;
  
  public static int NUM_READ_CHECK = 1000000;
  
  private String ref_name = null;
  
  public PolyPlotyping ha = null;
  
  public BAMfileReader(PolyPlotyping ha) {
    this.ha = ha;
   }
  
  public void init(VCFProcessor vcf_processor, String fragment_filename) 
  {
    this.vcf_processor = vcf_processor;
//    this.downsampler = new DownsampleSam();
  }
  
  public void init(VCFProcessor vcf_processor, String ref_name, String suffix_to_remove, String fragment_filename)
  {
    this.vcf_processor = vcf_processor;
    this.ref_name = ref_name;
//    this.suffix_to_remove = suffix_to_remove;
//    this.downsampler = new DownsampleSam();
  }  

  private boolean pass_qc(SAMRecord samRecord)
  {
    if (badRead(samRecord.getFlags())) {
     // if (!Utilities.is_primary(samRecord)) 
      if(!((samRecord.getFlags() & PRIMARYMASK) <= 0))
      {
        this.skipped_non_primary_alignments += 1;
      }
      this.number_reads_with_bad_flags += 1;
      return false; }
    if (!this.vcf_processor.vcf_reference_names.contains(samRecord.getReferenceName())) {
      this.number_reads_with_no_reference_alignment += 1;
      return false; }
    if ((this.ref_name != null) && (!samRecord.getReferenceName().equals(this.ref_name))) {
      this.number_reads_with_no_reference_alignment += 1;
      return false;
    }
    return true;
  }
  
  public ReducedFragment processRead(List<SAMRecord> samLines, int person_index) throws IOException
  {
    ReducedFragment rf = new ReducedFragment();
    List<List<SAMRecord>> superSamLines = new ArrayList<List<SAMRecord>>();
    
    String ref = ((SAMRecord)samLines.get(0)).getReferenceName();
    int this_not_matching = 0;
    for (int i = 1; i < samLines.size(); i++) {
      if (!ref.equals(((SAMRecord)samLines.get(i)).getReferenceName())) {
        this_not_matching++;
      }
    }
    if (this_not_matching > 0) {
      this.paired_reads_mapping_to_diff_chromosome += this_not_matching;
      
      for (SAMRecord rec : samLines) {
        boolean placed = false;
        for (List<SAMRecord> list : superSamLines) {
          if (((SAMRecord)list.get(0)).getReferenceName().equals(rec.getReferenceName())) {
            placed = true;
            list.add(rec);
          }
        }
        if (!placed) {
          superSamLines.add(new ArrayList<SAMRecord>());
          ((List<SAMRecord>)superSamLines.get(superSamLines.size() - 1)).add(rec);
        }
      }
    } else {
      superSamLines.add(samLines);
    }
    
    for (List<SAMRecord> record_group : superSamLines) {
      List<List<Integer>> variant_pointers = new ArrayList<List<Integer>>();
      List<List<Byte>> alleles = new ArrayList<List<Byte>>();
      List<List<String>> qualities = new ArrayList<List<String>>();
      List<String> refs = new ArrayList<String>();
      Set<Integer> positions_added = new HashSet<Integer>();
      variant_pointers.add(new ArrayList<Integer>());
      alleles.add(new ArrayList<Byte>());
      qualities.add(new ArrayList<String>());
      int goodVariants = 0;
      int it_number = 0;
      int added_this_read = 0;
      String read_chr = null;
      
      List<Integer> matchedVariantPtrs = new ArrayList<Integer>();
      
      for (SAMRecord samLine : record_group) {
        read_chr = samLine.getReferenceName();
        if (pass_qc(samLine)) //if the algined read pass some rules: chr, qualtiy, etc.
        {
          this.readsPassingFilters += 1;
          if (this.vcf_processor.vcf_reference_names.contains(read_chr)) {
//            int totalInsertionLength = 0;
 //           int totalDeletionLength = 0;
//            int start = samLine.getAlignmentStart();
            int current_pos_read = 1;
            int current_pos_ref = samLine.getAlignmentStart();
            String sequence = samLine.getReadString();
            List<CigarElement> cigarElements = samLine.getCigar().getCigarElements();
            for (int i = 0; i < cigarElements.size(); i++) {
              CigarElement e = (CigarElement)cigarElements.get(i);
              switch (e.getOperator().ordinal()) 
              {
              case 5: break;
              case 6: break;
              case 4: current_pos_read += e.getLength();break;
              case 3: current_pos_ref += e.getLength();break;
              
              case 2: 
//                totalDeletionLength += e.getLength();
                current_pos_ref += e.getLength();break;
              
              case 1: 
//                totalInsertionLength += e.getLength();
                current_pos_read += e.getLength();break;
              case 0: 
              case 7: 
              case 8: 
                matchedVariantPtrs = getAllVariants(current_pos_ref, e.getLength(), read_chr);
                for (Integer varptr : matchedVariantPtrs) {
                  if (positions_added.add(varptr)) {
                    if (i == cigarElements.size() - 1)
                    {
                      if (processSequenceMatchMismatch(
                    		  varptr.intValue(), 
                    		  sequence.substring(current_pos_read - 1, current_pos_read - 1 + e.getLength()), 
                    		  vcf_processor.ref_to_variant_list.get(read_chr).get(varptr.intValue()).genome_start_position - current_pos_ref, //- totalInsertionLength, 
                    		  variant_pointers.get(it_number), 
                    		  alleles.get(it_number), 
                    		  read_chr, 
                    		  qualities.get(it_number), 
                    		  samLine.getBaseQualityString())) 
                      {
                        goodVariants++;
                        added_this_read++;
                      } else {
                        positions_added.remove(varptr);
                      }
                      
                    }
                    else if (processSequenceMatchMismatch(varptr.intValue(), 
                    		sequence.substring(current_pos_read - 1, current_pos_read - 1 + e.getLength()), 
                    		vcf_processor.ref_to_variant_list.get(read_chr).get(varptr.intValue()).genome_start_position - current_pos_ref,  // - totalInsertionLength, 
                    		variant_pointers.get(it_number), 
                    		alleles.get(it_number), read_chr, 
                    		(CigarElement)cigarElements.get(i + 1), 
                    		current_pos_read + e.getLength(), 
                    		qualities.get(it_number), 
                    		samLine.getBaseQualityString())) 
                    {
                      goodVariants++;
                      added_this_read++;
                    } else {
                      positions_added.remove(varptr);
                    }
                  }
                }                
                current_pos_read += e.getLength();
                current_pos_ref += e.getLength();
                break;
              default:  throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
              }
              
            }
            int absolute_align_start = -1;
            int mateAlignmentStart = samLine.getMateAlignmentStart();
            if (mateAlignmentStart == 0)
              mateAlignmentStart = 2147483647;
            if (samLine.getAlignmentStart() < mateAlignmentStart)
              absolute_align_start = samLine.getAlignmentStart(); else
              absolute_align_start = samLine.getMateAlignmentStart();
            while ((ref_to_before_variant_ptrs.get(samLine.getReferenceName()).intValue() 
            		< vcf_processor.ref_to_variant_list.get(read_chr).size()) 
            		&& 
            (vcf_processor.ref_to_variant_list.get(read_chr).get(ref_to_before_variant_ptrs.get(samLine.getReferenceName()).intValue()).genome_start_position 
              < absolute_align_start)) 
            {
              ref_to_before_variant_ptrs.put(samLine.getReferenceName(), Integer.valueOf((ref_to_before_variant_ptrs.get(samLine.getReferenceName())).intValue() + 1));
            }
            ref_to_variant_ptrs.put(samLine.getReferenceName(), ref_to_before_variant_ptrs.get(samLine.getReferenceName()));
          } 
          else 
          { 
        	  readsPassingFiltersWrongReference += 1; 
         }
          if (added_this_read > 0) 
          {
            it_number++;
            refs.add(read_chr);
            variant_pointers.add(new ArrayList<Integer>());
            alleles.add(new ArrayList<Byte>());
            qualities.add(new ArrayList<String>());
          }
          added_this_read = 0;
        }
      }
      if (goodVariants > 0) 
      {
        variant_pointers.remove(variant_pointers.size() - 1);
        alleles.remove(alleles.size() - 1);
        qualities.remove(qualities.size() - 1);
        rf.addFragment(alleles.get(0), variant_pointers.get(0), qualities.get(0), samLines.get(0).getReferenceName());
      }
    }
    
    return rf;
  }
  
  private boolean badRead(int flags)
  {
    if ((flags & ERRORMASK) > 0)
      return true;
    return false;
  }
  

  
  private boolean processSequenceMatchMismatch(int var_ptr, String seq, int location_of_variant, List<Integer> variant_pointers, List<Byte> alleles, String reference, CigarElement c, int endOfCIGAR, List<String> qualities, String qual_seq)
  {
    int added_cnt = 0;
    Variant var = vcf_processor.ref_to_variant_list.get(reference).get(var_ptr);
    if ((location_of_variant + var.getLongestVariantLength() > seq.length()) || (location_of_variant < 0))
      return false;
    String allele_string = seq.substring(location_of_variant, location_of_variant + var.alleles[0].length());
    boolean matched = false;
    if (var.getType().intValue() == 0) {
      for (int i = 0; i < var.alleles.length; i++)
        if (allele_string.equals(var.alleles[i])) {
          matched = true;
          alleles.add(Byte.valueOf((byte)i));
          if (!qual_seq.equals("*"))
            qualities.add(qual_seq.substring(location_of_variant, location_of_variant + var.alleles[i].length())); else
            qualities.add("*");
          added_cnt++;
        }
    } else {
      if (var.getType().intValue() == 2147483647) {
        return false;
      }
      if (location_of_variant + var.getSharedSegmentLength() == endOfCIGAR)
      {
        if (c.getOperator().equals(CigarOperator.D))
        {
          if (c.getLength() >= var.length) {
            matched = true;
            alleles.add(Byte.valueOf((byte)1));
      //      qualities.add(StringUtils.leftPad("", var.length, '*'));
          }
        } else if (c.getOperator().equals(CigarOperator.I))
        {
          if (c.getLength() >= var.length) {
            matched = true;
            alleles.add(Byte.valueOf((byte)0));
            if (!qual_seq.equals("*"))
              qualities.add(qual_seq.substring(location_of_variant, location_of_variant + var.length)); else
              qualities.add("*");
          }
        }
      }
    }
    if (added_cnt > 1) {
//      Utilities.printErrorAndExit(
 //       "Found more than one allele matching reference sequence.\n" + seq + "\n" + var.toString(), new Exception(), logger);
    }
    
    if (matched) {
      vcf_processor.ref_to_variant_list.get(reference).get(var_ptr).readCount += 1;
      variant_pointers.add(Integer.valueOf(var_ptr));
      this.base_match_counter += 1;
      return true;
    }
    this.base_mismatch_counter += 1;
    return false;
  }
  
  private boolean processSequenceMatchMismatch(
		  int var_ptr, 
		  String seq, 
		  int location_of_variant, 
		  List<Integer> variant_pointers, 
		  List<Byte> alleles, 
		  String reference, 
		  List<String> qualities,
		  String qual_seq)
  {
    int added_cnt = 0;
    Variant var = vcf_processor.ref_to_variant_list.get(reference).get(var_ptr);
    if ((location_of_variant + var.getEffectiveLength().intValue() > seq.length()) || (location_of_variant < 0))
      return false;
    String allele_string = seq.substring(location_of_variant, location_of_variant + var.alleles[0].length());
    boolean matched = false;
    if (var.getType().intValue() == 0) {  //SNP
      for (int i = 0; i < var.alleles.length; i++) {
        if (allele_string.equals(var.alleles[i])) {
          matched = true;
          alleles.add(Byte.valueOf((byte)i));
          if (!qual_seq.equals("*"))
            qualities.add(qual_seq.substring(location_of_variant, location_of_variant + var.alleles[i].length())); 
          else
            qualities.add("*");
          added_cnt++;
        }
      }
    } else { 
      return false;
    }
    if (added_cnt > 1) {
//      Utilities.printErrorAndExit(
 //       "Found more than one allele matching reference sequence.\n" + seq + "\n" + var.toString(), new Exception(), logger);
    }
    
    if (matched) {
      vcf_processor.ref_to_variant_list.get(reference).get(var_ptr).readCount += 1;
      variant_pointers.add(Integer.valueOf(var_ptr));
      this.base_match_counter += 1;
      return true;
    }
    this.base_mismatch_counter += 1;
    return false;
  }
  
  
  private List<Integer> getAllVariants(int current_pos_ref, int length, String local_ref_name)
  {
    if ((this.ref_name != null) && (!local_ref_name.equals(this.ref_name))) {
      return Collections.emptyList();
    }     
   
    List<Variant> var_list = vcf_processor.ref_to_variant_list.get(local_ref_name);
    int n_variants = var_list.size();
    
    while (ref_to_variant_ptrs.get(local_ref_name) < n_variants) 
    {
    	if(var_list.get(ref_to_variant_ptrs.get(local_ref_name)).genome_start_position < current_pos_ref) 
    	{
    		int pre_order = ref_to_variant_ptrs.get(local_ref_name);
            this.ref_to_variant_ptrs.put(local_ref_name, pre_order + 1);
        }
    	else 
    		break;
    }
    
    if (ref_to_variant_ptrs.get(local_ref_name) == n_variants) {
      return Collections.emptyList();
    }
    
    List<Integer> vars = new ArrayList<Integer> ();
    
    int v_order = ref_to_variant_ptrs.get(local_ref_name);
        
    while ((v_order != n_variants) && (var_list.get(v_order).genome_start_position < current_pos_ref + length)) 
    {
       Variant var = var_list.get(v_order);      
      if ((var_list.get(v_order).genome_start_position >= current_pos_ref) && 
        (var_list.get(v_order).genome_start_position + var.getLongestVariantLength() <= current_pos_ref + length) && 
        (var_list.get(v_order).reference.equals(local_ref_name)))
      {
    	  vars.add(v_order);
      }
      v_order++;
      ref_to_variant_ptrs.put(local_ref_name, v_order);
    }
    return vars;
  } 
 
}
