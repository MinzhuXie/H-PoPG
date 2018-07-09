package dataStructure;

import java.util.Arrays;

public abstract class Variant
  implements Comparable<Variant>
{
  public int reduced_genome_start_position;
  public int genome_start_position;
  public String identifier;
  public String reference;
  public int readCount;
  public int length;
  public String[] alleles;
  public byte[] allele_encodings;
  
  public Variant(int reduced_genome_start_position, int genome_start_position, String identifier, String reference, int length, String[] alleles)
  {
    this.reduced_genome_start_position = reduced_genome_start_position;
    this.genome_start_position = genome_start_position;
    this.identifier = identifier;
    this.reference = reference;
    this.length = length;
    this.alleles = alleles;
    this.allele_encodings = new byte[alleles.length];
    for (byte i = 0; i < alleles.length; i = (byte)(i + 1)) {
      this.allele_encodings[i] = i;
    }
    this.readCount = 0;
  }
  
  public int compareTo(Variant o)
  {
    if (this.reference.compareTo(o.reference) < 0)
      return -1;
    if (this.reference.compareTo(o.reference) > 0) {
      return 1;
    }
    if (this.genome_start_position < o.genome_start_position)
      return -1;
    if (this.genome_start_position > o.genome_start_position)
      return 1;
    return 0;
  }
  
  public Integer getType()
  {
    if (this.alleles.length > 2)
      return Integer.valueOf(2147483647);
    if (this.alleles[0].length() > this.alleles[1].length())
      return Integer.valueOf(-1);
    if (this.alleles[0].length() < this.alleles[1].length()) {
      return Integer.valueOf(1);
    }
    return Integer.valueOf(0);
  }
  
  public Integer getEffectiveLength()
  {
    if (this.alleles[0].length() == this.alleles[1].length())
      return Integer.valueOf(1);
    if (this.alleles[0].length() > this.alleles[1].length())
    {
      return Integer.valueOf(this.alleles[0].length()); }
    if (this.alleles[0].length() < this.alleles[1].length())
    {
      return Integer.valueOf(1);
    }
    return null;
  }
  
  public int getSharedSegmentLength()
  {
    int minLen = 0;
    for (String allele : this.alleles) {
      if (allele.length() < minLen)
        minLen = allele.length();
    }
    return minLen;
  }
  
  public String toString()
  {
    StringBuilder builder = new StringBuilder();
    builder.append(this.reference);
    builder.append("\t");
    builder.append(this.genome_start_position);
    builder.append("\t");
    builder.append(this.identifier);
    builder.append("\t");
    builder.append(Arrays.toString(this.alleles));
    return builder.toString();
  }
  
  public int getLongestVariantLength() {
    int length = 0;
    for (String allele : this.alleles) {
      if (allele.length() > length)
        length = allele.length();
    }
    return length;
  }
}


