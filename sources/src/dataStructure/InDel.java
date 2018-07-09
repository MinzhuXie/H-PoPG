package dataStructure;

public class InDel extends Variant
{

    public InDel(int reduced_genome_start_position, int genome_start_position, String identifier, String reference, String alleles[])
    {
        super(reduced_genome_start_position, genome_start_position, identifier, reference, alleles[1].length() - alleles[0].length(), alleles);
    }

    public boolean isDeletion()
    {
        return alleles[0].length() > alleles[1].length();
    }
}