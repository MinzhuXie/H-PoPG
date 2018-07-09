package dataStructure;

public class ComplexVariant extends Variant
{

    public ComplexVariant(int reduced_genome_start_position, int genome_start_position, String identifier, String reference, String alleles[])
    {
        super(reduced_genome_start_position, genome_start_position, identifier, reference, -1, alleles);
    }
}