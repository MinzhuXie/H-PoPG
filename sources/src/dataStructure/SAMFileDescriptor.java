package dataStructure;

public class SAMFileDescriptor
{

    public SAMFileDescriptor(int person_index, String filename)
    {
        this.person_index = person_index;
        this.filename = filename;
        sample_prob = 1.0F;
    }

    public SAMFileDescriptor(int person_index, String filename, float sample_prob)
    {
        this.person_index = person_index;
        this.filename = filename;
        this.sample_prob = sample_prob;
    }

    public int getPerson_index()
    {
        return person_index;
    }

    public String getFilename()
    {
        return filename;
    }

    public float getSample_prob()
    {
        return sample_prob;
    }

    int person_index;
    String filename;
    float sample_prob;
}