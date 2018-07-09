package algorithms;

import io.VCFProcessor;
import dataStructure.Block;

public class HPBOP2Alg implements PolyplotypingAlgorithm {
	private int n_ploid; 
	private HPBOP2Partition par;
	private char[][] polyploid;	
	private int MAX_HEAP;	
	private double w_rate;
	VCFProcessor vcf_p;
	
	@Override
	public void setParameters(int n_p, int max_size, double weight_factor, VCFProcessor vcfp){
		n_ploid = n_p;
		MAX_HEAP = max_size;
		w_rate = weight_factor;
		vcf_p = vcfp;
	}	
	@Override
	public boolean buildHaplotype(Block b) {
		// TODO Auto-generated method stub
		HBOP2Builder builder = new HBOP2Builder(b, n_ploid, MAX_HEAP, w_rate, vcf_p);
		boolean ret_val;
		ret_val=builder.first_Row_Partitions();
		if(!ret_val)return false;
		try{
			ret_val = builder.iterate();
			if(!ret_val) return false;
		}
		catch (OutOfMemoryError e)
		{
			return (false);
		}
		builder.find_max_elem();
		par = builder.getPar();			
		b.object_opt_score = builder.getobj_score();		
		polyploid = builder.getPolyploid(b, par, n_ploid);
		b.setHaps(polyploid,n_ploid);	
		return true;		
	}
	
}
