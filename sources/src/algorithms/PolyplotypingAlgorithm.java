package algorithms;

import io.VCFProcessor;
import dataStructure.Block;

public interface PolyplotypingAlgorithm {
	public boolean buildHaplotype (Block b);
	public void setParameters(int n_poly, int p_int, double p_double, VCFProcessor vcfp);

}
