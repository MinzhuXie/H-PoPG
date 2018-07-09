package algorithms;

import java.util.List;

public class PolypartitionFunctions {
	int nreads, nboxes;
	private int idx_Stir;
	private boolean at_the_end; 
	List<Stirling2> nonEmptyPartitions;
	PolypartitionFunctions(int n, int m, List<Stirling2> nonEmpt)
	{
		nreads = n;
		nboxes = m;
		nonEmptyPartitions = nonEmpt;
		idx_Stir = 0;
		at_the_end = false;
	}
	void reset()
	{
		idx_Stir = 0;
		at_the_end = false;
	}
	boolean nextpartitionFunction(List<Integer> fun)
	{
		boolean ret = nonEmptyPartitions.get(idx_Stir).next_partitionFunction(fun);
		if(ret)
		{
			idx_Stir++;
			if(idx_Stir==nboxes || idx_Stir == nreads )
			{
				at_the_end = true;
			}			
		}
		return at_the_end;		
	}

	
	
	
	
	

}
