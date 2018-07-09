package algorithms;

import io.VCFProcessor;

import java.util.List;

import dataStructure.*;

public class  HBOP2Builder{
	private int n_ploid;
	private HPBOP2Partition par;
	private List<Fragment> fragments;
	private boolean useQualityScores;
	private HBOP2PartitionAssembly partitions;	
	private double obj_score;	
	private int MAX_SIZE_HEAP;
	private double w_rate;
	VCFProcessor vcf_p;	
	byte[] n_genotypes = null;
	int startpos;
	
	int init_rows = 0;
		
	public HBOP2Builder(Block b, int n_p, int max_size, double weight_factor, VCFProcessor vcf) {
		this(b, n_p, max_size, weight_factor, vcf, false);
	}
	/**
	 * PRE: Fragments are sorted
	 * @param fragments
	 */
	public HBOP2Builder(Block b, int n_p, int max_size, double weight_factor, VCFProcessor vcf, boolean useQual) {
		n_ploid = n_p;		
		//cut=null;
		fragments=b.getFragments();
		startpos = b.getFirstPos();
		MAX_SIZE_HEAP = max_size;
		w_rate = weight_factor;
		useQualityScores = useQual;
		vcf_p = vcf;
		//partitions=null;
		if(vcf!=null)
		{
			if(vcf.ref_to_population_genotypes!=null)
			{
				if(vcf.ref_to_population_genotypes.size()>0)
				{
					String chr = vcf.vcf_reference_names.get(0);
					if(vcf.ref_to_population_genotypes.get(chr)!=null)
					{
						if(vcf.ref_to_population_genotypes.get(chr).size()>0)
						{
							List<byte[]> genotypes = vcf.ref_to_population_genotypes.get(chr).get(0);
							int len = b.length();
							int p = b.getFirstPos(); 
							n_genotypes = new byte[len];
							for(int i = 0; i < len; i++ )
							{								
								byte[] gs = genotypes.get(p++);
								if(gs.length!=n_ploid)
								{
									System.out.println("The vcf file is not correct. Need "+ n_ploid + " alleles, but the genotypes provide by the file contains " + gs.length + " allels." );
									System.exit(-1);
								}
								n_genotypes[i] = 0;
								for(int j = 0; j<gs.length; j++)									
								{
									n_genotypes[i] += gs[j];
								}								
							}
						}
					}					
				}				
			}
		}		
	}	
	
	public boolean first_Row_Partitions(){
		partitions=new HBOP2PartitionAssembly(MAX_SIZE_HEAP, w_rate, n_ploid, n_genotypes, startpos);
		init_rows = partitions.init(0, fragments, useQualityScores);
		if(init_rows == -1)
			return false;
		return true;		
	}
	
	public boolean iterate() {		
		int rows = fragments.size();	
		HBOP2PartitionAssembly P_next;
		for(int i = init_rows;i<rows;i++) {			
			P_next = new HBOP2PartitionAssembly(MAX_SIZE_HEAP, w_rate, n_ploid, n_genotypes, startpos);
			boolean val_return = P_next.extend_from_prev(i, partitions);
			if(!val_return) return false;
			partitions = P_next;			
		}
		return true;
	}	
	
	public void find_max_elem(){		
		par = partitions.find_max_elem();		
		obj_score = partitions.get_max_score();
	}
	
	/**
	 * @return the cut
	 */
	public HPBOP2Partition getPar() {
		return par;
	}
	
	public double getobj_score() {
		return obj_score;
	}
	
	public char[][] getPolyploid(Block b, HPBOP2Partition par, int n_ploid)
    {
		boolean genotypeC = false;
		if(n_genotypes!=null) //n_geno[i] records the number of 1 at the i+startpos SNP locus
			genotypeC = true;		
		int rel = b.getFirstPos();
		int lastPos = b.getLastPos();
		int len = lastPos-rel+1;
		/*		
		if(rel != par.profiles.start_p || len != par.profiles.len)
		{
			System.out.println("There are something Error!");			
		}	
		*/
		if(len != par.profiles.len+par.profiles.start_p-rel)
		{
			System.out.println("There are something Error!");			
		}
		
		int cnt[][] = new int [2*n_ploid][len];		
		
		S_Profiles tmp_pfls = par.profiles;
		int tmp_counts[][]; 
		int relp, i, j, k, st_p, bl, n_g;		
		int map_1[] = new int [n_ploid];
		for(j = 0; j<n_ploid; j++){
			map_1[j] = j;			
		}					
		do{
			tmp_counts = tmp_pfls.counts;
			st_p = tmp_pfls.start_p;
			bl = tmp_pfls.len;
			n_g = tmp_pfls.n_groups;			
			for(i = 0; i < bl; i++)
			{
				relp = i + st_p - rel;
				for(j = 0; j < n_g; j++){					
					cnt[2*map_1[j]][relp] = tmp_counts[2*j][i];
					cnt[2*map_1[j]+1][relp] = tmp_counts[2*j+1][i];					
				}				
			}			
			if(tmp_pfls.prev!=null){
				if(tmp_pfls.map == null){
					tmp_pfls.map = new int[n_ploid];
					for(j = 0; j<n_ploid; j++){
						tmp_pfls.map[j] = j;			
					}					
				}
				for(j = 0; j<n_ploid; j++)
					tmp_pfls.map[j]=map_1[tmp_pfls.map[j]];
				map_1 = tmp_pfls.map;					
				tmp_pfls = tmp_pfls.prev;							
			}
			else
				break;			
		}while(true);	
		int start_position = tmp_pfls.start_p;
		if(rel != start_position )
		{
			System.out.println("There are something Error!");			
		}	
		
		
		char[][] poly_ploid = new char[len][n_ploid];		
		int votes[] = new int[n_ploid];
		int deltas[] = new int[n_ploid];
		int t_votes;	
		int neg_count;
		int pos_count;		
		for(i = 0; i < len; i++){			
			t_votes = 0;
			pos_count = 0;
			neg_count = 0;
			for(k = 0; k < n_ploid; k++)
			{				
				votes[k] = cnt[k*2][i] + cnt[k*2+1][i];
				if(votes[k]>0)
					t_votes ++;	
				deltas[k] =cnt[k*2+1][i] - cnt[k*2][i];	
				if(deltas[k]>0)						
					pos_count++;
				else if(deltas[k]<0)
					neg_count++;
				poly_ploid[i][k] = Fragment.ALLELE1CHAR;
			}								
			if(genotypeC)
			{
				int n_1 = n_genotypes[i]; 
				if( (t_votes >= n_ploid-1) || (pos_count >= n_1) || (neg_count >= n_ploid-n_1))
				{					
					if(n_1 < 1 || n_1>=n_ploid)
					{
						System.out.println("There should be at least one 0 and one 1 at the locus "+ i);
						System.exit(-1);
					}					
					for(k = 0; k < n_1; k++)
					{				
						int max_idx = 0;
						for(j = 1; j < n_ploid; j++ )
						{
							if(deltas[j] > deltas[max_idx])
								max_idx = j;
						}
						poly_ploid[i][max_idx] = Fragment.ALLELE2CHAR;				
						deltas[max_idx] = -100000; 
					}
				}				
				else
				{
					for(j = 0; j < n_ploid; j++ )
					{
						if(deltas[j] < 0){
							poly_ploid[i][j] = Fragment.ALLELE1CHAR;							
						}
						else if(deltas[j]>0){
							poly_ploid[i][j] = Fragment.ALLELE2CHAR;							
						}
						else					
							poly_ploid[i][j] = Fragment.NODATACHAR;
					}					
				}
			}
			else{				
				for(j = 0; j<n_ploid; j++)
				{
					if(votes[j] == 0)
						poly_ploid[i][j] = Fragment.NODATACHAR;
					else{
						if (deltas[j]>0) {
							poly_ploid[i][j] = Fragment.ALLELE2CHAR;
						} else if (deltas[j] <= 0){
							poly_ploid[i][j] = Fragment.ALLELE1CHAR;
						} 
//						else {
//							poly_ploid[i][j] = Fragment.NODATACHAR;
//						}						
					}							
				}	
			}			
		}
		return poly_ploid;
	}	
}
