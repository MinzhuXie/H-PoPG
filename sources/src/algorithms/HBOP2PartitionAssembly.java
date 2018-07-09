package algorithms;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import dataStructure.Fragment;

public class HBOP2PartitionAssembly{
	static final int split_len = 20;
	private int MAX_SIZE_HEAP;	
	private int current_row; //当前的片段的序号，从0开始
	private int n_rows_overlap; //储存与当前片段有重叠的前面的片段数 
	private int[] ids_rows_overlap; //储存与当前片段有重叠的前面的片段的序号 
	HBOPMinHeap ps;  //当前片段对应的分区的得分值
	private double max_score=0;
	private List<Fragment> fragments;
	private boolean useQualityScores;
	private double w_rate ;
	private int k_ploid;
	private byte[] n_geno;
	int startpos;
	
	/*
	public static final int [] Bits_set={0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 
		                                 0x100, 0x200, 0x400, 0x800,0x1000, 0x2000, 0x4000, 0x8000,
		                                 0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000,
		                                 0x1000000, 0x2000000, 0x4000000, 0x8000000,0x10000000, 0x20000000,
		                                 0x40000000, 0x80000000};
	public static final int [] Bits_clr={0xfffffffe, 0xfffffffd, 0xfffffffb, 0xfffffff7, 
										0xffffffef, 0xffffffdf, 0xffffffbf, 0xffffff7f,
										0xfffffeff, 0xfffffdff, 0xfffffbff, 0xfffff7ff,
										0xffffefff, 0xffffdfff, 0xffffbfff, 0xffff7fff,
										0xfffeffff, 0xfffdffff, 0xfffbffff, 0xfff7ffff,
										0xffefffff, 0xffdfffff, 0xffbfffff, 0xff7fffff,
										0xfeffffff, 0xfdffffff, 0xfbffffff, 0xf7ffffff,
										0xefffffff, 0xdfffffff, 0xbfffffff, 0x7fffffff};
	*/
	HBOP2PartitionAssembly(int max_size, double weight_factor, int n_ploid, byte[] genotypes, int st){
		MAX_SIZE_HEAP = max_size;
		w_rate = weight_factor;
		k_ploid = n_ploid;
		n_geno = genotypes;
		startpos = st;
	}
	
	
	
	//compute the haplotypes
	void comp_Hap(int[][] counts, char[][] haps, int maxg, int len)
	{		
		for(int i = 0; i < maxg; i++)
		{
			int ind1 = 2*i; 
			for(int j = 0; j < len; j++)
			{			
				if(counts[ind1][j] > counts[ind1+1][j])
				{
					haps[i][j] = Fragment.ALLELE1CHAR;					
				}
				else if(counts[ind1][j] < counts[ind1+1][j])
				{
					haps[i][j] = Fragment.ALLELE2CHAR;					
				}
				else            //if(counts[ind1][j] == counts[ind1+1][j])
				{
					if(counts[ind1][j] == 0)
						haps[i][j] = Fragment.NODATACHAR;
					else
					{
						haps[i][j] = Fragment.ALLELE1CHAR;					
					}					
				}				 
			}
		}		
	}
	
	//compute MEC
	int comp_MEC(int[][] counts, int maxg, int len)
	{
		int mec_score = 0;		
		for(int i = 0; i < maxg; i++)
		{
			int ind1 = 2*i; 
			for(int j = 0; j < len; j++)
			{			
				if(counts[ind1][j] >= counts[ind1+1][j])
				{					
					mec_score+=counts[ind1+1][j];
				}
				else if(counts[ind1][j] < counts[ind1+1][j])
				{					
					mec_score+=counts[ind1][j];
				}				 
			}
		}		
		return mec_score;
	}
	/*
	private int comp_genotype_Constrained_MEC(S_Profiles profiles, int startcol, int endcol)
	{
		int [][] counts = profiles.counts;
		int n_g = profiles.n_groups;
		int begin_p = profiles.start_p;
		int end_p = profiles.start_p+profiles.len-1;		
		if(startcol<begin_p){
			System.out.println("Error");
			startcol = begin_p;			
		}
		if(startcol>end_p){
			System.out.println("Error");
			return 0;			
		}
		if(endcol < startcol)
		{
			System.out.println("Error");
			return 0;			
		}
		if(endcol > end_p)
		{
			endcol = end_p;			
		}		
		endcol++;
		int[] deltas = new int[n_g];
		int mec = 0;
		for(int i = startcol; i < endcol; i++ )
		{
			int n_1 = n_geno[i-startpos];
			if(n_1<1 || n_1>=k_ploid)
			{
				System.out.println("There should be at least one 0 and one 1 at the locus "+ (startcol+i));
				System.exit(-1);
			}
			
			for(int j = 0; j < n_g; j++ )
			{
				deltas[j] =  counts[2*j+1][i-begin_p] - counts[2*j][i-begin_p];      // Bs[j][i] - As[j][i];
				mec += counts[2*j+1][i-begin_p];     //assume all 0
			}
			
			int max_n_1 = n_1, min_n_1 = n_1;
			if(n_g < k_ploid)   //how to deal with it?  proportion ? or min the MEC score?
			{
				//n_1 = (int)Math.rint(((float)(n_g*n_1))/k_ploid);	
				max_n_1 = n_g > n_1 ? n_1 : n_g;
				min_n_1 = n_1 + n_g - k_ploid;
				if(min_n_1<0)
					min_n_1 = 0;
			}
			
			int k = 0;
			for( ; k < min_n_1; k++)
			{				
				int max_idx = 0;
				for(int j = 1; j < n_g; j++ )
				{
					if(deltas[j] > deltas[max_idx])
						max_idx = j;
				}
				mec -= deltas[max_idx];				
				deltas[max_idx] = -100000; 
			}
			
			for( ; k < max_n_1; k++)
			{				
				int max_idx = 0;
				for(int j = 1; j < n_g; j++ )
				{
					if(deltas[j] > deltas[max_idx])
						max_idx = j;
				}
				if(deltas[max_idx]>0){
					mec -= deltas[max_idx];				
					deltas[max_idx] = -100000; 
				}								
				else
					break;	
			}
		}		
		return mec;		
	}	
	*/
	private int[] comp_genotype_Constrained_Scores(S_Profiles profiles, int startcol, int endcol)
	{
		int[] scores = new int[2];
		scores[0] = 0;
		scores[1] = 0;
		int [][] counts = profiles.counts;
		int n_g = profiles.n_groups;
		int begin_p = profiles.start_p;
		int end_p = profiles.start_p+profiles.len-1;		
		if(startcol<begin_p){
			System.out.println("Error");
			startcol = begin_p;			
		}
		if(startcol>end_p){
			System.out.println("Error");
			return scores;			
		}
		if(endcol < startcol)
		{
			System.out.println("Error");
			return scores;			
		}
		if(endcol > end_p)
		{
			endcol = end_p;			
		}		
		endcol++;
		int[] deltas = new int[n_g];
		char[] tmp_genotypes = new char[n_g];
		
		boolean[] undetermined = new boolean[n_g];
		
		int mec = 0;
		int distance = 0;
		for(int i = startcol; i < endcol; i++ )
		{
			int n_1 = n_geno[i-startpos];
			if(n_1<1 || n_1>=k_ploid)
			{
				System.out.println("There should be at least one 0 and one 1 at the locus "+ (startcol+i));
				System.exit(-1);
			}
			int idx = i - begin_p;
			for(int j = 0; j < n_g; j++ )
			{
				if(counts[2*j+1][idx]+counts[2*j][idx]>0)
					undetermined[j] = false;
				else
					undetermined[j] = true;		
				
				deltas[j] =  counts[2*j+1][idx] - counts[2*j][idx];      // Bs[j][i] - As[j][i];
				mec += counts[2*j+1][idx];     //assume all 0
		//		if(counts[2*j+1][idx] == 0 && counts[2*j+1][idx] == 0)
		//			tmp_genotypes[j] = '-';
		//		else
					tmp_genotypes[j] = '0';
			}
			
			int max_n_1 = n_1, min_n_1 = n_1;
			if(n_g < k_ploid)   //how to deal with it?  proportion ? or min the MEC score?
			{
				//n_1 = (int)Math.rint(((float)(n_g*n_1))/k_ploid);	
				max_n_1 = n_g > n_1 ? n_1 : n_g;
				min_n_1 = n_1 + n_g - k_ploid;
				if(min_n_1<0)
					min_n_1 = 0;
			}
			
			int k = 0;
			for( ; k < min_n_1; k++)
			{				
				int max_idx = 0;
				for(int j = 1; j < n_g; j++ )
				{
					if(deltas[j] > deltas[max_idx])
						max_idx = j;
				}
				mec -= deltas[max_idx];	
				tmp_genotypes[max_idx] = '1';
				deltas[max_idx] = -100000; 
			}
			
			for( ; k < max_n_1; k++)
			{				
				int max_idx = 0;
				for(int j = 1; j < n_g; j++ )
				{
					if(deltas[j] > deltas[max_idx])
						max_idx = j;
				}
				if(deltas[max_idx]>0){
					mec -= deltas[max_idx];	
					tmp_genotypes[max_idx] = '1';
					deltas[max_idx] = -100000; 
				}								
				else
					break;	
			}
			
			for(int j = 0; j < n_g; j++ )
			{
				if(undetermined[j])
					tmp_genotypes[j] = '-';
			}
			
			for(int j = 0; j < n_g; j++ ){
				if(tmp_genotypes[j]=='-')
					continue;
				for(k = j+1; k < n_g;  k++){					
					if(tmp_genotypes[k]!='-'&& tmp_genotypes[j]!=tmp_genotypes[k])
						distance++;	
				}				
			}			
		}		
		scores[0] = mec;
		scores[1] = distance;
		return scores;		
	}	
	
	
	//fragments  ids_rows_overlap  n_rows_overlap useQualityScores w_rate -->
	HPBOP2Partition calScore(List<Integer> fun)
	{		
		HPBOP2Partition tmp_P = new HPBOP2Partition();
		tmp_P.mec_score = 0 ;
		tmp_P.mix_score = 0;
		tmp_P.polypartition = new byte[n_rows_overlap-1];
		tmp_P.profiles = new S_Profiles();
		
		int i, j, maxg = 0, l_p, r_p;
		
		//the first row must be in group 0 and is not recorded
		for(i=1; i<n_rows_overlap; i++)
		{
			tmp_P.polypartition[i-1] = fun.get(i).byteValue();	
			if(maxg<tmp_P.polypartition[i-1])
				maxg = tmp_P.polypartition[i-1];			
		}
		maxg ++;
		Fragment tmp_f =  fragments.get(ids_rows_overlap[0]);
		
		l_p = tmp_f.getFirstPos();
		r_p = tmp_f.getLastPos();
		
		for(i=1; i<n_rows_overlap; i++)
		{
			tmp_f =  fragments.get(ids_rows_overlap[i]);
			if(r_p < tmp_f.getLastPos())
				r_p  = tmp_f.getLastPos();			
		}
		
		tmp_P.profiles.start_p = l_p;
		int len = r_p - l_p + 1; 
		tmp_P.profiles.len = len;
		int g_2 = 2*maxg;
		tmp_P.profiles.n_groups = maxg;
		tmp_P.profiles.counts = new int[g_2][len];	
		for(i = 0; i < g_2; i++)
		{
			for(j = 0; j < len; j++)
			{
				tmp_P.profiles.counts[i][j] = 0;
			}
		}		
		
		for(i=0; i<n_rows_overlap; i++)
		{
			tmp_f =  fragments.get(ids_rows_overlap[i]);
		    int p0 = tmp_f.getFirstPos();
		    int p1 = tmp_f.getLastPos();
		    int rel_p;
		    p1++;
		    int gr = fun.get(i);
		    for(j = p0; j<p1; j++)
		    {
		    	rel_p = j - l_p; 
		    	char ch = tmp_f.getCall(j); 
		    	if(ch == Fragment.ALLELE1CHAR)
		    	{
		    		tmp_P.profiles.counts[2*gr][rel_p]++;	    		 
		    	}
		    	else if(ch == Fragment.ALLELE2CHAR)
		    	{
		    		 tmp_P.profiles.counts[2*gr+1][rel_p]++;	    		 
		    	}		    	
		    }						
		}		
		
		//compute mec score and the haplotypes			
		int mec_score = 0;
		int d_score = 0;
		if(n_geno == null)	{
			mec_score = comp_MEC(tmp_P.profiles.counts, maxg, len);			
			//compute D(P)
			char[][] haps = new char[maxg][len];		
			comp_Hap(tmp_P.profiles.counts, haps, maxg, len);
			for(i = 0; i < maxg; i++)
			{
				for(j = i+1; j< maxg; j++)
				{
					for(int k= 0; k < len; k++)
					{
						if(haps[i][k]!=Fragment.NODATACHAR && haps[j][k]!=Fragment.NODATACHAR)
						{
							if(haps[i][k]!= haps[j][k])
								d_score++;
				////2015.4.23 comment it out
				//			else d_score--;							
						}						
					}
				}		
			}			
		}			
		else
		{					
			int startcol = tmp_P.profiles.start_p;
			int endcol = startcol + tmp_P.profiles.len - 1;	
			int[] scores = comp_genotype_Constrained_Scores(tmp_P.profiles, startcol, endcol);
//			mec_score = comp_genotype_Constrained_MEC(tmp_P.profiles, startcol, endcol);
			mec_score = scores[0];
			d_score = scores[1];
		}			
		tmp_P.mec_score = mec_score;		
		tmp_P.d_score = d_score;		
		tmp_P.mix_score = (1-w_rate)*d_score - w_rate* mec_score ;			
		return tmp_P;
	}	
	
	//动态规划的初始化，从第1行开始，序号为0
	public int init(int firstrow, List<Fragment> frags,boolean useQS)
	{
			if(firstrow!=0)
			{
				System.err.println("Error: the firstrow is "+firstrow+"(Must begin from row 0)");
				 return -1;
			}
			else{				
				int n_reads_s = 10; 
				if(k_ploid>2) 
					n_reads_s = 7; 
				if(n_reads_s>frags.size())
					n_reads_s=frags.size();
				for(int i=1; i<n_reads_s;i++)
					if(frags.get(i).getFirstPos()>frags.get(firstrow).getFirstPos())
					{
						n_reads_s = i;
						break;
					}						
				
				List<List<Stirling2>> arrayStir = new ArrayList<List<Stirling2>>();
				for(int i=0; i<n_reads_s;i++)
				{
					List<Stirling2> list = new ArrayList<Stirling2>() ;
					for(int j = 0; j < k_ploid; j++)
					{
						if(j > i)
							break;
						Stirling2 node = new Stirling2();
						node.init(i+1, j+1);
						if(j>0 && i>j) //non-leaf node				
						{
							node.setleftchild(arrayStir.get(i-1).get(j-1));
							node.setrightchildren(arrayStir.get(i-1).get(j));					
							
						}
						node.re_set_cur_child();				
						list.add(node);				
					}			
					arrayStir.add(list);
				}
				
				PolypartitionFunctions funs_all = new PolypartitionFunctions(n_reads_s, k_ploid, arrayStir.get(n_reads_s-1));
				
				current_row = n_reads_s-1;
				n_rows_overlap = n_reads_s;
				ids_rows_overlap = new int[n_reads_s];
				for(int i=0; i<n_reads_s;i++)
					ids_rows_overlap[i]=i;
				fragments=frags;
				useQualityScores=useQS;		
				ps=new HBOPMinHeap();	
				ps.init(MAX_SIZE_HEAP, false);			
				List<Integer> fun = new ArrayList<Integer>();
				while(true)
				{
					fun.clear();
					boolean ret = funs_all.nextpartitionFunction(fun); //ret is true when at the end
					HPBOP2Partition tmp_partition = calScore(fun);
					ps.add(tmp_partition);					
					if(ret)
				       break;					
				}				
				return n_reads_s;		
			}
	}
	
	//动态规划的递推实现 
	public boolean extend_from_prev(int next_row, HBOP2PartitionAssembly pAssembly1)
	{
			if(next_row!=pAssembly1.current_row+1)
			{
				System.err.println("The nex row " +next_row+" is not the one expected (the prevois row is "+ pAssembly1.current_row+")");
				return false;
			}
			
			current_row = next_row;
			fragments = pAssembly1.fragments;
			useQualityScores = pAssembly1.useQualityScores;		
			
			int pos_begin=fragments.get(current_row).getFirstPos();	
			
			//debug info
//			if(pos_begin > 19 && pos_begin < 23)
//			{
//				System.out.println("Current row: "+current_row);				
//			}
			
			int i; 
			int n1_prev_rows = pAssembly1.n_rows_overlap;	
			boolean extended = true;
			boolean [] del_rs=null; //存储pAssembly1.ids_rows_overlap中没有与当前行重叠的片段在pAssembly1.ids_rows_overlap中的序号
			
			int [] del_rs_id =null;
			int  n_del_rs = 0;
			int max_del_row_id = 0;
					
			if(fragments.get(current_row).getFirstPos() > fragments.get(pAssembly1.current_row).getFirstPos())
			{
				del_rs=new boolean[n1_prev_rows];
				del_rs_id=new int[n1_prev_rows];
				n_del_rs = 0;
				for(i=0;i<n1_prev_rows;i++){
					if(fragments.get(pAssembly1.ids_rows_overlap[i]).getLastPos()<pos_begin)
					{
						del_rs_id[n_del_rs++]=i;
						del_rs[i]=true;	
						if(max_del_row_id<i)
							max_del_row_id = i;
					}
					else
					{					
						del_rs[i]=false;	
					}				
				}			
				if(n_del_rs>0)
					extended = false;		
			}		
			if(extended)
			{
				//the first row in the ids_rows_previ is always partitioned into group 0			
				extendedPartitions_rev(pAssembly1);					
			}
			else
			{			    		
				//step 1. 投影
				projectPartitions(pAssembly1, del_rs_id, del_rs, n_del_rs);
				
				//step 2. if profiles is longer than 100, split it
				split_long_profiles(pAssembly1);
				
				//step 3. 扩展	
				extendedPartitions_rev(pAssembly1);				
			}
			return true;
	}
		
	private int projectPartitions(HBOP2PartitionAssembly pAssembly1, int[] del_rs_id, boolean[] del_rs, int n_del_rs)   //, int prev_l)
	{
			int i,j;
			int n1_prev_rows=pAssembly1.n_rows_overlap;
			int n_rows_left=n1_prev_rows-n_del_rs;			
			HBOPMinHeap tmp_ps_prev;
			tmp_ps_prev = pAssembly1.ps;
			int n_ps_prev = pAssembly1.ps.len;			
			byte[] cur_P;	
			
			//n_rows_lef必须大于0 	
			
			int heap_size = n_ps_prev;//MAX_SIZE_HEAP;				
			
			HBOPMinHeap tmp_ps_current = new HBOPMinHeap();
			
			if(n_rows_left<1)
			{
				System.out.println("Error! ");
			}
			
			tmp_ps_current.init(heap_size, true);	
			HPBOP2Partition tmp_elem;			
			
			//the first row should be group 0, and first row of group i should be come before than group i+1
			//sort groups such that equivalent partitions can be easily detected
			int []  map = new int[k_ploid]; 				
			for(i=0;i<n_ps_prev;i++)
			{
				    tmp_elem=tmp_ps_prev.elements[i];
				    
				    cur_P = new byte[n_rows_left-1];
				  
				  // bug found at April 29, 2016
				  // int tmp_n_g = tmp_elem.profiles.n_groups; 
				    int tmp_n_g = tmp_elem.profiles.counts.length/2; 
				    
					for(j= 0; j<k_ploid; j++)
						map[j] = -1;		
				    
				    int new_n_g = covert2(tmp_elem.polypartition, n1_prev_rows, del_rs, cur_P, map);				   
				   
				    int [][] tmp_count = new int[2*tmp_n_g][];
				   
				    for(j = 0; j<tmp_n_g; j++)
				    {				    	
				    	tmp_count[2*map[j]] = tmp_elem.profiles.counts[2*j];
				    	tmp_count[2*map[j]+1] = tmp_elem.profiles.counts[2*j+1];
				    }
				    
				    tmp_elem.profiles.counts = tmp_count;
			//	    if(tmp_elem.profiles.prev!=null)
			//	        tmp_elem.profiles.prev.map = map;
				    if( tmp_elem.profiles.prev!=null )
				    {
				    	if(tmp_elem.profiles.map==null)
				    	{
				    		tmp_elem.profiles.map = map;
				    		map = new int[k_ploid]; 				    		
				    	}
				    	else{
				    		for(j = 0; j<k_ploid; j++)
				    			tmp_elem.profiles.map[j]=map[tmp_elem.profiles.map[j]];				    		
				    	}				    	
				    }				
				  
				    tmp_elem.profiles.n_groups = new_n_g;
				   // bug ? 
				   // tmp_elem.profiles.n_groups = tmp_n_g;
				    tmp_elem.polypartition = cur_P;					    
				    tmp_ps_current.add(tmp_elem);				    
			}			
			int[] tmp_ids_rows = new int[n_rows_left]; 			
			int tmp_idx=0;
			for(j=0;j<n1_prev_rows;j++)
			{
				   if(del_rs[j])
							continue;
					tmp_ids_rows[tmp_idx++]=pAssembly1.ids_rows_overlap[j];					
			}			
			pAssembly1.ids_rows_overlap = tmp_ids_rows;
			pAssembly1.n_rows_overlap = n_rows_left;
			pAssembly1.ps = tmp_ps_current;
			return n_rows_left;
	}

	/*
	private int genotype_Constrained_MEC(int[][] As, int[][] Bs, int startcol, int endcol)
	{
		int l = endcol - startcol;
		int[] deltas = new int[k_ploid];
		int mec = 0;
		for(int i= 0; i < l; i++ )
		{
			int n_1 = n_geno[startcol+i-startpos];
			if(n_1<1 || n_1>=k_ploid)
			{
				System.out.println("There should be at least one 0 and one 1 at the locus "+ (startcol+i));
				System.exit(-1);
			}
			for(int j = 0; j < k_ploid; j++ )
			{
				deltas[j] = Bs[i][j] - As[i][j];
				mec += Bs[i][j];     //assume all 0
			}
			for(int k = 0; k < n_1; k++)
			{				
				int max_idx = 0;
				for(int j = 1; j < k_ploid; j++ )
				{
					if(deltas[j] > deltas[max_idx])
						max_idx = j;
				}
				mec -= deltas[max_idx];				
				deltas[max_idx] = -100000; 
			}			
		}		
		return mec;		
	}
	*/
	private void split_long_profiles(HBOP2PartitionAssembly pAssembly1){
		//if profiles of pAssembly1 is longer than 100, split it 
		int pos_begin = fragments.get(current_row).getFirstPos(); 
		if(pAssembly1.ps.len<1)
			return;
		int prev_begin = pAssembly1.ps.elements[0].profiles.start_p;
		if(pos_begin-prev_begin<split_len)
			return;
		HPBOP2Partition elm;
		int ps_len = pAssembly1.ps.len;
		for(int k=0; k<ps_len; k++){
			elm = pAssembly1.ps.elements[k];
			int old_g_2 = elm.profiles.counts.length;
			int len1 = pos_begin-prev_begin;
			int[][] cunts = new int[old_g_2][len1];
			for(int i = 0; i< old_g_2; i++){
				for(int j= 0; j<len1; j++){
					cunts[i][j] = elm.profiles.counts[i][j];
				}
			}
			S_Profiles prev_prof = new S_Profiles();
			prev_prof.counts = cunts;
			prev_prof.len = len1;
			prev_prof.n_groups = old_g_2/2;
			prev_prof.start_p = elm.profiles.start_p;
			prev_prof.prev = elm.profiles.prev;	
			prev_prof.map = elm.profiles.map;
			elm.profiles.map = null;
			
			int len2 = elm.profiles.len - len1;
			int n_g_2 = elm.profiles.n_groups*2;
			int[][] cunts_2 = new int[n_g_2][len2];
			for(int i = 0; i< n_g_2; i++){
				for(int j= 0; j<len2; j++){
					cunts_2[i][j] = elm.profiles.counts[i][j+len1];
				}
			}
			elm.profiles.counts = null;
			elm.profiles.counts = cunts_2;
			elm.profiles.start_p = pos_begin;
			elm.profiles.len = len2;
			elm.profiles.prev = prev_prof;		
		}	
	}
	
	private void extendedPartitions(HBOP2PartitionAssembly pAssembly1)
	{		
		int n1_prev_rows = pAssembly1.n_rows_overlap;	
		n_rows_overlap = n1_prev_rows + 1;			
		ids_rows_overlap = new int[n_rows_overlap];
		int i, j, k;	
		for(i=0;i<n1_prev_rows;i++)
			ids_rows_overlap[i]=pAssembly1.ids_rows_overlap[i];
		ids_rows_overlap[n1_prev_rows]=current_row;
		
		HBOPMinHeap prev_ps = pAssembly1.ps;
		int p1_len = prev_ps.len;		
		ps=new HBOPMinHeap();	
		
		int tmp_heap_size = k_ploid * p1_len;
		if(tmp_heap_size > MAX_SIZE_HEAP)
			tmp_heap_size =  MAX_SIZE_HEAP;	
		
		ps.init(tmp_heap_size,false);		
		HPBOP2Partition tmp_p, cur_p1;
		
		n1_prev_rows--; 
		
		for(i=0;i<p1_len;i++)
		{ 
			tmp_p = prev_ps.elements[i];	
			
			int cur_g = tmp_p.profiles.n_groups + 1;
			
//			int cur_g = 0;			
//			for(j = 0; j < n1_prev_rows; j++)
//			{
//				if(cur_g < tmp_p.polypartition[j])
//					cur_g = tmp_p.polypartition[j];				
//			}
//			cur_g = cur_g+2;
			if(cur_g > k_ploid)
				cur_g = k_ploid;			
			for(k = 0; k < cur_g; k++)
			{
				byte[] PKs = new byte[n_rows_overlap-1];	
				//for(j=0; j<tmp_p.polypartition.length;j++)
				for(j = 0; j < n1_prev_rows; j++)	
				{
					PKs[j] = tmp_p.polypartition[j];				
				}				
				PKs[j] = (byte) k;				
				cur_p1 = new HPBOP2Partition();
				cur_p1.polypartition = PKs;
				cur_p1.profiles = new S_Profiles();	
				if(tmp_p.profiles.map != null){
					cur_p1.profiles.map = new int[k_ploid];
					for(j = 0; j < k_ploid; j++)	
					{
						cur_p1.profiles.map[j] = tmp_p.profiles.map[j];				
					}					
				}				
				cur_p1.profiles.prev = tmp_p.profiles.prev;
				if(n_geno==null)
					calScore_n(cur_p1, tmp_p, k);
				else
					calScore_n_G(cur_p1, tmp_p, k);
				ps.add(cur_p1);								
			}				
		}
	}
	
	
	//the firstpos of f2 must not smaller than f1
	//if the overlapping parts are different, return -1; otherwise return the len of identical overlapping
	private int identical_len(int f1, int f2){
		Fragment f_r1 = fragments.get(f1);
		Fragment f_r2 = fragments.get(f2);
		int begin_ind = f_r2.getFirstPos();
		int end_ind = f_r1.getLastPos();
		if(end_ind > f_r2.getLastPos())
			end_ind = f_r2.getLastPos();
		end_ind++;
		int len_id = 0;
		for(int i = begin_ind; i < end_ind; i++){
			char ch1 = f_r1.getCall(i);
			char ch2 = f_r2.getCall(i);			
			if(ch1!=Fragment.NODATACHAR && ch2!=Fragment.NODATACHAR)
			{
				if(ch1 == ch2)
					len_id++;
				else				
					return 0;				
			}		
		}
		return len_id;		
	}
	
	private void extendedPartitions_rev(HBOP2PartitionAssembly pAssembly1)
	{		
		int n1_prev_rows = pAssembly1.n_rows_overlap;	
		n_rows_overlap = n1_prev_rows + 1;			
		ids_rows_overlap = new int[n_rows_overlap];
		int i, j, k;
		
		for(i=0;i<n1_prev_rows;i++)
			ids_rows_overlap[i]=pAssembly1.ids_rows_overlap[i];
		ids_rows_overlap[n1_prev_rows]=current_row;
		
		//check if there is a row before is identical to the current row 
		int len = 0;
		int pre_idx = -1;
		//for(i=n1_prev_rows-1;i>=0;i--)
		i = n1_prev_rows-1;
		if(i>=0)
		{
			int id_l = identical_len(ids_rows_overlap[i], current_row);
			if(len < id_l)
			{
				len = id_l;
				pre_idx = i;
			}
		}
		boolean extend_only_one = false;
		if(len > 0){
			if(fragments.get(current_row).getLastPos()==fragments.get(ids_rows_overlap[pre_idx]).getLastPos())
				if(fragments.get(ids_rows_overlap[pre_idx]).getFirstPos() == fragments.get(current_row).getFirstPos())
					extend_only_one = true;	
		}		
		
		HBOPMinHeap prev_ps = pAssembly1.ps;
		int p1_len = prev_ps.len;	
		
		ps=new HBOPMinHeap();	
		
		int tmp_heap_size;
		if(extend_only_one)
			tmp_heap_size = p1_len;
		else{
			tmp_heap_size = k_ploid * p1_len;
			if(tmp_heap_size > MAX_SIZE_HEAP)
				tmp_heap_size =  MAX_SIZE_HEAP;				
		}			
		
		ps.init(tmp_heap_size,false);		
		HPBOP2Partition tmp_p, cur_p1;
		
		n1_prev_rows--; 	
		
		if(extend_only_one)
		{
			pre_idx--;
			for(i=0;i<p1_len;i++)
			{ 
				tmp_p = prev_ps.elements[i];
				byte[] PKs = new byte[n_rows_overlap-1];				
				for(j = 0; j < n1_prev_rows; j++)	
				{
						PKs[j] = tmp_p.polypartition[j];				
				}				
				if(pre_idx < 0){
					PKs[j] = (byte) 0;
				}
				else
					PKs[j] = PKs[pre_idx];
				cur_p1 = new HPBOP2Partition();
				cur_p1.polypartition = PKs;
				cur_p1.profiles = new S_Profiles();	
				if(tmp_p.profiles.map != null){
					cur_p1.profiles.map = new int[k_ploid];
					for(j = 0; j < k_ploid; j++)	
					{
						cur_p1.profiles.map[j] = tmp_p.profiles.map[j];				
					}					
				}				
				cur_p1.profiles.prev = tmp_p.profiles.prev;
				if(n_geno==null)
						calScore_n(cur_p1, tmp_p, PKs[n1_prev_rows]);
				else
						calScore_n_G(cur_p1, tmp_p, PKs[n1_prev_rows]);
				ps.add(cur_p1);					
			}
		}
		else{			
			for(i=0;i<p1_len;i++)
			{ 
				tmp_p = prev_ps.elements[i];	
				
				int cur_g = tmp_p.profiles.n_groups + 1;			

				if(cur_g > k_ploid)
					cur_g = k_ploid;			
				for(k = 0; k < cur_g; k++)
				{
					byte[] PKs = new byte[n_rows_overlap-1];	
					//for(j=0; j<tmp_p.polypartition.length;j++)
					for(j = 0; j < n1_prev_rows; j++)	
					{
						PKs[j] = tmp_p.polypartition[j];				
					}				
					PKs[j] = (byte) k;				
					cur_p1 = new HPBOP2Partition();
					cur_p1.polypartition = PKs;
					cur_p1.profiles = new S_Profiles();	
					if(tmp_p.profiles.map != null){
						cur_p1.profiles.map = new int[k_ploid];
						for(j = 0; j < k_ploid; j++)	
						{
							cur_p1.profiles.map[j] = tmp_p.profiles.map[j];				
						}					
					}				
					cur_p1.profiles.prev = tmp_p.profiles.prev;
					if(n_geno==null)
						calScore_n(cur_p1, tmp_p, k);
					else
						calScore_n_G(cur_p1, tmp_p, k);
					ps.add(cur_p1);								
				}				
			}	
			
		}	
		
	
	}
	
	
	
	void calScore_n(HPBOP2Partition cupr_p1, HPBOP2Partition tmp_p, int k)
	{
		int deltaMec=0, deltaDscore=0;
		Fragment f = fragments.get(current_row);
		
		S_Profiles profiles = tmp_p.profiles;
			
		int groups = profiles.n_groups;
		if(k > groups)
		{
			System.out.println("ERROR!");
			System.exit(-1);
		}
		
		int s0 = f.getFirstPos() - profiles.start_p;
		if(s0 < 0)
		{
			System.out.println("ERROR!");
			System.exit(-1);
		}		
		
		int o_l = profiles.len;	
		int n_l = f.getLastPos() - profiles.start_p + 1; 
		if(n_l < o_l)
			n_l = o_l;
		
				
		int new_groups;
		if(k==groups)
			new_groups = groups+1;
		else
			new_groups = groups;		
		
		//There is a bug that found at April 29, 2016
		//int[][] cs = new int[2*new_groups][n_l];
		int old_cnt_size = profiles.counts.length;
		int new_cnt_size = 2*new_groups;
		if(old_cnt_size>new_cnt_size)
			new_cnt_size = old_cnt_size;		
		int[][] cs = new int[new_cnt_size][n_l];
		
		int i, j;
		//for(j = 0; j < 2*groups; j++)
		for(j = 0; j < old_cnt_size; j++)
		{			
			for(i = 0; i < o_l; i++ )
			{
				cs[j][i] = profiles.counts[j][i];					 
			}
			while(i<n_l)
			{
				cs[j][i] = 0;
				i++;
			}
		}
		//while(j<2*new_groups)
		while(j<new_cnt_size)
		{
			i = 0;
			while(i < n_l)
			{
				cs[j][i] = 0;
				i++;
			}
			j++;
		}
		
		cupr_p1.profiles.start_p = profiles.start_p;
		cupr_p1.profiles.len = n_l;		
		cupr_p1.profiles.n_groups = new_groups;
				
		int t0 = f.getLastPos() - profiles.start_p + 1;		
		for(i = s0; i < t0; i++ )		
		{		
			char ch = f.getCallbyIndex(i-s0);
			if(ch == Fragment.NODATACHAR)
				continue;
			if(ch == Fragment.ALLELE1CHAR)
			{
				if(cs[2*k][i]==0 && cs[2*k+1][i]==0) //the consensus hap at this pos was -
				{
					cs[2*k][i]++;   //now is 0
					for(j=0; j<groups; j++)
					{
						if(j==k)
							continue;
						if(cs[2*j][i]==0 && cs[2*j+1][i]==0) 						
							continue;  //'-'--'-'-->0--'-', don't change the distance
						if(cs[2*j][i]<cs[2*j+1][i])						
						{
							deltaDscore++;			 //'-'--1-->0--1 distance +1				
						}		
					}					
				}
				else  //the consensus hap at this pos was 0 or 1
				{
					cs[2*k][i]++;
					if(cs[2*k][i] < cs[2*k+1][i])  //the consensus hap at this pos doesn't change
					{
						deltaMec++;
					}
					else if(cs[2*k][i]==cs[2*k+1][i]) 
				    //the consensus hap at this pos changes from 1 to 0
					{					
						deltaMec++;
						for(j=0; j<groups; j++)
						{
							if(j==k)
								continue;
							if(cs[2*j][i]==0 && cs[2*j+1][i]==0) //0 ---'-', don't change the distance
								continue;
							if(cs[2*j][i]<cs[2*j+1][i])  //1--1 ---> 0--1, distance +1							
								deltaDscore++;							                    
							else				//1--0 ---> 0--0, distance -1				
								deltaDscore--;						
								
						}					
					}
								
				}		
					
			}
			else if(ch == Fragment.ALLELE2CHAR)
			{
				if(cs[2*k][i]==0 && cs[2*k+1][i]==0) //the consensus hap at this pos:  '-'-->'1'
				{					
					for(j=0; j<groups; j++)
					{
						if(j==k)
							continue;
						if(cs[2*j][i]==0 && cs[2*j+1][i]==0)
							continue;
						if(cs[2*j][i]>=cs[2*j+1][i])						
							deltaDscore++;           //'-'--0-->1--0 distance +1	
					}					
				}
				else if(cs[2*k+1][i]==cs[2*k][i]) //the consensus hap at this pos:  '0'-->'1'
				{					
					for(j=0; j<groups; j++)
					{
						if(j==k)
							continue;
						if(cs[2*j][i]==0 && cs[2*j+1][i]==0)
							continue;
						if(cs[2*j][i]<cs[2*j+1][i]) //'0'--'1'-->'1'--'1' distance -1			       
							deltaDscore--;						
						else		   //'0'--'0'-->'1'--'0' distance +1					
							deltaDscore++;						
					}					
				}
				else if(cs[2*k+1][i]<cs[2*k][i]){//the consensus hap at this pos doesn't change
					deltaMec++;				
				}
				cs[2*k+1][i]++;	
			}			
		}
		
//		if(n_geno != null)
//		{
//			int startcol = f.getFirstPos();
//			int endcol = f.getLastPos();				
//			int mec_score_1 = comp_genotype_Constrained_MEC(tmp_p.profiles, startcol, endcol);
//			int mec_score_2 = comp_genotype_Constrained_MEC(cupr_p1.profiles, startcol, endcol);
//			deltaMec = mec_score_2 - mec_score_1; 
//		}	
				
		cupr_p1.mec_score = tmp_p.mec_score + deltaMec;		
		cupr_p1.d_score = tmp_p.d_score + deltaDscore;
		cupr_p1.mix_score = (1-w_rate)*cupr_p1.d_score - w_rate*cupr_p1.mec_score;	
		cupr_p1.profiles.counts = cs;
	}
	
	void calScore_n_G(HPBOP2Partition cupr_p1, HPBOP2Partition tmp_p, int k)
	{
		int deltaMec=0, deltaDscore=0;
		Fragment f = fragments.get(current_row);
		
		S_Profiles profiles = tmp_p.profiles;
			
		int groups = profiles.n_groups;
		if(k > groups)
		{
			System.out.println("ERROR!");
			System.exit(-1);
		}
		
		int s0 = f.getFirstPos() - profiles.start_p;
		if(s0 < 0)
		{
			System.out.println("ERROR!");
			System.exit(-1);
		}
		
		int o_l = profiles.len;	
		int n_l = f.getLastPos() - profiles.start_p + 1; 
		if(n_l < o_l)
			n_l = o_l;	
				
		int new_groups;
		if(k==groups)
			new_groups = groups+1;
		else
			new_groups = groups;
		
		//There is a bug that found at April 29, 2016
		//int[][] cs = new int[2*new_groups][n_l];
		int old_cnt_size = profiles.counts.length;
		int new_cnt_size = 2*new_groups;
		if(old_cnt_size>new_cnt_size)
			new_cnt_size = old_cnt_size;		
		int[][] cs = new int[new_cnt_size][n_l];	
		
		int i, j;
		
		//for(j = 0; j < 2*groups; j++)
		for(j = 0; j < old_cnt_size; j++)	
		{			
			for(i = 0; i < o_l; i++ )
			{
				cs[j][i] = profiles.counts[j][i];					 
			}
			while(i<n_l)
			{
				cs[j][i] = 0;
				i++;
			}
		}
		//while(j<2*new_groups)
		while(j<new_cnt_size)		
		{
			i = 0;
			while(i < n_l)
			{
				cs[j][i] = 0;
				i++;
			}
			j++;
		}	
		
		cupr_p1.profiles.start_p = profiles.start_p;
		cupr_p1.profiles.len = n_l;
		cupr_p1.profiles.n_groups = new_groups;		

		int t0 = f.getLastPos() - profiles.start_p + 1;		
		for(i = s0; i < t0; i++ )		
		{		
			char ch = f.getCallbyIndex(i-s0);					
			if(ch == Fragment.ALLELE1CHAR)			
				cs[2*k][i]++;		
			else if(ch == Fragment.ALLELE2CHAR)							
				cs[2*k+1][i]++;
			else if(ch == Fragment.NODATACHAR)
				continue;
			else{
				System.out.print("There is an illegal char in Fragment: "+ f.getId() + " "+ch );				
				System.exit(-1);				
			}						
		}
		
		cupr_p1.profiles.counts = cs;
		
//		if(n_geno != null)
//		{
			int startcol = f.getFirstPos();
			int endcol = f.getLastPos();
			int[] scores_1 = comp_genotype_Constrained_Scores(tmp_p.profiles, startcol, endcol);
			int[] scores_2 = comp_genotype_Constrained_Scores(cupr_p1.profiles, startcol, endcol);
			deltaMec = scores_2[0] - scores_1[0]; 
			deltaDscore = scores_2[1] - scores_1[1]; 			
//		}	
		
		cupr_p1.mec_score = tmp_p.mec_score + deltaMec;		
		cupr_p1.d_score = tmp_p.d_score + deltaDscore;
//		if(cupr_p1.d_score == 0){
//			System.out.print("o");
//		}
		
		cupr_p1.mix_score = (1-w_rate)*cupr_p1.d_score - w_rate*cupr_p1.mec_score;		
	}
	
	
	
	public HPBOP2Partition find_max_elem()
	{
		int i;	
		
		List<HPBOP2Partition> max_s = new ArrayList<HPBOP2Partition>();
		HPBOP2Partition ps_max = ps.elements[0];
		max_s.add(ps_max);
		int len_ps=ps.len;
		for(i=1;i<len_ps;i++)
		{			
			if(ps_max.mix_score < ps.elements[i].mix_score || (ps_max.mix_score == ps.elements[i].mix_score && ps_max.mec_score > ps.elements[i].mec_score))
			{
				max_s.clear();
				ps_max=ps.elements[i];
				max_s.add(ps_max);
			}
			else if(ps_max.mix_score == ps.elements[i].mix_score && ps_max.mec_score == ps.elements[i].mec_score)
			{				
				max_s.add(ps.elements[i]);
			}			
		}
		
		max_score= ps_max.mix_score;
		
		Random r= new Random(System.currentTimeMillis());
		int rg = max_s.size();
		ps_max = max_s.get(r.nextInt(rg));
		
//		System.out.println(" " + ps_max.mec_score + " "+ ps_max.mix_score);
		
		return ps_max;
	}	
	
	
	private int covert2(byte[] p, int n_rows, boolean[] del_rs, byte[] par_ret, int [] map)
	{				
		boolean first_row;
		int next_new_g_id = 0;
		if(del_rs[0])
		{
			first_row = true;			
		}
		else{
			first_row = false;
			map[0] = 0;
			next_new_g_id = 1;
		}	
		int idx = 1, old_idx = 0, new_idx = 0; //第0行没有编码
		int tmp_g ;
		while(idx < n_rows)
		{				
			if(!del_rs[idx])
			{
				tmp_g =  p[old_idx];
				if(first_row)
				{
					map[tmp_g] = 0;
					next_new_g_id = 1;
					first_row = false;
				}
				else
				{
					if(map[tmp_g] == -1)
					{
						map[tmp_g] =  next_new_g_id;
						next_new_g_id++;
					}
					par_ret[new_idx] = (byte)map[tmp_g];
					new_idx++;			
				}					
			}
			idx++; 
			old_idx++;						
		}		
		int new_g = next_new_g_id;
		for(int i= 0; i<k_ploid; i++)
			if(map[i] == -1)
				map[i] = next_new_g_id++;
		return new_g;
	}	
	
	public double get_max_score()
	{
		return max_score;
	}
	
	
}