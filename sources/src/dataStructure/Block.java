package dataStructure;

import java.util.ArrayList;
import java.util.List;

public class Block {
	private List<Fragment> fragments;
	private int firstPos;
	private int lastPos;
//	private String haplotype;
	private int MEC=0;
	private int phased=0;
	private int calls=0;
//	private double hapScores [];
//	public int n_holes = 0;
//	public int n_undetermined = 0;
		
	int[] positions;
	char[][] phased_hs;	
	int n_ploid;
	String[] strHaps;
	
	public int coverage[], co_rows[]; //
		
	public double object_opt_score=0;	
	
	private void sum_calls(){
		for(Fragment f:fragments) {			
			calls+=f.getnCalls();
			}		
	}	
	
	public Block(List<Fragment> fragments) {
		super();
		this.fragments = fragments;
		this.firstPos = fragments.get(0).getFirstPos();
		this.lastPos = calculateLastPos(fragments); 
		//xie added 2013/7/30
		sum_calls();
	}
	private int calculateLastPos(List<Fragment> fragments) {
		int answer =0;
		for(Fragment f:fragments) {
			if(answer < f.getLastPos()) {
				answer = f.getLastPos();
			}
		}
		return answer;
	}
	
		
	public void setHaps(char [][] hs, int k_ploid )
	{
		phased_hs = hs;
		n_ploid = k_ploid;
		strHaps = new String[n_ploid];
		int len_h = lastPos - firstPos + 1;			
		for(int j = 0; j < n_ploid; j++)
		{
			StringBuilder hap = new StringBuilder("");
			for(int i=0;i<len_h;i++) {
				hap.append(phased_hs[i][j]);
			}
			strHaps[j]= hap.toString();					
		}		
		updateSt();
		
	}
	
	private void updateSt() 
	{			
			phased=0;		
			int len_h = lastPos - firstPos + 1;	
			
			for(int i=0;i<len_h;i++) {
				for(int j = 0; j < n_ploid; j++)
				{
					if(phased_hs[i][j]!=Fragment.NODATACHAR)
					{
						phased++;
					}
				}
			}	
					
			MEC = 0;
			for(Fragment f:fragments) {
				if(f.getLastPos()<=firstPos)
					continue;
				if(f.getFirstPos()>=lastPos)
					break; 
				int tmp_mec = 100000;
				for(int j = 0; j < n_ploid; j++)
				{					
					int d1 = f.getHammingDistance(strHaps[j],firstPos);
					if(tmp_mec > d1)
						tmp_mec = d1;
					if(tmp_mec == 0)
						break;
				}
//				if(tmp_mec > 0)
//				{
//					System.out.println(" " + f.getFirstPos()+"---"+f.getExtendedCalls());
//					for(int j = 0; j < n_ploid; j++)
//					{	
//						System.out.println(strHaps[j]);					
//					}					
//				}
				MEC += tmp_mec; 		
			}	
			
	}
	

	public char[][] getPolyH() {
		return phased_hs;
	}
	
	
	public int getMEC() {
		return MEC;
	}
	public int getPhased() {
		return phased;
	}
	
	public int getCalls() {
		return calls;
	}
	
	public List<Fragment> getFragments() {
		return fragments;
	}
	
	public int getFirstPos() {
		return firstPos;
	}
	
	public int getLastPos() {
		return lastPos;
	}
	public int length() {
		return lastPos-firstPos+1;
	}
	public List<Fragment> getFragments(int first,int last) {
		List<Fragment> answer = new ArrayList<Fragment>();
		for(Fragment f:fragments) {
			if(f.getFirstPos()<=last && f.getLastPos()>=first) {
				answer.add(f);
			}
		}
		return answer;
	}	
	
}

