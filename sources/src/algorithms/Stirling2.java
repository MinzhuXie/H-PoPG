package algorithms;

import java.util.ArrayList;
import java.util.List;

public class Stirling2 {
	int n_reads, m_boxes;
	Stirling2 leftchild, rightchildren;
	int next_child; 
	void init(int n1, int m1)
	{
		n_reads = n1; m_boxes = m1;			
	}
	void re_set_cur_child()
	{
		if(m_boxes == 1 || n_reads == m_boxes)
		   next_child = -1; //-1 if leaf, else 0
		else
			next_child = 0;
		
	}
    void setleftchild(Stirling2 child)
    {
    	leftchild = child;
    }
    void setrightchildren(Stirling2 child)
    {
    	rightchildren = child;
    }
    void setchild_ind(int idx)
    {
    	next_child = idx;
    }
    boolean next_partitionFunction(List<Integer> fun)
    //return true if it is the last function, false else
    {    	
    	boolean ret = true; 
    	if(next_child < 0) //leaf node next_child == -1
    	{    		
    		if(m_boxes > 1) //n_reads == m_boxs >1
    		{
    			for(int i = 0; i < m_boxes; i++)
    				fun.add(i);  			
    			
    		}
    		else //m_boxs == 1
    		{
    			for(int i = 0; i < n_reads; i++)
    				fun.add(0);     			
    		}    		
    	}    	
    	else  //non-leaf node   next_child = 0, ..., m_boxes
    	{  
    		if(next_child == 0)  //left branch    	
	    	{
	    		ret = leftchild.next_partitionFunction(fun);
	    		fun.add(m_boxes-1);    		
	    	}
    		else  //right branch
    		{
    			ret = rightchildren.next_partitionFunction(fun);
    			fun.add(next_child-1);    			
    		}
    		if(ret)
    		{
    			next_child ++;
    			if(next_child > m_boxes)
    			{
    				next_child = 0;    				
    			}
    			else
    				ret = false;   			
    		}
    	}
    	return ret;    	
    }
    public static void main(String[] args){
    	
    	List<List<Stirling2>> arrayStir = new ArrayList<List<Stirling2>>();
		for(int i=0; i<6;i++)
		{
			List<Stirling2> list = new ArrayList<Stirling2>() ;
			for(int j = 0; j < 5; j++)
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
		List<Integer> fun = new ArrayList<Integer>();
		PolypartitionFunctions funs_all = new PolypartitionFunctions(6, 5, arrayStir.get(5));
		while(true)
		{
			fun.clear();
			boolean ret = funs_all.nextpartitionFunction(fun); //ret is true when at the end
			
				
			for(int i= 0; i<fun.size(); i++)
			{
				System.out.print(" "+fun.get(i));
			}
			System.out.println();
			
			if(ret)
			       break;
			
		}				
    	
    	
    	
    	
    }

}
