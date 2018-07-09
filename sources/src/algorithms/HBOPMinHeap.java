package algorithms;

public class HBOPMinHeap {
	private int max_size;  //the maximum number of elements that the heap can store
	int len; //the number of elements in the heap
	HPBOP2Partition[] elements;
	private boolean use_map;  //if maximum_root is true, the maximum element is elements[0]
//	private HashMap map;
	public void init(int size, boolean need_map)
	{
		//max_size = size+1;
		max_size = size;		
		elements = new HPBOP2Partition[max_size];
		use_map = need_map;
		len = 0;
	//	if(use_map)
	//		map = new HashMap<Long, HBOPPartition>((int)((max_size+1)*1.3),0.8f);		
	}
	public int is_exist(byte[] p_code){   //-1 means there is no element in the heap whose partition equals to p_code, else the return integer is the index
		byte[] tmp_p;
		int i, j, n_bytes;
		for(i=0;i<len;i++)
		{
			tmp_p=elements[i].polypartition;
			n_bytes=tmp_p.length;
			if(n_bytes!=p_code.length)
			{
				System.out.print("Two partitions cannot be compared!");
				System.exit(-1);
			}
			for(j=0;j<n_bytes;j++)
			{
				if(tmp_p[j]!=p_code[j])
					break;
			}
			if(j==n_bytes)
				return i;		
		}
		return -1;	
	}
	
	//if use_map is true, check whether the key is in the heap	
	public boolean add(HPBOP2Partition elem){
		byte[] p_code=elem.polypartition;
		HPBOP2Partition tmp_elem;
		
		if(use_map)
		{
			int tmp_idx = is_exist(p_code);
			if(tmp_idx==-1) //p_code is unique
			{
	//			if(len<max_size)
				{
					elements[len++]=elem;					
	//				if(len==max_size)
	//					build_heap();					
				}
				/*
				else
				{
					//if(elem.mec_score<elements[0].mec_score||(elem.mec_score==elements[0].mec_score && elem.score>elements[0].score))
					if(elem.mix_score > elements[0].mix_score)
					{
						fixHeap(0, len-1, elem);					
					}
					else
					{
						return false;						
					}
				}
				*/	
				
			}
			else{              //p_code is duplicate
				tmp_elem = elements[tmp_idx];
				if(tmp_elem.mix_score < elem.mix_score || (tmp_elem.mix_score == elem.mix_score && tmp_elem.mec_score > elem.mec_score))
				{			
					elements[tmp_idx] = elem;
			//		if(len==max_size)
			//			build_heap();					
				}
				else
					return false;
			}		
		}
		else
		{
			if(len<max_size){
				elements[len++]=elem;				
				if(len==max_size)
					build_heap();					
			}
			else{
				//if(elem.mec_score<elements[0].mec_score||((elements[0].mec_score==elem.mec_score)&&(elements[0].score<elem.score)))
				if(elem.mix_score > elements[0].mix_score || (elem.mix_score == elements[0].mix_score && elem.mec_score < elements[0].mec_score))
				{					
					fixHeap(0, len-1, elem);
				}
				else
				{
					return false;						
				}
			}			
		}
		return true;		
	}
	
	
	
	// Establish the heap property.
	 public void build_heap()
     {
   	  int i, end = len-1;
   	  for (i = (end-1) / 2; i >= 0; i-- )
             fixHeap(i, end, elements[i]); 
     }	
	/**
	       Assuming that the partial order tree
	       property holds for all descendants of
	       the element at the root, make the
	       property hold for the root also.
	 
	       @param root the index of the root of the current subtree
	       @param end  the highest index of the heap
	 */
    private void fixHeap( int root, int end, HPBOP2Partition key )
    {
       int child = 2*root + 1; // left child    
       // Find the smaller child.
      // if ( child < end && ((elements[child].mec_score < elements[child + 1].mec_score) || (elements[child].mec_score == elements[child + 1].mec_score && elements[child].score > elements[child + 1].score)))
       if ( child < end && elements[child].mix_score > elements[child + 1].mix_score) 
             child++;  // right child is smaller
       // If the smaller child is smaller than the
       // element at the root, move the smaller child
       // to the root and filter the former root 
       // element down into the "smaller" subtree.
      // if ( child <= end && ( (key.mec_score < elements[child].mec_score) || (key.mec_score == elements[child].mec_score && key.score > elements[child].score)))
       if ( child <= end && key.mix_score > elements[child].mix_score) 
       {
            elements[root] = elements[child];
            fixHeap( child, end, key );
       }
       else
       {
    	   elements[root] = key;
//    	   if( key.mec_score == 0)//&&  key.profiles.len > 3
//    	   {
//    		   System.out.println("Pause! "+ root + " "+ key.profiles.start_p +"---" + key.profiles.len +":"+ key.mec_score);
//    		   
//    	   }
       }
    }   
    
    /**
         Swaps two entries of the array.
         @param i the first position to swap
         @param j the second position to swap
       */
      private void swap(int i, int j)
      {
    	  HPBOP2Partition temp = elements[i];
    	  elements[i] = elements[j];
    	  elements[j] = temp;
      }        
      public void sort()
      {
    	  int i, end=len-1;
    	  
         // Establish the heap property.
    	  build_heap();  
    	  
         // Now place the smallest element last,
         // 2nd smallest 2nd last, etc.
         for (i=end; i>0; i--){
        	 
         // elements[0] is the smallest element.
              swap(0, i);    
              
         // Heap shrinks by 1 element.
              fixHeap(0, i-1, elements[0]);
         }
      }
}

class HPBOP2Partition{

	S_Profiles profiles;  
	
	//public byte[] prev_partition; //储存在当前片段之前，且与当前片段没有重叠的片段的当前的分组

	//public double score; //储存在该分区对应的片段割的得分，分区由该分区在数组中的序号决定
	
	public double mix_score;  //(1-w_rate)*d_score-w_rate*mec_score
	
	public int mec_score;
	
	public int d_score;
	
	public byte[] polypartition; //the first row must be in group 0
}

