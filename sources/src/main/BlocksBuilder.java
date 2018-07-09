package main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import dataStructure.Block;
import dataStructure.Fragment;

public class BlocksBuilder {
	private int [] assignments;
	private int lastBlock=0;
	private List<Join>joins=new ArrayList<Join>(); 
	public List<Block> buildBlocks(Block bigBlock) {
		fillAssignments(bigBlock);
		joinBlocks();
		return distributeFragments(bigBlock.getFragments());
	}
	private void fillAssignments(Block bigBlock) {
		assignments = new int [bigBlock.getLastPos()-bigBlock.getFirstPos()+1];
		Arrays.fill(assignments, 0);
		List<Fragment> fragments = bigBlock.getFragments();
		int firstPos = bigBlock.getFirstPos();
		for(Fragment f:fragments) {
			int fragBlock = 0;
			String calls = f.getExtendedCalls();
			int relPos = f.getFirstPos()-firstPos;
			for(int k=0;k<calls.length();k++) {
				if(calls.charAt(k)!=Fragment.NODATACHAR) {
					int block = assignments[relPos+k];
					if(block > 0) {
						if(fragBlock == 0) {
							fragBlock = block;
						} else if(block > fragBlock) {
							joins.add(new Join(block,fragBlock));
						} else if(fragBlock > block) {
							joins.add(new Join(fragBlock,block));
							fragBlock = block;
						}
					}
				}
			}
			if(fragBlock==0) {
				lastBlock++;
				fragBlock = lastBlock;
			}
			for(int k=0;k<calls.length();k++) {
				if(calls.charAt(k)!=Fragment.NODATACHAR) {
					assignments[relPos+k] = fragBlock;
				}
			}
		}
		
	}

	private void joinBlocks() {
		if(joins.size()<1)
			return;
		Collections.sort(joins,new JoinComparator());		
		int lastGroupId = joins.get(0).group1;
		int size_m = lastGroupId+1;
		int[] parents = new int[size_m];
		for(int i=0;i<size_m;i++){
			parents[i] = i;
		}
		for(Join join:joins ) {
			int root1 =join.group1;
			int root2 =join.group2;
			//join two sets, the root should be the smallest number in the joined set
			while(root1!=parents[root1])
				root1 = parents[root1];
			while(root2!=parents[root2])
				root2 = parents[root2];
			if(root1<root2)
				parents[root2]=root1;
			else if(root2<root1)
				parents[root1]=root2;
		}
		for(int j=0;j<assignments.length;j++) {
			int r = assignments[j];
			if(r>lastGroupId)
				continue;
			while (parents[r]!=r)
				r = parents[r];			
			assignments[j] = r;
		}			
	}
		
	private List<Block> distributeFragments(List<Fragment> fragments) {
		List<List<Fragment>> distFragments = new ArrayList<List<Fragment>>(lastBlock);
		for(int i=0;i<lastBlock;i++) {
			distFragments.add(new ArrayList<Fragment>());
		}
		int firstPos = fragments.get(0).getFirstPos();
		for(Fragment f:fragments) {
			int pos = assignments[f.getFirstPos()-firstPos]-1;
			distFragments.get(pos).add(f);
		}
		List<Block> blocks = new ArrayList<Block>();
		for(List<Fragment> l:distFragments) {
			if(l.size()>0) {
				blocks.add(new Block(l));
			}
		}
		return blocks;
	}
}
class Join {
	int group1;
	int group2;
	public Join(int group1, int group2) {
		super();
		this.group1 = group1;
		this.group2 = group2;
	}
	public int getGroup1() {
		return group1;
	}
	public int getGroup2() {
		return group2;
	}
	
	
}
class JoinComparator implements Comparator<Join> {

	@Override
	public int compare(Join j1, Join j2) {
		if(j1.group1!=j2.group1) {
			return j2.group1-j1.group1;
		}
		return j2.group2-j1.group2;
	}
	
}
