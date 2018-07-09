package io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Comparator;

import main.PolyPlotyping;


import dataStructure.Fragment;

public class FragsfileReader {
	PolyPlotyping alg;
	
	
	public FragsfileReader(PolyPlotyping pB){
		alg = pB; 
	}
	public List<Fragment> loadFragments(List<String> fs) 
	throws IOException {
		List<Fragment> fragments = new ArrayList<Fragment>();
		for(String filename: fs)
		{
			FileInputStream fis = new FileInputStream(filename);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			String line=in.readLine();
			while(line!=null) {
				String [] items = line.split(" |\t");
				if(items.length > 2) {
					try {
						StringBuilder calls = new StringBuilder();
						String qualityScores = null;
						int startPos = Integer.parseInt(items[2])-1;
						int pos = startPos;
						for(int j=2;j<items.length;j+=2) {
							//Assumes that start positions are one based
							if(j<items.length-1) {
								int nextStart = Integer.parseInt(items[j])-1;
								while(pos<nextStart) {
									calls.append(Fragment.NODATACHAR);
									pos++;
								}
								calls.append(items[j+1]);
								pos+=items[j+1].length();
							} else {
								qualityScores = items[j];
							}
						}					
						//Xie add Mar. 22. 2014
						if(calls.length()>1)
						{					
						//Mar. 22					
							int k=0;
							double [] extendedProbabilities = null;
							Fragment f;
							if(qualityScores!=null) {
								extendedProbabilities = new double[calls.length()];
								Arrays.fill(extendedProbabilities, 0);
								for(int j=0;j<calls.length();j++) {
									if(calls.charAt(j)!=Fragment.NODATACHAR) {
										extendedProbabilities [j] = calculateProbability(qualityScores.charAt(k));
										k++;
									}
								}
								
								f = new Fragment(items[1],startPos, calls.toString(),extendedProbabilities);
							} else {
								f = new Fragment(items[1],startPos, calls.toString());
							}
							fragments.add(f);					
						}		
						
					} catch(Exception e) {
						in.close();
						throw new IOException("Error reading line: "+ line+ " of file "+filename,e);
					}
				}				
				line=in.readLine();
			}
			in.close();
			fis.close();			
		}		
		Collections.sort(fragments,new FragmentsComparator());
		return fragments;
	}
	private double calculateProbability(char qual) {
		double score = (double)(qual -33);
		return 1 - Math.pow(10, -score /10);
	}
	public List<Integer> loadVariantPositions(String filename, int posCol) throws IOException {
		List<Integer> variantPos = new ArrayList<Integer>();
		FileInputStream fis = new FileInputStream(filename);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line=in.readLine();
		while(line!=null) {
			String [] items = line.split(" |\t");
			try {
				variantPos.add(Integer.parseInt(items[posCol]));
			} catch(NumberFormatException e) {
				System.err.println("WARN: Could not read variant id in line: "+line );
			}
			line=in.readLine();
		}
		fis.close();
		return variantPos;
	}	
}

class FragmentsComparator implements Comparator<Fragment> {
	@Override
	public int compare(Fragment f1, Fragment f2) {		
		if(f1.getFirstPos()==f2.getFirstPos())
		{
			if(f1.getLastPos()==f2.getLastPos())
			{
				return f2.getExtendedCalls().compareTo(f1.getExtendedCalls());				
			}
			else
			   return f1.getLastPos()-f2.getLastPos();
			
		}			
		else 
			return f1.getFirstPos()-f2.getFirstPos();
	}
}
