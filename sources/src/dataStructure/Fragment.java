package dataStructure;

public class Fragment {
	private String id;
	private int firstPos=0;
	private int lastPos=0;
	private int nCalls=0;
	private String extendedCalls = null;
	
	private double [] extendedProbabilities;
	public static final char ALLELE1CHAR = '0';
	public static final char ALLELE2CHAR = '1';
	public static final char NODATACHAR = '-';
	public Fragment(String id, int firstPos, String extendedCalls) {
		this (id,firstPos,extendedCalls,null);
	}
	public Fragment(String id, int firstPos, String extendedCalls,double [] extendedProbabilities) {
		super();
		this.id = id;
		this.firstPos = firstPos;
		this.lastPos = firstPos +extendedCalls.length()-1;
		this.extendedCalls = extendedCalls;
		nCalls = 0;
		for(int i=0;i<extendedCalls.length();i++) {
			if(extendedCalls.charAt(i)!=NODATACHAR) {
				nCalls++;
			}
		}
		this.extendedProbabilities = extendedProbabilities;
	}
	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}
	
	/**
	 * @return the firstPos
	 */
	public int getFirstPos() {
		return firstPos;
	}
	/**
	 * @return the lastPos
	 */
	public int getLastPos() {
		return lastPos;
	}
	
	public int getnCalls() {
		return nCalls;
	}
	public int length() {
		return lastPos - firstPos + 1;
	}
	
	public int getDistanceWithQuals (Fragment f2) {
		int relPos1 = Math.max(0, f2.getFirstPos()-firstPos);
		int relPos2 = Math.max(0, firstPos-f2.getFirstPos());
		String calls2 = f2.getExtendedCalls();
		double [] extProbs2 = f2.getExtendedProbabilities();
		double distance =0;
		while(relPos1 <extendedCalls.length() && relPos2<calls2.length()) {
			char c1 = extendedCalls.charAt(relPos1);
			char c2 = calls2.charAt(relPos2);
			if(c1!=NODATACHAR && c2!=NODATACHAR) {
				double p1 = extendedProbabilities[relPos1];
				double p2 = extProbs2[relPos2];
				if(c1 != c2) {
					distance+=(p1+p2)/2;
				} else {
					distance-=(p1+p2)/2;
				}
			}
			relPos1++;
			relPos2++;
		}
		return (int)Math.round(distance);
	}
	
	public int getHamming2(String sequence,int start) {
		int relPos1 = Math.max(0, start-firstPos);
		int relPos2 = Math.max(0, firstPos-start);
		String calls = getExtendedCalls();
		int disagree =0;
//		int total =0;
		while(relPos1 <calls.length() && relPos2<sequence.length()) {
			char c1 = calls.charAt(relPos1);
			char c2 = sequence.charAt(relPos2);
			if(c1!=NODATACHAR && c2!=NODATACHAR) {
				if(c1 != c2) {
					disagree++;
				} 
//				total++;
			}
			relPos1++;
			relPos2++;
		}
//		if(total > 0) {
//			return 2*disagree-total;			
//		}
//		return 0;
		return disagree;
	}
	public int getHamming2(Fragment f2) {
		return getHamming2(f2.getExtendedCalls(),f2.getFirstPos());
	}
	public int getHammingDistance(String sequence,int start) {
		int relPos1 = Math.max(0, start-firstPos);
		int relPos2 = Math.max(0, firstPos-start);
		String calls = getExtendedCalls();
		int disagree =0;
		while(relPos1 <calls.length() && relPos2<sequence.length()) {
			char c1 = calls.charAt(relPos1);
			char c2 = sequence.charAt(relPos2);
			if(c1!=NODATACHAR && c2!=NODATACHAR) {
				if(c1 != c2) {
					disagree++;
				} 
			}
			relPos1++;
			relPos2++;
		}
		return disagree;
	}
	public int getHammingDistance(Fragment f2) {
		return getHammingDistance(f2.getExtendedCalls(), f2.getFirstPos());
	}
	public int getOverlappingCount(Fragment f2) {
		return getOverlappingCount(f2.getExtendedCalls(),f2.getFirstPos());
	}
	public int getOverlappingCount(String sequence,int start) {
		int relPos1 = Math.max(0, start-firstPos);
		int relPos2 = Math.max(0, firstPos-start);
		String calls = getExtendedCalls();
		int overlap =0;
		while(relPos1 <calls.length() && relPos2<sequence.length()) {
			char c1 = calls.charAt(relPos1);
			char c2 = sequence.charAt(relPos2);
			if(c1!=NODATACHAR && c2!=NODATACHAR) {
				overlap++; 
			}
			relPos1++;
			relPos2++;
		}
		return overlap;
	}
	public String getExtendedCalls() {
		return extendedCalls;
	}
	
	
	public void setExtendedCalls(String eCalls) {
		extendedCalls = eCalls;
	}
	
	
	/**
	 * @return the extendedProbabilities
	 */
	public double[] getExtendedProbabilities() {
		return extendedProbabilities;
	}
	/**
	 * @param extendedProbabilities the extendedProbabilities to set
	 */
	public void setExtendedProbabilities(double[] extendedProbabilities) {
		this.extendedProbabilities = extendedProbabilities;
	}
	public char getCall(int pos) {
		int relPos = pos-firstPos;
		if(relPos>=0 && relPos < extendedCalls.length()) {
			return extendedCalls.charAt(relPos);
		}
		return NODATACHAR;
	}
	public char getCallbyIndex(int relPos) {		
		if(relPos>=0 && relPos < extendedCalls.length()) {
			return extendedCalls.charAt(relPos);
		}
		return NODATACHAR;
	}
	
}
