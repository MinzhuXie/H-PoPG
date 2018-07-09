package main;

import io.BAMfileReader;
import io.FragsfileReader;
import io.VCFProcessor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import algorithms.PolyplotypingAlgorithm;

import dataStructure.Block;
import dataStructure.Fragment;
import dataStructure.SAMFileDescriptor;

import shares.Tools;

public class PolyPlotyping {	
	int n_ploid = 3;
	String bamFileName = null;
	String version = "0.1.1"; //added by Xie July 8, 2018
	List<String> fragsFileNames = new ArrayList<String>();
	String vcfFileName = null;
	String outPhasedfileName = "phased_polyplotyping_BoP.txt"; 
	VCFProcessor vcf_processor;
	BAMfileReader sam_processor;
	FragmentComputer frag_computer;
	FragsfileReader frag_processor;
	public List<SAMFileDescriptor> sam_file_descriptors = new ArrayList<SAMFileDescriptor>();
	public List<Fragment> f;
	int max_size = 0; 
	double weight_factor = 0.9;
	private String algorithmClassName = "algorithms.HPBOP2Alg";
	private PolyplotypingAlgorithm algClass = null;
	
	public String getAlgorithmClassName() {
		return algorithmClassName;
	}
	
	public void setAlgorithmClassName(String algorithmClassName) {
		this.algorithmClassName = algorithmClassName;
		try {
			algClass = (PolyplotypingAlgorithm)Class.forName(algorithmClassName).newInstance();
		} catch (InstantiationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			System.out.println("Access denied: "+algorithmClassName);
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			System.out.println("the class not found: "+algorithmClassName);
			e.printStackTrace();
		}
	}
	
	private void printUsage()
	{
		System.out.println("H-PopG "+version+". Author: Xie Minzhu. Email: xieminzhu@hotmail.com \n"); //added by Xie July 8, 2018
		System.out.println("Usage: \n \t java -jar H-PoPG [-p ploidy (default is 3) -m maxsize (default is ploidy*ploidy*10) -w weight (default is 0.9)] -v[cf] vcffilename -b[am] bamfilenme -o outputPhasedFilename");
		System.out.println("or \t java -jar  H-PoPG [-p ploidy (default is 3) -m maxsize (default is ploidy*ploidy*10) -w weight (default is 0.9)] -f[rag] alignedFragsFileName [-v[cf] vcffilename] -o outputPhasedFilename");
		System.out.println("If a BAM File has been provide, a VCF file is needed and the Frag File is ignored!");
		System.out.println("\t java -jar H-PoPG -h to print this help information");
		System.exit(-1);
	}
	private void checkParameters(){
		if(bamFileName==null && fragsFileNames.size() == 0)
		{
			System.out.println("There are no input file.");
			printUsage();
		}
		if(bamFileName!=null  && vcfFileName == null && fragsFileNames.size() == 0)
		{
			System.out.println("Need vcfFile to parse bamFile");
			printUsage();			
		}
		if(max_size == 0)
			max_size = n_ploid * n_ploid * 10;
	}
	
	public void go(){
		long time = System.currentTimeMillis();
		if(vcfFileName!=null){
			vcf_processor = new VCFProcessor(n_ploid, 1, this);
			 try
	         {
	                vcf_processor.loadVCF(vcfFileName);
	         }	
			 catch(IOException e)
	         {
	                System.err.println("Error while reading VCF file "+ vcfFileName + e.getMessage());   
	                System.exit(-1);
	         }
		}
		if(bamFileName!=null) //read the bam file
		{
			String[] fNames = new String[1];
			fNames[0] = bamFileName;					        	 
	       	List<File> files =Tools.getFiles(fNames, ".bam");
	       	for (int j = 0; j < files.size(); j++) 
	       	{
	            sam_file_descriptors.add(new SAMFileDescriptor(0, ((File)files.get(j)).getAbsolutePath()));
	        }
					
	       	String frag_output = "xieminzhuBAM2frag.txt";
	         
	        String ref_name = vcf_processor.vcf_reference_names.get(0);
	         
	         sam_processor = new BAMfileReader(this);
	         
	         sam_processor.init(vcf_processor, ref_name, null, frag_output);
		         
		     frag_processor = new FragsfileReader(this);
		     
		     try
	         {
	                frag_computer = new FragmentComputer(frag_output, null, this);
	         }
		     catch(IOException e)
	         {
	                System.err.println("Could not create fragment file " +frag_output + e.getMessage());
	                System.exit(-1);	
	         }
		     try
	         {
	                int added = frag_computer.go();
	                if(added > 0)
	                	fragsFileNames.add(frag_output);
	         }
	         catch(IOException e)
		     {
		                System.err.println("Error while processing fragment file " +frag_output + e.getMessage());
		                System.exit(-1);
		
		     }
		}		 
		
		//read the frag file 
		FragsfileReader loader = new FragsfileReader(this);
		try{	
			f = loader.loadFragments(fragsFileNames); 
		}
		catch(IOException e)
	    {
	                System.err.println("Error while reading all the fragment files " + e.getMessage());
	                System.exit(-1);
	
	    }	
		
		BlocksBuilder blocksBuilder = new BlocksBuilder();
		List<Block> blocks = blocksBuilder.buildBlocks(new Block(f));
		// blocks �а�����ÿ���飬ÿ�����а����˸ÿ��е�Ƭ����Ϣ
		
		for(Block b:blocks) {
			//��ָ���㷨��������			
			if(!Polyphasing (b))
			{
				System.out.println("Failure to run the algorithm: " + algorithmClassName);
				System.out.println("The block length: "+ b.length() + " \t The number of fragments: " + b.getFragments().size());
				return; 
			}
		}
		//�����㷨ִ�����
		System.out.println("Number of blocks:" +blocks.size());
		
		PrintStream out = null;
		try {
			out = new PrintStream(outPhasedfileName);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			System.err.println("Can not create the output file: " + outPhasedfileName);
			e.printStackTrace();
		}
		int totalPhased=0;
		int totalCalls=0;
		int totalMEC =0;		
			
		for(Block b:blocks) {
			//	System.out.println("Offset: "+b.getFirstPos()+" Phased: "+b.getPhased()+ " Calls: "+b.getCalls()+ " MEC: "+b.getMEC());
				totalPhased += b.getPhased();
				totalCalls += b.getCalls();
				totalMEC += b.getMEC();
				printBlock(b, out);
		}
		
		out.close();
		double diff = System.currentTimeMillis() - time;
		diff/=1000;
		System.out.println("Phased: "+totalPhased+ " Calls: "+totalCalls+ " MEC: "+totalMEC+" Time(s): "+diff);	
		
	}
	
	private void printBlock(Block b, PrintStream out) {
		if(b.getPhased()<2)
			return; 
		int len = b.length();
		out.println("BLOCK: offset: "+(b.getFirstPos()+1)+" len: "+b.length()+" phased: "+b.getPhased());
		char[][] hap = b.getPolyH();
	//	boolean printReverse = r.nextBoolean();
		for(int i=0;i<len;i++) 
		{		
			int pos = b.getFirstPos()+i;			
			pos++;
			out.print(""+pos);			
			for(int j = 0; j<n_ploid; j++)
				out.print("\t"+hap[i][j]);	
			out.print("\n");
		}
		out.println("********");
		out.flush();		
	}
	
	
	public boolean Polyphasing (Block b) 
	{		
		if(algClass == null) {
			System.out.println("the class for the algorithm is null.");
			return false;
		}
		algClass.setParameters(n_ploid, max_size, weight_factor, vcf_processor);	
		return(algClass.buildHaplotype(b));
	}
	
	public static void main(String[] args) {		
		
		PolyPlotyping alg = new PolyPlotyping(); 
		String algorithmName = "HPBOP2";   
		int i = 0;
		//get the command line parameters
		while(i<args.length && args[i].startsWith("-")) {
			if("-p".equals(args[i])) 
			{
				i++;
				alg.n_ploid = Integer.parseInt(args[i]); 				
			} 
//			else if("-a".equals(args[i])) {
//				i++;
//				algorithmName = args[i];  //-a �㷨��
//			}
			else if("-m".equalsIgnoreCase(args[i])) {   //-m  max number of itermediate solutions 
				i++;
				alg.max_size = Integer.parseInt(args[i]);
			}else if("-w".equalsIgnoreCase(args[i])) {    //-w  weight-rate
				i++;
				alg.weight_factor = Double.parseDouble(args[i]); 
			}
			else if("-vcf".equals(args[i]) || "-v".equals(args[i])) 
			{
				i++;
				alg.vcfFileName = args[i]; 
			} 
			else if("-bam".equals(args[i]) || "-b".equals(args[i]))
			{
				i++;
				alg.bamFileName = args[i]; 
			}
			else if("-frag".equals(args[i]) || "-f".equals(args[i]))
			{
				i++;
				alg.fragsFileNames.add(args[i]); 
			}
			else if("-o".equals(args[i]))
			{
				i++;
				alg.outPhasedfileName = args[i]; 
			}			
			else if("-h".equals(args[i]))
			{
				alg.printUsage();
			}
			i++;
		}
		alg.checkParameters();
		
		alg.setAlgorithmClassName("algorithms."+algorithmName+"Alg");
		
		alg.go();
	}
}
