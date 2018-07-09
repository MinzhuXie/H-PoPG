package shares;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMRecord;

import org.apache.commons.io.filefilter.WildcardFileFilter;


public class Tools {
	private static int PRIMARYMASK = 256;
	
	public static List<File> getFiles(String filenames[], String dir_filter)
    {
        File samFiles = new File(filenames[0]);
        List<File> files = new ArrayList<File>();
        if(samFiles.isDirectory())
        {
            FileFilter fileFilter = new WildcardFileFilter(dir_filter);
            File dir = samFiles;
            File tempFileList[] = dir.listFiles(fileFilter);
            for(int j = 0; j < tempFileList.length; j++)
                files.add(tempFileList[j]);

        } else
        {
            for(int j = 0; j < filenames.length; j++)
            {
                String filename = filenames[j];
                int indexOf = filename.lastIndexOf(File.separator);
                if(indexOf == -1)
                    indexOf = filename.lastIndexOf("\\");
                if(indexOf == -1)
                    indexOf = filename.lastIndexOf("/");
                String base = filename.substring(0, indexOf + 1);
                if(base.equals(""))
                    base = ".";
                String rFilter = filename.substring(indexOf + 1);
                if(filename.contains("*"))
                {
                    File dir = new File(base);
                    FileFilter fileFilter = new WildcardFileFilter(rFilter);
                    File files_to_process[] = dir.listFiles(fileFilter);
                    for(int k = 0; k < files_to_process.length; k++)
                        files.add(files_to_process[k]);

                } else
                {
                    files.add(new File(filename));
                }
            }

        }
        return files;
    }
	
	 public static String getReadName(SAMRecord samLine, String suffix_to_remove)
	 {
	        if(suffix_to_remove != null)
	        {
	            int index = samLine.getReadName().lastIndexOf(suffix_to_remove);
	            if(index >= 0)
	                return new String(samLine.getReadName().substring(0, samLine.getReadName().lastIndexOf(suffix_to_remove)));
	            else
	                return new String(samLine.getReadName());
	        } else
	        {
	            return new String(samLine.getReadName());
	        }
	  }
	 
	  public static boolean is_primary(SAMRecord samRecord)
	  {
	        return (samRecord.getFlags() & PRIMARYMASK) <= 0;
	  }	  

}
