package diJest;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
/**
 * Class computing a R script for several files
 * @see RloopInterface
 * @author Clément DELESTRE
 * @version 1.0
 *
 */
public class Rloop implements RloopInterface{
	/**
	 * The R script
	 */
	protected File Rscript;
	/**
	 * File(s)
	 */
	protected ArrayList<String> files;
	/**
	 * Constructor with only a script
	 * @param rscript
	 */
	public Rloop(File rscript) {
		Rscript = rscript;
		files = new ArrayList<String>();
	}
	/**
	 * Constructor with script and file(s)
	 * @param rscript
	 * @param files
	 */
	public Rloop(File rscript,ArrayList<String> files) {
		Rscript = rscript;
		this.files = files;
	}
	@Override
	public void setScript(File script){
		this.Rscript=script;
	}
	@Override
	public void setFiles(ArrayList<String> files){
		this.files=files;
	}
	@Override
	public void lauchScript() {
		for (String file : files){
			System.out.println("Computing R script for "+file+"...");
			Runtime rt = Runtime.getRuntime();
			Process pr = null;
			try {
				pr = rt.exec("Rscript "+Rscript+" "+file);		
			}
			catch (Exception e) {				
				e.printStackTrace();			
			}
			try {
				pr.waitFor();
			} 
			catch (InterruptedException e1) {
				e1.printStackTrace();
			}
			String msg = "";
			if (pr.exitValue()!=0){
				String line = "";
				BufferedReader reader = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
				try {
					while((line = reader.readLine()) != null) {
						msg +=line;	
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
				System.err.println("***** ERROR with R script \n"+msg);
			}
		}
	}
}
