package diJest;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
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
	protected Path Rscript;
	/**
	 * File(s)
	 */
	protected ArrayList<Path> files;
	/**
	 * Constructor with only a script
	 * @param rscript
	 */
	public Rloop(Path rscript) {
		Rscript = rscript;
		files = new ArrayList<Path>();
	}
	/**
	 * Constructor with script and file(s)
	 * @param rscript
	 * @param files
	 */
	public Rloop(Path rscript,ArrayList<Path> files) {
		Rscript = rscript;
		this.files = files;
	}
	@Override
	public void setScript(Path script){
		this.Rscript=script;
	}
	@Override
	public void setFiles(ArrayList<Path> files){
		this.files=files;
	}
	@Override
	public void lauchScript() {
		for (Path p : files){
			System.out.println("Computing R script for "+p+"...");
			Runtime rt = Runtime.getRuntime();
			Process pr = null;
			try {
				pr = rt.exec("Rscript "+Rscript+" "+p.toString());		
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
