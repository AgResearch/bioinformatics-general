package diJest;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
/**
 * Class looping on restrict command (belonging to EMBOSS library) for several enzymes. Output files are in format input_Enzyme_restrict.out
 * @author Clément DELESTRE
 * @version 1.0
 * @see LoopingRestrictInterface
 */
public class LoopingRestrict implements LoopingRestrictInterface {
	/**
	 * The input file
	 */
	protected Path sequenceFile;
	/**
	 * The enzymes file
	 */
	protected Path enzymesFile;
	/**
	 * list of enzymes
	 */
	protected ArrayList<String> listofEnzymes;
	/**
	 * list of outputfiles
	 */
	protected ArrayList<String> outputFiles;
	/**
	 * Creating a LoopingRestrict object
	 * @param sequenceFile
	 * @param enzymesFile
	 */
	public LoopingRestrict(Path inputFile, Path enzymesFile) {
		this.sequenceFile=inputFile;
		this.enzymesFile=enzymesFile;
		outputFiles=new ArrayList<String>();
	}
	@Override
	public void restrictRun() throws IOException {
		computeEnzymesFile();
		for (String enzyme : listofEnzymes){
			Runtime rt = Runtime.getRuntime();
			Process pr = null;		
			File file = new File("");
			String outputFileName = file.getAbsolutePath()+File.separator+sequenceFile.getFileName()+"_"+enzyme+"_restrict.out";
			try {
				pr = rt.exec("restrict -sequence "+sequenceFile+" -outfile "+outputFileName+" -sitelen 2 -enzymes "+enzyme+" -fragments Y");		
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
				while((line = reader.readLine()) != null) {
					msg +=line;	
				}
				reader.close();
				System.out.println("***** ERROR : \n"+msg);
			}
			outputFiles.add(outputFileName);
			System.out.println("Restrict command done for enzyme "+enzyme);
		}
		System.out.println("Restrict command finish.");
	}
	/**
	 *  Read enzyme(s) file and add each enzyme to the listofEnzymes arrayList
	 * @throws IOException
	 */
	private void computeEnzymesFile() throws IOException {
		listofEnzymes=new ArrayList<String>();
		BufferedReader br =  Files.newBufferedReader(enzymesFile,StandardCharsets.UTF_8);
		String line;
		while ((line=br.readLine())!=null){
			listofEnzymes.add(line);
		}
		br.close();
	}
	@Override
	public ArrayList<String> getEnzymes() {
		return listofEnzymes;
	}
	@Override
	public ArrayList<String> getOutputFiles() {
		return outputFiles;
	}
}
