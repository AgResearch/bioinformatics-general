package diJest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
/**
 * Class parsing restrict output files and extract fragments sequence thanks to extractseq command (belonging to the EMBOSS package)
 * @author Clément DELESTRE
 * @version 1.0
 * @see LoopingRestrict
 * @see ExtractSequences
 */

public class Extractor implements ExtractSequences {

	/**
	 * Regex to get the sequence name file
	 */
	protected Pattern sequenceLine = Pattern.compile("^# Sequence: (\\w+)(.*)from: (\\d+)(.*)to: (\\d+)");
	/**
	 * Regex to get the number of hits
	 */
	protected Pattern hitLine = Pattern.compile("^# HitCount: (\\d+)(.*)");
	/**
	 * Regex to get the fragments line
	 */
	protected Pattern fragmentsLine = Pattern.compile("^(\\s*)(\\d+)(\\s+)(\\d+)(\\s+)(\\+|-)(\\s+)(\\w+)(\\s+)(\\w+)(\\s+)(\\d+)(\\s+)(\\d+)(\\s+)(.*)");
	/**
	 * list of enzymes
	 */
	protected ArrayList<String> listofEnzymes;

	/**
	 * Sequence file name
	 */
	protected  String seqFileName;

	/**
	 * List of inpute files
	 */
	protected ArrayList<String> inputFiles;
	/**
	 * Maximum fragments length
	 */
	protected int maxLength=150;
	/**
	 * Minimum fragments length
	 */
	protected int minLength=100;
	/**
	 * Intermediate File
	 */
	protected Path intermediateFile=Paths.get("temp_extract");
	/**
	 * OutputFile
	 */
	protected File outputFile;
	/**
	 * Writer
	 */
	protected PrintWriter writer;
	/**
	 * First part of extract command
	 */
	protected String beginExtractSeq; 
	/**
	 * Second part of extract command
	 */
	protected String middleExtractSeq="";
	/**
	 * Third part of extract command
	 */
	protected String endExtractSeq="-outseq "+intermediateFile+" -separate Y";
	/**
	 * Beginning of extraction
	 */
	protected int begin;
	/**
	 * End of extraction
	 */
	protected int end;
	/**
	 * Number of enzymes cut
	 */
	protected int nbCut;
	/**
	 * Number of fragments
	 */
	protected int nbFrags;
	/**
	 * Use extract command each N times
	 */
	protected int moduloExtFrags=42;
	/**
	 * Current line when parsing the input file
	 */
	protected int numberlines;
	/**
	 * Runtime to exec command
	 */
	protected Runtime rt;
	/**
	 * Process to exec command
	 */
	protected Process pr;
	/**
	 * Name of output files based on inputfiles name
	 */
	protected HashMap<String,String> outputfiles;
	/**
	 * The current sequence's name
	 */
	protected String seqName;
	/**
	 * Extractor need one or several inputFile to parse and the sequence file
	 * @param inputFiles
	 * @param seqFileName
	 */
	public Extractor( ArrayList<String> inputFiles,String seqFileName) {
		this.inputFiles = inputFiles;
		outputfiles=new HashMap<String,String>();
		this.seqFileName=seqFileName;
		beginExtractSeq="extractseq -sequence "+seqFileName+" -regions";
	}

	@Override
	public void setInputFiles(ArrayList<String> files){
		inputFiles=files;
	}

	@Override
	public void extractSequences() {
		if (checkTempFile()){
			try {	
				run();
			} catch (IOException e) {
				System.err.println("*** ERROR WITH FILE : ");
				e.printStackTrace();
			}
		}
		else{
			System.err.println("A file nammed like the temp file use for extract sequence already exist (nbFrags.e : "+intermediateFile+" ). Please delete it or move it.\nSystem will exit.");
			System.exit(1);	
		}

		boolean b = intermediateFile.toFile().delete();
		if (!b)
			System.err.println("Error deleting temporary file "+intermediateFile);
	}


	/**
	 * Main method
	 * @throws IOException
	 */
	private void run() throws IOException {
		computeOutputFiles();
		for (String s : inputFiles){
			outputFile = new File(outputfiles.get(s));
			System.out.println("creation du fichier : "+outputfiles.get(s));
			System.out.println(outputfiles.get(s));
			System.out.println("S : "+s);
			try {
				outputFile.createNewFile();
			} catch (IOException e) {
				System.out.println("Can't create the output file "+outputFile+ " system will exit");
				System.exit(1);
				e.printStackTrace();
			}
			extractSequences(Paths.get(s));
		}
		System.out.println("Extraction for all fragments done.");
	}

	/**
	 * Method that parse file and execute restrict command
	 * @param input
	 * @throws IOException
	 */
	private void extractSequences(Path input) throws IOException {

		BufferedReader br =  Files.newBufferedReader(input,StandardCharsets.UTF_8);
		// current line
		String line;
		// check if we are at the first line or not
		boolean first=true;
		// number of current line
		numberlines=1;
		// total number of fragments
		nbFrags=0;
		// where was the last fragment
		int lastnb=0;

		while ((line=br.readLine())!=null){
			Matcher m = sequenceLine.matcher(line);
			if ( m.matches()){ 
				seqName= m.group(1);
				begin= Integer.parseInt(m.group(3));
				end=Integer.parseInt(m.group(5));
			}
			m = hitLine.matcher(line);
			if ( m.matches()){ 
				nbCut=Integer.parseInt(m.group(1));
			}
			m = fragmentsLine.matcher(line);
			if ( m.matches()){ 
				int nb=Integer.parseInt(m.group(12));
				// we have to be carefull when it s the 1st fragment
				if (first){
					first=false;
					if (isLengthOk(begin,end)){
						middleExtractSeq=begin+","+lengthMax(begin,nb)+",";
					}
					lastnb=nb;
				}
				else {
					numberlines++;
					lastnb++;
					if (isLengthOk(lastnb,nb)){
						middleExtractSeq+= ""+lastnb+","+lengthMax(lastnb,nb);
						nbFrags++;
						if (nbFrags%moduloExtFrags==0){
							execExtract();
							middleExtractSeq="";
						}
						else {
							middleExtractSeq= middleExtractSeq+",";
						}
					}
					lastnb=nb;
					// If it s the last one, we have to execute the extract command.
					if (numberlines==nbCut){
						lastnb++;
						if (isLengthOk(lastnb,end)){
							middleExtractSeq+=""+lastnb+","+lengthMax(lastnb,end);
							execExtract();
						}
					}
				}

			}
		}
	}
	/**
	 * Execute extract command
	 * @throws IOException
	 */
	private void execExtract() throws IOException {
		float loading = (float)numberlines/nbCut*100;
		DecimalFormat df = new DecimalFormat("###.##");
		System.out.print("Loading... "+df.format(loading)+"%\r");
		rt = Runtime.getRuntime();
		pr = null;
		try {
			pr = rt.exec(beginExtractSeq+" "+middleExtractSeq+" "+endExtractSeq);
		}
		catch (IOException e) {
			System.err.println("IO Exception here.");
			e.printStackTrace();
		}

		try {
			pr.waitFor();
		} 
		catch (InterruptedException e1) {
			System.err.println("Interrupted Exception here.");
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
			System.err.println("Failed to use the command extract: \n"+msg);
		}

		pr.destroy();
		writeOutputFile();
	}
	private void writeOutputFile() throws IOException {
		BufferedReader br =  Files.newBufferedReader(intermediateFile,StandardCharsets.UTF_8);
		String line;
		String toWrite="";
		writer = new PrintWriter ( new BufferedWriter(new FileWriter(outputFile, true)));
		while ((line=br.readLine())!=null){
			toWrite+="\n"+line;
		}
		writer.write(toWrite+"\n");
		writer.close();
	}
	/**
	 * Check the length, if fragment is too big, return little+maxLength if not return big
	 * @param little
	 * @param big
	 * @return maximum possible length
	 */
	private int lengthMax(int little,int big){
		if (  big-little >  maxLength ){
			return little+maxLength;
		}
		return big;
	}
	/**
	 * Check the length
	 * @param little
	 * @param big
	 * @return boolean
	 */
	private boolean isLengthOk(int little, int big) {
		if (  big-little <  minLength ){
			return false;
		}
		return true;
	}
	/**
	 * Compute outputfiles name
	 */
	private void computeOutputFiles() {
		for(String s : inputFiles){
			outputfiles.put(s, s+"_extractFragmentsSequences.fasta");
		}
	}
	/**
	 * Check if the temp file already exist or not
	 * @return true if exist
	 */
	private boolean checkTempFile() {
		if  (intermediateFile.toFile().exists())
			return false;
		else 
			return true;
	}
	@Override
	public void setMaxLength(int maxLength) {
		this.maxLength=maxLength;
	}
	@Override
	public void setMinLength(int minLength) {
		this.minLength=minLength;
	}
}
