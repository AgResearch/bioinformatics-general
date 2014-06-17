package diJest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
/**
 * Class parsing restrict results then write the fragments length in file(s) nammed in the following format : input_Enzyme_Sequence.out.
 * Each fragment length is separated by a line break.
 * @author Clément DELESTRE
 * @version 1.0
 * @see LoopingRestrict
 * @see ParsingResult
 */
public class ParsingRestrictResult implements ParsingResult {


	/**
	 * Do we need to parsing or not
	 */
	protected boolean parsing;

	/**
	 * The input file
	 */
	protected Path inputFile;
	/**
	 * The enzymes file
	 */
	protected Path enzymesFile;

	/**
	 * list of enzymes
	 */
	protected ArrayList<String> listofEnzymes;

	/**
	 * list of sequenceFile
	 */
	protected ArrayList<String> listofInputFile;
	/**
	 * Regex to get the sequence name
	 */
	protected Pattern sequenceLine = Pattern.compile("^#(.*)Sequence:(\\W+)(\\w+)(.*)");
	/**
	 * Regex to get the 1st line of a fragment part
	 */
	protected Pattern fragmentFirstLine = Pattern.compile("^#(.*)Fragment(.*)lengths:$");
	/**
	 * Regex to get the last line of a fragment part
	 */
	protected Pattern fragmentLastLine = Pattern.compile("^#-+");
	/**
	 *  Regex to get the fragment part
	 */
	protected Pattern fragmentLine = Pattern.compile("^#(\\s+)(\\d+)$");

	/**
	 * All fragments' size
	 */
	protected ArrayList<Integer> fragments;

	/**
	 * All the output files
	 */
	protected ArrayList<Path> outputFiles;
	/**
	 * ParsingRestrictResult object
	 * @param sequenceFile
	 * @param enzymesFile
	 */
	public ParsingRestrictResult(Path inputFile,Path enzymesFile){
		this.inputFile=inputFile;
		this.enzymesFile=enzymesFile;
		parsing=false;
		fragments=new ArrayList<Integer>();
		outputFiles=new ArrayList<Path>();
	}


	@Override
	public void parsingRun() throws IOException {
		computeEnzymesFile();
		computeInputFile();
		parsingRestrictResults();
	}


	private void parsingRestrictResults() throws IOException {
		for (String file : listofInputFile){
			BufferedReader br =  Files.newBufferedReader(Paths.get(file),StandardCharsets.UTF_8);
			String line;
			System.out.println("Let's parse the file : "+file);
			boolean toWrite=false;
			String seqname=new String();
			while ((line=br.readLine())!=null){
				Matcher m = sequenceLine.matcher(line);
				if ( m.matches()){ 
					seqname=m.group(3);
					fragments.clear();
				}
				m = fragmentFirstLine.matcher(line);
				if ( m.matches()){ 
					System.out.println("Begin parsing fragments for "+seqname+" sequence...");
					parsing=true;

				}
				m = fragmentLastLine.matcher(line);
				if ( ( m.matches() ) && (parsing) ){ 
					System.out.println("Finish parsing fragments for "+seqname+" sequence...");
					parsing=false;
					toWrite=true;
				}
				m = fragmentLine.matcher(line);
				if ( ( m.matches() ) && (parsing) ){ // parsing fragments
					fragments.add(Integer.parseInt(m.group(2)));
				}

				if (toWrite){
					runWrite(seqname,file);
					toWrite=false;
				}
			}
			br.close();
		}
	}



	private void runWrite(String seqname,String file) throws IOException {
		int i=0;
		String toWrite = new String();
		for (Integer fragment : fragments){
			if (i==0){
				toWrite+=fragment;
			}
			else  {
				toWrite+="\n"+fragment;
			}
			i++;
		}

		String newFileName = file.replaceAll("_restrict.out", "_"+seqname+".out");
		outputFiles.add(Paths.get(newFileName));
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(newFileName), StandardCharsets.UTF_8);
		writer.write(toWrite);
		writer.close();
	}

	@Override
	public ArrayList<String> getEnzymes() {
		return listofEnzymes;
	}
	@Override
	public ArrayList<String> getInputFilesName(){
		return listofInputFile;
	}

	private void computeEnzymesFile() throws IOException {
		listofEnzymes=new ArrayList<String>();
		BufferedReader br =  Files.newBufferedReader(enzymesFile,StandardCharsets.UTF_8);
		String line;
		while ((line=br.readLine())!=null){
			listofEnzymes.add(line);
		}
	}

	private void computeInputFile() {
		listofInputFile=new ArrayList<String>();
		for (String enzyme : listofEnzymes){
			listofInputFile.add(inputFile.getFileName()+"_"+enzyme+"_restrict.out");
		}
	}

	@Override
	public ArrayList<Path> getOutputFiles() {
		return outputFiles;
	}
}