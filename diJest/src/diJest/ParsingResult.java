package diJest;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
/**
 * Interface to parsing result files
 * @author Clément DELESTRE
 * @version 1.0
 */
public interface ParsingResult {

	
	/**
	 * Parsing file
	 * @throws IOException 
	 */
	public void parsingRun() throws IOException;
	
	/**
	 * Get all the enzymes
	 * @return list of enzymes
	 */
	public ArrayList<String> getEnzymes();
	
	/**
	 * Get all the input files
	 * @return list of input files
	 */
	public ArrayList<String> getInputFilesName();

	/**
	 * Get all the output files
	 * @return list of input files
	 */
	public ArrayList<Path> getOutputFiles();
	
	
}
