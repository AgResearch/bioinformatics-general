package diJest;

import java.io.IOException;
import java.util.ArrayList;
/**
 * Loop on restrict command (belonging to the EMBOSS package)
 * @author Clément DELESTRE
 * @version 1.0
 */
public interface LoopingRestrictInterface {

	/**
	 * Launch the restrict command
	 * @throws IOException 
	 */
	public void restrictRun() throws IOException;
	/**
	 * Get all the enzymes
	 * @return list of enzymes
	 */
	public ArrayList<String> getEnzymes();
	/**
	 * Get all the outputfile names
	 * @return list of outputfile names
	 */
	public ArrayList<String> getOutputFiles();
}
