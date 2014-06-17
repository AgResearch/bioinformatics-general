package diJest;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
/**
 * Interface that execute a R script for several files
 * @author Clément DELESTRE
 * @version 1.0
 */
public interface RloopInterface {
	/**
	 * Set the files
	 * @param files
	 */
	public void setFiles(ArrayList<Path> files);
	/**
	 * Launch the script
	 */
	public void lauchScript();
	/**
	 * Set the script.
	 * @param R script
	 */
	public void setScript(Path script);
}
