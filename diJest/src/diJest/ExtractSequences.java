package diJest;

import java.nio.file.Path;
import java.util.ArrayList;
/**
 * Interface to extracting sequence between a maximum and a minimum given size
 * @author Clément DELESTRE
 * @version 1.0
 *
 */
public interface ExtractSequences {
	/**
	 * Extract sequences
	 */
	public void extractSequences();
	/**
	 * Set the maximum length of sequence to extract
	 * @param maxLength
	 */
	public void setMaxLength(int maxLength);
	/**
	 * Set the minimum length of sequence to extract
	 * @param minLength
	 */
	public void setMinLength(int minLength);
	/**
	 * Set files to parse
	 * @param files
	 */
	public void setInputFiles(ArrayList<Path> files);
}
