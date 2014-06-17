package jast;

import java.io.File;
import java.io.FilenameFilter;
import java.util.regex.Pattern;
/**
 * Class that filter files.
 * @author Clément DELESTRE
 * @version 1.0
 * @since 1.0
 */

public class Filter implements FilenameFilter{
	/**
	 * The Pattern
	 * @since 1.0
	 */
	protected Pattern pattern;
	/**
	 * Create a filter.
	 * @param regex
	 */
	public Filter(String regex){
		pattern=Pattern.compile(regex);
	}

	@Override
	public boolean accept(File dir, String name) {
		return pattern.matcher(name).matches();
	}

}
