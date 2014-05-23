package diJest;

import java.io.File;
import java.util.ArrayList;
/**
 * Class deleting several files
 * @author Clément DELESTRE
 * @version 1.0
 * @see DeleteI
 */
public class Deletor implements DeleteI {
	protected ArrayList<String> files; 
/**
 * Create Deletor object
 * @param files to delete
 */
	public Deletor( ArrayList<String> files){
		this.files=files;
	}

	@Override
	public  void deleteFiles() { 
		for (String file : files){
			boolean b = new File(file).delete();
			if (!b)
				System.err.println("Clean failed for file : "+file);
		}
	}
}
