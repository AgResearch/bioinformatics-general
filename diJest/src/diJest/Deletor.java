package diJest;

import java.nio.file.Path;
import java.util.ArrayList;
/**
 * Class deleting several files
 * @author Clément DELESTRE
 * @version 1.0
 * @see DeleteI
 */
public class Deletor implements DeleteI {
	protected ArrayList<Path> files; 
	/**
	 * Create Deletor object
	 * @param files to delete
	 */
	public Deletor( ArrayList<Path> files){
		this.files=files;
	}

	@Override
	public  void deleteFiles() { 
		for (Path file : files){
			boolean b = file.toFile().delete();
			if (!b)
				System.err.println("Clean failed for file : "+file);
		}
	}
}
