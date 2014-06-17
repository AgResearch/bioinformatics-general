package jast;

/**
 * Exception throwed when user specify a forbidden option in a config file.
 * @author Clément DELESTRE
 * @version 1.0
 * @since 1.0
 */
public class ForbiddenOptionsException extends Exception{
	public ForbiddenOptionsException(String option,String command){
		super("The following option : "+option+" is forbidden for the command "+command);
	}
}
