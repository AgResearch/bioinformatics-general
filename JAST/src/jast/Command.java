package jast;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

/**
 * Abstract class that execute commands
 * @author Clément DELESTRE
 * @version 1.0
 * @since 1.0
 */

public abstract class Command {
	/**
	 * Forbidden Options for one command
	 */
	protected  String[] forbiddenOptions=null;
	/**
	 * The command
	 */
	protected String command;
	/**
	 * The command and its option(s)
	 */
	protected ArrayList<String> totalCommand;
	/**
	 * The config file
	 */
	protected Path config; 
	/**
	 * The output file
	 */
	protected Path outputFile;
	/**
	 * Do we need to use a config file or not ? Sometime is not usefull
	 */
	protected boolean useConfig;

	/**
	 * Runtime to exec command
	 */
	protected static Runtime rt;
	/**
	 * Process to exec command
	 */
	protected static Process pr;
	/**
	 * Create a new command
	 * @param config file
	 * @param array Of Forbidden Options
	 */
	public Command(Path config,String[] arrayOfForbiddenOptions){
		this.config=config;
		useConfig=true;
		totalCommand = new ArrayList<String>();
		forbiddenOptions=arrayOfForbiddenOptions;
	}
	/**
	 * Get the output file
	 * @return outputfile
	 */
	public Path getOutputFile(){
		return outputFile;
	}
	/**
	 * Main method that execut the command
	 */
	public void exec(){
		totalCommand.add(0,command);
		if (useConfig && config!=null){
			try {
				readConfig();
			} catch (IOException e2) {
				System.err.println("[JAST_Error] Error with file "+config.toAbsolutePath()+" please check it.");
				e2.printStackTrace();
				System.err.println("[JAST_Error] System will exit.");
			}
			catch (ForbiddenOptionsException e2) {
				System.err.println("[JAST_Error] Error with file "+config.toAbsolutePath()+" : it contains an option already use by JAST. Please check it.");
				e2.printStackTrace();
				System.err.println("[JAST_Error] System will exit.");
				System.exit(1);
			}
		}
		System.out.println("[JAST] Running command : "+command);
		rt = Runtime.getRuntime();
		pr = null;
		try {
			//Note : exec(String[]) method >>>> exec(String) method
			String temp="";
			for (String s : totalCommand){
				temp+=" "+s+" ";
			}
			System.out.println("[JAST] Trying to execute : "+temp);
			pr = rt.exec(totalCommand.toArray(new String[totalCommand.size()]));
		}
		catch (IOException e) {
			System.err.println("[JAST_Error] IO Exception here.");
			e.printStackTrace();
		}
		try {
			pr.waitFor();
		} 
		catch (InterruptedException e1) {
			System.err.println("[JAST_Error] Interrupted Exception here.");
			e1.printStackTrace();
		}
		System.out.println("exit value : "+pr.exitValue());
		if (pr.exitValue()!=0 ){ 
			printError();
		}
		else {
			printOutputCommand();
			System.out.println("[JAST] Output is in file : "+outputFile.toString());
		}
		pr.destroy();
	}


	private void printOutputCommand() {
		System.out.println("[JAST] "+command+" output  : ");
		BufferedReader reader = new BufferedReader(new InputStreamReader(pr.getInputStream()));
		String line;
		try {
			while((line = reader.readLine()) != null) {
				System.out.println(line);
			}
		} catch (IOException e) {
			System.err.println("[JAST_Error] IO excpetion reader.");
			e.printStackTrace();
		}
		try {
			reader.close();
		} catch (IOException e) {
			System.err.println("[JAST_Error] IO excpetion closing reader.");
			e.printStackTrace();
		}
	}

	private void printError() {
		String msg = "";
		String line = "";
		BufferedReader reader = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
		try {
			while((line = reader.readLine()) != null) {
				msg +=line+"\n";	
			}
		} catch (IOException e) {
			System.err.println("[JAST_Error] IO excpetion reader.");
			e.printStackTrace();
		}
		try {
			reader.close();
		} catch (IOException e) {
			System.err.println("[JAST_Error] IO excpetion closing reader.");
			e.printStackTrace();
		}
		System.err.println("[JAST_Error] Failed to use the command : "+command+" :\n"+msg+"\n[JAST_Error] System will exit.");
		System.exit(1);

	}

	private void readConfig() throws IOException, ForbiddenOptionsException {
		BufferedReader br =  Files.newBufferedReader(config,StandardCharsets.UTF_8);
		String linetemp;
		while ((linetemp=br.readLine())!=null){
			if (checkOption(linetemp))
				totalCommand.add(linetemp);
			else 
				throw new ForbiddenOptionsException(linetemp, command);
		}
		br.close();
	}

	protected boolean checkOption(String optionToCheck) {
		if (forbiddenOptions==null)
			return true;
		for (String option : forbiddenOptions){
			if (optionToCheck.equals(option)){
				return false;
			}
		}
		return true;
	}
}
