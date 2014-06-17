package jast;


import jast.commands.A5Command;
import jast.commands.BowtieIndexCommand;
import jast.commands.BowtieMapCommand;
import jast.commands.ColombusCommand;
import jast.commands.FlexbarCommand;
import jast.commands.SSPACEcommand;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;

/**
 * JAST (Java Assembling and Scaffolding Tool) main class
 * @author Clément DELESTRE
 * @version 1.0
 * @since 1.0
 */
public class JASTmain {

	public static final String appliName ="JAST (Java Assembling and Scaffolding Tool)";
	public static final float version=(float)1.0;
	public static final String author ="Clément DELESTRE";
	/**
	 * The SSPACE's file created will finish by this extension
	 */
	public static final String libSSPACEext="_JAST";
	/**
	 * Forbidden options for SSPACE command
	 */
	private static String[] forbiddenSSPACE={"-l","-s","-b"};
	/**
	 * Forbidden options for Colombus command
	 */
	private static String[] forbiddenColombus={"-f","-p"};
	/**
	 * Forbidden options for Flexbar command
	 */
	private static String[] forbiddenFlexbar={"-t","-r","-p"};
	/**
	 * Forbidden options for A5 command
	 */
	private static String[] forbiddenA5=null;
	/**
	 * Forbidden options for Bowtie build command
	 */
	private static String[] forbiddenBowtieIndex=null;
	/**
	 * Forbidden options for Bowtie-x command
	 */
	private static String[] forbiddenBowtieMap={"-x","-q","-1","-2","-S"};


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		JSAP jsap = new JSAP();
		FileStringParser fsp =  FileStringParser.getParser(); // call the factory
		fsp.setMustExist(true); // the file must exist

		// option for display help
		Switch help = new Switch("help")
		.setShortFlag('h')
		.setLongFlag("help");
		help.setHelp("Print help and exit.");

		//first input file
		FlaggedOption firstFile = new FlaggedOption("First input reads")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setShortFlag('r')
		.setLongFlag("reads");
		firstFile.setHelp("First input file of paired reads (fastq format)");

		//second input file
		FlaggedOption secondFile = new FlaggedOption("Second input reads")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setShortFlag('p')
		.setLongFlag("paired");
		secondFile.setHelp("Second input file of paired reads (fastq format)");

		//final output file
		FlaggedOption output = new FlaggedOption("Output file")
		.setStringParser(JSAP.STRING_PARSER) 
		.setRequired(true) 
		.setShortFlag('o')
		.setLongFlag("output");
		output.setHelp("Final output name");

		//Reference
		FlaggedOption reference = new FlaggedOption("Reference file")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setLongFlag("ref");
		reference.setHelp("Referencec sequence (fasta format)");

	
		//Bowtie File Build
		FlaggedOption bowtieConfigBuild = new FlaggedOption("Bowtie config Build")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setLongFlag("bowbuild");
		bowtieConfigBuild.setHelp("Config file use for bowtie-build command.");

		//Bowtie File Map
		FlaggedOption bowtieConfigMap = new FlaggedOption("Bowtie config Map")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setLongFlag("bowmap");
		bowtieConfigMap.setHelp("Config file use for bowtie-x command. All mandatory options must be specified except : "+printForbiddenOptions(forbiddenBowtieMap));


		//SSPACE library file
		FlaggedOption sspacelib = new FlaggedOption("SSPACE library")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setShortFlag('l')
		.setLongFlag("sspacelib");
		sspacelib.setHelp("SSPACE library file (see SSPACE doc, equivalent to file you specify with '-l' option).\nWarning : first line must contains only the 3 last columns : the other information (name and fastq file) will be writen in a file nammed with _JAST extension and will be nammed JASTlib and the flexbar results files will be added.");

		// Flexbar
		FlaggedOption flexbarFile = new FlaggedOption("flexbar config file")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setLongFlag("flexbar");
		flexbarFile.setHelp("Flexbar config file. All mandatory options must be specified except : "+printForbiddenOptions(forbiddenFlexbar));

		// Colombus
		FlaggedOption colombusFile = new FlaggedOption("Colombus config file")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setLongFlag("colombus");
		colombusFile.setHelp("Colombus config file. All mandatory options must be specified except : "+printForbiddenOptions(forbiddenColombus));

		// SSPACE
		FlaggedOption sspaceFile = new FlaggedOption("SSPACE config file")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setLongFlag("sspace");
		sspaceFile.setHelp("SSPACE config file. All mandatory options must be specified except : "+printForbiddenOptions(forbiddenSSPACE));


		try {
			jsap.registerParameter(help);
			jsap.registerParameter(firstFile);
			jsap.registerParameter(secondFile);
			jsap.registerParameter(output);
			jsap.registerParameter(bowtieConfigBuild);
			jsap.registerParameter(bowtieConfigMap);
			jsap.registerParameter(flexbarFile);
			jsap.registerParameter(colombusFile);
			jsap.registerParameter(reference);
			jsap.registerParameter(sspaceFile);
			jsap.registerParameter(sspacelib);
		} catch (JSAPException e) {
			System.err.println("[JAST_Error] with JSAP Parameter :");
			e.printStackTrace();
		}
		JSAPResult config = jsap.parse(args);
		if (config.getBoolean("help"))
			displayUsage(jsap);
		System.out.println("[JAST] Starting "+appliName+"...");
		if (config.success()) {
			FlexbarCommand flexbar = new FlexbarCommand(config.getFile("flexbar config file").toPath(),config.getFile("First input reads").toPath(),config.getFile("Second input reads").toPath(),forbiddenFlexbar);
			flexbar.exec();
			A5Command a5 = new A5Command(null,flexbar.getOutputFile(),config.getString("Output file"),forbiddenA5);
			a5.exec();
			//We must create a new file for SSPACE
			File NewSSPACElib = new File("");
			try {
				NewSSPACElib=computeSSPACELibrairie(config.getFile("SSPACE library").toPath(), flexbar.getOutputFile());
			} catch (IOException e) {
				System.err.println("[JAST_Error] I/O excpetion with file "+config.getFile("SSPACE library").getPath()+" and/or "+flexbar.getOutputFile());
				e.printStackTrace();
				System.err.println("[JAST_Error] System will exit.");
				System.exit(1);
			}



			
			BowtieIndexCommand bowtieIndex=new BowtieIndexCommand(config.getFile("Bowtie config Build").toPath(), config.getFile("Reference file").toPath(),forbiddenBowtieIndex);
			bowtieIndex.exec();
			BowtieMapCommand bowtieMap=new BowtieMapCommand(config.getFile("Bowtie config Map").toPath(), config.getFile("Reference file").toPath(), bowtieIndex.getOutputFile() , config.getFile("First input reads").toPath(), config.getFile("Second input reads").toPath(),forbiddenBowtieMap);
			bowtieMap.exec();
			ColombusCommand colombus = new ColombusCommand(config.getFile("Colombus config file").toPath(), a5.getOutputFile(), config.getFile("Reference file").toPath(), bowtieMap.getOutputFile(), flexbar.getOutputFile(),config.getString("Output file"),forbiddenColombus);

			colombus.exec();
			File[] files=new File(".").listFiles(new Filter(colombus.getOutputFile().toString()));
			File contigsAfterColombus=new File("");
			if (files.length!=1){
				System.err.println("[JAST_Error] Only one directory should match with this pattern : "+colombus.getOutputFile()+" (instead of "+files.length+"))\nIf directories exist beacause of early use, please either change the output name or delete the directory.\nSystem will exit");
				System.exit(1);
			}
			else {
				contigsAfterColombus = new File(files[0]+"/contigs.fa");
			}

			SSPACEcommand sspaceFinal = new SSPACEcommand(config.getFile("SSPACE config file").toPath(), contigsAfterColombus.toPath(), 2,config.getString("Output file"),NewSSPACElib.toPath(),forbiddenSSPACE);
			sspaceFinal.exec();
			System.out.println("[JAST] Finsih ! output file : "+sspaceFinal.getOutputFile());

		}
		else {
			for (Iterator errs = config.getErrorMessageIterator();
					errs.hasNext();) {
				System.err.println("[JAST_Error] : " + errs.next());
			}
		}

	}

/**
 * Get string containing forbidden options. 
 * @param forbiddenOptions
 * @return forbiddenOptions
 */
	private static String printForbiddenOptions(String[] forbiddenOptions) {
		String toReturn="\n";
		for (String s : forbiddenOptions){
			toReturn+=s+"\n";
		}
		return toReturn;
	}


	public static void displayUsage(JSAP jsap){
		System.out.println(appliName+" version "+version+"\t"+author+"\n");
		System.out.println("Desctiption : ");
		System.out.println("Program that do assembling and scaffolding from paired-end Illumina files.\n"+appliName+" use following softwares :\n\t*Flexbar\n\t*A5\n\t*Bowtie\n\t*Colombus (belonging to Velvet)\n\t*SSPACE\nPlease check they all are installed on your computer\n\n");
		System.out.println("Usage: java -jar jast.jar"); 
		System.out.println("\n\t\t" + jsap.getUsage()+"\n");
		System.out.println(jsap.getHelp());
		//System.out.println(booleanRules());
		System.exit(1);
	}
	/**
	 * Could be usefull to inform the user of what is "false" and what is "true".
	 * @return rules
	 */
	private static String booleanRules() {
		return " \nThe following arguments are interpreted as TRUE: 1 t true y yes (case-insensitive)\nThe following arguments are interpreted as FALSE: 0 f false n no (case-insensitive)\n\n";
	}

	private static File computeSSPACELibrairie(Path library,Path flexbaRoot) throws IOException {
		BufferedReader br =  Files.newBufferedReader(library,StandardCharsets.UTF_8);
		File newLib = new File(library+libSSPACEext);
		BufferedWriter writer = Files.newBufferedWriter(newLib.toPath(), StandardCharsets.UTF_8);
		String linetemp;
		boolean first=true;
		while ((linetemp=br.readLine())!=null){
			if (first){ // We only change the first line
				String newLine = "JASTlib "+flexbaRoot+"_1.fastq " +flexbaRoot+"_2.fastq "+linetemp;
				writer.write(newLine);
				first=false;
			}
			else {
				writer.write("\n"+linetemp);
			}
		}
		br.close();
		writer.close();
		return newLib;
	}
}
