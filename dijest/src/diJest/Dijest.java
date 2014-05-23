package diJest;


import java.io.IOException;
import java.util.Iterator;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;
/**
 * Main class of dijest project
 * @author Clément DELESTRE
 * @version 1.0
 *
 */
public class Dijest {
	public static final String appliName ="dijest";
	public static final float version=(float)1.0;
	public static final String author ="Clément DELESTRE";
	public static void main(String[] args){
		JSAP jsap = new JSAP();

		FileStringParser fsp =  FileStringParser.getParser(); // call the factory
		fsp.setMustExist(true); // the file must exist

		// option for display help
		Switch help = new Switch("help")
		.setShortFlag('h')
		.setLongFlag("help");
		help.setHelp("Print help and exit.");

		// option for input file
		FlaggedOption input = new FlaggedOption("Sequence File")
		.setStringParser(fsp) 
		.setRequired(true) 
		.setShortFlag('s') 
		.setLongFlag("sequence");
		input.setHelp("Input genome sequence (fasta format).");
		
		//option for enzyme file
		FlaggedOption enzyme = new FlaggedOption("Enzyme File")
		.setStringParser(fsp) 
		.setRequired(true)
		.setShortFlag('e')
		.setLongFlag("enz");
		enzyme.setHelp("File with several enzymes name separated by a line break.");
		
		//option for clean
		Switch clean = new Switch("clean")
		.setShortFlag('c') 
		.setLongFlag("clean");
		clean.setHelp("Delete output files created by "+appliName+" (not restrict output files).");

		//option for rscript
		FlaggedOption r = new FlaggedOption("rscript")
		.setStringParser(fsp) 
		.setShortFlag('r');
		r.setHelp("R script launch with output file containing fragments length as argument");

		//option for extract sequence
		Switch extract = new Switch("extract fragments sequences")
		.setShortFlag('x') 
		.setLongFlag("ext");
		extract.setHelp("Extract fragments sequences (thanks to extractseq command) after parsing results and record it in file in the following format : input_enzyme_restrict.out_extractFragmentsSequences.fasta\n/!\\ Warning : extractseq command doesn't support multi-FASTA format");

		//compute all options
		try {
			// mandatory option
			jsap.registerParameter(input);
			jsap.registerParameter(enzyme);
			// optionnal options
			jsap.registerParameter(clean);
			jsap.registerParameter(r);
			jsap.registerParameter(extract);
			// help
			jsap.registerParameter(help);
		} catch (JSAPException e) {
			System.err.println("ERROR ****  ");
			e.printStackTrace();
		}
		JSAPResult config = jsap.parse(args);    
		if (config.getBoolean("help"))
			displayUsage(jsap);
		if (config.success()) {
			LoopingRestrictInterface lr = new LoopingRestrict(config.getFile("Sequence File").toPath(),config.getFile("Enzyme File").toPath());
			ParsingResult pr = new ParsingRestrictResult(config.getFile("Sequence File").toPath(),config.getFile("Enzyme File").toPath());
			// compute the 2 main class
			try {
				lr.restrictRun();
				pr.parsingRun();
			} catch (IOException e) {
				System.err.println("ERROR WITH FILE");
				e.printStackTrace();
				System.err.println("*********\n\n System will exit *********\n\n");
				System.exit(1);
			}

			// compute R script
			if  (config.getFile("rscript")!=null){
				RloopInterface rl = new Rloop(config.getFile("rscript"),pr.getOutputFiles());
				rl.lauchScript();
			}
			
			// extract sequence from output restrict file
			if  (config.getBoolean("extract fragments sequences")){
				ExtractSequences es = new Extractor(lr.getOutputFiles(),config.getFile("Sequence File").toString());
				es.extractSequences();
			}

			// now clean
			if (config.getBoolean("clean")){
				DeleteI deletor = new Deletor(pr.getOutputFiles());
				deletor.deleteFiles();
			}

		}
		else {
			for (Iterator errs = config.getErrorMessageIterator();
					errs.hasNext();) {
				System.err.println("Error: " + errs.next());
			}
		}
	}

	public static void displayUsage(JSAP jsap){
		System.out.println(appliName+" version "+version+"\t"+author+"\n");
		System.out.println("Desctiption : ");
		System.out.println("Program that loop on restrict command (belonging to EMBOSS library) for several enzymes then parse the results to store only the fragments length.\n"+appliName+" can also compute a R script and/or extract fragments sequence.\nOutput files containing fragments length separated by a line break are in format input_Enzyme_Sequence.out and can be delete with -c option. Restrict output files are in format : input_enzyme_restrict.out.\n\n");
		System.out.println("Usage: java -jar "+appliName+".jar"); 
		System.out.println("\n\t\t" + jsap.getUsage()+"\n");
		System.out.println(jsap.getHelp());
		//System.out.println(booleanRules());
		System.exit(1);
	}

	private static String booleanRules() {
		return " \nThe following arguments are interpreted as TRUE: 1 t true y yes (case-insensitive)\nThe following arguments are interpreted as FALSE: 0 f false n no (case-insensitive)\n\n";
	}
}
