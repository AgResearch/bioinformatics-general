package jast.commands;

import jast.Command;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class SSPACEcommand extends Command {
	public SSPACEcommand(Path config, Path input,int numTime,String output,Path lib,String [] arrayOfForbbidenOptions) {
		super(config,arrayOfForbbidenOptions);
		String sspaceExtenstion=".final.scaffolds.fasta";
		command="SSPACE_Basic_v2.0.pl";
		totalCommand.add("-l");
		totalCommand.add(lib.toString());
		totalCommand.add("-s");
		totalCommand.add(input.toString());
		totalCommand.add("-b");
		if (output==null){
			totalCommand.add("SSPACE"+numTime+"stTimeJAST");
			outputFile=Paths.get( "SSPACE"+numTime+"stTimeJAST"+sspaceExtenstion);
		}
		else{
			totalCommand.add(output);
			outputFile=Paths.get( output+sspaceExtenstion);
		}
	}

}

