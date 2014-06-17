package jast.commands;

import jast.Command;

import java.nio.file.Path;
import java.nio.file.Paths;

public class BowtieMapCommand extends Command {

	public BowtieMapCommand (Path config, Path genome,Path indexFile,Path readOne,Path readTwo,String [] arrayOfForbbidenOptions){
		super(config,arrayOfForbbidenOptions);
		command="bowtie2";
		String outputName = genome.toString()+"_VSreads.sam";
		totalCommand.add("-x");
		totalCommand.add(indexFile.toString());
		totalCommand.add("-q");
		totalCommand.add("-1");
		totalCommand.add(readOne.toString());
		totalCommand.add("-2");
		totalCommand.add(readTwo.toString());
		totalCommand.add("-S");
		totalCommand.add(outputName);
		outputFile=Paths.get(outputName);
	}

}
