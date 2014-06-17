package jast.commands;

import jast.Command;

import java.nio.file.Path;
import java.nio.file.Paths;

public class BowtieIndexCommand  extends Command  {

	public BowtieIndexCommand (Path config, Path ref,String[] arrayOfForbbidenOptions){
		super(config,arrayOfForbbidenOptions);
		String outputName = ref.toString()+"_BowtieIndex";
		command="bowtie2-build";
		totalCommand.add(ref.toString());
		totalCommand.add(outputName);
		outputFile=Paths.get(outputName);
	}

}
