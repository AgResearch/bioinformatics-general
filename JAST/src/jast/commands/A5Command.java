package jast.commands;

import jast.Command;
import java.nio.file.Path;
import java.nio.file.Paths;

public class A5Command extends Command{
	public A5Command(Path config,Path input,String outputroot,String [] arrayOfForbbidenOptions) {
		super(config,arrayOfForbbidenOptions);
		command="a5_pipeline.pl";
		totalCommand.add(input.toString()+"_1.fastq");
		totalCommand.add(input.toString()+"_2.fastq");
		totalCommand.add(outputroot+"_scaffoldsAfterA5");
		outputFile=Paths.get( outputroot+"_scaffoldsAfterA5.final.scaffolds.fasta");
	}
}
