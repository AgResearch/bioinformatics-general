package jast.commands;

import java.nio.file.Path;
import java.nio.file.Paths;

import jast.Command;

public class ColombusCommand  extends Command {
	public ColombusCommand(Path config, Path contigs,Path reference,Path sam,Path flexbarOutput,String output,String [] arrayOfForbbidenOptions) {
		super(config,arrayOfForbbidenOptions);
		command="VelvetOptimiser.pl";
		totalCommand.add("-f");
		totalCommand.add("-reference "+reference+" -shortPaired -sam "+sam+" -separate -fastq "+flexbarOutput+"_1.fastq "+flexbarOutput+"_2.fastq -long -fasta "+contigs);
		totalCommand.add("-p");
		totalCommand.add("JASTColombusDir_"+output);
		outputFile=Paths.get( "JASTColombusDir_"+output+"_data_(\\d+)");
	}
}
