package jast.commands;

import java.nio.file.Path;
import java.nio.file.Paths;

import jast.Command;

public class FlexbarCommand extends Command {

	public FlexbarCommand(Path config, Path readsOne, Path readsTwo,String [] arrayOfForbbidenOptions) {
		super(config,arrayOfForbbidenOptions);
		command="flexbar";
		String outputFileName = readsOne.toString().replace(".fastq", "");
		outputFileName = outputFileName.replace(".fq","");
		outputFileName = outputFileName.replace("_1","");
		outputFileName+="_afterFlexbar"; // results will be outputFileName_1.fastq and outputFileName_2.fastq
		outputFile = Paths.get(outputFileName);
		totalCommand.add("-t");
		totalCommand.add(outputFileName);
		totalCommand.add("-r");
		totalCommand.add(readsOne.toString());
		totalCommand.add("-p");
		totalCommand.add(readsTwo.toString());
	}
}
