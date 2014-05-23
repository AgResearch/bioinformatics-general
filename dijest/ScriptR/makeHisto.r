args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0){
	print("makeHisto should have one input file with a name in format : Organism_Enzyme_Sequence.out")
	print("Example : Rscript makeHisto.r EColi_ApekI_s128.out")
}else  {

	nameInputFile<-args[1] #Take the input file
	#Input file name is in format :  OrganismInput_Enzyme_Sequence.out 
	name<-gsub("\\.*.out$","", nameInputFile) 
	vectorName<-strsplit(name,"_")[[1]]  
	size<-length(vectorName)
	organismInput<-vectorName[size-2]
	enzymeName<-vectorName[size-1]
	seqName<-vectorName[size]
	#Load data
	data<-read.table(nameInputFile,header=FALSE,sep="\n")
	realdata<-c(data[,1])

#Compute percentage
	
	oneH<-0
	twoH<-0
	threeH<-0
	fourH<-0
	fiveH<-0
	sixH<-0
	sevenH<-0
	eightH<-0
	nineH<-0
	thousand<-0
	moreThousand<-0

	
	 for (n in 1:length(realdata)){ 
	 	if (realdata[n]<100){ 
	 		oneH<-oneH+1 
	 	} else if (realdata[n]<200){ 
	 		twoH<-twoH+1
	 	}
	 	else if (realdata[n]<300){ 
	 		threeH<-threeH+1
	 	}
	 	else if (realdata[n]<400){ 
	 		fourH<-fourH+1
	 	}
	 
	 	else if (realdata[n]<500){ 
	 		fiveH<-fiveH+1
	 	}
	 	else if (realdata[n]<600){ 
	 		sixH<-sixH+1
	 	}
	 	else if (realdata[n]<700){ 
	 		sevenH<-sevenH+1
	 	}
	 	else if (realdata[n]<800){ 
	 		eightH<-eightH+1
	 	}
	 	else if (realdata[n]<900){ 
	 		nineH<-nineH+1
	 	}
	 	else if (realdata[n]<=1000){ 
	 		thousand<-thousand+1
	 	}
	 	else{
	 		moreThousand<-moreThousand+1
	 	}
	 	
	 	
	}
	#abs=c(oneH,twoH,threeH,fourH,fiveH,sixH,sevenH,eightH,nineH,thousand,moreThousand)

	totalsize<-length(realdata)
	percent=c(oneH/totalsize*100,twoH/totalsize*100,threeH/totalsize*100,fourH/totalsize*100,fiveH/totalsize*100,sixH/totalsize*100,sevenH/totalsize*100,eightH/totalsize*100,nineH/totalsize*100,thousand/totalsize*100,moreThousand/totalsize*100)

	#Just Checking
	#print(paste(percent))	

	#Name and title for PNG
	title<-paste("Percentage of fragments between 0 and 1000pb for sequence ",seqName," digest by ",enzymeName,sep="")
	namePNG<-paste(name,".png",sep="")
	png(namePNG,800,1200)
	
	#Do the picture
	barplot(percent,names=c("0-100","100-200","200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000",">1000"),xlab="Fragments size",ylab="% of total fragments",main=title)
	graphics.off()
	#Advertising user
	print(paste("Done. Picture saved on ",namePNG))	
}

