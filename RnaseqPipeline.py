#!/Users/vv39/anaconda2/bin/python #change path as necessary to run python

'''
This is a script that will help run rna-seq Analysis.

Dependencies:
	>Python 2.7
	For pre-processing:
		>FastQC
		>Trimmomatic
	For fast assmebly:
		>Kallisto
	For general assembly:
		>Bwa
		>eXpress
	For Novel Transcript Detection:
		>Bowtie2
		>TopHat
		>Cufflinks
'''

import sys,os
import argparse
import logging
import re
import string
import warnings
import string
import collections
import math
import sets
from time import strftime
import subprocess

from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

from qcmodule import SAM

from subprocess import call

#Global Variable:
post = "False"

def Quality_Assessment(read,kmer,logger,out,type):
	#Function to run Fastqc
	logger.info('Running FastQC')
	call(['fastqc','--kmers',str(kmer),'--extract',read,'--outdir',out+"/"+type+"FastqcMetrics/"])  
	
	#Parsing FastQC Results
	print "Basic Statistics:"
	print "--------------------------------"
	with open(out+"/"+type+"FastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/fastqc_data.txt",'r') as f:
		for line in f:
			lines = line.strip().split("\t")
			if "Sequence length" in line:
				print lines[0] +"\t" + lines[1]
			elif "%GC" in line:
				print lines[0] +"\t" + lines[1]
			elif "Sequences flagged as poor quality" in line:
				print lines[0] +"\t" + lines[1]
			elif "Total Sequences" in line:
				print lines[0] +"\t" + lines[1]
			elif ">>END_MODULE" in line:
				break

	print "\nOther Statistics Summary:"
	print "--------------------------------"
	with open(out+"/"+type+"FastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/summary.txt",'r') as f:	 	
		print "Metrics \t\t Result"
		for line in f:
			line = line.split('\t')
			print line[1]+"\t"+line[0]
			if "Per base sequence content" in line:
				if line[0] == "FAIL":
					print "\t\t\tThis could've failed due to fragmentation bias in the first 13bp. Checking..."
					with open(out+"/"+type+"FastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/fastqc_data.txt",'r') as f:
						for i, line in enumerate(f, 1):
							if "Per base sequence content" in line:
								start = i
							elif "Per sequence GC content" in line:
								end = i

					with open(out+"/"+type+"FastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/fastqc_data.txt",'r') as f:
						lines = f.readlines()
					i = start+1
					bp = 1
					check = 1
					while i < end-2:
						if bp > 13:
							values = lines[i].strip().split('\t')
							if abs(float(values[1]) - float(values[4])) > 20 or abs(float(values[3])-float(values[4])) > 20:
								print "\t\tThe error isn't just due to fragment bias. Please check the reads to remove possible adapter contamination"
								check = 0
								break
						i+=1
						bp+=1
					if check == 1:
						print "\t\t Only fragment bias in the first 13bp found."

def TrimmingPE(read1,read2,thread,phred,lead,trail,crop,minlen,window,qual,out,logger):
	#Function to trim reads using Trimmomatic
	logger.info('Trimming the reads')
	#call(["trimmomatic","PE","-trimlog","Trimmomatic.log","-threads",str(thread),"-"+phred,read1,read2,out+"/"+out+"_paired1.fq",out+"/"+out+"_paired2.fq",out+"/"+out+"_unpaired1.fq",out+"/"+out+"_unpaired2.fq","HEADCROP:"+str(crop),"LEADING:"+str(lead),"TRAILING:"+str(trail),"MINLEN:"+str(minlen),"SLIDINGWINDOW:"+str(window)+":"+str(qual)])
	call(["trimmomatic","PE","-trimlog","Trimmomatic.log","-threads",str(thread),"-"+phred,read1,read2,out+"/"+out+"_paired1.fq",out+"/"+out+"_paired2.fq",out+"/"+out+"_unpaired1.fq",out+"/"+out+"_unpaired2.fq","HEADCROP:"+str(crop),"LEADING:"+str(lead),"TRAILING:"+str(trail)])	

def khmer(read1,read2,kmer,logger):
	#function to normalize reads based on given kmer length
	logger.info("Normalizing the given reads...")
	call(["normalize-by-median.py","-k",kmer,read1,read2])

def Kallisto(read1,read2,file,kmer,boost,thread,logger,out):
	#Function to run Kallisto
	logger.info('Starting kallisto Indexing step...')
	call(["kallisto","index","-i","kallisto_indexfile.idx","-k",str(kmer),file])
	
	logger.info("Kallisto Indexing over. Starting Quantification...")
	call(["kallisto","quant","-t",str(thread),"-b",str(boost),"-i","kallisto_indexfile.idx","-o",out+"/Kallisto_Output",read1,read2])

def Bwa(read1,read2,file,algo,extra,logger,out):
	global post
	#Function to run bwa
	logger.info("Starting bwa indexing step...")
	call(["bwa","index",file])

	logger.info("Bwa Indexing over. Starting Quantification...")
	cmd = "bwa "+algo+" "+file+" "+read1+" "+read2+" > "+out+"/Bwa_output.sam"
	os.system(cmd)
	logger.info("Sorting sam file...")
	cmd = "sort -k 1 "+out+"/Bwa_output.sam > "+out+"/Bwa.hits.sam.sorted"
	os.system(cmd)

	cmd = "echo 'Post Mapping Metrics' > "+out+"/PostMappingMetrics.txt"
	os.system(cmd)
	PostMapping(out+"/Bwa.hits.sam.sorted",logger,out)
	if post=="True":
		PostQuality(out+"/Bwa.hits.sam.sorted",out,logger)
	eXpress(file, out+"/Bwa.hits.sam.sorted",extra,logger)

def eXpress(file,sam,extra,logger):
	#Function to run express analysis
	logger.info("Running analysis on output files...")

	cmd = "express "+file+" " +sam+" "+extra
	os.system(cmd)
	logger.info("Results are written into results.xprs")

def BowtieIndex(file,out,logger,libtype,bowalgo,read1,read2,thread,gtf,multi):
	#Function to build bowtie index
	logger.info("Starting Bowtie...")

	cmd = "bowtie2-build "+file+" "+file
	os.system(cmd)
	TopHat(read1,read2,file,out,libtype,bowalgo,thread,logger,gtf,multi)

def TopHat(read1,read2,bowref,out,libtype,bowalgo,thread,logger,gtf,multi,mask):
	#Function to run tophat
	logger.info("Running TopHat...")
	cmd = "tophat "+libtype+" "+bowalgo+" -o "+out+"/TopHat "+"-p "+str(thread)+" "+bowref+" "+read1+" "+read2
	os.system(cmd)
	cmd = "cat "+ out+"/TopHat/align_summary.txt > "+out+"/PostMappingMetrics.txt"
	os.system(cmd)
	cmd = "echo ' ' >> "+out+"/PostMappingMetrics.txt"
        os.system(cmd)
        cmd = "echo ' ' >> "+out+"/PostMappingMetrics.txt"
        os.system(cmd)
	PostMapping(out+"/TopHat/accepted_hits.bam",logger,out)
	if post =="True":
		PostQuality(out+"/TopHat/accepted_hits.bam",out,logger)
	Cufflinks(out+"/TopHat/accepted_hits.bam",gtf,thread,libtype,multi,out,logger)

def Cufflinks(file,gtf,thread,libtype,multi,out,logger,mask):
	logger.info("Running Cufflinks...")
	if gtf == "Na":
		logger.warning("Running cufflinks without GTF file! This is not recommended")
		gtf = ''
	else:
		gtf = "--GTF-guide "+gtf

	if mask == "Na":
                logger.warning("Running cufflinks without mask file! This is NOT recommended")
                mask = ''
        else:
                mask = " --mask-file  "+mask

	call(["cufflinks",file,"-o",out+"/Cufflinks","-u",multi,"--library-type",libtype,"-p",str(thread),gtf,mask])
	logger.info("Outputs can be accessed from "+out+"/Cufflinks/isoforms.fpkm_tracking")

def PostMapping(file,logger,out):
	logger.info("Doing Post Mapping analysis...")
	cmd = "echo ' ' >> "+out+"/PostMappingMetrics.txt"
	os.system(cmd)
	cmd = "echo ' ' >> "+out+"/PostMappingMetrics.txt"
        os.system(cmd)
	cmd = "samtools flagstat "+file+" >> "+out+"/PostMappingMetrics.txt"
	os.system(cmd)
	logger.info("Post Mapping Metrics can be accessed from "+out+"/PostMappingMetrics.txt")

def PostQuality(file,out,logger):
	out = out+"/"+out+"_"
	#Function to reimplement Read_quality.py and read_duplication.py from RSeQC-2.6.3 package written by Liguo Wang
	logger.info("Plotting post quality graphs...")
	if (os.path.exists(file)):
                obj = SAM.ParseBAM(file)
                obj.readsQual_boxplot(outfile=out, q_cut = 30, shrink = 1000)
                try:
                        subprocess.call("Rscript " + out+ ".qual.r",shell=True)
                except:
                        pass

		obj.readDupRate(outfile=out,up_bound=500, q_cut = 30)
		try:
                        subprocess.call("Rscript " + out +  ".DupRate_plot.r", shell=True)
                except:
                        pass
        else:
                logger.error("Failed to plot graphs")
                sys.exit(0)	

def main():
	global post
	parser = argparse.ArgumentParser(prog='rnaseq',description="Program to help run RNA-seq Analysis")
	required = parser.add_argument_group("Required Options")
	required.add_argument('-t','--type',dest='type',help="Enter the type of analysis you want to do (F=Fast, S=Slow, N=Novel) [default=F]",choices=['F','S','N'],required=True, default="F")
	required.add_argument('-r1','--read1',dest='read1',help="Enter the location of read 1", required=True,metavar="\b",default="Na")
	required.add_argument('-r2','--read2',dest='read2',help="Enter the location of read 2", required=True,metavar="\b",default="Na")
	required.add_argument('-o','--output',dest="out",help="Enter the prefix to your output files", required=True, metavar="\b",default="out")	

	fastqc = parser.add_argument_group("FastQC options")
	fastqc.add_argument('-fk','--fastqck',dest='fastqck',help="Set the kmer size for running fastqc [Default=7]",default=7,type=int,metavar="kmer-size(int)")
	
	trimmomatic = parser.add_argument_group("Trimmomatic options")
	parser.add_argument('--threads',dest="thread",help="Set the number of threads to use for trimmomatic [Default=4]",default=4,type=int,metavar="number(int)")
	trimmomatic.add_argument('--phred33|phred64',dest="phred",default="phred64",choices=["phred33","phred64"])
	trimmomatic.add_argument('--leading',dest="trimlead",help="Specify the minimum quality to keep a base towards the start of the sequence [Default = 20]", type=int, default=20,metavar="quality(int)")
	trimmomatic.add_argument('--trailing',dest="trimtrail",help="Specify the minimum quality to keep a base towards the end of the sequence [Default = 20]",type=int,default=20,metavar="quality(int)")
	trimmomatic.add_argument('--headcrop',dest="trimcrop",help="Specify the number of bases to remove, from the start of the read [Default = 0]",type=int,default=0,metavar="length(int)")	
	trimmomatic.add_argument('--minlen',dest="trimlen",help="Specify the minimum lengths of read to keep [Default = 0]",type=int,default=0,metavar="length(int)")
	trimmomatic.add_argument('--window',dest="trimwindow",help="Specify the number of bases to average across for sliding window [Default=0]",type=int,default=0,metavar="size(int)")
	trimmomatic.add_argument('--quality',dest="trimq",help="Specify the average quality required [Default = 0]", type=int, default=0, metavar="quality(int)")

	kallisto = parser.add_argument_group("Fast Alginment Option - Kallisto")
	kallisto.add_argument('--kindex',dest='kindex',help="Enter the fasta file to be used to construct index",metavar="<path to file>",default="Na")
	kallisto.add_argument('--kkmer',dest='kkmer',help="Specify the kmer length to be used [Default = 31]", type=int, default=31,metavar="kmer-size(int)")
	kallisto.add_argument('--kboost',dest='kboost',help="Specify the number of bootstraps to perform [Default = 0]", type=int, default=0, metavar="num(int)")	

	slow = parser.add_argument_group("Slow Allignment Option - BWA + eXpress")
	slow.add_argument('--bindex',dest='bindex',help="Enter the fasta file to be used to construct index", metavar="<path to file>", default="Na")
	slow.add_argument('--algo',dest='algo',help="Specify the algorithm to be used for algorithm [Default = mem]", choices=["mem","bwasw","sampe"], default="mem") 
	slow.add_argument('--strand',dest='strand',help="Specify Strand option to run eXpress", choices=['fr-stranded','rf-stranded','f-stranded','r-stranded'], default="None")	

	novel = parser.add_argument_group("Novel Transcription Detection - TopHat + Cufflinks")
	novel.add_argument('--bowref',dest="bowref",help="Enter the file to be used to construct bowtie index", metavar="<path to file>", default="Na")
	novel.add_argument('--bowindex',dest="bowindex", help="Enter the bowtie index", metavar="<path to file>", default="Na")
	novel.add_argument('--libtype',dest="libtype",help="Specify the library type", choices=['fr-unstranded','fr-firststrand','fr-secondstrand','Na'], default="Na")
	novel.add_argument('--bowmode',dest="bowalgo",help="Specify the bowtie mode", choices=['b2-very-fast','b2-fast','b2-sensitive','b2-very-sensitive'], default="Na")		
	novel.add_argument('-g','--gtf',dest="gtf",help="Enter the GTF file to be used for annotations",metavar="<path to file>",default="Na")	
	novel.add_argument('-u','--multi-read-correct',dest="multi",help="Use 'rescue method' for multi-reads(more accurate)",choices=['True','False'],default="False")
	novel.add_argument('-m','--mask-file',dest="mask",help="Ignore all alignment within transcripts in this file",metavar="<path to file>",default="Na")

	post = parser.add_argument_group("Post Mapping Metrics Graph")
	post.add_argument('-p','--post',dest="post",help="Use this option if you want post mapping metrics graph [Default=False]", action="store_const",const="True",default="False")
	

	parser.add_argument("--skippre",dest="skippre",help="Use this option if you want to skip pre processing steps [Default=False]",action="store_const",const="True",default="False")
	parser.add_argument("--normalize", dest="norm",help="Use this option if you want to normalize the input data [Default=False]", action="store_const",const="True",default="False")
	parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	post = args.post

	logger = logging.getLogger('Main')
	logger.setLevel(logging.DEBUG)
        stream = logging.StreamHandler()
        stream.setLevel(logging.DEBUG)
        formats = logging.Formatter('[%(asctime)s %(name)s %(levelname)s]: %(message)s')
        stream.setFormatter(formats)
        logger.addHandler(stream)

	#Check the existance of files
	if(args.read1 == 'Na' or not os.path.exists(args.read1)):
		logger.error("The read 1 file is not readable")
		parser.print_help()
		sys.exit()

	if(args.read2 == 'Na' or not os.path.exists(args.read2)):
                logger.error("The read 2 file is not readable")
                parser.print_help()
                sys.exit()

	if(args.type == "F"):
		if (args.kindex == 'Na' or not os.path.exists(args.kindex)):
			logger.error("The fasta file for building index is not readable [--kindex]")
			parser.print_help()
			sys.exit()
	
	if(args.type == "S"):
                if (args.bindex == 'Na' or not os.path.exists(args.bindex)):
                        logger.error("The fasta file for building index is not readable [--bindex]")
                        parser.print_help()
                        sys.exit()

	if(args.type == "N"):
		if (args.bowindex == "Na" and args.bowref == "Na"):
			logger.error("Please enter either a reference file to build bowtie reference [--bowref] or a bowtie2 index (.ebwt) [--bowindex]")
			sys.exit()

                elif (args.bowindex == 'Na' and not os.path.exists(args.bowref)):
                        logger.error("The file for building index is not readable")
                        sys.exit()

		elif (args.bowref == 'Na' and not os.path.exists(args.bowindex+".1.bt2")):
			logger.error("Trouble reading bowtie files (.1.bt2)")
			sys.exit()
		
		elif (args.bowref == 'Na' and not os.path.exists(args.bowindex+".rev.1.bt2")):
                        logger.error("Trouble reading bowtie reference files (.rev.1.bt2)")
                        sys.exit()

		elif (args.bowref == 'Na' and not os.path.exists(args.bowindex+".fa")):
                        logger.error("Trouble reading bowtie reference file (.fa)")
                        sys.exit()

		if (args.gtf == "Na"):
			logger.warning("GTF file for annotation not specified. Might lead to inaccurate analysis.")
		else:
			if(not os.path.exists(args.gtf)):
				logger.error("GTF file not readable")
				sys.exit()

		if (args.mask == "Na"):
                        logger.warning("Mask file for ignoring rRNA contamination not specified. This is a recommended option.")
                else:
                        if(not os.path.exists(args.mask)):
                                logger.error("Mask file not readable")
                                sys.exit()

	read1 = args.read1
	read2 = args.read2
	
	#Creating Directory Structure
	call(['mkdir',args.out])

	#Pre Processing Data
	if (args.skippre=="False"):
		log = logging.getLogger('Pre Processing')
        	log.setLevel(logging.DEBUG)
        	log.addHandler(stream)
		
		call(['mkdir',args.out+'/preFastqcMetrics'])
		Quality_Assessment(args.read1,args.fastqck,log,args.out,"pre")
		Quality_Assessment(args.read2,args.fastqck,log,args.out,"pre")

		TrimmingPE(args.read1,args.read2,args.thread,args.phred,args.trimlead,args.trimtrail,args.trimcrop,args.trimlen,args.trimwindow,args.trimq,args.out,logger)
		read1 = args.out+"/"+args.out+"_paired1.fq"
		read2 = args.out+"/"+args.out+"_paired2.fq"
	
		call(['mkdir',args.out+'/postFastqcMetrics'])
                Quality_Assessment(args.read1,args.fastqck,log,args.out,"post")
                Quality_Assessment(args.read2,args.fastqck,log,args.out,"post)
	else:
		logger.info("User opted to skip Pre Processing step")

	if (args.norm=="True"):
		khmer(read1,read2,args.fastqck,logger)
		read1 = read1+".keep"
		read2 = read2+".keep"

	#Running RNA-seq according to what the user selected:
	log = logging.getLogger('RNA Seq')
	log.setLevel(logging.DEBUG)
	log.addHandler(stream)

	if args.type=="F":
		log.info("User opted for fast analysis, Prepping for Kallisto...")
		Kallisto(read1,read2,args.kindex,args.kkmer,args.kboost,args.thread,logger,args.out)

	elif args.type=="S":
		log.info("User opted for slow analysis, Prepping for bwa and eXpress...")
		extra = ""
		if args.strand != "None":
			extra = "--"+args.strand
		Bwa(read1,read2,args.bindex,args.algo,extra,logger,args.out)
	else:
		log.info("User opted to perform analysis to detect novel transcripts, Prepping for Tophat, Cufflinks and eXpress...")	
		libtype = ""
		bowalgo = ""
                if args.libtype != "Na":
                	libtype = "--library-type "+args.libtype
                if args.bowalgo != "Na":
                	bowalgo = "--"+args.bowalgo

		if args.bowref != "Na":
			log.info("Building Bowtie2 index before proceeding with TopHat")
			BowtieIndex(args.bowref,args.out,logger,libtype,bowalgo,read1,read2,args.thread,args.gtf,args.multi)
		else:
			print args.bowindex
			TopHat(read1,read2,args.bowindex,args.out,libtype,bowalgo,args.thread,logger,args.gtf,args.multi,args.mask)

if __name__ == "__main__":
	main()
