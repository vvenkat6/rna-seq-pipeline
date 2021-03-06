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
		>samtools sort
		>eXpress

	For Novel Transcript Detection:
		>Bowtie2
		>TopHat
		>Cufflinks

	For report generation:
		>pylatex
		>pdflatex (or maclatex)
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
import glob
from time import strftime
import subprocess
from subprocess import call
import numpy as np
from cStringIO import StringIO

from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
import pysam

#modules from RSEQC
from qcmodule import SAM
from qcmodule import BED
from qcmodule import bam_cigar

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Package, Matrix, Command, Itemize
from pylatex.utils import italic, NoEscape

#Global Variable:
post = "False"

def DeterminePhred(read,out,type,logger):
	logger.info("Determining read quality encoding...")
	with open(out+"/"+type+"FastqcMetrics/"+"triM"+read.split('/')[-1].split('.')[0]+"_fastqc/fastqc_data.txt",'r') as f:
		for line in f:
                        lines = line.strip().split("\t")
			if "Encoding" in line:
				print lines[1].split(" ")[-1]
				encoding = float(lines[1].split(" ")[-1])
			elif ">>END_MODULE" in line:
                                break
		
	if encoding >= 1.8:
		return("--phred33")
	else:
		return("--phred64")
	

def Quality_Assessment(read,kmer,logger,out,type):
	#Function to run Fastqc
	logger.info('Running FastQC')
	call(['fastqc','--kmers',str(kmer),'--extract',read,'--outdir',out+"/"+type+"FastqcMetrics/"])
	output = out+"/"+read.split('/')[-1].split('.')[0]+"_"+type+"process_Summary.txt"
	fout = open(output,'w')  
	
	#Parsing FastQC Results
	fout.write("\nBasic Statistics:\n")
	fout.write("--------------------------------\n\n")
	with open(out+"/"+type+"FastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/fastqc_data.txt",'r') as f:
		for line in f:
			lines = line.strip().split("\t")
			if "Sequence length" in line:
				fout.write(lines[0] +"\t" + lines[1]+"\n")
			elif "%GC" in line:
				fout.write(lines[0] +"\t" + lines[1]+"\n")
			elif "Sequences flagged as poor quality" in line:
				fout.write(lines[0] +"\t" + lines[1]+"\n")
			elif "Total Sequences" in line:
				fout.write(lines[0] +"\t" + lines[1]+"\n")
			elif ">>END_MODULE" in line:
				break

	fout.write("\nOther Statistics Summary:\n")
	fout.write("--------------------------------\n")
	with open(out+"/"+type+"FastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/summary.txt",'r') as f:	 	
		fout.write("Metrics \t\t Result")
		for line in f:
			line = line.split('\t')
			fout.write(line[1]+"\t"+line[0]+"\n")
			if "Per base sequence content" in line:
				if line[0] == "FAIL":
					fout.write("\t\t\tThis could've failed due to fragmentation bias in the first 13bp. Checking...\n")
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
								fout.write("\t\tThe error isn't just due to fragment bias. Please check the reads to remove possible adapter contamination\n")
								check = 0
								break
						i+=1
						bp+=1
					if check == 1:
						fout.write("\t\t Only fragment bias in the first 13bp found.")
	fout.close()

def TrimmingPE(read1,read2,thread,lead,trail,crop,minlen,window,qual,out,adapter,logger):
	#Function to trim paired end reads using Trimmomatic
	logger.info('Trimming the reads')
	call(["trimmomatic","PE","-trimlog","Trimmomatic.log","-threads",str(thread),read1,read2,out+"/triM"+read1.split('/')[-1],out+"/"+out+"_unpaired1.fq.gz",out+"/triM"+read2.split('/')[-1],out+"/"+out+"_unpaired2.fq.gz","ILLUMINACLIP:"+adapter+":2:30:10","HEADCROP:"+str(crop),"LEADING:"+str(lead),"TRAILING:"+str(trail),"MINLEN:"+str(minlen),"SLIDINGWINDOW:"+str(window)+":"+str(qual)])	

def TrimmingSE(read,thread,lead,trail,crop,minlen,window,qual,out,adapter,logger):
        #Function to trim single ended using Trimmomatic
        logger.info('Trimming the reads')
        call(["trimmomatic","SE","-trimlog","Trimmomatic.log","-threads",str(thread),read,out+"/triM"+read.split('/')[-1],"ILLUMINACLIP:"+adapter+":2:30:10","HEADCROP:"+str(crop),"LEADING:"+str(lead),"TRAILING:"+str(trail),"MINLEN:"+str(minlen),"SLIDINGWINDOW:"+str(window)+":"+str(qual)]) 

def khmer(read1,read2,kmer,logger):
	#function to normalize reads based on given kmer length
	logger.info("Normalizing the given reads...")
	call(["normalize-by-median.py","-k",kmer,read1,read2])

def Kallisto(read1,read2,read,index,ref,kmer,boost,thread,logger,out,length,sd):
	#Function to run Kallisto
	if ref != 'Na':
		logger.info('Starting kallisto Indexing step...')
		call(["kallisto","index","-i",out+"/kallisto_indexfile.idx","-k",str(kmer),ref])
		index = out+"/kallisto_indexfile.idx"
	
	logger.info("Starting Kallisto Quantification...")
	if read == "Na":
		call(["kallisto","quant","-t",str(thread),"-b",str(boost),"-i",index,"-o",out+"/Kallisto_Output",read1,read2])
	else:
		call(["kallisto","quant","-t",str(thread),"-b",str(boost),"-i",index,"-o",out+"/Kallisto_Output","--single",read,"-l",str(length),"-s",str(sd)])

def Bwa(read1,read2,read,index,ref,algo,extra,logger,out,thread,ribo):
	global post
	#Function to run bwa
	if ref != 'Na':
		logger.info("Starting bwa indexing step...")
		call(["bwa","index",ref])
		index = ref

	logger.info("Starting BWA Quantification...")
	if read == 'Na':
		cmd = "bwa "+algo+" "+"-t "+str(thread)+" "+index+" "+read1+" "+read2+" > "+out+"/Bwa_output.sam"
	else:
		cmd = "bwa "+algo+" "+"-t "+str(thread)+" "+index+" "+read+" > "+out+"/Bwa_output.sam"
	os.system(cmd)

	logger.info("Sorting sam file...")
	cmd = "samtools sort -n --threads "+str(thread)+" "+out+"/Bwa_output.sam > "+out+"/Bwa.hits.sam.sorted"
	os.system(cmd)

	call(["rm",out+"/Bwa_output.sam"])

	cmd = "echo 'Post Mapping Metrics' > "+out+"/PostMappingMetrics.txt"
	os.system(cmd)

	if ribo != 'Na':
		RibosomalStat(out+"/Bwa.hits.sam.sorted",ribo,out,logger)
		cmd = "cat "+out+"/RibosomalStat.txt > "+out+"/PostMappingMetrics.txt"
                os.system(cmd)

	#PostMapping(out+"/Bwa.hits.sam.sorted",logger,out)
	if post=="True":
		print ""
		PostQuality(out+"/Bwa.hits.sam.sorted",out,logger)
	eXpress(index, out+"/Bwa.hits.sam.sorted",extra,logger,out)

def eXpress(file,sam,extra,logger,out):
	#Function to run express analysis
	logger.info("Running analysis on output files...")

	cmd = "express "+file+" " +sam+" -o "+out+"/express "+extra
	os.system(cmd)
	logger.info("Results are written into results.xprs")

def BowtieIndex(file,out,logger,libtype,bowalgo,read1,read2,read,thread,gtf,multi,mask,ribo):
	#Function to build bowtie index
	logger.info("Starting Bowtie indexing...")

	cmd = "bowtie2-build "+file+" "+file
	os.system(cmd)
	TopHat(read1,read2,read,file,out,libtype,bowalgo,thread,logger,gtf,multi,mask,ribo)

def BowtieAlign(extra,read1,read2,read,index,ref,logger,out,thread,ribo,phred):
	#Function to run bowtie2
	logger.info("Starting Bowtie2...")

	if read == 'Na':
		call(["bowtie2",phred,"-p",str(thread),"-x",index,"-1",read1,"-2",read2,"-S", out+"/Bowtie2.sam","--very-sensitive-local","--met-file",out+"/Bowtie2.metrics"])
	else:
		call(["bowtie2",phred,"-p",str(thread),"-x",index,"-U",read,"-S", out+"/Bowtie2.sam","--very-sensitive-local","--met-file",out+"/Bowtie2.metrics"])

	logger.info("Sorting sam file...")
        cmd = "samtools sort -n --threads "+str(thread)+" "+out+"/Bowtie2.sam > "+out+"/Bowtie2.sam.sorted"
        os.system(cmd)
	
	call(["rm",out+"/Bowtie2.sam"])

        cmd = "echo 'Post Mapping Metrics' > "+out+"/PostMappingMetrics.txt"
        os.system(cmd)

	if ribo != 'Na':
                RibosomalStat(out+"/Bowtie2.sam.sorted",ribo,out,logger)
                cmd = "cat "+out+"/RibosomalStat.txt > "+out+"/PostMappingMetrics.txt"
                os.system(cmd)

        if post=="True":
                PostQuality(out+"/Bowtie2.sam.sorted",out,logger)

	eXpress(index, out+"/Bowtie2.sorted",extra,logger,out)


def TopHat(read1,read2,read,bowref,out,libtype,bowalgo,thread,logger,gtf,multi,mask,ribo):
	#Function to run tophat
	logger.info("Running TopHat...")
	if read == 'Na':
		cmd = "tophat "+libtype+" "+bowalgo+" -o "+out+"/TopHat "+"-p "+str(thread)+" "+bowref+" "+read1+" "+read2
	else:
                cmd = "tophat "+libtype+" "+bowalgo+" -o "+out+"/TopHat "+"-p "+str(thread)+" "+bowref+" "+read
		
	os.system(cmd)
	cmd = "cat "+ out+"/TopHat/align_summary.txt > "+out+"/PostMappingMetrics.txt"
	os.system(cmd)
	cmd = "echo ' ' >> "+out+"/PostMappingMetrics.txt"
        os.system(cmd)
        cmd = "echo ' ' >> "+out+"/PostMappingMetrics.txt"
        os.system(cmd)

	if ribo != 'Na':
		RibosomalStat(out+"/TopHat/accepted_hits.bam",ribo,out,logger)
		cmd = "cat "+out+"/RibosomalStat.txt > "+out+"/PostMappingMetrics.txt"
		os.system(cmd)

	#PostMapping(out+"/TopHat/accepted_hits.bam",logger,out)
	if post =="True":
		PostQuality(out+"/TopHat/accepted_hits.bam",out,logger)
	Cufflinks(out+"/TopHat/accepted_hits.bam",gtf,thread,libtype,multi,out,logger,mask)

def Cufflinks(file,gtf,thread,libtype,multi,out,logger,mask):
	logger.info("Running Cufflinks...")
	if gtf == "Na":
		logger.warning("Running cufflinks without GTF file! This is not recommended")
		gtf = ''
	else:
		gtf = " --GTF-guide "+gtf

	if mask == "Na":
                logger.warning("Running cufflinks without mask file! This is NOT recommended")
                mask = ''
        else:
                mask = " --mask-file  "+mask

	if libtype!="":
		libtype = " --library-type "+libtype

	call(["cufflinks",file,"-o",out+"/Cufflinks","-u",multi,"--library-type",libtype,"-p",str(thread),gtf,mask])
	cmd = "cufflinks "+file+" -o "+out+"/Cufflinks -u "+multi+" -p "+str(thread)+" "+gtf+" "+mask+" "+libtype
	os.system(cmd)
	logger.info("Outputs can be accessed from "+out+"/Cufflinks/isoforms.fpkm_tracking")	
	TPMCalculator(out,logger)

def TPMCalculator(out,logger):
	#Function to add TPM metrics to Cufflinks output file
	logger.info("Adding TPM metrics to Cufflinks output...")

	with open(out+"/Cufflinks/isoforms.fpkm_tracking") as f:
		sum = 0
		for line in f:
			line = line.strip().split('\t')
			if line[0] != 'tracking_id':
				sum += float(line[9])

	fout = open(out+"/Cufflinks/Final.isoforms.fpkm_tracking",'w')
	with open(out+"/Cufflinks/isoforms.fpkm_tracking") as f:
		for line in f:
                        lines = line.strip().split('\t')
			if lines[0] != 'tracking_id':
				tpm = (float(lines[9])/sum)*(10**6)
				fout.write(line.strip()+"\t"+str(tpm)+"\n")
			else:
				fout.write(line.strip()+"\tTPM\n")

	logger.info("Outputs can be accessed from "+out+"/Cufflinks/Final.isoforms.fpkm_tracking")
			

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
	#Function to reimplement Read_quality.py and read_duplication.py from RSeQC-2.6.3 package written by Liguo Wang
	if (os.path.exists(file)):
		logger.info("Reading in alignment file to compile post alignment stats...")
                obj = SAM.ParseBAM(file)

		#bam_stats.py implementation from RSeqC
		fout = open(out+"/tempStat.txt",'w')
		backup = sys.stdout
        	sys.stdout = StringIO()
        	obj.stat(q_cut = 30)
        	op = sys.stdout.getvalue()
        	sys.stdout.close()  
        	sys.stdout = backup
		logger.info("Restored backup")
        	fout.write(op)
        	fout.close()
	 		
		cmd = "cat "+out+"/tempStat.txt >> "+out+"/PostMappingMetrics.txt"
		print cmd
		os.system(cmd)
		logger.info("Post Mapping Metrics can be accessed from "+out+"/PostMappingMetrics.txt")
				
		logger.info("Plotting post mapping quality graph...")
		obj2 = SAM.ParseBAM(file)
                obj2.readsQual_boxplot(outfile=out+"/"+out, q_cut = 30, shrink = 1000)
                try:
                        subprocess.call("Rscript " + out+"/"+out+ ".qual.r",shell=True)
                except:
                        pass
	
		logger.info("Plotting duplication rate graph")
		
		backobj = SAM.ParseBAM(file)
		backobj.readDupRate(outfile=out+"/"+out,up_bound = 500, q_cut = 30)
		try:
                        subprocess.call("Rscript " + out +"/"+out+  ".DupRate_plot.r", shell=True)
                except:
                        pass
        else:
                logger.error("Failed to plot graphs")
                sys.exit(0)

###Functions from RSEQC###	
def searchit(exon_range, exon_list):
        '''return 1 if find, return 0 if cannot find'''
        for chrom, st, end in exon_list:
                if chrom.upper() not in exon_range:
                        return 0
                elif len(exon_range[chrom].find(st,end)) >=1:
                        return 1
        return 0

def build_bitsets(list):
        '''build intevalTree from list'''
        ranges={}
        for l in list:
                chrom =l[0].upper()
                st = int(l[1])
                end = int(l[2])
                if chrom not in ranges:
                        ranges[chrom] = Intersecter()
                ranges[chrom].add_interval( Interval( st, end ) )
        return ranges
###Functions from RSEQC###

def RibosomalStat(bam,ribo,out,logger):
	#Funcetion to call split_bam.py to calculate ribosomal stats
	logger.info("Reading Ribosome bed file to compute contamination...")
	obj = BED.ParseBED(ribo)
        exons = obj.getExon()
        exon_ranges = build_bitsets(exons)
	
	samfile = pysam.Samfile(bam,'rb')
	fout = open(out+"/RibosomalStat.txt",'w') 

	total_alignment = 0
        in_alignment = 0
        ex_alignment = 0
        bad_alignment = 0
	
	logger.info("Computing Contamination...")
	try:
                while(1):
                        aligned_read = samfile.next()
                        total_alignment += 1
                        
                        if aligned_read.is_qcfail:
                                bad_alignment +=1
                                continue
                        if aligned_read.is_unmapped:
                                bad_alignment +=1
                                continue
                        
                        chrom = samfile.getrname(aligned_read.tid)
                        chrom=chrom.upper()     
                        read_start = aligned_read.pos
                        mate_start = aligned_read.mpos
                                
                        #read_exons = bam_cigar.fetch_exon(chrom, aligned_read.pos, aligned_read.cigar)
                        if aligned_read.mate_is_unmapped:       #only one end mapped
                                if chrom not in exon_ranges:
                                        ex_alignment += 1
                                        continue                
                                else:           
                                        if len(exon_ranges[chrom].find(read_start, read_start +1)) >= 1:
                                                in_alignment += 1
                                                continue
                                        elif len(exon_ranges[chrom].find(read_start, read_start +1)) == 0:
                                                ex_alignment += 1
                                                continue
                        else:                                                   #both end mapped
                                if chrom not in exon_ranges:
                                        ex_alignment += 1
                                        continue
                                else:
                                        if (len(exon_ranges[chrom].find(read_start, read_start +1)) >= 1) or (len(exon_ranges[chrom].find(mate_start, mate_start +1)) >= 1):
                                                in_alignment += 1
                                        else:
                                                ex_alignment += 1
	except StopIteration:
                logger.info("Done")

	fout.write("\nRibosomal RNA Content Stats:\n")
	fout.write("-----------------------------------\n")
	fout.write("\nTotal Records: "+total_alignment)
	fout.write("\nReads consumed by Ribosomal input file: "+in_alignment)
	fout.write("\nReads not consumed by Ribosomal input file: "+ex_alignment)
	fout.write("\nQC Failed/Unmapped reads: "+ad_alignment)
	fout.close()

def GeneratePDF(pre1,pre2,pre,read1,read2,read,out,logger):
	#Function to generate pdf report summary of all the stats
	logger.info("Generating PDF report of summary data")
	doc = Document()
	doc.packages.append(Package('geometry', options=['tmargin=1cm','lmargin=1cm','rmargin=1cm']))

	doc.preamble.append(Command('title', 'Combined Summary Report'))
	doc.append(NoEscape(r'\maketitle'))

	with doc.create(Section('Pre-Processing Data')):
		doc.append(italic("This section contains the input read data before processing."))

        with doc.create(Subsection('Base Quality Graph')):
		doc.append(italic("This graph shows an overview of the range of quality values across all bases at each position in the FastQ file."))
		if read == 'Na':
			with doc.create(Figure(position='h!')) as pic:
                		pic.add_image("preFastqcMetrics/"+read1.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                		pic.add_caption(pre1.split('/')[-1].split('.')[0])
		
			with doc.create(Figure(position='h')) as pic2:
                        	pic2.add_image("preFastqcMetrics/"+read2.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                        	pic2.add_caption(pre2.split('/')[-1].split('.')[0])
		else:
			with doc.create(Figure(position='h')) as pic2:
                                pic2.add_image("preFastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                                pic2.add_caption(pre.split('/')[-1].split('.')[0])

	with doc.create(Subsection('Per Base Sequence Content Graph')):
                doc.append(italic("This graph plots out the proportion of each base position in a file for which each of the four normal DNA bases has been called."))
                if read == 'Na':
                        with doc.create(Figure(position='h!')) as pic:
                                pic.add_image("preFastqcMetrics/"+read1.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                pic.add_caption(pre1.split('/')[-1].split('.')[0])

                        with doc.create(Figure(position='h')) as pic2:
                                pic2.add_image("preFastqcMetrics/"+read2.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                pic2.add_caption(pre2.split('/')[-1].split('.')[0])
                else:
                        with doc.create(Figure(position='h')) as pic2:
                                pic2.add_image("preFastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                pic2.add_caption(pre.split('/')[-1].split('.')[0])

	with doc.create(Subsection('Data Summary')):
		with doc.create(Itemize()) as itemize:
			if read == 'Na':
				itemize.add_item(pre1.split('/')[-1].split('.')[0])
				with open(out+"/"+read1.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt") as f:
					for line in f:
						doc.append(line)

				itemize.add_item(pre2.split('/')[-1].split('.')[0])
                        	with open(out+"/"+read2.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt") as f:
                                	for line in f:
                                        	doc.append(line)
			else:
				itemize.add_item(pre.split('/')[-1].split('.')[0])
                                with open(out+"/"+read.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt") as f:
                                        for line in f:
                                                doc.append(line)
	
	with doc.create(Section('Post-Processing Data')):
                doc.append("")

        with doc.create(Subsection('Base Quality Graph')):
		doc.append(italic("This graph shows an overview of the range of quality values across all bases at each position in the FastQ file."))
                if read == "Na":
			with doc.create(Figure(position='h!')) as pic3:
                        	pic3.add_image("postFastqcMetrics/"+read1.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                        	pic3.add_caption(pre1.split('/')[-1].split('.')[0])

                	with doc.create(Figure(position='h')) as pic4:
                        	pic4.add_image("postFastqcMetrics/"+read2.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                        	pic4.add_caption(pre2.split('/')[-1].split('.')[0])
		else:
			with doc.create(Figure(position='h')) as pic4:
                                pic4.add_image("postFastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                                pic4.add_caption(pre.split('/')[-1].split('.')[0])
	
	with doc.create(Subsection('Per Base Sequence Content Graph')):
                doc.append(italic("This graph plots out the proportion of each base position in a file for which each of the four normal DNA bases has been called."))
                if read == 'Na':
                        with doc.create(Figure(position='h!')) as pic:
                                pic.add_image("postFastqcMetrics/"+read1.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                pic.add_caption(pre1.split('/')[-1].split('.')[0])

                        with doc.create(Figure(position='h')) as pic2:
                                pic2.add_image("postFastqcMetrics/"+read2.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                pic2.add_caption(pre2.split('/')[-1].split('.')[0])
                else:
                        with doc.create(Figure(position='h')) as pic2:
                                pic2.add_image("postFastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                pic2.add_caption(pre.split('/')[-1].split('.')[0])

        with doc.create(Subsection('Data Summary')):
		with doc.create(Itemize()) as itemize:
			if read == 'Na':
                		itemize.add_item(pre1.split('/')[-1].split('.')[0])
                        	with open(out+"/"+read1.split('/')[-1].split('.')[0]+"_postprocess_Summary.txt") as f:
                                	for line in f:
                                        	doc.append(line)

                		itemize.add_item(pre2.split('/')[-1].split('.')[0])
                        	with open(out+"/"+read2.split('/')[-1].split('.')[0]+"_postprocess_Summary.txt") as f:
                                	for line in f:
                                        	doc.append(line)
			else:
				itemize.add_item(pre.split('/')[-1].split('.')[0])
                                with open(out+"/"+read.split('/')[-1].split('.')[0]+"_postprocess_Summary.txt") as f:
                                        for line in f:
                                                doc.append(line)

	if os.path.exists(out+"/PostMappingMetrics.txt"):
		with doc.create(Section('Post Mapping Data')):
			doc.append("")

		with doc.create(Subsection('Data Summary')):
			with open(out+"/PostMappingMetrics.txt") as f:
				for line in f:
					doc.append(line)
		
		if os.path.exists(out+"/"+out+".qual.heatmap.pdf"):
			cmd = "sips -s format png "+out+"/"+out+".qual.heatmap.pdf --out "+out+"/"+out+".qual.heatmap.png"
			os.system(cmd)
		
			cmd = "sips -s format png "+out+"/"+out+".DupRate_plot.pdf --out "+out+"/"+out+".Duprate_plot.png"
			os.system(cmd)

			with doc.create(Subsection('Visual Summary')):
				doc.append(italic("This graph shows an overview of the range of quality values across all bases at each position in the alignment file."))
				with doc.create(Figure(position='h!')) as pic:
					pic.add_image(out+".qual.heatmap.png", width='160px')
					pic.add_caption("Alignment Quality")

				doc.append(italic("This graph shows read duplication rate."))
				with doc.create(Figure(position='h')) as pic2:
					pic2.add_image(out+".DupRate_plot.png", width='160px')
					pic2.add_caption("Duplication Rate")
				with doc.create(Itemize()) as itemize:
					itemize.add_item(italic("Sequence based: reads with identical sequence are regarded as duplicated reads."))
					itemize.add_item(italic("Mapping based: reads mapped to the exactly same genomic location are regarded as duplicated reads."))
		
	else:
		logger.info("No Post Mapping Data available. The rest of the informtation can be accessed from 'Summary_Report.pdf'")	
	
	doc.generate_pdf(out+'/Summary_Report', clean=False)
	

def main():
	global post
	parser = argparse.ArgumentParser(prog='rnaseq',description="Program to help run RNA-seq Analysis")
	required = parser.add_argument_group("Required Options")
	required.add_argument('-t','--type',dest='type',help="Enter the type of analysis you want to do (F=Fast, S=Slow, N=Novel) [default=F]",choices=['F','S','N'],required=True, default="F")
	required.add_argument('-r1','--read1',dest='read1',help="Enter the location of read 1", metavar="\b",default="Na")
	required.add_argument('-r2','--read2',dest='read2',help="Enter the location of read 2", metavar="\b",default="Na")
	required.add_argument('-r','--read',dest='read',help="Enter the location of single ended read",metavar="\b",default="Na")
	required.add_argument('-o','--output',dest="out",help="Enter the prefix to your output files", required=True, metavar="\b",default="out")	

	fastqc = parser.add_argument_group("FastQC options")
	fastqc.add_argument('-fk','--fastqck',dest='fastqck',help="Set the kmer size for running fastqc [Default=7]",default=7,type=int,metavar="kmer-size(int)")
	
	trimmomatic = parser.add_argument_group("Trimmomatic options")
	parser.add_argument('--threads',dest="thread",help="Set the number of threads to use for trimmomatic [Default=4]",default=4,type=int,metavar="number(int)")
	trimmomatic.add_argument('--leading',dest="trimlead",help="Specify the minimum quality to keep a base towards the start of the sequence [Default = 5]",type=int,default=5,metavar="quality(int)")
	trimmomatic.add_argument('--trailing',dest="trimtrail",help="Specify the minimum quality to keep a base towards the end of the sequence [Default = 5]",type=int,default=5,metavar="quality(int)")
	trimmomatic.add_argument('--headcrop',dest="trimcrop",help="Specify the number of bases to remove, from the start of the read [Default = 0]",type=int,default=0,metavar="length(int)")	
	trimmomatic.add_argument('--minlen',dest="trimlen",help="Specify the minimum lengths of read to keep [Default = 0]",type=int,default=0,metavar="length(int)")
	trimmomatic.add_argument('--window',dest="trimwindow",help="Specify the number of bases to average across for sliding window [Default=0]",type=int,default=0,metavar="size(int)")
	trimmomatic.add_argument('--quality',dest="trimq",help="Specify the average quality required [Default = 0]", type=int, default=0, metavar="quality(int)")
	trimmomatic.add_argument('--adapter',dest="adapter",help="Specify the adapter file to use to clip from reads", default="adapters/TruSeq2-PE.fa", metavar="<path to adapter file>")

	kallisto = parser.add_argument_group("Fast Alginment Option - Kallisto")
	kallisto.add_argument('--kref',dest='kref',help="Enter the fasta file to be used to construct index",metavar="<path to file>",default="Na")
	kallisto.add_argument('--kindex',dest='kindex',help="Enter the kallisto index to be used for quantification",metavar="<path to file>",default="Na")
	kallisto.add_argument('--kkmer',dest='kkmer',help="Specify the kmer length to be used [Default = 31]", type=int, default=31,metavar="kmer-size(int)")
	kallisto.add_argument('--kboost',dest='kboost',help="Specify the number of bootstraps to perform [Default = 0]", type=int, default=0, metavar="num(int)")		
	kallisto.add_argument('-l','--fragment-length-mean',dest='length',help="Specify the fragment length mean of the read",type=int, metavar="num(int)",default=0)
	kallisto.add_argument('-s','--std-dev',dest='sd',help="Specify the standard deviation of the read",type=int,metavar="num(int)",default=-1)

	slow = parser.add_argument_group("Slow Allignment Option - BWA + eXpress")
	slow.add_argument('--bref',dest='bref',help="Enter the fasta file to be used to construct the index", metavar="<path to .fa file>",default="Na")
	slow.add_argument('--bindex',dest='bindex',help="Enter the bwa index file to be used for quantification", metavar="<path to bwa index file>", default="Na")
	slow.add_argument('--algo',dest='algo',help="Specify the algorithm to be used for algorithm [Default = mem]", choices=["mem","bwasw","sampe"], default="mem") 
	slow.add_argument('--strand',dest='strand',help="Specify Strand option to run eXpress", choices=['fr-stranded','rf-stranded','f-stranded','r-stranded'], default="None")	
	slow.add_argument('--bwa',dest='useBwa',help="Use this option if you want the pipeline to force use BWA [Default=False]", action="store_const",const="True")

	novel = parser.add_argument_group("Novel Transcription Detection - TopHat + Cufflinks")
	novel.add_argument('--bowref',dest="bowref",help="Enter the file to be used to construct bowtie index", metavar="<path to file>", default="Na")
	novel.add_argument('--bowindex',dest="bowindex", help="Enter the bowtie index", metavar="<path to file>", default="Na")
	novel.add_argument('--libtype',dest="libtype",help="Specify the library type", choices=['fr-unstranded','fr-firststrand','fr-secondstrand','Na'], default="Na")
	novel.add_argument('--bowmode',dest="bowalgo",help="Specify the bowtie mode", choices=['b2-very-fast','b2-fast','b2-sensitive','b2-very-sensitive'], default="Na")		
	novel.add_argument('-g','--gtf',dest="gtf",help="Enter the GTF file to be used for annotations",metavar="<path to file>",default="Na")	
	novel.add_argument('-u','--multi-read-correct',dest="multi",help="Use 'rescue method' for multi-reads(more accurate)",choices=['True','False'],default="False")
	novel.add_argument('-m','--mask-file',dest="mask",help="Ignore all alignment within transcripts in this file",metavar="<path to file>",default="Na")

	post = parser.add_argument_group("Post Mapping Metrics Graph")
	post.add_argument('-p','--post',dest="post",help="Use this option if you want post mapping metrics graph [Default=True]", action="store_const",const="True",default="True")
	post.add_argument('--ribosome',dest='rrna',help="Specifiy the ribosomal bed file in order to calculate statistics of ribosomal content", default='Na', metavar="<path to file>")

	parser.add_argument("--skippre",dest="skippre",help="Use this option if you want to skip pre processing steps [Default=False]",action="store_const",const="True",default="False")
	parser.add_argument("--normalize", dest="norm",help="Use this option if you want to normalize the input data [Default=False]", action="store_const",const="True",default="False")
	parser.add_argument('-v','--version', action='version', version='%(prog)s 2.1.1')
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
	if(args.read == 'Na' and args.read1 == 'Na' and args.read2 == 'Na'):
                logger.error("No reads given")
                parser.print_help()
                sys.exit()

	if(args.read != 'Na' and args.read1 != 'Na' and args.read2 != 'Na'):
                logger.error("Please choose between single end read and paired end read")
                parser.print_help()
                sys.exit()

	if(args.read != 'Na' and not os.path.exists(args.read)):
                logger.error("Single read not readable")
                sys.exit()

	if(args.read1 != 'Na' and not os.path.exists(args.read1) and args.read == 'Na'):
		logger.error("The read 1 file is not readable")
		parser.print_help()
		sys.exit()

	if(args.read2 != 'Na' and not os.path.exists(args.read2) and args.read == 'Na'):
                logger.error("The read 2 file is not readable")
                parser.print_help()
                sys.exit()

	adapter = args.adapter

	if(not os.path.exists(args.adapter)):
		logger.warning("The adapter file to clip reads is not specified")
		adapter = 'Na'

	if(args.type == "F"):
		if (args.kindex == "Na" and args.kref == "Na"):
                        logger.error("Please enter either a reference file to build kallisto reference [--kref] or a kallisto index [--kindex]")
			sys.exit()
			
                elif (args.kindex == 'Na' and not os.path.exists(args.kref)):
                        logger.error("The file for building index is not readable")
                        sys.exit()

		elif (args.kref == 'Na' and not os.path.exists(args.kindex)):
                        logger.error("Trouble reading kallisto index file")
                        sys.exit()
		
		if args.read != "Na":
			if args.length == 0 or args.sd == -1:
				logger.error("Fragment length mean and sd must be supplied for single-end reads using -l and -s")
				sys.exit()
	
	if(args.type == "S"):
		if (args.bindex == "Na" and args.bref == "Na"):
                        logger.error("Please enter either a reference file to build reference [--bref] or an index file [--bindex]")
                        sys.exit()

                elif (args.bindex == 'Na' and not os.path.exists(args.bref)):
                        logger.error("The file for building index is not readable")
                        sys.exit()

                elif (args.bref == 'Na' and not os.path.exists(args.bindex)):
                        logger.error("Trouble reading bwa index file")
                        sys.exit()	

		if args.useBwa != "True":
			if (args.bref == 'Na' and not os.path.exists(args.bindex+".1.bt2")):
                        	logger.error("Trouble reading bowtie files (.1.bt2)")

                	if (args.bref == 'Na' and not os.path.exists(args.bindex+".rev.1.bt2")):
                        	logger.error("Trouble reading bowtie reference files (.rev.1.bt2)")

	if(args.type == "N"):
		if (args.bowindex == "Na" and args.bowref == "Na"):
			logger.error("Please enter either a reference file to build bowtie reference [--bowref] or a bowtie2 index (.ebwt) [--bowindex]")
			sys.exit()

                elif (args.bowindex == 'Na' and not os.path.exists(args.bowref)):
                        logger.error("The file for building index is not readable")
                        sys.exit()

		elif (args.bowref == 'Na' and not os.path.exists(args.bowindex+".1.bt2")):
			logger.error("Trouble reading bowtie files (.1.bt2)")
		
		elif (args.bowref == 'Na' and not os.path.exists(args.bowindex+".rev.1.bt2")):
                        logger.error("Trouble reading bowtie reference files (.rev.1.bt2)")

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
	read = args.read
	
	#Creating Directory Structure
	call(['mkdir',args.out])
		
	#Pre Processing Data
	if (args.skippre=="False"):
		log = logging.getLogger('Pre Processing')
                log.setLevel(logging.DEBUG)
                log.addHandler(stream)

		if (args.read == 'Na'):			
			call(['mkdir',args.out+'/preFastqcMetrics'])
			Quality_Assessment(args.read1,args.fastqck,log,args.out,"pre")
			Quality_Assessment(args.read2,args.fastqck,log,args.out,"pre")

			TrimmingPE(args.read1,args.read2,args.thread,args.trimlead,args.trimtrail,args.trimcrop,args.trimlen,args.trimwindow,args.trimq,args.out,adapter,logger)
			read1 = args.out+"/triM"+args.read1.split('/')[-1]
			read2 = args.out+"/triM"+args.read2.split('/')[-1]

			if (args.norm=="True"):
                		khmer(read1,read2,args.fastqck,logger)
                		read1 = read1+".keep"
                		read2 = read2+".keep"
			
			#Renaming preProcessing Summary files for reads
			call(['mv',args.out+"/preFastqcMetrics/"+args.read1.split('/')[-1].split('.')[0]+"_fastqc",args.out+"/preFastqcMetrics/"+read1.split('/')[-1].split('.')[0]+"_fastqc"])
			call(['mv',args.out+"/preFastqcMetrics/"+args.read2.split('/')[-1].split('.')[0]+"_fastqc",args.out+"/preFastqcMetrics/"+read2.split('/')[-1].split('.')[0]+"_fastqc"])
			call(['mv',args.out+"/"+args.read1.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt",args.out+"/"+read1.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt"])
			call(['mv',args.out+"/"+args.read2.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt",args.out+"/"+read2.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt"])

			call(['mkdir',args.out+'/postFastqcMetrics'])
                	Quality_Assessment(read1,args.fastqck,log,args.out,"post")
                	Quality_Assessment(read2,args.fastqck,log,args.out,"post")
		else:
			logger.info("Running Pre Processing for single end reads")
			call(['mkdir',args.out+'/preFastqcMetrics'])
			Quality_Assessment(args.read,args.fastqck,log,args.out,"pre")

			TrimmingSE(args.read,args.thread,args.trimlead,args.trimtrail,args.trimcrop,args.trimlen,args.trimwindow,args.trimq,args.out,adapter,logger)
                        read = args.out+"/triM"+args.read.split('/')[-1]

			#Renaming preProcessing Summary files for reads
                        call(['mv',args.out+"/preFastqcMetrics/"+args.read.split('/')[-1].split('.')[0]+"_fastqc",args.out+"/preFastqcMetrics/"+read.split('/')[-1].split('.')[0]+"_fastqc"])
			call(['mv',args.out+"/"+args.read.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt",args.out+"/"+read.split('/')[-1].split('.')[0]+"_preprocess_Summary.txt"])

			call(['mkdir',args.out+'/postFastqcMetrics'])
                        Quality_Assessment(read,args.fastqck,log,args.out,"post")
	else:
		logger.info("User opted to skip Pre Processing step")
	
	#Running RNA-seq according to what the user selected:
	log = logging.getLogger('RNA Seq')
	log.setLevel(logging.DEBUG)
	log.addHandler(stream)

	if args.type=="F":
		log.info("User opted for fast analysis, Prepping for Kallisto...")
		Kallisto(read1,read2,read,args.kindex,args.kref,args.kkmer,args.kboost,args.thread,logger,args.out,args.length,args.sd)

	elif args.type=="S":
		extra = ""

		if args.useBwa=="True":
			log.info("User opted for slow analysis, Prepping for bwa and eXpress...")
			if args.strand != "None":
				extra = "--"+args.strand
			Bwa(read1,read2,read,args.bindex,args.bref,args.algo,extra,logger,args.out,args.thread,args.rrna)
		else:
			log.info("User opted for slow analysis, Prepping for bowtie2 and eXpress...")
	
			#determine phred score of reads
			if read != "Na":
				phred = DeterminePhred(read,args.out,"pre",logger)
			else:
				phred = DeterminePhred(read1,args.out,"pre",logger)

			#calling bowtie2
			BowtieAlign(extra,read1,read2,read,args.bindex,args.bref,logger,args.out,args.thread,args.rrna,phred)
			
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
			BowtieIndex(args.bowref,args.out,logger,libtype,bowalgo,read1,read2,read,args.thread,args.gtf,args.multi,args.mask,args.rrna)
		else:
			print args.bowindex
			TopHat(read1,read2,read,args.bowindex,args.out,libtype,bowalgo,args.thread,logger,args.gtf,args.multi,args.mask,args.rrna)
	
	GeneratePDF(args.read1,args.read2,args.read,read1,read2,read,args.out,logger)	

if __name__ == "__main__":
	main()
