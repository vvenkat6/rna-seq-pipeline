#!/Users/vv39/anaconda2/bin/python #change path as necessary to run python

#This is a wrapper script for rnaseq.py to run runseq.py for all the single ended samples in a directory

import sys,os
import argparse
import logging
from os import listdir
import glob
from subprocess import call

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Package, Matrix, Command, Itemize, SubFigure
from pylatex.utils import italic, NoEscape

def main():
	parser = argparse.ArgumentParser(prog='rnaseq',description="Program to help run RNA-seq Analysis")
        required = parser.add_argument_group("Required Options")
	
	required.add_argument('-t','--type',dest='type',help="Enter the type of analysis you want to do (F=Fast, S=Slow, N=Novel) [default=F]",choices=['F','S','N'],required=True, default="F")
        required.add_argument('-f','--folder',dest='folder',help="Enter the location of the folder containing the reads", metavar="\b",default="Na",required=True)
	required.add_argument('-i','--index',dest='index',help="Enter the location of the index file for the run u choose",default="Na",required=True)

	parser.add_argument('-g','--gtf',dest='gtf',help="Enter the location of the gtf file", default='Na')
	parser.add_argument('-m','--ribo',dest='ribo',help="Enter the location of the ribosomal file", default='Na')
	args = parser.parse_args()

	logger = logging.getLogger('Main')
        logger.setLevel(logging.DEBUG)
        stream = logging.StreamHandler()
        stream.setLevel(logging.DEBUG)
        formats = logging.Formatter('[%(asctime)s %(name)s %(levelname)s]: %(message)s')
        stream.setFormatter(formats)
        logger.addHandler(stream)
	
	if not os.path.exists(args.index):
		logger.error("The index file is not readable")
		sys.exit()

	reads =  os.listdir(args.folder)	
	extension = ''	
	
	#Running the pipeline
	for read in reads:
		logger.info("Processing "+read+" ...")
		out = read.split('/')[-1].split('.')[0]
		read = args.folder+"/"+read
		
		if args.type == "F" and not os.path.exists(out+"_Kallisto"):
			call(["time","./rnaseq.py","-t","F","-r",read,"--kindex",args.index,"-l","54","-s","10","-o",out+"_Kallisto"])
			logger.info("Ouputs written into "+out+"_Kallisto")
			extension = "_Kallisto"

		elif args.type == "S" and not os.path.exists(out+"_BWA"):
			extra = ''
			if args.ribo != 'Na':
				extra = "-m "+args.ribo
			call(["time","./rnaseq.py","-t","S","-r",read,"--bindex",args.index,"-o",out+"_BWA",extra])
			logger.info("Ouputs written into "+out+"_BWA")
			extension = "_BWA"

		elif args.type == "N" and not os.path.exists(out+"_TopHat"):
			extra = '' 
			if args.gtf != 'Na':
				extra = "-g "+args.gtf
			if args.ribo != 'Na':
				extra = extra+" -m "+args.ribo
			call(["time","./rnaseq.py","-t","N","-r",read,"--bowindex",args.index,"-o",out+"_TopHat",extra])
			extension = "_TopHat"

		else:
			logger.info("Skipping this read as folder already exists...")
	
	#Converging Report
        doc = Document()

        doc.packages.append(Package('geometry', options=['tmargin=1cm','lmargin=1cm','rmargin=1cm']))

	with doc.create(Section('Base Quality Graph')):
                doc.append(italic("This graph shows an overview of the range of quality values across all bases at each position in the FastQ file."))

	extension = "_BWA"
	for read in reads:
                logger.info("Converging Pre Graphs Reports "+read+" ...")
                out = read.split('/')[-1].split('.')[0]+extension
			
		with doc.create(Figure(position='h!')) as baseQuality:
			with doc.create(SubFigure(position='b',width=NoEscape(r'0.45\linewidth'))) as prepic:
                		prepic.add_image(out+"/preFastqcMetrics/triM"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                       		prepic.add_caption("Pre "+read.split('/')[-1].split('.')[0])

 			with doc.create(SubFigure(position='b',width=NoEscape(r'0.45\linewidth'))) as postpic:
                        	postpic.add_image(out+"/postFastqcMetrics/triM"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_quality.png", width='160px')
                        	postpic.add_caption("Post "+read.split('/')[-1].split('.')[0])

	with doc.create(Section('Per Base Sequence Content Graph')):
                doc.append(italic("This graph plots out the proportion of each base position in a file for which each of the four normal DNA bases has been called."))

        for read in reads:
                out = read.split('/')[-1].split('.')[0]+extension

                with doc.create(Figure(position='h!')) as baseQuality:
                        with doc.create(SubFigure(position='b',width=NoEscape(r'0.45\linewidth'))) as prepic:
                                prepic.add_image(out+"/preFastqcMetrics/triM"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                prepic.add_caption("Pre "+read.split('/')[-1].split('.')[0])

                        with doc.create(SubFigure(position='b',width=NoEscape(r'0.45\linewidth'))) as postpic:
                                postpic.add_image(out+"/postFastqcMetrics/triM"+read.split('/')[-1].split('.')[0]+"_fastqc/Images/per_base_sequence_content.png", width='160px')
                                postpic.add_caption("Post "+read.split('/')[-1].split('.')[0])	
		
	with doc.create(Section('Data Summary')):
                doc.append(italic("This section contains post trimming FastQC summary results."))
	
	for read in reads:
                out = read.split('/')[-1].split('.')[0]+extension
		
		with doc.create(Subsection(read)):
			doc.append("")
			with open(out+"/triM"+read.split('/')[-1].split('.')[0]+"_postprocess_Summary.txt") as f:
                		for line in f:
                        		doc.append(line)

	if os.path.exists(out+"/PostMappingMetrics.txt"):
        	with doc.create(Section('Post Mapping Data')):
                	doc.append("")

		for read in reads:
                	out = read.split('/')[-1].split('.')[0]+extension

                	with doc.create(Subsection(read)):
                        	with open(out+"/PostMappingMetrics.txt") as f:
                                	for line in f:
                                        	doc.append(line)
	
			if os.path.exists(out+"/"+out+".qual.heatmap.pdf"):
                        	cmd = "sips -s format png "+out+"/"+out+".qual.heatmap.pdf --out "+out+"/"+out+".qual.heatmap.png"
                        	os.system(cmd)
                
                        	cmd = "sips -s format png "+out+"/"+out+".DupRate_plot.pdf --out "+out+"/"+out+".Duprate_plot.png"
                        	os.system(cmd)
		
		with doc.create(Section('Post Mapping Visual Data')):
                        doc.append("")
				
			for read in reads:
                		out = read.split('/')[-1].split('.')[0]+extension
				
				if os.path.exists(out+"/"+out+".qual.heatmap.png"):	
                			with doc.create(Figure(position='h!')) as VisualSummary:
                        			with doc.create(SubFigure(position='b',width=NoEscape(r'0.45\linewidth'))) as pic:
							pic.add_image(out+"/"+out+".qual.heatmap.png", width='160px')
                                        		pic.add_caption("Alignment Quality")
                        			with doc.create(SubFigure(position='b',width=NoEscape(r'0.45\linewidth'))) as pic2:
							pic2.add_image(out+"/"+out+".DupRate_plot.png", width='160px')
                                        		pic2.add_caption("Duplication Rate")	
						VisualSummary.add_caption(read)

	doc.generate_pdf('Summary_Report', clean=False)
        	

if __name__ == "__main__":
        main()
