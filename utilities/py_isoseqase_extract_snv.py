#!/usr/bin/env python
import sys,re,time,argparse
from multiprocessing import cpu_count,Pool

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	extract_snv_from_vcf(args.input,args.output,args.sample_id)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()


def extract_snv_from_vcf(input_vcf,output_snv_txt,sample_id):
	base4 = ["A","T","C","G"] 
	for line in input_vcf:
		if line.startswith("##"): continue # meta-information lines
		if line.startswith("#CHROM"): # header line
			sample_list = line.strip().split("\t")[9:]
			flag = 0
			for sample in sample_list:
				if sample == sample_id:
					sample_idx = sample_list.index(sample)
					flag = 1
					break
			if flag == 0:
				sys.stderr.write("Error: cannot find the sample ID in VCF file\n")
				sys.exit()
		
		chrom,pos,id,ref,alt_set,qual,filter,info,format_set = line.strip().split("\t")[:9]
		if ref not in base4: continue # single base for ref
		if len(alt_set.split(",")) > 2: continue # only handle diploid currently
		if (alt_set.split(",")[0] not in base4) or (alt_set.split(",")[-1] not in base4): continue # single base for alt

		snv =  (ref + "," + alt_set).split(",")
		format_index = format_set.split(":").index("GT")
		idv_gt = line.rstrip().split("\t")[9:][sample_idx].split(":")[format_index]
		if "|" in idv_gt:
			if snv[int(idv_gt.split("|")[0])] == snv[int(idv_gt.split("|")[1])]: continue # homologous
			print >>output_snv_txt, "\t".join([chrom,pos,snv[int(idv_gt.split("|")[0])],snv[int(idv_gt.split("|")[1])]])
		else:
			if snv[int(idv_gt.split("/")[0])] == snv[int(idv_gt.split("/")[1])]: continue # homologous
			print >>output_snv_txt, "\t".join([chrom,pos,snv[int(idv_gt.split("/")[0])],snv[int(idv_gt.split("/")[1])]])

def do_inputs():
	output_txt_format = '''
1. chromosome id
2. position
3. allele1
4. allele2'''

	parser = argparse.ArgumentParser(description="Function: extract SNV (single nucleotide variant) information from VCF file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: vcf file with header line showing the sample ID of your interest")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: txt file")
	parser.add_argument('-s','--sample_id',type=str,required=True,help="Sample id shown after 9th field")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
