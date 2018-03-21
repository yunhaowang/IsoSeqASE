#!/usr/bin/env python
import sys,re,time,argparse,collections
from multiprocessing import cpu_count,Pool
import motility

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	dic_chr_pos_snv = parse_snv(args.input_snv)
	dic_lr_info = parse_lr_seq(args.input_lr)
	output_gpd = args.output
	p = Pool(processes=args.cpu)
	csize = 10
	results = p.imap(func=phase,iterable=generate_tx(args.input_iso),chunksize=csize)
	for res in results:
		if not res: continue
		output_gpd.write(res+"\n")
	output_gpd.close()
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_snv(input_snv_txt): # parse SNV information from snv TXT file, TXT file must be 'sort-k1,1 -k2,2n -u'
	global dic_chr_pos_snv
	dic_chr_pos_snv = {}
	for line in input_snv_txt:
		chr,pos,snv1,snv2 = line.strip().split("\t")[:4]
		if chr not in dic_chr_pos_snv.keys():
			dic_chr_pos_snv[chr] = {}
			dic_chr_pos_snv[chr][int(pos)] = [snv1,snv2]
			dic_chr_pos_snv[chr]["pos_list"] = []
			dic_chr_pos_snv[chr]["pos_list"].append(int(pos))
		else:
			dic_chr_pos_snv[chr][int(pos)] = [snv1,snv2]
			dic_chr_pos_snv[chr]["pos_list"].append(int(pos))
	input_snv_txt.close()
	return dic_chr_pos_snv

def parse_lr_seq(input_polish_lr_gpd): #  parse long read information from polished long read GPD file
	global dic_lr_info
	dic_lr_info = {}
	for line in input_polish_lr_gpd:
		read_id,read_id,chrom,strand,tss,tts,mapq,sf,exon_number,exon_start,exon_end,sam_flag,seq = line.rstrip("\n").split("\t")[:13]
		dic_lr_info[read_id] = [tss,tts,exon_number,exon_start,exon_end,seq]
	input_polish_lr_gpd.close()
	return dic_lr_info


def get_exon_snv_pos(pos,exon_number,tss,tts,exon_start,exon_end): # determine the SNV in the exon region
	new_pos = ""
	if int(exon_number) > 1:
		exon_start_list = [int(i) for i in exon_start.split(",")[:-1]]
		exon_end_list = [int(i) for i in exon_end.split(",")[:-1]]
		for i in range(0,int(exon_number)):
			if pos > exon_start_list[i] and pos <= exon_end_list[i]:
				new_pos = pos
				break
	else:
		if pos > int(tss) and pos <= int(tts):
			new_pos = pos
	return new_pos


def get_snv_arrary(dic_chr_pos_snv,chr,pos_list,read_info_list): # for each long read, make the SNV string
	tss,tts,exon_number,exon_start,exon_end,seq = read_info_list
	snv_list = []
	if int(exon_number) == 1:
		exon_seq = seq
		for pos in pos_list:
			if pos > int(tss) and pos <= int(tts): # in exon region
				if exon_seq[pos-int(tss)-1] in dic_chr_pos_snv[chr][pos]: # nucleotide in vcf file
					snv = exon_seq[pos-int(tss)-1]
				else:
					snv = "-"
			else:
				snv = "-"
			snv_list.append(snv)
		snv_set = "".join(snv_list)

	else:
		exon_start_list = [int(i) for i in exon_start.split(",")[:-1]]
		exon_end_list = [int(i) for i in exon_end.split(",")[:-1]]
		for pos in pos_list:
			new_seq = seq
			if pos > int(tss) and pos <= int(tts): # in isoform region
				for i in range(0,int(exon_number)):
					exon_seq = new_seq[:(exon_end_list[i]-exon_start_list[i])]
					if pos > exon_start_list[i] and pos <= exon_end_list[i]: # in exon region
						if exon_seq[pos-exon_start_list[i]-1] in dic_chr_pos_snv[chr][pos]: # nucleotide in vcf file
							snv = exon_seq[pos-exon_start_list[i]-1]
						else:
							snv = "-"
						break
					else:
						snv = "-"
					new_seq = new_seq[(exon_end_list[i]-exon_start_list[i]):]
			else:
				snv = "-"
			snv_list.append(snv)
		snv_set = "".join(snv_list)
	return snv_set

def get_consensus_snv_arrary(dic_chr_pos_snv,chr,snv_pos_list,snv_arrary_list): # get consensus SNV (also called 'phasing') using Position-Weight Matrices (commonly used in DNA/RNA/Protein motif finding)
	# first step: filter some SNV positions with many missing data
	retain_snv_idx = []
	retain_snv_pos_list = []
	for i in range(len(snv_pos_list)):
		ret_snv_count = 0
		for snv_arrary in snv_arrary_list:
			if snv_arrary[i] != "-":
				ret_snv_count += 1
		if float(ret_snv_count)/len(snv_arrary_list) >= 0.5: # require >50% long read have nucleotide
			retain_snv_idx.append(i)
			retain_snv_pos_list.append(str(snv_pos_list[i]))
	retain_snv_pos = ",".join(retain_snv_pos_list)
	
	# second step: get the consensus (phased major allele)
	if retain_snv_idx != []: # have valid SNV position for constructing new SNV arrary
		# obtain new SNV arrary
		new_snv_arrary_list = []
		for snv_arrary in snv_arrary_list:
			new_snv_arrary = ""
			for idx in retain_snv_idx:
				new_snv_arrary += snv_arrary[idx]
			new_snv_arrary_list.append(new_snv_arrary)
		
		# obtain SNV arrary for constructing consensus 
		consus_snv_arrary_list = []
		for new_snv_arrary in new_snv_arrary_list:
			if "-" not in new_snv_arrary:
				consus_snv_arrary_list.append(new_snv_arrary)

		if consus_snv_arrary_list != []: # have valid SNV arrary for constructing consensus
			# get major allele by most common element in the list
			counter_cm_ele_dic = collections.Counter(consus_snv_arrary_list)
			most_cm_ele = counter_cm_ele_dic.most_common(2)
			if most_cm_ele[0][1] >= 5: # if the snv_arrary is support by >= 5 long reads
				major_allele = most_cm_ele[0][0]
			else:
				# get major allele by PWM
				pwm = motility.make_pwm(consus_snv_arrary_list) # PWM calculate
				dic_snv_uniq_score = {}
				for snv_uniq in set(consus_snv_arrary_list):
					snv_score = float(pwm.calc_score(snv_uniq))
					if snv_score not in dic_snv_uniq_score.keys():
						dic_snv_uniq_score[snv_score] = []
						dic_snv_uniq_score[snv_score].append(snv_uniq)
					else:
						dic_snv_uniq_score[snv_score].append(snv_uniq)
					max_consus_snv = dic_snv_uniq_score[max(dic_snv_uniq_score.keys())][0] # get SNV arrary with maximal PWM score as major allele 
				major_allele = max_consus_snv

			# get minor allele
			minor_phase_list = []
			for i in range(len(retain_snv_pos_list)):
				if dic_chr_pos_snv[chr][int(retain_snv_pos_list[i])][0] != major_allele[i]:
					minor_phase_list.append(dic_chr_pos_snv[chr][int(retain_snv_pos_list[i])][0])
				else:
					minor_phase_list.append(dic_chr_pos_snv[chr][int(retain_snv_pos_list[i])][1])
			minor_allele = "".join(minor_phase_list)

			# quantify allele-specific read count for both alleles
			major_allele_c,minor_allele_c,uncertain_allele_c = 0,0,0
			for new_snv_arrary in new_snv_arrary_list:
				common_count_major = sum(1 for a,b in zip(major_allele,new_snv_arrary) if a==b)
				common_count_minor = sum(1 for a,b in zip(minor_allele,new_snv_arrary) if a==b)
				if float(common_count_major)/len(major_allele) > 0.5:
					major_allele_c += 1
				elif float(common_count_minor)/len(minor_allele) > 0.5:
					minor_allele_c += 1
				else:
					uncertain_allele_c += 1
			
			output_res = "\t".join([retain_snv_pos,major_allele,minor_allele,str(major_allele_c),str(minor_allele_c),str(uncertain_allele_c)])
		else:
			output_res = "\t".join(["*","*","*","*","*","*"])
	else:
		output_res = "\t".join(["*","*","*","*","*","*"])
	return output_res

def generate_tx(input_iso_gpd):
	z = 0
	for line in input_iso_gpd:
		z += 1
		yield (line,z)
	input_iso_gpd.close()

def phase(inputs):
	(line,z) = inputs
	gene_id,iso_id,chr,strand,tss,tts,fl_lr,tt_lr,exon_number,exon_start,exon_end = line.strip().split("\t")[:11]
	read_id_set = line.strip().split("\t")[-1].split(",")

	if chr not in dic_chr_pos_snv.keys():
		output_ase_gpd = "\t".join(["\t".join(line.strip().split("\t")[:-1]),"*","*","*","*","*","*","*","*"])
	else:
		if len(dic_chr_pos_snv[chr].keys()) == 1:
			output_ase_gpd = "\t".join(["\t".join(line.strip().split("\t")[:-1]),"*","*","*","*","*","*","*","*"])
		else:
			new_pos_list = []
			for pos in dic_chr_pos_snv[chr]["pos_list"]:
				new_pos = get_exon_snv_pos(pos,exon_number,tss,tts,exon_start,exon_end)
				if pos < int(tss):
					continue
				elif pos > int(tts):
					break
				else:
					new_pos = get_exon_snv_pos(pos,exon_number,tss,tts,exon_start,exon_end)
					if new_pos != "":
						new_pos_list.append(new_pos)

			if new_pos_list == []:
				output_ase_gpd = "\t".join(["\t".join(line.strip().split("\t")[:-1]),"*","*","*","*","*","*","*","*"])
			else:
				snv_arrary_list = []
				for read_id in read_id_set:
					read_info_list = dic_lr_info[read_id]
					snv_arrary = get_snv_arrary(dic_chr_pos_snv,chr,new_pos_list,read_info_list)
					snv_arrary_list.append(snv_arrary)
				snv_arrary_set = ",".join(snv_arrary_list)
				phase_res = get_consensus_snv_arrary(dic_chr_pos_snv,chr,new_pos_list,snv_arrary_list)
				output_ase_gpd = "\t".join(["\t".join(line.strip().split("\t")[:-1]),",".join([str(i) for i in new_pos_list]),snv_arrary_set,phase_res])
	return output_ase_gpd

def do_inputs():
	output_txt_format = '''

1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS (+)
6. TTS (+)
7. number of support full-length long reads
8. number of support total long reads
9. exon count
10. exon start set
11. exon end set
12. For novel isoform, derived genic locus
13. For novel isoform, overlap percentage with derived genic locus
14. For novel singleton isoform, if it is located at the last exon of any known isoform. If yes, isoform ID otherwise '-'
15. For novel singleton isoform, the overlap percentage with the the last exon
16. For novel multi-exon isoform, number of splice sites are detected by anno and/or short reads; and the total number of splice sites
17. For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of known multi-exon isoform, isoform ID if yes otherwise '-'
18. For novel isoform, maximal length of polyA track in defined region
19. For novel isoform, maximal percentage of nucleotide A in defined region
20. Positions with SNV, split by ','
21. SNV pattern for all long reads, split by ','
22. Positions with SNV after filtering some positions with missing value (the nucleotide at the position of long reads is not annotated by SNV TXT file, 50%)
23. Major allele pattern (assumed)
24. Minor allele pattern (assumed)
25. Long read count for major allele
26. Long read count for minor allele
27. Long read count for uncertain allele'''

	parser = argparse.ArgumentParser(description="Function: haplotype and quantify ASE",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input_iso',type=argparse.FileType('r'),required=True,help="Input: constructed isoform gpd file generated by 'py_isoseqase_generate_output.py")
	parser.add_argument('-l','--input_lr',type=argparse.FileType('r'),required=True,help="Input: polished long read gpd file generated by 'py_isoseqase_polish.py'")
	parser.add_argument('-s','--input_snv',type=argparse.FileType('r'),required=True,help="Input: snv txt file generated by 'py_isoseqase_extract_snv.py' and then sort by 'sort -k1,1 -k2,2n -u'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed ASE-specific isoform gpd file")
	parser.add_argument('-p','--cpu',type=int,default=cpu_count(),help="Number of threads")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
