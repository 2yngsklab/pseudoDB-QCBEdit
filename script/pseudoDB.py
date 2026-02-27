import re
import os
import sys
import shutil
import argparse
import subprocess

def put_file(src, dest, softlink = False):
	"""
	Function to either copy files or create softlinks.
	If softlink creation fails, it will fall back to copying.

	Arguments:
	- src: Path to source file
	- dest: Path to destination
	- softlink: Create softlink if true. Otherwise, copy the file.
	"""
	if os.path.abspath(src) == os.path.abspath(dest):
		print(f"Source and destination are the same file: {src}, {dest}. Skipping.")
		return
	
	if softlink:
		try:
			os.symlink(src, dest)
		except FileExistsError:
			print(f"Symlink or file already exists: {dest}. Skipping.")
		except Exception as e:
			print(f"Failed to create symlink: {e}")
			print(f"Attempting to copy file instead...")
			return put_file(src, dest, softlink = False)
		else:
			print(f"Created symlink: {dest} -> {src}")
	else:
		try:
			shutil.copy2(src, dest)
		except Exception as e:
			raise type(e)(f"Failed to copy file {src} to {dest}: {e}") from e
		else:
			print(f"Copied file {dest} from {src}")

def process_sample_list(sample_paths_file):
	"""
	Reads the sample_paths_file, returns a dictionary of {sample_acc: (sample_r1_path, sample_r2_path)}.

	Arguments:
	- sample_paths_file: Path to file containing a list of paths to samples. 1 sample per line, with pattern _(1|2)\.(fastq|fq)(\.gz)?$
	
	Returns:
	- sample_dict: A dictionary of samples with sample name as key and (r1_path, r2_path) as values.
	"""

	# Read and check existence of every file
	sample_file_paths = []
	with open(sample_paths_file, "r") as f:
		sample_file_paths = [file_path.strip() for file_path in f if file_path.strip()]

	for sample_file_path in sample_file_paths:
		if not os.path.exists(sample_file_path):
			raise FileNotFoundError(f"Sample file not found or is not accessible: {sample_file_path}")
		
	
	# Pair sample files
	sample_dict = {}
	for sample_file_path in sample_file_paths:
		file_name_components = re.match(r'^(.+)_(1|2)\.(fastq|fq)(\.gz)?$', os.path.basename(sample_file_path))
		if file_name_components:
			sample = file_name_components.group(1)
			read_num = int(file_name_components.group(2))
			
			read_idx = read_num - 1
			if sample in sample_dict.keys() and sample_dict[sample][read_idx]:
				print(f"Duplicate R{read_num} detected for sample: {sample}. The first one will be used.")
			else:
				read_files = sample_dict.get(sample, ["", ""])
				read_files[read_idx] = sample_file_path
				sample_dict[sample] = read_files

	# Check pairs
	has_missing = False
	for sample, reads in sample_dict.items():
		r1, r2 = reads
		if not r1:
			print(f"Missing R1 for sample {sample} (R2: {r2})")
			has_missing = True
		if not r2:
			print(f"Missing R2 for sample {sample} (R1: {r1})")
			has_missing = True

	if has_missing:
		raise FileNotFoundError(f"One or more pairs were not found")

	return sample_dict

def set_wd(output_dir_path, sample_paths_file, src_fasta_path, src_database_path = None, softlink = False):
	"""
	Set output directory of the following structure:
		output_dir_path
		|-- data
		|   |-- db -> Location of VCF database files
		|   |-- fastq -> Location of FASTQ files
		|   `-- ref -> Location of reference FASTA file and index files
		`-- module
			|-- align
			|-- error
			|-- machine
			|-- model
			`-- variants
	This function also copies (or softlinks) source input files to appropriate destinations.

	Arguments:
	- output_dir_path: Path for root output directory.
	- sample_paths_file: Path to file containing a list of paths to samples. 1 sample per line, with pattern _(1|2)\.(fastq|fq)(\.gz)?$
	- src_fasta_path: Path to input FASTA file.
	- src_database_path: Path to input database VCF file (optional).

	Returns:
	- sample_dict: A dictionary of samples with sample name as key and (r1_path, r2_path) as values.
	- fasta_path: Path to the usable FASTA file (in data/ref).
	- db_path: Path to the usable database VCF file (in data/db).
	"""

	print(f"Setting up output directory: {output_dir_path}")

	# Create folders
	for data_dir in ["db", "fastq", "ref"]:
		os.makedirs(os.path.join(output_dir_path, "data", data_dir), exist_ok = True)
	
	for module_dir in ["align", "error", "machine", "model", "variants"]:
		os.makedirs(os.path.join(output_dir_path, "module", module_dir), exist_ok = True)

	# Copy/softlink files
	fasta_path = os.path.join(output_dir_path, "data", "ref", os.path.basename(src_fasta_path))
	put_file(src_fasta_path, fasta_path, softlink = softlink)
	
	db_path = None
	if src_database_path is not None:
		db_path = os.path.join(output_dir_path, "data", "db", os.path.basename(src_database_path))
		put_file(src_database_path, db_path, softlink = softlink)

	src_sample_dict = process_sample_list(sample_paths_file)
	sample_dict = {}

	for sample, files in src_sample_dict.items():
		src_r1, src_r2 = files
		
		dest_r1 = os.path.join(output_dir_path, "data", "fastq", os.path.basename(src_r1))
		dest_r2 = os.path.join(output_dir_path, "data", "fastq", os.path.basename(src_r2))

		put_file(src_r1, dest_r1, softlink = softlink)
		put_file(src_r2, dest_r2, softlink = softlink)
		
		sample_dict[sample] = [dest_r1, dest_r2]

	print("All directories and input files are ready.")

	return sample_dict, fasta_path, db_path

def pre_align(fasta_path):
	"""
	Generate index files for reference FASTA.

	Arguments:
	- fasta_path: Path to the usable FASTA file (in data/ref).
	"""
	
	ref_folder = os.path.dirname(fasta_path)

	fasta_dict_path = re.sub(r'\.(fa|fna|fasta)(.gz)?$', '.dict', fasta_path)
	file_exts = [
		f"{fasta_path}.0123",
		f"{fasta_path}.amb",
		f"{fasta_path}.ann",
		f"{fasta_path}.fai",
		f"{fasta_path}.bwt.2bit.64",
		f"{fasta_path}.pac",
		fasta_dict_path
	]

	# Check for missing reference index files
	missing_files = [f for f in file_exts if not os.path.exists(f)]
	if len(missing_files) == 0:
		print(f"All reference and index files are ready in {ref_folder}.")
		return

	# Start indexing if there are missing reference index files
	print(f"Missing reference index files:")
	for f in missing_files:
		print(f" - {f}")
	print(f"Start indexing...")

	# 1. delete existing dict file 
	if os.path.exists(fasta_dict_path):
		os.remove(fasta_dict_path)

	# 2. bwa-mem2 index
	subprocess.run(["bwa-mem2", "index", fasta_path], check=True)

	# 3. samtools faidx
	subprocess.run(["samtools", "faidx", fasta_path], check=True)

	# 4. picard CreateSequenceDictionary
	subprocess.run([
		"picard", "CreateSequenceDictionary",
		f"R={fasta_path}",
		f"O={fasta_dict_path}"
	], check=True)

	print("Reference preprocessing completed successfully.")

def align_fastq(sample_dict, fasta_path, output_dir_path, n_thread):
	"""
	Align FASTQ file of single samples to the reference.

	Arguments:
	- sample_dict: A dictionary of samples with sample name as key and (r1_path, r2_path) as values.
	- fasta_path: Path to FASTA file in data/ref.
	- output_dir_path: Path for root output directory.
	- n_thread: Number of threads passed to BWA MEM.
	"""

	module_align_dir = os.path.join(output_dir_path, "module", "align")

	for sample_name, reads_set in sample_dict.items():
		if os.path.exists(os.path.join(module_align_dir, f"{sample_name}_aligned.bam")):
			continue

		print(f"--- Processing Sample: {sample_name} ---")
	
		# mapping to reference
		r1, r2 = reads_set
		rg_header = f"@RG\\tID:{sample_name}\\tLB:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

		tmp_dir = os.path.join(module_align_dir, 'temp')
		init_sam_path = os.path.join(module_align_dir, f'{sample_name}_init.sam')
		sorted_sam_path = os.path.join(module_align_dir, f'{sample_name}_sorted.sam')
		aligned_bam_path = os.path.join(module_align_dir, f'{sample_name}_aligned.bam')
		metrics_txt_path = os.path.join(module_align_dir, f'{sample_name}_metrics.txt')

		cmd_align=(f"bwa-mem2 mem -M -t {n_thread} -R '{rg_header}' {fasta_path} {r1} {r2} > {init_sam_path}") 
		subprocess.run(cmd_align, shell=True, check=True)

		# Mark Duplicate and Sort
		os.makedirs(tmp_dir, exist_ok = True)
		subprocess.run(["picard","SortSam",f"I={init_sam_path}", f"TMP_DIR={tmp_dir}", \
				f"O={sorted_sam_path}", "SORT_ORDER=coordinate"],check=True)
	
		subprocess.run(["picard","MarkDuplicates", f"I={sorted_sam_path}", \
				f"O={aligned_bam_path}", \
				f"M={metrics_txt_path}","MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000","CREATE_INDEX=true"],check=True)
		
		if os.path.exists(init_sam_path):
			os.remove(init_sam_path)
	
		if os.path.exists(sorted_sam_path):
			os.remove(sorted_sam_path)

		if os.path.exists(metrics_txt_path):
			os.remove(metrics_txt_path)

def pseudo_db(species_name, sample_list, fasta_path, output_dir_path):
	"""
	Construct a pseudo database by using all samples in the align directory.

	Arguments:
	- species_name: Name of species that will be used in output file names.
	- sample_list: List of sample names.
	- fasta_path: Path to FASTA file in data/ref.
	- output_dir_path: Path for root output directory.

	Return:
	- output_vcf_path: Path to output pseudoDB VCF file.
	"""

	module_align_dir = os.path.join(output_dir_path, "module", "align")
	output_vcf_path = os.path.join(output_dir_path, "data", "db", f"{species_name}_pseudoDB.vcf.gz")

	sample_list_strs = []
	missing_sample_strs = []
	for sample in sample_list:
		sample_aligned_bam = os.path.join(module_align_dir, f"{sample}_aligned.bam")

		if os.path.exists(sample_aligned_bam):
			sample_list_strs.append(f"-I {sample_aligned_bam}")
		else:
			missing_sample_strs.append(f" - Aligned file for sample {sample} not found or is not accessible in path {sample_aligned_bam}")

	if len(missing_sample_strs) > 0:
		print(f"Missing BAM files:")
		for ms_str in missing_sample_strs:
			print(ms_str)
		raise FileNotFoundError(f"One or more aligned BAM files for pseudoDB construction are missing.")

	if len(sample_list_strs) == 0:
		raise ValueError("No samples to process.")
	
	# UnifiedGenotyper caller
	print(f"Start creating a pseudoDB with {len(sample_list_strs)} samples")
	vcf_cmd = f"gatk -T UnifiedGenotyper -R {fasta_path} {' '.join(sample_list_strs)} -o {output_vcf_path} --genotype_likelihoods_model BOTH"
	subprocess.run(vcf_cmd, shell=True, check=True)

	return output_vcf_path

def qs_recal(sample_list, fasta_path, db_name, db_path, output_dir_path):
	"""
	Recalibrate base quality score from samples.

	Arguments:
	- sample_list: List of sample names.
	- fasta_path: Path to FASTA file in data/ref.
	- db_name: Name of database VCF that will be used in output file names.
	- db_path: Path to the database VCF.
	- output_dir_path: Path for root output directory.
	"""

	module_align_dir = os.path.join(output_dir_path, "module", "align")
	module_machine_dir = os.path.join(output_dir_path, "module", "machine")

	if not os.path.exists(db_path) :
		raise FileNotFoundError(f"Database VCF file not found or is not accessible in path {db_path}")

	for sample in sample_list:
		aligned_bam_path = os.path.join(module_align_dir, f"{sample}_aligned.bam")
		
		recal_table_path = os.path.join(module_machine_dir, f"{sample}_{db_name}_recal.table")
		recalibrated_bam_path = os.path.join(module_machine_dir, f"{sample}_{db_name}_recalibrated.bam")

		if not os.path.exists(recalibrated_bam_path):
			if not os.path.exists(aligned_bam_path):
				print(f"Warning: Aligned BAM file missing for {sample} (expected path: {aligned_bam_path}). Skipping.")
			else:
				# BaseRecalibrator
				subprocess.run(f"gatk -T BaseRecalibrator -R {fasta_path} -I {aligned_bam_path}  -knownSites {db_path} -o {recal_table_path} &> {recal_table_path}.log", shell=True, check=True)

				# PrintReads
				subprocess.run(f"gatk -T PrintReads -R {fasta_path} -I {aligned_bam_path}  -BQSR {recal_table_path} -o {recalibrated_bam_path} &> {recalibrated_bam_path}.log", shell=True, check=True)

				# delete file
				os.remove(f"{recal_table_path}")
				os.remove(f"{recal_table_path}.log")
				os.remove(f"{recalibrated_bam_path}.log")

def variant_call(species_name, sample_list, fasta_path, db_name, output_dir_path):
	"""
	Call genetic variants - step 1 : Call variants from each samples.

	Arguments:
	- species_name: Name of species that will be used in output file names.
	- sample_list: list of sample names
	- fasta_path: Path to FASTA file in data/ref.
	- db_name: Name of database VCF that will be used in output file names.
	- output_dir_path: Path for root output directory.

	Return:
	- output_vcf: Path to resulting VCF file from variant calling.
	"""

	module_machine_dir = os.path.join(output_dir_path, "module", "machine")
	module_variants_dir = os.path.join(output_dir_path, "module", "variants")

	sample_list_strs = []
	missing_sample_strs = []
	for sample in sample_list:
		recalibrated_bam_path = os.path.join(module_machine_dir, f"{sample}_{db_name}_recalibrated.bam")

		if os.path.exists(recalibrated_bam_path) :
			sample_list_strs.append(f"-I {recalibrated_bam_path}")
		else:
			missing_sample_strs.append(f" - Recalibrated BAM file for sample {sample} not found or is not accessible in path {recalibrated_bam_path}")
		

	if len(missing_sample_strs) > 0:
		print(f"Missing recalibrated BAM files:")
		for ms_str in missing_sample_strs:
			print(ms_str)
		raise FileNotFoundError(f"One or more recalibrated BAM files are missing.")
	
	# UnifiedGenotyper caller
	output_vcf = os.path.join(module_variants_dir, f"{species_name}_{db_name}_variant_calling.vcf.gz")

	print(f"Start genetic variants calling with {len(sample_list_strs)} samples")
	vcf_cmd = f"gatk -T UnifiedGenotyper -R {fasta_path} {' '.join(sample_list_strs)} -o {output_vcf} --genotype_likelihoods_model BOTH &> {output_vcf}.log"
	subprocess.run(vcf_cmd, shell=True, check=True)
	
	# delete file
	os.remove(f"{output_vcf}.log")

	return output_vcf

def error_rate(species_name, sample_list, fasta_path, db_name, db_path, output_dir_path):
	"""
	Estimate sample error rate.

	Arguments:
	- species_name: Name of species that will be used in output file names.
	- sample_list: list of sample names.
	- fasta_path: Path to FASTA file in data/ref.
	- db_name: Name of database VCF that will be used in output file names.
	- db_path: Path to the database VCF.
	- output_dir_path: Path for root output directory.

	Return:
	- module_error_dir: Path to directory containing error rate calculation result.
	"""
	
	data_db_dir = os.path.join(output_dir_path, "data", "db")
	module_align_dir = os.path.join(output_dir_path, "module", "align")
	module_error_dir = os.path.join(output_dir_path, "module", "error")

	## database check
	if not os.path.exists(db_path):
		raise FileNotFoundError(f"Database VCF file not found or is not accessible in path {db_path}")
	if os.path.splitext(db_path)[-1] == ".gz":
		db_path = os.path.splitext(db_path)[0]
		subprocess.run(f"zcat {db_path}.gz > {db_path}", shell=True, check=True)

	db_uniq_check = None
	for sample in sample_list:
		aligned_bam_path = os.path.join(module_align_dir, f"{sample}_aligned.bam")

		if not os.path.exists(aligned_bam_path) :
			print(f"Warning: Aligned BAM file missing for {sample} (expected location: {aligned_bam_path})")
			continue
		
		sample_error_path = os.path.join(module_error_dir, f"{sample}_error")
		sample_error_analysis_path = os.path.join(module_error_dir, f"{sample}_error_analysis")

		os.system(f"samtools mpileup -Bf {fasta_path} {aligned_bam_path} > {sample_error_path}")

		infile=open(sample_error_path, "r") # mpileup output file load	
		outfile=open(sample_error_analysis_path, "w")

		line=infile.readline()
		line_list=line.strip().split("\t")

		while line !="" :
			if line_list[3]!="0" :
				d=line_list[4].find("^")   # start of read segment 
				while d !=-1 :
					line_list[4]=line_list[4].replace(line_list[4][d:d+2],"")
					d=line_list[4].find("^")

				line_list[4]=line_list[4].replace("$","")   # end of a read segment
				line_list[4]=line_list[4].replace("*","")   #
				line_list[4]=line_list[4].replace(".","")   # match to the refernece base on the forward strand
				line_list[4]=line_list[4].replace(",","")   # match to the reference base on the reverse strand

				if line_list[4]!="" :
					indelnum=0
					indelnum=indelnum+line_list[4].count("+")   # insertion from the reference
					indelnum=indelnum+line_list[4].count("-")   # deletion from the reference
					tmpgeno=line_list[4]
					i=tmpgeno.find("+")
					while i!=-1 :
						if tmpgeno[i+1:i+3].isdigit()==True :
							n=int(tmpgeno[i+1:i+3])
							tmpgeno=tmpgeno.replace(tmpgeno[i:i+3+n],"")
						else :
							n=int(tmpgeno[i+1:i+2])
							tmpgeno=tmpgeno.replace(tmpgeno[i:i+2+n],"")
						i=tmpgeno.find("+")
					i=tmpgeno.find("-")
					while i!=-1 :
						if tmpgeno[i+1:i+3].isdigit()==True :
							n=int(tmpgeno[i+1:i+3])
							tmpgeno=tmpgeno.replace(tmpgeno[i:i+3+n],"")         
						else :
							n=int(tmpgeno[i+1:i+2])
							tmpgeno=tmpgeno.replace(tmpgeno[i:i+2+n],"")
						i=tmpgeno.find("-")           
					mnum=len(tmpgeno)+indelnum
					outfile.write(f"{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{mnum}\t{line_list[4]}\n")
			line=infile.readline()
			line_list=line.strip().split("\t")

		os.remove(sample_error_path)
		infile.close()
		outfile.close()

		## database unique position check
		db_uniq_check = os.path.join(data_db_dir, f"{species_name}_{db_name}_uniq_pos")
		sample_uniq_check = os.path.join(module_error_dir, f"{sample}_error_analysis_uniq_pos")
		sample_db_analysis = os.path.join(module_error_dir, f"{sample}_{db_name}_analysis")
		sample_db_common = os.path.join(module_error_dir, f"{sample}_{db_name}_common")
		sample_db_variant_pos = os.path.join(module_error_dir, f"{sample}_{db_name}_variant_pos")
		error_rate_file = os.path.join(module_error_dir, f"{sample}_{db_name}_erate")
		mismatch_name = os.path.join(module_error_dir, f"{sample}_mismatch") 

		if not os.path.exists(db_uniq_check):
			snp_extract=f'grep -v "^#" {db_path}  | cut -f1,2 | uniq > {db_uniq_check}' # database uniq position search
			os.system(snp_extract)

		if not os.path.exists(sample_uniq_check):
			sample_extract=f"cut -f1,2 {sample_error_analysis_path} > {sample_uniq_check}" # sample uniq position search
			os.system(sample_extract)

		sdiff_exe=f"sdiff {db_uniq_check} {sample_uniq_check} > {sample_db_analysis}" ## database and sample analysis file 
		os.system(sdiff_exe)

		os.remove(sample_uniq_check)

		awk_cmd="awk '{if(NF==4) print $0;}'"
		sdiff_extract=f"{awk_cmd} {sample_db_analysis} > {sample_db_common}"
		os.system(sdiff_extract)

		os.remove(sample_db_analysis)

		eff_variant=f"cut -f1,2 {sample_db_common} > {sample_db_variant_pos}"
		os.system(eff_variant)
   
		os.remove(sample_db_common)
	
		eff_name=sample_db_variant_pos

		sample_infile=open(sample_error_analysis_path,"r")
		eff_infile=open(eff_name,"r")

		error_rate_handle=open(error_rate_file,"w")

		mismatch_str="awk '{ sum+=$5} END { print sum;}'"
		mismatch_cmd=f"{mismatch_str} {sample_error_analysis_path} > {mismatch_name}"
		os.system(mismatch_cmd)

		mismatch_infile=open(mismatch_name,"r")
		mismatch_num=int(mismatch_infile.readline())

		eff_num=0
		while True :
			eff_base=eff_infile.readline()

			if eff_base=="" :
				break

			eff_list=eff_base.strip().split("\t")
		
			while True :
				base_sample=sample_infile.readline()

				if base_sample=="" :
					break

				base_list=base_sample.split('\t') 
				if eff_list[0]==base_list[0] and eff_list[1]==base_list[1] :
					eff_num=eff_num+int(base_list[4])
					break

		error_rate_handle.write(f"{sample}\t{(mismatch_num-eff_num)/mismatch_num}")

		os.remove(sample_db_variant_pos)
	
		os.remove(mismatch_name)
	
		os.remove(sample_error_analysis_path)
	
		mismatch_infile.close()
		sample_infile.close()
		eff_infile.close()
		error_rate_handle.close()

	if db_uniq_check is not None:
		os.remove(db_uniq_check)

	# os.remove(f"{species}/data/db/{species_name}_{dbtype}.vcf") # Maybe don't delete the VCF file?

	return module_error_dir

def qs_model(sample_list, db_name, output_dir_path):
	"""
	Estimate model-adjusted base quality score.

	Arguments:
	- sample_list: list of sample names
	- db_name: Name of database VCF that will be used in output file names.
	- output_dir_path: Path for root output directory.

	Return:
	- module_model_dir: Path to directory containing estimation of model-adjusted base quality score.
	"""

	module_machine_dir = os.path.join(output_dir_path, "module", "machine")
	module_model_dir = os.path.join(output_dir_path, "module", "model")

	for sample in sample_list:
		recalibrated_bam_path = os.path.join(module_machine_dir, f"{sample}_{db_name}_recalibrated.bam")
		recalibrated_sam_path = os.path.join(module_model_dir, f"{sample}_{db_name}_recalibrated.sam")

		if not os.path.exists(recalibrated_bam_path) :
			print(f"Warning: Recalibrated BAM file missing for {sample}. Skipping.")
			continue

		os.system(f"samtools view -h {recalibrated_bam_path} > {recalibrated_sam_path}")

		sample_infile=open(recalibrated_sam_path,"r")
	
		q_count=[]
		for i in range(100) :
			q_count.append(0)

		line=sample_infile.readline()

		while line[0]=="@" :
			line=sample_infile.readline()

		while line!="" :
			line_list=line.strip().split("\t")
			i=0
			while i < len(line_list[10]) :
				qscore=ord(line_list[10][i])-33
				q_count[qscore]=q_count[qscore]+1
				i=i+1
			line=sample_infile.readline()
   
	
		sample_infile.close()
	
		sample_outname = os.path.join(module_model_dir, f"{sample}_{db_name}_qs")
		sample_outfile=open(sample_outname,"w")

		hap=0
		hhap=0
  
		for i in range(len(q_count)) :
			hap=hap+q_count[i]
			hhap=hhap+i*q_count[i]
	
		sample_outfile.write(f"{sample}\t{hhap/hap}")
		sample_outfile.close()

		os.remove(recalibrated_sam_path)

	return module_model_dir

def main(species_name, src_fasta_path, sample_paths_file, output_dir_path, src_database_path = None, db_name = None, n_thread = 1, softlink = False):
	sample_dict, fasta_path, db_path = set_wd(
											output_dir_path = output_dir_path,
											sample_paths_file = sample_paths_file,
											src_fasta_path = src_fasta_path,
											src_database_path = src_database_path,
											softlink = softlink
										)
	sample_list = list(sample_dict.keys())

	pre_align(fasta_path = fasta_path)
	align_fastq(
		sample_dict = sample_dict,
		fasta_path = fasta_path,
		output_dir_path = output_dir_path,
		n_thread = n_thread
	)
	
	if src_database_path is None:
		pseudo_db(
			species_name = species_name,
			sample_list = sample_list,
			fasta_path = fasta_path,
			output_dir_path = output_dir_path
		)
	else:
		db_name = db_name if db_name is not None else re.sub(r'\.(vcf|bcf)(\.gz)?$', '', os.path.basename(src_database_path))

		qs_recal(
			sample_list = sample_list,
			fasta_path = fasta_path,
			db_name = db_name,
			db_path = db_path,
			output_dir_path = output_dir_path
		)

		vc_vcf_path = variant_call(
			species_name = species_name,
			sample_list = sample_list,
			fasta_path = fasta_path,
			db_name = db_name,
			output_dir_path = output_dir_path
		)
		print(f"VCF result of variant-calling saved: {vc_vcf_path}")

		err_rate_dir = error_rate(
			species_name = species_name,
			sample_list = sample_list,
			fasta_path = fasta_path,
			db_name = db_name,
			db_path = db_path,
			output_dir_path = output_dir_path
		)
		print(f"Error rate saved in directory: {err_rate_dir}")

		qs_model_dir = qs_model(
			sample_list = sample_list,
			db_name = db_name,
			output_dir_path = output_dir_path
		)
		print(f"Model-adjusted quality scores saved in directory: {qs_model_dir}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-sp", "--species",			required=True,		type=str,	help=f"Target species name for output file naming.")
	parser.add_argument("-fa", "--fasta",			required=True,		type=str,	help=f"Path to reference FASTA file.")
	parser.add_argument("-s",  "--sample-list",		required=True,		type=str,	help=f"Path to file containing sample FASTQ paths.")
	parser.add_argument("-o",  "--output-dir",		required=True,		type=str,	help=f"Path to output directory.")
	parser.add_argument("-db", "--database",		default=None,		type=str,	help=f"Path to database VCF to use.")
	parser.add_argument("-dn", "--database-name",	default=None,		type=str,	help=f"Name of database for output file naming.")
	parser.add_argument("-t",  "--threads",			default=1,			type=int,	help=f"Number of CPU threads to use [1].")
	parser.add_argument("-sl", "--softlink",		action="store_true",			help=f"Create softlinks instead of copying source input files into output directory.")
	args = parser.parse_args()

	main(
		species_name = args.species,
		src_fasta_path = args.fasta,
		sample_paths_file = args.sample_list,
		output_dir_path = args.output_dir,
		src_database_path = args.database,
		db_name = args.database_name,
		n_thread = args.threads,
		softlink = args.softlink
	)