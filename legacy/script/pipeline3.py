import os
import sys
import subprocess
import argparse

# working directory 
def set_wd(species) :
    path_dir=os.path.join(species,"module")
    # module
    # align : result of aligning FASTQ to reference resulting BAM
    # machine : result of recalibrating maching-provided base quality score
    # error : result of estimating sample error rate
    # model : result of estimating model-adjusted base quality score
    # variants :  result of genetic variants
    sub_dirs=["align","machine","error","model","variants"]
    
    print(f"\nSetting up directories for {species}...")
    for sub in sub_dirs :
        full_path=os.path.join(path_dir,sub)
        os.makedirs(full_path, exist_ok=True)

    print("All directories are ready.")

# end of set_wd()


# reference section
def pre_align(species, reference_file) :
    path_dir=f"{species}/data/ref"
    file_list=os.listdir(path_dir)

    file_exts=[f"{reference_file}.0123",f"{reference_file}.amb",f"{reference_file}.ann",f"{reference_file}.fai",\
            f"{reference_file}.amb",f"{reference_file}.bwt.2bit.64",f"{reference_file}.pac",f"{reference_file[:reference_file.find('.f')]}.dict"]

    missing_files=[]

    for f in file_exts :
        full_path=os.path.join(path_dir,f)
        if not os.path.exists(full_path) :
            missing_files.append(f)

    if not missing_files :
        print(f"\nAll reference files for {species} are ready.")
        return

    print(f"\nMissing {missing_files}: Starting indexing...")

    ref_path=f"{species}/data/ref/{reference_file}"
    dict_file = f"{species}/data/ref/{reference_file[:reference_file.find('.f')]}.dict"
    
    # 1. delete existing dict file 
    if os.path.exists(dict_file):
        os.remove(dict_file)

    # 2. bwa-mem2 index
    subprocess.run(["bwa-mem2", "index", ref_path], check=True)

    # 3. samtools faidx
    subprocess.run(["samtools", "faidx", ref_path], check=True)

    # 4. picard CreateSequenceDictionary
    subprocess.run([
        "picard", "CreateSequenceDictionary",
        f"R={ref_path}",
        f"O={dict_file}"
    ], check=True)

    print("Reference preprocessing completed successfully.")
   
# end of pre_align()


# Align FASTQ file of single samples to the reference
def align_fastq(species, reference_file, n_thread, file_list) :
    
    with open(file_list, "r") as f:
        sample_list = [line.strip() for line in f if line.strip()]

    exists_samples=[]
    missing_samples=[]

    for sample in sample_list :
        bam_path=f"{species}/module/align/{sample}_aligned.bam"

        if  os.path.exists(bam_path) :
            exists_samples.append(sample)
        else :
            missing_samples.append(sample)

    for sample_name in missing_samples :         
        print(f"\n--- Processing Sample: {sample_name} ---")
    
        # mapping to reference
        r1 = f"{species}/data/fastq/{sample_name}_1.fastq.gz"
        r2 = f"{species}/data/fastq/{sample_name}_2.fastq.gz" 
        ref_path=f"{species}/data/ref/{reference_file}"
        rg_header = f"@RG\\tID:{sample_name}\\tLB:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

        cmd_align=(f"bwa-mem2 mem -M -t {n_thread} -R '{rg_header}' {ref_path} {r1} {r2} > {species}/module/align/{sample_name}_init.sam")
        subprocess.run(cmd_align, shell=True, check=True)

        # Mark Duplicate and Sort
        subprocess.run(["picard","SortSam",f"I={species}/module/align/{sample_name}_init.sam", "TMP_DIR=temp",  \
                f"O={species}/module/align/{sample_name}_sorted.sam", "SORT_ORDER=coordinate"],check=True)
    
        if os.path.exists(f"{species}/module/align/{sample_name}_init.sam"):
            os.remove(f"{species}/module/align/{sample_name}_init.sam")
    
        subprocess.run(["picard","MarkDuplicates", f"I={species}/module/align/{sample_name}_sorted.sam", \
                f"O={species}/module/align/{sample_name}_aligned.bam",  \
                f"M={species}/module/align/{sample_name}_metrics.txt","MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000","CREATE_INDEX=true"],check=True)
    
        if os.path.exists(f"{species}/module/align/{sample_name}_sorted.sam"):
            os.remove(f"{species}/module/align/{sample_name}_sorted.sam")

        if os.path.exists(f"{species}/module/align/{sample_name}_metrics.txt"):
            os.remove(f"{species}/module/align/{sample_name}_metrics.txt")


# end of align_fastq()


# Construct a pseudo database by using all samples in the align directory
def pseudo_db(species, reference_file, file_list):

    with open(file_list, "r") as f:
        sample_name = [line.strip() for line in f if line.strip()]
    
    missing_samples=[]
    sample_list_str=""

    for sample in sample_name :
        bam_path=f"{species}/module/align/{sample}_aligned.bam"

        if os.path.exists(bam_path) :
            sample_list_str +=f"-I {bam_path} "
        else :
            missing_samples.append(sample)

    if missing_samples :
        print(f"\nWarning: Bam files for the following samples are missing: {missing_samples}")
        return

    # UnifiedGenotyper caller
    species_name=species.split('/')[-1] 
    output_vcf=f"{species}/data/db/{species_name}_pseudoDB.vcf.gz"
    vcf_cmd=f"gatk -T UnifiedGenotyper -R {species}/data/ref/{reference_file} {sample_list_str} -o {output_vcf} --genotype_likelihoods_model BOTH"
    print(f"Start creating a psesudoDB with {len(sample_name)} samples")
    os.system(vcf_cmd)

# end of pseudo_db()


# Recalibrate base quality score from samples
def qs_recal(species, reference_file, dbtype, file_list) :

    with open(file_list, "r") as f:
        sample_list = [line.strip() for line in f if line.strip()]
    
    exists_samples=[]
    missing_samples=[]

    for sample in sample_list :
        bam_path=f"{species}/module/machine/{sample}_{dbtype}_recalibrated.bam"

        if os.path.exists(bam_path) :
            exists_samples.append(sample)
        else :
            missing_samples.append(sample)
 
    species_name=species.split('/')[-1] 
    db_path=f"{species}/data/db/{species_name}_{dbtype}.vcf.gz"
    if not os.path.exists(db_path) :
        sys.exit(f"Not found database : {db_path}")		
    

    # run each sample
    for sample_name in missing_samples :
        bam_path=f"{species}/module/align/{sample_name}_aligned.bam"

        if not os.path.exists(bam_path) :
            print(f"\nWarning: Bam files missing for {sample_name}")
            continue
    
        machine_dir = f"{species}/module/machine"
        ref_path = f"{species}/data/ref/{reference_file}"
        recal_table = f"{machine_dir}/{sample_name}_{dbtype}_recal.table"
        output_bam = f"{machine_dir}/{sample_name}_{dbtype}_recalibrated.bam" 
  
        # BaseRecalibrator
        os.system(f"gatk -T BaseRecalibrator -R {ref_path} -I {bam_path}  -knownSites {db_path} -o {recal_table} &> {recal_table}.log")

        # PrintReads
        os.system(f"gatk -T PrintReads -R {ref_path} -I {bam_path}  -BQSR {recal_table} -o {output_bam} &> {output_bam}.log")

        # delete file
        os.remove(f"{recal_table}")
        os.remove(f"{recal_table}.log")
        os.remove(f"{output_bam}.log")

# end of qs_recal()

# Call genetic variants - step 1 : Call variants from each samples. 
def variant_call(species, reference_file, dbtype, file_list):
    
    with open(file_list, "r") as f:
        sample_name = [line.strip() for line in f if line.strip()]

    missing_samples=[]
    sample_list_str=""

    for sample in sample_name :
        bam_path=f"{species}/module/machine/{sample}_{dbtype}_recalibrated.bam"

        if os.path.exists(bam_path) :
            sample_list_str +=f"-I {bam_path} "
        else :
            missing_samples.append(f"{sample}_{dbtype}_recalibrated.bam")

    if missing_samples :
        print(f"\nWarning: Bam files for the following samples are missing: {missing_samples}")
        return

    # UnifiedGenotyper caller
    species_name=species.split('/')[-1] 
    output_vcf=f"{species}/module/variants/{species_name}_{dbtype}_variant_calling.vcf.gz"
    ref_path = f"{species}/data/ref/{reference_file}"

    vcf_cmd=f"gatk -T UnifiedGenotyper -R {ref_path} {sample_list_str} -o {output_vcf} --genotype_likelihoods_model BOTH &> {output_vcf}.log"

    print(f"\nStart genetic variants calling  with {len(sample_name)} samples")
    os.system(vcf_cmd)

    # delete file
    os.remove(f"{output_vcf}.log")

# end of variant_call()


# Estimate sample error rate
def error_rate(species, reference_file, dbtype, file_list) :
    
    with open(file_list, "r") as f:
        sample_list1 = [line.strip() for line in f if line.strip()]

    species_name=species.split('/')[-1]
    database=f"{species_name}_{dbtype}.vcf.gz"

    ## database check
 
    db_path=f"{species}/data/db/{species_name}_{dbtype}.vcf.gz"
    if not os.path.exists(db_path) :
        sys.exit(f"Not found database : {db_path}")		
    
   
    if database in db_path and ".gz" in database :
        os.system(f"zcat {db_path} > {species}/data/db/{species_name}_{dbtype}.vcf")
   
    database=database[:database.find(".gz")]

    ref_path = f"{species}/data/ref/{reference_file}"
    for sample in sample_list1 :
        bam_path=f"{species}/module/align/{sample}_aligned.bam"

        if not os.path.exists(bam_path) :
            print(f"Warning: Bam files missing for {sample}")
            continue
    
        os.system(f"samtools mpileup -Bf {ref_path} {bam_path} > {species}/module/error/{sample}_error\n")
    
        infile_name=f"{species}/module/error/{sample}_error"  # mileup output file load
        infile=open(infile_name,"r")

    
        outfile_name=f"{species}/module/error/{sample}_error_analysis"
        outfile=open(outfile_name,"w")

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

        os.remove(f"{species}/module/error/{sample}_error")
        infile.close()
        outfile.close()
    
	
   
        ## database unique position check
        db_name=f"{species}/data/db/{database}"
        db_uniq_check=f"{species}/data/db/{species_name}_{dbtype}_uniq_pos"

        db_dir=f"{species}/data/db"
        db_list=os.listdir(db_dir)
        if db_uniq_check not in db_list :
            snp_extract=f'grep -v "^#" {db_name}  | cut -f1,2 | uniq > {species}/data/db/{species_name}_{dbtype}_uniq_pos'    # database uniq position search
            os.system(snp_extract)
    
        sample_uniq_check=f"{species}/module/error/{sample}_error_analysis_uniq_pos"
        sample_dir=f"{species}/module/error"
        sample_list=os.listdir(sample_dir)

        if sample_uniq_check not in sample_list :
            sample_name=f"{species}/module/error/{sample}_error_analysis"
            sample_extract=f"cut -f1,2 {sample_name} > {species}/module/error/{sample}_error_analysis_uniq_pos"     # sample uniq position search
            os.system(sample_extract)

        sdiff_exe=f"sdiff {species}/data/db/{species_name}_{dbtype}_uniq_pos  {species}/module/error/{sample}_error_analysis_uniq_pos \
                > {species}/module/error/{sample}_{dbtype}_analysis"   ## database and sample analysis file 
        os.system(sdiff_exe)

        os.remove(f"{species}/module/error/{sample}_error_analysis_uniq_pos")

        awk_cmd="awk '{if(NF==4) print $0;}'"
        sdiff_extract=f"{awk_cmd} {species}/module/error/{sample}_{dbtype}_analysis > {species}/module/error/{sample}_{dbtype}_common"
        os.system(sdiff_extract)

        os.remove(f"{species}/module/error/{sample}_{dbtype}_analysis")

        eff_variant=f"cut -f1,2 {species}/module/error/{sample}_{dbtype}_common  > {species}/module/error/{sample}_{dbtype}_variant_pos"
        os.system(eff_variant)
   
        os.remove(f"{species}/module/error/{sample}_{dbtype}_common")
    
        sample_name=f"{species}/module/error/{sample}_error_analysis"
        eff_name=f"{species}/module/error/{sample}_{dbtype}_variant_pos"

        sample_infile=open(sample_name,"r")
        eff_infile=open(eff_name,"r")

        error_rate_file=f"{species}/module/error/{sample}_{dbtype}_erate"
        error_rate=open(error_rate_file,"w")

        mismatch_str="awk '{ sum+=$5} END { print sum;}'"
        mismatch_cmd=f"{mismatch_str} {sample_name} > {species}/module/error/{sample}_mismatch"
        os.system(mismatch_cmd)

        mismatch_name=f"{species}/module/error/{sample}_mismatch"

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

        error_rate.write(f"{sample}\t{(mismatch_num-eff_num)/mismatch_num}")

        os.remove(f"{species}/module/error/{sample}_{dbtype}_variant_pos")
    
        os.remove(f"{species}/module/error/{sample}_mismatch")
    
        os.remove(f"{species}/module/error/{sample}_error_analysis")
    
	
        sample_infile.close()
        eff_infile.close()
        error_rate.close()

    os.remove(f"{species}/data/db/{species_name}_{dbtype}_uniq_pos")

    os.remove(f"{species}/data/db/{species_name}_{dbtype}.vcf")
	
# end of error_rate()

# Estimate model-adjusted base quality score.
def qs_model(species, dbtype, file_list) :

    with open(file_list, "r") as f:
        sample_list = [line.strip() for line in f if line.strip()]

    for sample_name in sample_list :
        bam_path=f"{species}/module/machine/{sample_name}_{dbtype}_recalibrated.bam"
        sam_path=f"{species}/module/model/{sample_name}_{dbtype}_recalibrated.sam"
        if not os.path.exists(bam_path) :
            print(f"\nWarning: Bam files missing for {sample_name}")
            continue

        os.system(f"samtools view -h {bam_path} > {sam_path}")
     
        sample_file=f"{sam_path}"

        sample_infile=open(sample_file,"r")
    
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
	
        sample_outname=f"{species}/module/model/{sample_name}_{dbtype}_qs"
        sample_outfile=open(sample_outname,"w")

        hap=0
        hhap=0
  
        for i in range(len(q_count)) :
            hap=hap+q_count[i]
            hhap=hhap+i*q_count[i]
    
        sample_outfile.write(f"{sample_name}\t{hhap/hap}")
        sample_outfile.close()

        os.remove(f"{sam_path}")

# end of qs_model()


def main() :
    #os.system("curl -L -O https://raw.githubusercontent.com/infoLab204/pseudoDB/refs/heads/main/pipeline/pipeline3.py")

    """
    parser = argparse.ArgumentParser(usage='python %(prog)s [options]', prog='pipeline3.py' )

    parser.add_argument('--species', '-sp', help='the target speices name', required=True)
    parser.add_argument('--reference', '-ref', help='the filename of the reference sequence file', required=True)
    parser.add_argument('--sample', '-s', help='the filename of the sample data', required=True)
    parser.add_argument('--database', '-db', help='the specific database of known variants to be used', required=True)
    parser.add_argument('--thread', '-nt', type=int, help='the number of CPU threads to allocate for the process', required=True)

    args = parser.parse_args()

    species=args.species # the target species name
    ref=args.reference  # the filename of the reference sequence file
    sample=args.sample  # the filename of the FASTQ data
    dbtype=args.database  # the specific database of known variants to be used
    thread=args.thread  # the number of CPU threads to allocate for the process
    """
    if len(sys.argv) < 6 : 
        print("\n[Error] Insufficient arguments.")
        print("python pipeline3.py <species> <ref> <threads> <db_type> <list_of_sample_files>")
        print("=" * 100)
        print("""
        Case 1: Generating a New Pseudo-Database
        python pipeline3.py human GRCh38_full_analysis_set_plus_decoy_hla.fa 16 null list_human  
        Note : You will need to prepare a text file titled list_human. 
               Inside this file, list the sample names you wish to process, separated by enterkey (e.g., HG00096\n HG00097\n HG00098\n).
        Note : Once the process is complete, you can find your output VCF files in the db/ directory.   

        Case 2: Variant Calling with a dbSNP
        python pipeline3.py human GRCh38_full_analysis_set_plus_decoy_hla.fa 16 dbSNP list_human 
        Note : Once the process is complete, you can find your output VCF files in the variants/ directory.   
       
        Case 3: Variant Calling with a Pseudo-Database
        python pipeline3.py human GRCh38_full_analysis_set_plus_decoy_hla.fa 16 pseudoDB list_human
        Note : Once the process is complete, you can find your output VCF files in the variants/ directory.    
        """)
        sys.exit(1)

    #working_dict=sys.argv[1] # working directory
    species=sys.argv[1] # the target species name
    ref=sys.argv[2]  # the filename of the reference sequence file
    thread=int(sys.argv[3]) # the number of CPU threads to allocate for the process
    dbtype=sys.argv[4]  # the specific database of known variants to be used

    if dbtype.upper()=="NULL" :
        pseudo_file=sys.argv[5]  # the filename of construct a pseudo-database using  pseudo_sample list
    else :
        analysis_file=sys.argv[5]  # the filename of construct a pseudo-database using  pseudo_sample list

    species=os.getcwd()
    os.chdir(species)
     
    # Create subdirectories under directory "module".
    set_wd(species)

    # Create file names for the alignment under directory "ref".
    pre_align(species, ref)
    
    if dbtype.upper()=="NULL" :  # Generating a New Pseudo-Database
        # Align FASTQ file of single samples to the reference.
        align_fastq(species, ref, thread, pseudo_file)
        
        # Construct a pseudo-database using pseudo_sample list
        pseudo_db(species, ref, pseudo_file ) 
    else :  
        # Align FASTQ file of single samples to the reference.
        align_fastq(species, ref, thread, analysis_file)

        # Recalibrate base quality score from sample.
        qs_recal(species, ref, dbtype, analysis_file)
        
        # Call genetic variants.
        variant_call(species, ref, dbtype, analysis_file)

        # Estimate sample error rate.
        error_rate(species,ref, dbtype, analysis_file)

        # Estimate model-adjusted base quality score.
        qs_model(species, dbtype, analysis_file)
    
    
# end of main()

if __name__ == "__main__":
    main()

