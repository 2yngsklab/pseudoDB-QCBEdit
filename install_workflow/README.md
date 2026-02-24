# Installation and Workflow Guide
This tutorial is organized into two sections: tool installation and the genetic variant calling workflow. The first section guides you through using Conda to install essential software, such as BWA, Samtools, Picard, and GATK. The second section provides how to prepare the list of parameters for the workflow script. The resulting directory structure is highlighted in Fig. 1. In the following examples:    
-	Lines starting with # are comments providing context.
-	Commands inside code blocks are executable.
  
Note: Ensure all required data has been fully downloaded before proceeding.

<img width="800" height="350" alt="image" src="https://github.com/user-attachments/assets/7ebf5eb5-5dc2-4506-ab0a-8ed463fdefb2" />  <br>
*Fig. 1 : The overall structure of the directories.*
<br><br>
## Section1: Installation Procedure
This section provides a step-by-step guide for installing the necessary bioinformatics tools, including BWA, Samtools, Picard, and GATK v3 (UnifiedGenotyper) using the Conda. To find the specific version of Conda compatible with your system, visit the Anaconda Archive at https://repo.anaconda.com/archive/

Note: To find the specific version of Conda compatible with your system, visit the Anaconda Archive at https://repo.anaconda.com/archive/

### 1. Install Conda and Update Environment
Download the installer, run the installation script, and refresh your shell environment variables:    

- #### Download the Anaconda installer    
```   
wget https://repo.anaconda.com/archive/Anaconda3-2025.06-0-Linux-x86_64.sh     
```

- #### Execute the installation script    
```   
bash Anaconda3-2025.06-0-Linux-x86_64.sh    
```

- #### Reload environmental variables    
```   
source ~/.bashrc
```

### 2. Tool Installation
This workflow is for pipelines requiring the GATK3 UnifiedGenotyper. Note that GATK3 requires manual registration due to licensing.    
   
```   
conda create -n gatk3        # Create a dedicated virtual environment    
```
```   
conda activate gatk3         # Activate the environment    
```
```   
conda install bioconda::bwa-mem2    
```
```   
conda install bioconda::samtools=1.13    
```
```   
conda install bioconda::picard=2.0    
```
```   
conda install bioconda::gatk=3.8    
```   
Manually download and register GATK 3.8    
```  
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2     
```
```  
tar -jxvf  GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 
```
```  
mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk3.8    
```
```  
gatk-register gatk3.8/GenomeAnalysisTK.jar    
```
```  
rm -rf  GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2  gatk3.8 
```
<br><br>
## Section2: Variant Calling Workflow    
This section details the execution of the variant calling pipeline. The workflow requires the following input parameters:  
-	**species**: The target species name (refer to the species list in the dataset directory).
-	**ref**: The filename of the reference sequence.   
-	**thread**: The number of CPU threads to allocate for the process.
- **db_type**: The specific database of known variants to be used (see "Note" below).     
-	**sample**: The filename of the FASTQ data.    <br><br>

Note on db_type Values:    
    1. **null** : Use this when constructing a new pseudo-database.    
    2. **\<dbSNP Name\>**: Use this when calling variants using an existing dbSNP.    
    3. **\<pseudoDB Name\>**: Use this when calling variants using a previously generated pseudo-database.    <br><br>
Note: Certain species may lack established dbSNP resources.

### GATK3 Workflow
Copy the pipeline3.py module from the **script/** repository to your **default** working directory (e.g., human).         
 
```   
curl -L -O https://raw.githubusercontent.com/infoLab204/pseudoDB/refs/heads/main/script/pipeline3.py  # download "pipeline3.py" module   
```   

```   
conda activate gatk3         # Activate the environment    
```

### Running the Pipeline
Navigate to your default working directory (e.g., human) and execute the variant calling functions by following this command structure:        
#### python pipeline3.py \<species\> \<ref\> \<threads\> \<db_type\> \<list_of_sample_files>   <br><br>

We have provided three example use cases based on the human dataset to guide you.        
  
- ##### Case 1: Generating a new Pseudo-Database    
    ```  
    python pipeline3.py human GRCh38_full_analysis_set_plus_decoy_hla.fa 16 null list_human    
    ```
 Note: You will need to prepare a text file titled list_human. Inside this file, list the sample names you wish to process, separated by spaces (e.g., HG00096 HG00097 HG00098).    
 Note: Once the process is complete, you can find your output VCF files in the **db/** directory.    
 
- ##### Case 2: Variant Calling with an dbSNP    
    ```  
    python pipeline3.py human GRCh38_full_analysis_set_plus_decoy_hla.fa 16 dbSNP list_human      
    ```
     	Note: Once the process is complete, you can find your output VCF files in the **variants/** directory.   
- ##### Case 3: Variant Calling with a Pseudo-Database     
    ```  
    python pipeline3.py human GRCh38_full_analysis_set_plus_decoy_hla.fa 16 pseudoDB list_human       
    ```  
    Note: Once the process is complete, you can find your output VCF files in the **variants/** directory.    





