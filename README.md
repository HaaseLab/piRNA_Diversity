[View Haase Lab HERE](https://www.niddk.nih.gov/research-funding/at-niddk/labs-branches/laboratory-cell-molecular-biology/rna-biology-section). 

<img src="https://user-images.githubusercontent.com/11409899/109845179-d802fb00-7c1a-11eb-8e45-43b0dbf50de6.png" width="600">

#
### __*‘Cellular abundance shapes function in piRNA-guided genome defense’*__. 

__Authors__  
*"Pavol Genzor\*, Parthena Konstantinidou\*, Daniel Stoyko\**     *\*equal contributors*  
*Amirhossein Manzourolajdad, Celine Marlin Andrews, Alexandra R. Elchert, Constantinos Stathopoulos, Astrid D. Haase*  
  
  
This vignette describes the computational materials & methods associated with this manuscript. Please visit [**HaaseLab/piRNA_Diversity github repository**](https://github.com/HaaseLab/piRNA_Diversity) to download functions for the script. Please reffer to the GEO dataset [**GSE156058**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156058) associated with this study for adapter sequences and raw data. The analysis in this vignette was not performed with full data sets, but only subset of the data to demonstrate the materials and methods.

####
__About small RNA Libraries__  
To be able to account for the PCR duplication when quantifying individueal piRNA sequences, cloning procedure utilized adapter sequences containing multiple random nucleotides. Each ligated read contains 10 N's or **U**nique **M**olecular **I**dentifiers **(UMIs)**. There are 8N at the 5-prime and 2N at the 3-prime of the small RNA that allows for 4^10 or 1,048,576 possible combinations. 

__Library structure__  
5-FivePrimeAdapter--*NNNNNNNN*--**smallRNA**--*NN*--ThreePrimeAdapter-3

####
__Pre-requisites__  
* Acquire the raw sequencing data from your facility or GEO at NCBI
* Ensure that you have the appropriate 5-prime and 3-prime adapter sequences
* Ensure you have access to computing cluster and a computer with R environment
* Ensure that you have all the necessary software installed and running

####
__Vignette Content__  
  
A. DATA PREPARATION

  1. Pre-process fastq files
  2. Generate unique sequence fasta files
  3. Remove structural contaminants and align to the genome
  4. Load and process files in R.  

B. FIGURE-RELATED SCRIPTS. 

  * Figure 1
  * Figure 2
  * Figure S2 & 3B
  * Figure 3 & 4
  * Additinoal figures
  
C. *__bash__* Scripts  

