Created: 4/24 - N.Q.
Last Modified: 4/24 11:15 AM - N.Q.

Genotype Likelihood Consists of a multistep processs: 

The first step is utilizing the reference genome and our sequencing files (merged/filtered BAM files) to create the vcf files -> utilizing GLIMPSE 
    NOTE: WE ARE USING THE UPDATED GLIMPSE2 for all analysis as it is the newer version and it is tailored for low-coverage sequencing analysis 

The files that we require for this step is: 
1. GLIMPSE
2. GLIMPSE_chunk bash script (can be copied from my directory or acquired on GLIMPSE github page for executable)
3. gl_calculation.sh <- my script to call glimpse to do the genotype likelihood calculations>
4. Making the .txt file with all the bam file paths to input into the gls calculation script 


**THESE SHOULD ALL BE RUN ON AMAREL**, not on here. I am dropping the scripts here for reference only 
