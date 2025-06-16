# request all output stays in one .o file
# request 16 cores to run job
# request job only runs on host "colony-02"

#$ -j y
#$ -pe smp 16
#$ -l h=colony-02

echo $HOSTNAME
# sleep 10s
# SRR_list=("SRR13704012" "SRR13704098" "SRR13704283" "SRR13704108" "SRR13574198")
#(cancer transcript #1, healthy transcript #1, cancer transcipt #2, healthy transcript #2)
# 1<60 and 2>60
declare -A unhealthy_tissue_SRA=([LC_C1]=ERR164550  [LC_C7]=ERR164556 [LC_S20]=ERR164602 [LC_S14]=ERR164594)
declare -A healthy_tissue_SRA=([LC_C1]="ERR164473" [LC_C7]="ERR164477" [LC_S20]="ERR164517" [LC_S14]="ERR164511")
declare -A patient_age=([LC_C1]=54 [LC_C7]=82 [LC_S20]=58 [LC_S14]=79)
ref=~/data/transcriptome/reference/GCF_000001405.40_GRCh38.p14_rna.fna

#mapping database
# db=/home/chris/Teaching/BIOL4225/data/SARS-CoV2/SARS-CoV-2_Wuhan-Hu-1_genome.fasta 
# #set minid
# minid=.92

#get unhealthy tissue reads for each patient
for patient in ${!unhealthy_tissue_SRA[@]}
do
	#set value for unhealthy tissue ERR value for patient
	SRR=${unhealthy_tissue_SRA[$patient]}
	#echo "$SRR"

	mkdir -p unhealthy_tissue/${patient}/untrimmed
	output_dir=~/project/unhealthy_tissue/${patient}/untrimmed
	file_path=~/project/unhealthy_tissue/${patient}/untrimmed/${SRR}_1.fastq

	#grab reads using accession number
	#If the SRA file for patient does not exist, download data 
	if [ ! -f ${output_dir}/${SRR}_1.fastq ]; then
		echo "downloading reads from SRA for $patient"
    	fasterq-dump ${SRR} -O ${output_dir}
    #else, print message indicating that reads for patients already exists
    else
    	echo "unhealthy transcriiptome SRA file for $patient already exists"
	fi



	trimmmed_dir=~/project/unhealthy_tissue/${patient}/trimmed
	mkdir -p ${trimmmed_dir}

	# trim reads
	#If the read files have not been trimmed, trim reads
	if [ ! -f ${trimmmed_dir}/${SRR}_1.fastq ]; then

		echo "trimming unhealthy tissue reads for patient $patient"
    	cutadapt --cores 16 -q 10 --minimum-length 80 -o ${trimmmed_dir}/${SRR}_1.fastq -p ${trimmmed_dir}/${SRR}_2.fastq ${output_dir}/${SRR}_1.fastq ${output_dir}/${SRR}_2.fastq

    #else, print message indicating that reads for patients already exists
    else

    	echo "unhealthy tissue read files have already been trimmed for patient $patient"

	fi

	mapped_dir=~/project/unhealthy_tissue/${patient}/mapped
	mkdir -p ${mapped_dir} 
	#map reads
	
	#map reads if they have not been mapped yet
	if [ ! -f ${mapped_dir}/${patient}_mapped_80.bam ]; then

		echo "aligning unhealthy tissue reads for patient $patient"
		bbmap.sh threads=8 ref=~/data/transcriptome/reference/GCF_000001405.40_GRCh38.p14_rna.fna in=${trimmmed_dir}/${SRR}_1.fastq in2=${trimmmed_dir}/${SRR}_2.fastq trimreaddescriptions=t killbadpairs=t pairedonly=t minid=0.80 outm=${mapped_dir}/${patient}_mapped_80.bam

	else

		echo "unhealthy tissue reads have already been mapped for patient $patient"

	fi
	
	#sort and index mapped .bam files
	if [ ! -f ${mapped_dir}/${patient}_mapped_80_sorted.bam ]; then

		echo "sorting mapped unhealthy tissue .bam file for patient $patient"
		samtools sort ${mapped_dir}/${patient}_mapped_80.bam  > ${mapped_dir}/${patient}_mapped_80_sorted.bam 
		samtools index ${mapped_dir}/${patient}_mapped_80_sorted.bam 

	else

		echo "mapped unhealthy tissue .bam file has already been sorted for patient $patient"

	fi

	#call SNPs 

	#SNP call directory
	SNP_dir=~/project/unhealthy_tissue/${patient}/SNP_calls
	mkdir -p $SNP_dir

	#call variants
	if [ ! -f ${SNP_dir}/${patient}_unhealthy_calls_80.bcf ]; then

		echo "Unhealthy variant calling for patient $patient"
		bcftools mpileup --threads 16 -Ou -f ~/data/transcriptome/reference/GCF_000001405.40_GRCh38.p14_rna.fna ${mapped_dir}/${patient}_mapped_80_sorted.bam | bcftools call --threads 16 -mv -Ob -o ${SNP_dir}/${patient}_unhealthy_calls_80.bcf


	else

		echo "Unhealthy variants have already been called for patient $patient"

	fi

	
done

#get healthy tissue reads for each patient
for patient in ${!healthy_tissue_SRA[@]}
do
	#set value for healthy tissue ERR value for patient
	SRR=${healthy_tissue_SRA[$patient]}
	
	#make sure directory for healthy untrimmed reads exists
	mkdir -p healthy_tissue/${patient}/untrimmed
	#set output directory variable
	output_dir=~/project/healthy_tissue/${patient}/untrimmed
	#set file_path variable
	file_path=~/project/healthy_tissue/${patient}/untrimmed/${SRR}_1.fastq

	#grab reads using accession number
	#If the SRA file for patient does not exist, download data 
	if [ ! -f ${file_path} ]; then
		echo "downloading reads from SRA for $patient"
    	fasterq-dump ${SRR} -O ${output_dir}
    #else, print message indicating that reads for patients already exists
    else
    	echo "healthy transcriiptome SRA file for $patient already exists"
	fi

	trimmmed_dir=~/project/healthy_tissue/${patient}/trimmed
	mkdir -p ${trimmmed_dir}

	# trim reads
	#If the read files have not been trimmed, trim reads
	if [ ! -f ${trimmmed_dir}/${SRR}_1.fastq ]; then

		echo "trimming healthy tissue reads for patient $patient"
    	cutadapt --cores 16 -q 10 --minimum-length 80 -o ${trimmmed_dir}/${SRR}_1.fastq -p ${trimmmed_dir}/${SRR}_2.fastq ${output_dir}/${SRR}_1.fastq ${output_dir}/${SRR}_2.fastq

    #else, print message indicating that reads for patients already exists
    else

    	echo "healthy tissue read files have already been trimmed for patient $patient"

	fi

	mapped_dir=~/project/healthy_tissue/${patient}/mapped
	mkdir -p ${mapped_dir} 
	#map reads
	
	#map reads if they have not been mapped yet
	if [ ! -f ${mapped_dir}/${patient}_mapped_80.bam ]; then

		echo "aligning healthy tissue reads for patient $patient"
		bbmap.sh threads=16 ref=~/data/transcriptome/reference/GCF_000001405.40_GRCh38.p14_rna.fna in=${trimmmed_dir}/${SRR}_1.fastq in2=${trimmmed_dir}/${SRR}_2.fastq trimreaddescriptions=t killbadpairs=t pairedonly=t minid=0.80 outm=${mapped_dir}/${patient}_mapped_80.bam

	else

		echo "healthy tissue reads have already been mapped for patient $patient"

	fi

	#sort and index mapped .bam files
	if [ ! -f ${mapped_dir}/${patient}_mapped_80_sorted.bam ]; then

		echo "sorting mapped healthy tissue .bam file for patient $patient"
		samtools sort ${mapped_dir}/${patient}_mapped_80.bam  > ${mapped_dir}/${patient}_mapped_80_sorted.bam 
		samtools index ${mapped_dir}/${patient}_mapped_80_sorted.bam 

	else

		echo "mapped healthy tissue .bam file has already been sorted for patient $patient"

	fi

	#SNP call directory
	SNP_dir=~/project/healthy_tissue/${patient}/SNP_calls
	mkdir -p $SNP_dir

	#call variants
	if [ ! -f ${SNP_dir}/${patient}_healthy_calls_80.bcf ]; then

		echo "Healthy variant calling for patient $patient"
		bcftools mpileup --threads 16 -Ou -f ~/data/transcriptome/reference/GCF_000001405.40_GRCh38.p14_rna.fna ${mapped_dir}/${patient}_mapped_80_sorted.bam | bcftools call --threads 16 -mv -Ob -o ${SNP_dir}/${patient}_healthy_calls_80.bcf


	else

		echo "Healthy variants have already been called for patient $patient"

	fi


	unhealthy_calls=~/project/unhealthy_tissue/${patient}/SNP_calls/${patient}_unhealthy_calls_80.bcf

	#filter out non-somatic SNPs
	if [ ! -f ~/project/healthy_tissue/${patient}/cancer_only_calls_80.vcf ]; then

		echo "filtering out non-somatic SNPs for patient $patient"
		bcftools index ${unhealthy_calls}
		bcftools index ${SNP_dir}/${patient}_healthy_calls_80.bcf
		bcftools isec --threads 16 --complement ${unhealthy_calls} ${SNP_dir}/${patient}_healthy_calls_80.bcf -w 1 -o ~/project/healthy_tissue/${patient}/cancer_only_calls_80.vcf
		# bcftools index 	~/project/healthy_tissue/${patient}/cancer_only_calls_80.vcf


	else

		echo "unhealthy SNPs have already been filtered out for patient $patient"

	fi

python3 variant_analysis.py

done

