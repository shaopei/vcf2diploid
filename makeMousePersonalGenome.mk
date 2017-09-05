USAGE := "make -f makePersonalGenome.mk DATA_DIR=/path/to/your/VCFdir OUTPUT_SAMPLE_NAME=NA12878 MAT_SAMPLE_NAME=mat AT_SAMPLE_NAME=pat FILE_NAME_BAM=filename.in.DATA_DIR.bam FILE_NAME_VCF=filename.in.DATA_DIR.vcf"

##
## Define fixed, experiment-wide analysis parameters
##
# location of fasta file (combined)
FASTA_PATH := /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/GRCm38_68.fa

##
## Machine-specific variables
##
N_THREADS := 32
MAX_RAM := 100000000000

## Cluster-specific variables
EXE_DIR := /workdir/sc2457/tools/vcf2diploid_v0.2.6a
JAVA_EXE := java

##
## Run-specific variables
##
OUTPUT_SAMPLE_NAME := NULL
MAT_SAMPLE_NAME := NULL
PAT_SAMPLE_NAME := NULL
FILE_NAME_BAM := NULL
FILE_NAME_VCF := NULL
DATA_DIR := NULL
OUTPUT_DIR := $(DATA_DIR)/PersonalGenome_$(OUTPUT_SAMPLE_NAME)
VCF2SNP_SNP_INDEL := 1 ## only SNPs for CNV calc; to include indels, set this to 0

VCF_sampleID := $(OUTPUT_SAMPLE_NAME)

all: processSample



##
## Make output directoy

$(OUTPUT_DIR): 
	@echo -e "$(USAGE)"
	mkdir $(OUTPUT_DIR)
	mkdir $(OUTPUT_DIR)/AltRefMother
	mkdir $(OUTPUT_DIR)/AltRefFather


##
## VCF to diploid
##
$(OUTPUT_DIR)/maternal.chain:$(OUTPUT_DIR) 
	$(JAVA_EXE) -Xmx$(MAX_RAM) -jar $(EXE_DIR)/vcf2diploid.jar -id $(VCF_sampleID) -pass -chr $(FASTA_PATH) -vcf $(DATA_DIR)/$(FILE_NAME_VCF) -outDir $(OUTPUT_DIR) >& $(OUTPUT_DIR)/vcf2diploid.log
	cat $(OUTPUT_DIR)/maternal.chain | awk ' ($$1 == "chain") {print $$1,$$2,$$8,$$9,$$10,$$11,$$12,$$3,$$4,$$5,$$6,$$7,$$13}; ($$1 != "chain"){print $$1, $$3, $$2 }'>  $(OUTPUT_DIR)/mat2ref.chain
	cat $(OUTPUT_DIR)/paternal.chain | awk ' ($$1 == "chain") {print $$1,$$2,$$8,$$9,$$10,$$11,$$12,$$3,$$4,$$5,$$6,$$7,$$13}; ($$1 != "chain"){print $$1, $$3, $$2 }'>  $(OUTPUT_DIR)/pat2ref.chain

##
## Build Maternal genome
##
$(OUTPUT_DIR)/bowtie_build.maternal.log: $(OUTPUT_DIR)/maternal.chain
	bowtie-build --offrate 2 $(OUTPUT_DIR)/1_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/2_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/3_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/4_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/5_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/6_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/7_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/8_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/9_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/10_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/11_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/12_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/13_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/14_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/15_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/16_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/17_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/18_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/19_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/X_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/Y_$(VCF_sampleID)_maternal.fa,$(OUTPUT_DIR)/MT_$(VCF_sampleID)_maternal.fa $(OUTPUT_DIR)/AltRefMother/AltRefMother > $(OUTPUT_DIR)/bowtie_build.maternal.log


##
## Build Paternal genome
##
$(OUTPUT_DIR)/bowtie_build.paternal.log: $(OUTPUT_DIR)/bowtie_build.maternal.log
	bowtie-build --offrate 2 $(OUTPUT_DIR)/1_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/2_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/3_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/4_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/5_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/6_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/7_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/8_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/9_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/10_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/11_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/12_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/13_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/14_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/15_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/16_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/17_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/18_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/19_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/X_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/Y_$(VCF_sampleID)_paternal.fa,$(OUTPUT_DIR)/MT_$(VCF_sampleID)_paternal.fa $(OUTPUT_DIR)/AltRefFather/AltRefFather > $(OUTPUT_DIR)/bowtie_build.paternal.log


##
## Count reads in each snp region
##
$(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.counts: $(OUTPUT_DIR)/bowtie_build.paternal.log
	$(EXE_DIR)/vcf2snp -c $(VCF_sampleID) -m $(MAT_SAMPLE_NAME) -d $(PAT_SAMPLE_NAME) -p 1 -r 1 -s $(VCF2SNP_SNP_INDEL) $(DATA_DIR)/$(FILE_NAME_VCF) > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.unfiltered
	python $(EXE_DIR)/Remove_snps_skipped_in_vcf2diploid.py $(OUTPUT_DIR)/vcf2diploid.log $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.unfiltered $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.skipped
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp | awk '{print $$1"\t"($$2-1000)"\t"($$2+1000)"\t"$$1"_"($$2-1000)"_"($$2+1000)}' | grep -v "-" > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.unsort
	LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/  --parallel=10 -k1,1 -k2,2n $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.unsort > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed
	intersectBed -sorted -wo -bed -a $(DATA_DIR)/$(FILE_NAME_BAM) -b $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed | awk '{print $$16}' | LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/  --parallel=10 | uniq -c > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.counts


##
## Compute SNV from read coverage and create .snp and .cnv files that will be used as input to alleleseq
##
$(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleSeqInput.snp: $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.counts
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.counts | awk '{ total += $$1; count++ } END { print total/count }' > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.meancount
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.counts | LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/  --parallel=10 -n -k 1 | awk '{ lines[NR]=$$0; } END { print lines[int(NR/2)+1] }' | awk '{print $$1}' > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.mediancount
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.counts | awk '{getline avg<"$(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp.bed.mediancount"; print $$2"\t"$$1/avg }' > $(OUTPUT_DIR)/tmp.cnv
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snp | awk '{print $$1"_"($$2-1000)"_"($$2+1000)"\t"$$0}' | LC_ALL=C sort --temporary-directory=/workdir/sc2457/tmp/  --parallel=10 -k 1 > $(OUTPUT_DIR)/tmp.snp
	## Combine the snp data and the cnv data by region ID
	awk 'NR==FNR {h[$$1] = $$0; next} {print h[$$1]"\t"$$0}' $(OUTPUT_DIR)/tmp.snp $(OUTPUT_DIR)/tmp.cnv > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snpANDcnv
	## Create the snp and cnv input files for alleleseq
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snpANDcnv | awk '{print $$2"\t"$$3"\t"$$4"\t"$$5"\t"$$6"\t"$$7"\t"$$8}' > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleSeqInput.snp
	echo -e "chrm\tsnppos\trd" > $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleSeqInput.cnv
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).snpANDcnv | awk '{print $$2"\t"$$3"\t"$$10}' >> $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleSeqInput.cnv
	cat $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleSeqInput.snp | awk '{OFS="\t"}{FS="\t"}{print "chr"$$1, $$2-1, $$2, $$6 }' >  $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleDBInput.snp.bed
	rm $(OUTPUT_DIR)/tmp.cnv
	rm $(OUTPUT_DIR)/tmp.snp


processSample: $(OUTPUT_DIR)/$(OUTPUT_SAMPLE_NAME).alleleSeqInput.snp

