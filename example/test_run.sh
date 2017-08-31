#!/bin/sh

java -jar ../vcf2diploid.jar \
          -id NA12878 \
          -chr human_chr20_hg18.fa \
	  -vcf CEU.trio.chr20.2010_03.genotypes.vcf