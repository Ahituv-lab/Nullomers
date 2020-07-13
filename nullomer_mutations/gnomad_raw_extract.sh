#!/bin/bash
#

filename="path/to/gnomad/file/eg/gnomad.genomes.r3.0.sites.vcf.gz/it/is/assumed/to/be/in/gz/format"

zcat $filename | grep -v -e '#' | grep -e 'AF_raw' | awk '{if(length($4)==1 && length($5)==1){ print $0; }}' | awk -F\; '{af=0; OFS="\t";
                        for(i=1;i<=NF;i++){ 
                                if($i ~ /AF_raw=/){
                                        split($i,afraw,"=");
                                        af = afraw[2];
                                        break; 
                                }
                        }
                        print $0, af; }' $1 | cut -f 1-5,9 > /path/to/directory/snps_only_gnomad_afraw.vcf
