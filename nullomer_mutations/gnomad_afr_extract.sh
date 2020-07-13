#!/bin/bash
#

filename="/path/to/gnomad/file/eg/gnomad.genomes.r3.0.sites.vcf.gz"

zcat $filename | grep -v -e '#' | grep -e 'AF_afr' | awk -F\; '{success=1; a=0.05; b=0.01;;
                        for(i=1;i<=NF;i++){ 
                                if($i ~ /AF_afr=/){
                                        split($i,afafr,"="); 
                                        if(afafr[2]<a){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_fin=/){
                                        split($i,affin,"="); 
                                        if(affin[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_asj=/){
                                        split($i,afasj,"="); 
                                        if(afasj[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_eas=/){
                                        split($i,afeas,"="); 
                                        if(afeas[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_nfe=/){
                                        split($i,afnfe,"="); 
                                        if(afnfe[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_ami=/){
                                        split($i,afami,"="); 
                                        if(afami[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_amr=/){
                                        split($i,afamr,"="); 
                                        if(afamr[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_sas=/){
                                        split($i,afsas,"="); 
                                        if(afsas[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_raw=/){
                                        split($i,afraw,"="); 
                                        if(afraw[2]>b){
                                                success=0
                                                break
                                        }
                                }
                                else if($i ~ /AF_oth=/){
                                        split($i,afoth,"="); 
                                        if(afoth[2]>b){
                                                success=0
                                                break
                                        }
                                }
                        }
                        if(success==1){ print $0; }}' $1 > /path/to/directory/gnomad_afr.vcf
