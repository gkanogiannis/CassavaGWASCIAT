#
# README.txt.sh
# Script part of CassavaGWASCIAT project
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

basename=Cassava_RAD_CIAT_GWAS.snps ;
input_vcf=/data_A/Cassava/Cassava_RAD/Cassava_RAD_CIAT/vcf/Cassava_RAD_CIAT.snps.filter_info_DP3.vcf.gz ;
pop_groups=`pwd`/Cassava_RAD_CIAT.LRLAC_MAF5.K7.sample_group.txt ;

# Keep GWAS 2020 samples from Cassava_RAD_CIAT
if [[ ! -f ${basename}.vcf.gz ]]; then
vcftools \
	--gzvcf input_vcf \
	--recode --recode-INFO-all \
	--keep samples.GWAS.txt \
	--stdout | bgzip -c > ${basename}.vcf.gz ;
fi ;

pushd geno ;
# Filtering and relatedness analysis
bash /home/agkanogiannis/CIAT/scripts/initialFilteringAndRelatedness.sh ${basename}.vcf.gz 0.2 0.5 0.03 1e-10 0.85 ${pop_groups} ;

# Look geno/${basename}.*.clones.tsv for related (MLE best) and look pheno.all.tsv for clones with low data
# Save a list for removal to ${basename}.remove_clones.txt
# Save any other sample to be removed (ie from het) to ${basename}.remove.txt 
if [[ ! -f ${basename}.remove_clones.txt ]]; then
	exit 0 ;
fi ;
vcftools --vcf ${basename}.initialFiltered.vcf --recode --recode-INFO-all --remove <(cat ${basename}.remove_clones.txt ${basename}.remove.txt) --stdout > ${basename}.filtered.vcf ;
vcftools --vcf ${basename}.initialFiltered.ld.vcf --recode --recode-INFO-all --remove <(cat ${basename}.remove_clones.txt ${basename}.remove.txt) --stdout > ${basename}.filtered.ld.vcf ;

# Calculate kinship (in ld thinned snps)
vcftools --vcf ${basename}.filtered.ld.vcf --relatedness2 --out ${basename}.filtered.ld ;
rm -f ${basename}.filtered.ld.log ;

# Impute 
java \
	-jar /home/agkanogiannis/software/beagle/beagle.12Jul19.0df.jar \
	gt=${basename}.filtered.vcf \
	out=${basename}.filtered.imp ;

# Convert to HapMap
/home/agkanogiannis/software/TASSEL5/run_pipeline.pl \
	-Xms64G -Xmx64G -fork1 \
	-vcf ${basename}.filtered.imp.vcf.gz \
	-export ${basename}.filtered.imp \
	-exportType HapmapDiploid ;

popd ;

# Perform GAPIT GWAS
pushd GAPIT ;
ln -sfn ../geno/${basename}.filtered.imp.hmp.txt ${basename}.hmp.txt ;
ln -sfn ../geno/${basename}.filtered.ld.relatedness2 ${basename}.kinship.txt ;


exit 0;
