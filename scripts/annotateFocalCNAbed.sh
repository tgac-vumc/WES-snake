##########################################################################################
#!/bin/bash
# date: 25 august 2017	, December 2017 - changed to work in snakemake pipeline Tjitske
# author: Matias Mendeville

# annotate bed file with all focal CNAs per sample...
# Function to annotate the genes (dbSNP, Cosmic) that are located in a region of a CNA

# Input: BED file
# Output: BED file with genes in column ..

# TODO: stricter selection of genes to annotate
##########################################################################################

input=$1
sname=$2
outdir=$3
output=$4

tmp_all_genes=$outdir$sname
tmp_geneAnn=$tmp_all_genes".geneAnn.bed"
tmp_cosmic_file=$tmp_all_genes".geneAnn_CosmicCensus.bed"


	# Add 'normal' genes on the selected regions (input bed file):
	perl scripts/regions2genes_loc.pl -o $tmp_all_genes -bed $input -p

	# Get number of columns in file (the column where the genes are listed)
	genecol=`awk '{print NF}' $tmp_geneAnn | sort -nu | tail -n 1`

	# Add the genes listed in COSMIC database: (from the genes listed in previous function)
	# NB: -col refers to the column where the genes are listed.
	perl scripts/addCosmicCencus.pl -bed $tmp_geneAnn -col $genecol

    # remove de column containing all non-Cosmic genes
	# TODO: adjust code so that een-na-laatste column wordt verwijderd
	#awk '!($$genecol="")' CVU17-183.geneAnn_CosmicCensus.bed > file w only cosmic genes
	awk '!($8="")' $tmp_cosmic_file > $output


	# remove unnecessary files:
	rm $tmp_geneAnn $tmp_cosmic_file
