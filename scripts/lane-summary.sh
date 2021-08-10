#!/bin/bash

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8

title=$1
bamfolder=$2
samples=$bamfolder*.bam



echo -e '<?xml version="1.0" encoding="UTF-8"?>'
echo -e '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN"'
echo -e '\t"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
echo -e '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'
echo -e '\t<head>'
echo -e "\t\t<title>$title</title>"
echo -e '\t\t<style type="text/css">'
echo -e '\t\t\ttable {border-collapse: collapse; font-size: smaller;}'
echo -e '\t\t\tth {border: 1px solid gray; padding: 5px;}'
echo -e '\t\t\ttd {border: 1px solid gray; padding: 5px; text-align: right;}'
echo -e '\t\t</style>'
#echo -e '\t\t<link rel="stylesheet" href="http://ccagc-gen01.vumc.nl/js/lightbox2-master/dist/css/lightbox.css">'
echo -e '\t\t<link rel="stylesheet" href="../QDNAseq.snakemake/lb2/css/lightbox.css">'
echo -e '\t</head>'
echo -e '\t<body>'
echo -e "\t\t<h1>$title</h1>"
#echo -e '\t\t<p><a href="../index.html">Back to overview</a></p>'
echo -e '\t\t<table>'
echo -e '\t\t\t<tr>'
echo -e '\t\t\t\t<th>sample</th>'
echo -e '\t\t\t\t<th>total</th>'
echo -e '\t\t\t\t<th colspan="2">aligned</th>'
echo -e '\t\t\t\t<th colspan="2">unique</th>'
echo -e '\t\t\t\t<th colspan="2">q1</th>'
echo -e '\t\t\t\t<th colspan="2">q37</th>'
echo -e '\t\t\t\t<th>qc</th>'
echo -e '\t\t\t\t<th><a href="profiles/corrected/index.html">corrected</a></th>'
echo -e '\t\t\t\t<th><a href="profiles/dewaved/index.html">dewaved</a></th>'
echo -e '\t\t\t\t<th><a href="profiles/segmented/index.html">segmented</a></th>'
echo -e '\t\t\t\t<th><a href="profiles/called/index.html">called</a></th>'
echo -e '\t\t\t\t<th><a href="profiles/reCalled/index.html">reCalled</a></th>'
echo -e '\t\t\t</tr>'

sumtotal=0
sumaligned=0
sumunique=0
sumq1=0
sumq37=0
n=0

for sample in $samples
do
  n=`echo $n + 1 | bc -l`
  sample=`basename $sample .bam`
  sample3=`echo ../qc-fastq/${sample}*_fastqc.html`
  sample2=`echo $sample | sed -r 's/^[0-9]{6}_[A-Z0-9]{9}_L[1-8]{1}_//'`   #TOTO - remove line and remove sample lines
  total=`cat "../stats/${sample}.reads.all"`
  aligned=`cat "../stats/${sample}.reads.aligned"`
  unique=`cat "../stats/${sample}.reads.unique"`
  q1=`cat "../stats/${sample}.reads.q1"`
  q37=`cat "../stats/${sample}.reads.q37"`
  sumtotal=`echo $sumtotal + $total | bc -l`
  sumaligned=`echo $sumaligned + $aligned | bc -l`
  sumunique=`echo $sumunique + $unique | bc -l`
  sumq1=`echo $sumq1 + $q1 | bc -l`
  sumq37=`echo $sumq37 + $q37 | bc -l`
  echo -e '\t\t\t<tr>'
  echo -e '\t\t\t\t<td style="text-align: left;">'$sample2'</td>'
  printf "\t\t\t\t<td>%'i</td>\n" $total
  printf "\t\t\t\t<td>%'i</td>\n" $aligned
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $aligned/$total*100 | bc -l`
  printf "\t\t\t\t<td>%'i</td>\n" $unique
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $unique/$aligned*100 | bc -l`
  printf "\t\t\t\t<td>%'i</td>\n" $q1
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $q1/$unique*100 | bc -l`
  printf "\t\t\t\t<td>%'i</td>\n" $q37
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $q37/$unique*100 | bc -l`
  echo -en '\t\t\t\t<td><a href='$sample3'>fastq</a>, '
  echo -e '<a href="../qc-bam/'$sample'_fastqc.html">bam</a></td>'
  echo -e '\t\t\t\t<td><a href="profiles/corrected/'$sample'.png">corrected</a></td>'
  echo -e '\t\t\t\t<td><a href="profiles/dewaved/'$sample'.png">dewaved</a></td>'
  echo -e '\t\t\t\t<td><a href="profiles/segmented/'$sample'.png">segmented</a></td>'
  echo -e '\t\t\t\t<td><a href="profiles/called/'$sample'.png">called</a></td>'
  echo -e '\t\t\t\t<td><a href="profiles/reCalled'$sample'.png">reCalled</a></td>'  
  echo -e '\t\t\t</tr>'
done
echo -e '\t\t\t<tr style="border-top: double;">'
echo -e '\t\t\t\t<td style="text-align: left;">Total</td>'
printf "\t\t\t\t<td>%'i</td>\n" $sumtotal
printf "\t\t\t\t<td>%'i</td>\n" $sumaligned
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumaligned/$sumtotal*100 | bc -l`
printf "\t\t\t\t<td>%'i</td>\n" $sumunique
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumunique/$sumaligned*100 | bc -l`
printf "\t\t\t\t<td>%'i</td>\n" $sumq1
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumq1/$sumunique*100 | bc -l`
printf "\t\t\t\t<td>%'i</td>\n" $sumq37
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumq37/$sumunique*100 | bc -l`
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t</tr>'
echo -e '\t\t\t<tr>'
echo -e '\t\t\t\t<td style="text-align: left;">Average</td>'
printf "\t\t\t\t<td>%'.0f</td>\n" `echo $sumtotal / $n | bc -l`
printf "\t\t\t\t<td>%'.0f</td>\n" `echo $sumaligned / $n | bc -l`
echo -e '\t\t\t\t<td>&nbsp;</td>'
printf "\t\t\t\t<td>%'.0f</td>\n" `echo $sumunique / $n | bc -l`
echo -e '\t\t\t\t<td>&nbsp;</td>'
printf "\t\t\t\t<td>%'.0f</td>\n" `echo $sumq1 / $n | bc -l`
echo -e '\t\t\t\t<td>&nbsp;</td>'
printf "\t\t\t\t<td>%'.0f</td>\n" `echo $sumq37 / $n | bc -l`
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t</tr>'

echo -e '\t\t</table>'
#echo -e '\t<script src="http://ccagc-gen01.vumc.nl/js/lightbox2-master/dist/js/lightbox-plus-jquery.min.js"></script>'
echo -e '\t<script src=""../QDNAseq.snakemake/lb2/js/lightbox-plus-jquery.min.js"></script>'
echo -e '\t</body>'
echo -e '</html>'
