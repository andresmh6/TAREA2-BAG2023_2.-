#!/bin/bash

#SBATCH --export=ALL

WD=$1
i=$2
PEAK=$3
NUMREP=$4
INSDIR=$5
EXP=$6
CHR=$7
TSSUP=$8
TSSDOWN=$9
GENOME=${10}

echo ""
echo "===================="
echo "PROCESSING SAMPLE $i"
echo "===================="
echo ""

cd $WD/$EXP/samples/replica_$i

cd chip

if [ -f sample_chip_{$i}_2.fq.gz ]
then
    fastqc sample_chip_${i}_1.fq.gz
    fastqc sample_chip_${i}_2.fq.gz
    bowtie2 -x ../../../genome/index -1 sample_chip_{$i}_1.fq.gz -2 sample_chip_{$i}_2.fq.gz -S chip_$i.sam

else
    fastqc sample_chip_$i.fq.gz
    bowtie2 -x ../../../genome/index -U sample_chip_$i.fq.gz -S chip_$i.sam
fi

cd ../input

#Se realiza el análisis de calidad de las muestras y la generación de los .sam en función de si las muestras son PAIRED o UNPAIRED

if [ -f sample_input_{$i}_2.fq.gz ]
then
    fastqc sample_input_${i}_1.fq.gz
    fastqc sample_input_${i}_2.fq.gz
    bowtie2 -x ../../../genome/index -1 sample_input_${i}_1.fq.gz -2 sample_input_${i}_2.fq.gz -S input_$i.sam
else
    fastqc sample_input_$i.fq.gz
    bowtie2 -x ../../../genome/index -U sample_input_$i.fq.gz -S input_$i.sam
fi

#Generamos los archivos .bam a partir de los .sam y eliminamos estos últimos.

samtools sort -o input_$i.bam input_$i.sam
rm input_$i.sam
cd ../chip
samtools sort -o chip_$i.bam chip_$i.sam
rm chip_$i.sam

#Generamos los índices de los archivos .bam para una posible visualización en IGV.
samtools index input_$i.bam
samtools index chip_$i.bam

cd ../replica_results

#Tras comprobar el tipo de pico se realiza la llamada con MACS2.

if [ $PEAK -eq 1 ]
then
    macs2 callpeak -t ../chip/chip_$i.bam -c ../input/input_$i.bam -f BAM -n $i

elif [ $PEAK -eq 2 ]
then
    macs2 callpeak --broad -t ../chip/chip_$i.bam -c ../input/input_$i.bam -f BAM -n $i
fi

echo "Peak calling $i done!" >> ../../../results/merged_list.txt

cd ../../..

echo ""
echo "=================="
echo "REPLICA $i CREATED"
echo "=================="
echo ""

#Generamos variable NUMPROC que contiene el número de muestras que se han procesado
NUMPROC=$(wc -l results/merged_list.txt | awk '{print($1)}')

#Cuando NUMPROC indica que todas las muestras son procesadas, estas son enviadas al siguiente script: peak_call.sh
if [ $NUMPROC -eq $NUMREP ]
then
    sbatch --job-name=peak_call_$i --output=peaks_$i --error=err_peak_$i $INSDIR/peak_call.sh $WD $EXP $PEAK $NUMREP $INSDIR $CHR $TSSUP $TSSDOWN $WD/$EXP/genome/genome.fa 
fi

