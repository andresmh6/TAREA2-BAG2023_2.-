#!/bin/bash

if [ $# -ne 1 ]
then
    echo "The number of arguments is: $#"
    echo "Usage: rnapipe.sh <params.file>"
    echo ""
    echo "params.file: Input file with arguments"
    echo "An example of params.file can be found in the test folder"
    exit
fi

PARAMS=$1

echo ""
echo "====================="
echo "LOADING PARAMETERS"
echo "==================="
echo ""

INSDIR=$(grep installation $PARAMS | awk '{print($2)}')
echo "Installation directory: $INSDIR"

WD=$(grep working $PARAMS | awk '{print($2)}')
echo "Working directory: $WD"

EXP=$(grep experiment $PARAMS | awk '{print($2)}')
echo "Experiment name: $EXP"

NUMREP=$(grep number_replicas $PARAMS | awk '{print($2)}')
echo "Number of replica: $NUMREP"

GENOME=$(grep genome $PARAMS | awk '{print($2)}')
echo "Genome path: $GENOME"

ANNOTATION=$(grep annotation $PARAMS | awk '{print($2)}')
echo "Genome path: $ANNOTATION"

CHR=$(grep chromosomes $PARAMS | awk '{print($2)}')
echo "Universe of chromosomes: $CHR"

PEAK=$(grep peak $PARAMS | awk '{print($2)}')
echo "Peak type: $PEAK"

ENDTYPE=$(grep single $PARAMS | awk '{print($2)}')
echo "Single or paired: $ENDTYPE"

TSSUP=$(grep upstream $PARAMS | awk '{print($2)}')
echo "TSS upstream region: $TSSUP"

TSSDOWN=$(grep downstream $PARAMS | awk '{print($2)}')
echo "TSS downstream region: $TSSDOWN"

#Generamos 2 arrays que contienen muestras ChIP y muestras control tipo INPUT por separado.

CHIPS=()
INPUTS=()
i=0

if [ $ENDTYPE -eq 1 ]
then
    while [ $i -lt $NUMREP ]
    do
        j=$(($i + 1))
        CHIPS[$i]=$(grep path_sample_chip_$j $PARAMS | awk '{print($2)}')
        INPUTS[$i]=$(grep path_sample_input_$j $PARAMS | awk '{print($2)}')
        ((i++))
    done

elif [ $ENDTYPE -eq 2 ]
then
    while [ $i -lt $NUMREP ]
    do
        j=$(($i + 1))
        k=$(($i * 2))
        l=$(($k + 1))
        CHIPS[$k]=$(grep path_sample_chip_$j: $PARAMS | awk '{print($2)}')
        CHIPS[$l]=$(grep path_sample_chip_$j: $PARAMS | awk '{print($3)}')
        INPUTS[$k]=$(grep path_sample_input_$j: $PARAMS | awk '{print($2)}')
        INPUTS[$l]=$(grep path_sample_input_$j: $PARAMS | awk '{print($3)}')
        ((i++))
    done
else
    echo "Not permitted input to determine if reads are single or paired"
fi

echo "Samples:"
echo "${CHIPS[@]}"
echo "${INPUTS[@]}"

#Creamos el espacio de trabajo

echo ""
echo "====================="
echo "GENERATING WORK SPACE"
echo "====================="
echo ""

cd $WD
mkdir $EXP
cd $EXP
mkdir genome annotation results samples scripts

cp $GENOME genome/genome.fa
cp $ANNOTATION annotation/annotation.gtf
cd samples

if [ $ENDTYPE -eq 1 ]
then
    i=1
    while [ $i -le $NUMREP ]
    do
        mkdir replica_$i
        cd replica_$i
        mkdir chip input replica_results
        j=$(($i-1))
        cp ${CHIPS[$j]} chip/sample_chip_$i.fq.gz
        cp ${INPUTS[$j]} input/sample_input_$i.fq.gz
        cd ..
        ((i++))
    done

elif [ $ENDTYPE -eq 2 ]
then
    i=1
    while [ $i -le $NUMREP ]
    do
        mkdir replica_$i
        cd replica_$i
        mkdir chip input replica_results
        cd chip
        j=$(($i - 1))
        k=$(($j * 2))
        l=$(($k + 1))
        cp ${CHIPS[$k]} sample_chip_${i}_1.fq.gz
        cp ${CHIPS[$l]} sample_chip_${i}_2.fq.gz
        cd ..
        cd input
        cp ${INPUTS[$k]} sample_input_${i}_1.fq.gz
        cp ${INPUTS[$l]} sample_input_${i}_2.fq.gz
        cd ../..
        ((i++))
    done
else
    echo "Not permitted input to determine if reads are single or paired"
fi

echo ""
echo "================="
echo "WORKPLACE CREATED"
echo "================="
echo ""

#Generación del índice del genoma

echo ""
echo "=============="
echo "CREATING INDEX"
echo "=============="
echo ""

cd ../genome

bowtie2-build genome.fa index

echo ""
echo "============="
echo "INDEX CREATED"
echo "============="
echo ""

cd ../results

#Enviamos cada muestra al procesamiento del script sample_proc.sh

i=1
while [ $i -le $NUMREP ]
do
    sbatch --job-name=proc_sam_$i --output=sam_$i --error=err_sam_$i $INSDIR/sample_proc.sh $WD $i $PEAK $NUMREP $INSDIR $EXP $CHR $TSSUP $TSSDOWN $GENOME
    ((i++))
done
