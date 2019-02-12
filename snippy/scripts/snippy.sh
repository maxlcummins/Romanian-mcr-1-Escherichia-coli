#make timestamped user directory in scratch
SCRATCH="/scratch/${USER}_$(date +%Y%m%d%H%M%S)"
mkdir ${SCRATCH}
mkdir ${SCRATCH}/reads
# Change to the PBS working directory where qsub was started from.
# Copy your input files from there to scratch.
cd ${PBS_O_WORKDIR}
cp ${1}/*.fastq.gz ${SCRATCH}/reads
cp ${2} ${SCRATCH}


###############
# Start the Job
###############

# Change to the scratch directory where you copied your 
# input files to before you start. 
cd ${SCRATCH}

#activate environment
source activate snippy

#define array of input fastq files
array=($(ls ${SCRATCH}/reads/*.fastq.gz))

for ((i=0; i<${#array[@]}; i+=2));
do snippy --outdir ${array[i]}.out --ref $2 --pe1 ${array[i]} --pe2 ${array[i+1]};
done

#move snippy output to read directory
mv $SCRATCH/reads/*.out $1

#remove scratch directory
rm -rf $SCRATCH
