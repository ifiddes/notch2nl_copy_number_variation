batchSystem="gridEngine"
maxThreads="15"
defaultMemory="8589934592"
jobTree=".jobTree/"
log="log.txt"
maxCpus="250"


export PYTHONPATH=./:${PYTHONPATH}
export PATH=./sonLib/bin:./submodules/jobTree/bin:${PATH}


if [ -d ${jobTree} ]; then 
    rm -rf ${jobTree}
fi

if [ $1 == "fastq" ]; then
    python src/fastqPipeline.py --fastq $2 --maxThreads=${maxThreads} --batchSystem=${batchSystem} --maxCpus ${maxCpus} --defaultMemory=${defaultMemory} --jobTree ${jobTree} --logLevel DEBUG --stats &> ${log}
elif [ $1 == "fastq_list" ]; then
    python src/fastqPipeline.py --fastq_list $2 --maxThreads=${maxThreads} --batchSystem=${batchSystem} --maxCpus ${maxCpus} --defaultMemory=${defaultMemory} --jobTree ${jobTree} --logLevel DEBUG --stats &> ${log}
fi