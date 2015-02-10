batchSystem = gridEngine
maxThreads = 15
defaultMemory = 8589934592
jobTree = .jobTree
log = log.txt
maxCpus = 250

ifeq ($(strip $(queries)),)
queries = queries/queries.pickle
endif

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:${PATH}


all :
	cd sonLib && make
	cd jobTree && make

run : all
	if [ -d ${jobTree} ]; then rm -rf ${jobTree}; fi
	python src/bamSlicerPipeline.py -q ${queries} --maxThreads=${maxThreads} --batchSystem=${batchSystem} --maxCpus ${maxCpus} --defaultMemory=${defaultMemory} --jobTree ${jobTree} --logLevel DEBUG --stats &> ${log}
