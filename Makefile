batchSystem = singleMachine
maxThreads = 15
defaultMemory = 8589934592
jobTree = .jobTree
log = log.txt

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:${PATH}

all :
	cd sonLib && make
	cd jobTree && make

run : all
	if [ -d ${jobTree} ]; then rm -rf ${jobTree}; fi
	python src/bamSlicerPipeline.py --maxThreads=${maxThreads} --batchSystem=${batchSystem} \
	--defaultMemory=${defaultMemory} --jobTree ${jobTree} --logLevel DEBUG &> ${log}
