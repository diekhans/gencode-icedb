from __future__ import print_function
import sys
import os
from pycbio.sys import fileOps
from gencode_icedb.rnaSeqData import setDatabase, RnaSeqData
from gencode_icedb import pipelineOps, config
from peewee import SqliteDatabase
from toil.job import Job


def starGenerateGenome(refGeneomeName, refGeneomeFa, readLength, geneAnnotationSetName,
                       geneAnnotationSetGtf, numThreads):
    """Generate STAR genome given reference genome, read, and annotation set.
    Create a flag file to indicate it is complete.
    """
    pass
    # cmd = ["starGenerateGenome", "--numThreads={}".format(numThreads), "refGeneomeFa",
    #        geneAnnotationSetGtf, readLength, genomeDir]
    # NOT IMPLEMENTED YET


def starSpliceJunctionMapFn(job, genomeDir, readsFile, readsFile2, numThreads, sjOut):
    """Run STAR to map reads to splice junctions.
    """
    cmd = ["starSpliceJunctionMap", "--numThreads={}".format(numThreads),
           genomeDir, readsFile]
    if readsFile2 is not None:
        cmd.append("--readsFile2={}".format(+readsFile2))
    sjOutTmp = fileOps.atomicTmpFile(sjOut)
    pipelineOps.runCmd(cmd)
    fileOps.atomicInstall(sjOutTmp, sjOut)


def checkGenomeExists(pathConfig, readLength):
    "make sure the genome done flag is there"
    doneFlag = pathConfig.starGenomeDoneFlag(readLength)
    if not os.path.exists(doneFlag):
        raise Exception("STAR genome done flag not found: {}".format(doneFlag))


def maybeMakeSplitJunctionJob(job, pathConfig, rnaSeqData, numThreads):
    "check if output exists, if not add job to generate it"
    checkGenomeExists(pathConfig, rnaSeqData.readlength)
    sjOut = pathConfig.rnaSeqSetAnalysisSjOut(rnaSeqData.setname, rnaSeqData.runname)
    if not os.path.exists(sjOut):
        readsPath = pathConfig.rnaSeqSetDataFile(rnaSeqData.setname, rnaSeqData.readsfile)
        readsPath2 = pathConfig.rnaSeqSetDataFile(rnaSeqData.setname, rnaSeqData.readsfile2) if rnaSeqData.readsfile2 is not None else None
        job.addChildJobFn(starSpliceJunctionMapFn,
                          pathConfig.starGenomeDir(rnaSeqData.readlength),
                          readsPath, readsPath2, numThreads, sjOut)


def masterJobFn(job):
    pass  # nothing yet


def createJobs(opts):
    pathConfig = config.PathConfig(opts.rootDir, opts.organism, opts.assembly, opts.geneSet)
    database = SqliteDatabase(pathConfig.dataDatabase())
    setDatabase(database)
    masterJob = Job.wrapJobFn(masterJobFn)
    for rnaSeqData in RnaSeqData.select():
        maybeMakeSplitJunctionJob(masterJob, pathConfig, rnaSeqData, opts.numThreads)
        print("ONLY ONE JOB", file=sys.stderr)
        break
    return masterJob
