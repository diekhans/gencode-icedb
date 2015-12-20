""""functions to operate on read data"""
import os
from pycbio.sys import PycbioException
from pycbio.sys import fileOps

def isFastq(readsFile):
    """does file appear to be a fastq by naming convention"""
    base = fileOps.compressBaseName(readsFile)
    return base.endswith(".fastq")

def getReadsCatCommand(readsFile):
    """Get command to cat reads file based on extension from either a bam or a fastq, possibly compressed"""
    base = fileOps.compressBaseName(readsFile)
    if base.endswith(".bam") or base.endswith(".sam"):
        return ["samtools", "bam2fq", readsFile]
    else:
        return [fileOps.compressCmd(readsFile), readsFile]


