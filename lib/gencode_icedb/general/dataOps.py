""""general functions to operate on file or memory data"""
import os
import socket
import re
from pycbio.sys import fileOps
import pipettor
import logging


def ensureList(strOrList):
    if isinstance(strOrList, str):
        return [strOrList]
    else:
        return strOrList


def isChrYPar(annotTrans):
    """Is this a PAR gene on Y"""
    # UCSC db puts PAR tag on both chrX and chrY, since it tracks only by transcript id
    return ("PAR" in annotTrans.attrs.tags) and (annotTrans.chrom.name in ("chrY", "Y"))


def ensemblIdSplit(ensIds):
    """split into gene and transcript ids"""
    # split by matching pattern
    ensIds = ensureList(ensIds)
    geneIds = []
    transIds = []
    for ensId in ensIds:
        if re.match("ENS[A-Z]*G[0-9]+\\.[0-9]+$", ensId):
            geneIds.append(ensId)
        elif re.match("ENS[A-Z]*T[0-9]+\\.[0-9]+$", ensId):
            transIds.append(ensId)
        else:
            raise Exception("not a valid Ensembl id: {}".format(ensId))
    return geneIds, transIds


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
        return [fileOps.decompressCmd(readsFile), readsFile]


class TmpUncompress(object):
    """Wrapper for possibly compressed files to pass to programs that can't read
    them."""
    def __init__(self, dataFile, tmpDir=None):
        self._dataFile = dataFile
        self._tmpUncompressed = None
        if fileOps.isCompressed(dataFile):
            self._tmpUncompressed = self._getTmpFile(dataFile, tmpDir)
            pipettor.run([fileOps.decompressCmd(dataFile), dataFile], stdout=self._tmpUncompressed, logger=logging.getLogger())

    def _getTmpFile(self, dataFile, tmpDir):
        # drop compressed extension, keep next extension
        uncompName = os.path.basename(os.path.splitext(dataFile)[0])
        uncompBase, uncompExt = os.path.splitext(uncompName)
        uncompExt = uncompExt[0:-1] if len(uncompExt) > 0 else None  # drop dot
        return fileOps.tmpFileGet(uncompBase, uncompExt, tmpDir=tmpDir)

    @property
    def path(self):
        "get path to use to access file"
        return self._tmpUncompressed if self._tmpUncompressed is not None else self._dataFile

    def __str__(self):
        return self.path

    def __del__(self):
        self.finish()

    def finish(self):
        """done with file, purge tmp file if create"""
        if self._tmpUncompressed is not None:
            if os.path.exists(self._tmpUncompressed):
                os.unlink(self._tmpUncompressed)
            self._tmpUncompressed = None


def getNewTmpDir(tmpDir):
    """need a temporary directory that doesn't exist, since STAR wants
     the directory not to exist"""
    if tmpDir is not None:
        fileOps.ensureDir(tmpDir)
    cnt = 0
    while True:
        path = os.path.join(fileOps.findTmpDir(tmpDir),
                            "star.{}.{}.{}.tmp".format(socket.gethostname(), os.getpid(), cnt))
        if not os.path.exists(path):
            return path
        cnt += 1
