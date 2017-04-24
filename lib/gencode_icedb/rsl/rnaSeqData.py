"""
Data associate with RNA seq database, store in an sqlite database
"""
import os
from peewee import Proxy, Model, IntegerField, CharField
from pycbio.sys.symEnum import SymEnum
import pipettor

database_proxy = Proxy()


def rslEstimateReadLength(readsPath):
    "run rslEstimateReadLength program to get read length from read data"
    return int(pipettor.runout(["rslEstimateReadLength", readsPath, "/dev/stdout"]))


def setDatabase(database):
    "bind the proxy to a database"
    database_proxy.initialize(database)


class RnaSeqData(Model):
    """specification of a single set of RNA-Seq reads.
    name is a symbolic, unique and is used to create output directories.  Reads can be
    SAM, BAM, or fastq and can be compressed.  readsFile path is relative
    to data root."""
    id = IntegerField(primary_key=True)
    setname = CharField(index=True)
    runname = CharField(unique=True)  # important this is unique
    organism = CharField(index=True)
    readsfileurl = CharField(unique=True)
    readsfile = CharField()
    readsfile2url = CharField(unique=True, null=True)
    readsfile2 = CharField(null=True)
    readlength = IntegerField()
    description = CharField()
    tissue = CharField()

    class Meta:
        database = database_proxy


class RnaSeqDataLoader(object):
    "load into RnaSeqData table"

    class Status(SymEnum):
        "status of added a record"
        added = 1
        skipped = 2

    def __init__(self, database, pathConfig, rnaSeqSet, organism, skipExisting):
        setDatabase(database)
        self.pathConfig = pathConfig
        self.rnaSeqSet = rnaSeqSet
        self.organism = organism
        self.skipExisting = skipExisting

    def __validateReadsFilePath(self, readsFile):
        """validate that the name of reads file is valid and file
        is in the right location"""
        if len(os.path.split(readsFile)[0]) > 0:
            raise Exception("RNASeq file name should be a base file name without directory: " + readsFile)
        readsPath = self.pathConfig.rnaSeqSetDataFile(self.rnaSeqSet, readsFile)
        if not os.path.exists(readsPath):
            raise Exception("RNASeq file does not exist in expected directory: " + readsPath)
        return readsPath

    def __addRecord(self, rnaSeqRun,
                    description, tissue,
                    readsFileUrl, readsFile,
                    readsFile2Url, readsFile2,
                    readLength):
        rec = RnaSeqData(
            setname=self.rnaSeqSet,
            runname=rnaSeqRun,
            organism=self.organism,
            readsfileurl=readsFileUrl,
            readsfile=readsFile,
            readsfile2url=readsFile2Url,
            readsfile2=readsFile2,
            readlength=readLength,
            description=description,
            tissue=tissue)
        rec.save()

    def __checkForExisting(self, rnaSeqRun):
        # FIXME: should validate?
        return RnaSeqData.select().where(RnaSeqData.setname == self.rnaSeqSet,
                                         RnaSeqData.runname == rnaSeqRun).count() > 0

    def registerRun(self, rnaSeqRun,
                    description, tissue,
                    readsFileUrl, readsFile,
                    readsFile2Url, readsFile2):
        readsPath = self.__validateReadsFilePath(readsFile)
        if readsFile2 is not None:
            self.__validateReadsFilePath(readsFile2)
        if not RnaSeqData.table_exists():
            RnaSeqData.create_table()
        exists = self.__checkForExisting(rnaSeqRun)
        if self.skipExisting and exists:
            return RnaSeqDataLoader.Status.skipped
        else:
            self.__addRecord(rnaSeqRun,
                             description, tissue,
                             readsFileUrl, readsFile,
                             readsFile2Url, readsFile2,
                             rslEstimateReadLength(readsPath))
            return RnaSeqDataLoader.Status.added
