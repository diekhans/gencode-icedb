"""
Objects use to implement caches of database model objects
"""
from __future__ import print_function
from gencode_icedb.rsl.dbModel import PutativeIntron


class PutativeIntronCache(object):
    """Cache of PutativeIntron objects, indeed by (chrom, start, end, strand)"""
    def __init__(self, dbconnLoadAll=None):
        "load all if dbconnLoadAll is specified"
        self.byLoc = dict()
        if dbconnLoadAll is not None:
            self.loadAll(dbconnLoadAll)

    def __len__(self):
        return len(self.byLoc)

    @staticmethod
    def mkLocKey(chrom, start, end, strand):
        return (chrom, start, end, strand)

    def add(self, pintron):
        "add a new intron to the cache"
        key = self.mkLocKey(pintron.chrom, pintron.start, pintron.end, pintron.strand)
        if key in self.byLoc:
            raise Exception("duplicate intron {}".format(key))
        self.byLoc[key] = pintron

    def loadAll(self, dbconn):
        "load all introns currently in the database"
        for pintron in PutativeIntron.select().iterator():
            self.add(pintron)
