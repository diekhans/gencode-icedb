"""
Objects use to implement caches of database model objects
"""
from __future__ import print_function
from gencode_icedb.rsl.rslModel import PutativeIntron, MappingParameters


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

    def fetchByLoc(self, chrom, start, end, strand):
        "return record or error if it doesn't exist."
        # currently assume all have been loaded
        key = self.mkLocKey(chrom, start, end, strand)
        pintron = self.byLoc.get(key)
        if pintron is None:
            raise Exception("intron not found in database: {}".format(key))
        return pintron

    def loadAll(self, dbconn):
        "load all introns currently in the database"
        for pintron in PutativeIntron.select().iterator():
            self.add(pintron)


class MappingParametersCache(object):
    """Cache of MappingParameters objects, indexed by """
    def __init__(self, dbconnLoadAll=None):
        "load all if dbconnLoadAll is specified"
        self.bySymId = dict()

    def getBySymId(self, mapping_param_symid):
        """get by symbolic id, caching result, or None if not in database"""
        if mapping_param_symid in self.bySymId:
            return self.bySymId[mapping_param_symid]
        recs = list(MappingParameters.select().where(mapping_param_symid == mapping_param_symid))
        if len(recs) == 0:
            return None
        else:
            self.bySymId[mapping_param_symid] = recs[0]
            return recs[0]

    def fetchBySymId(self, mapping_param_symid):
        """get by symbolic id, caching result, or error if not in database"""
        mappingParams = self.getBySymId(mapping_param_symid)
        if mappingParams is None:
            raise Exception("mapping parameters not found: {}".format(mapping_param_symid))
        return mappingParams

    def create(self, mapping_param_symid, assembly, gene_set, commands, comments):
        rec = MappingParameters(mapping_param_symid=mapping_param_symid, assembly=assembly, gene_set=gene_set, commands=commands, comments=comments)
        rec.save()
        self.bySymId[mapping_param_symid] = rec
        return rec
