"""
Objects use to implement caches of database model objects
"""
from __future__ import print_function
from gencode_icedb.rsl.rslModel import MappingParameters


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
