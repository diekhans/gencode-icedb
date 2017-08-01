"""
TSV file used to track result files.  Paths are relative to directory containing
the TSV.
"""
from __future__ import print_function
import os
from pycbio.tsv import TsvReader


class StarMappingParameters(list):
    """load TSV with mapping parameters"""
    headers = ("mapping_param_symid", "assembly", "gene_set", "commands", "comments")

    def __init__(self, mappingParamsTsv):
        self.byMappingParamSymId = dict()
        for row in TsvReader(mappingParamsTsv):
            self.__loadRow(row)

    def __loadRow(self, row):
        if row.mapping_param_symid in self.byMappingParamSymId:
            raise Exception("supplicated mapping_param_symid: {}".format(row.mapping_param_symid))
        self.byMappingParamSymId[row.mapping_param_symid] = row
        self.append(row)

    def updateDatabase(self, mappingParamsCache):
        """update the database from the TSV via a MappingParametersCache
        object."""
        for row in self:
            self.__updateDatabaseEntry(mappingParamsCache, row)

    def __updateDatabaseEntry(self, mappingParamsCache, row):
        mappingParams = mappingParamsCache.getBySymId(row.mapping_param_symid)
        if mappingParams is None:
            mappingParamsCache.create(row.mapping_param_symid, row.assembly, row.gene_set, row.commands, row.comments)

    def getBySymId(self, mapping_param_symid):
        return self.byMappingParamSymId[mapping_param_symid]

class StarResultsDir(list):
    """load results directory TSV, add sjoutPath field for abs path to sjout file"""
    headers = ("run_acc", "mapping_param_symid", "mapping_symid", "sjout",)
    typeMap = {}

    def __init__(self, starResultsTsv):
        self.rootDir = os.path.abspath(os.path.dirname(starResultsTsv))
        for row in TsvReader(starResultsTsv, typeMap=self.typeMap):
            self.__loadRow(row)
        if len(self) == 0:
            raise Exception("no data in STAR results TSV: {}".format(starResultsTsv))

    def __loadRow(self, row):
        setattr(row, "sjoutPath", os.path.join(self.rootDir, row.sjout))
        self.append(row)
