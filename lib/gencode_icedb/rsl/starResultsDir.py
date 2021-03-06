"""
TSV file used to track result files.  Paths are relative to directory containing
the TSV.
"""
import os
from pycbio.tsv import TsvReader


class StarMappingParameters(list):
    """load TSV with mapping parameters"""
    headers = ("mapping_param_symid", "assembly", "gene_set", "commands", "comments")

    def __init__(self, mappingParamsTsv):
        self.byMappingParamSymId = dict()
        for row in TsvReader(mappingParamsTsv):
            self._loadRow(row)

    def _loadRow(self, row):
        if row.mapping_param_symid in self.byMappingParamSymId:
            raise Exception("duplicate mapping_param_symid: {}".format(row.mapping_param_symid))
        self.byMappingParamSymId[row.mapping_param_symid] = row
        self.append(row)

    def getBySymId(self, mapping_param_symid):
        return self.byMappingParamSymId[mapping_param_symid]


class StarResultsDir(list):
    """load results directory TSV, add sjoutPath field for abs path to sjout file"""
    headers = ("run_acc", "mapping_param_symid", "mapping_symid", "sjout",)
    typeMap = {}

    def __init__(self, starResultsTsv):
        self.rootDir = os.path.abspath(os.path.dirname(starResultsTsv))
        for row in TsvReader(starResultsTsv, typeMap=self.typeMap):
            self._loadRow(row)
        if len(self) == 0:
            raise Exception("no data in STAR results TSV: {}".format(starResultsTsv))

    def _loadRow(self, row):
        setattr(row, "sjoutPath", os.path.join(self.rootDir, row.sjout))
        self.append(row)
