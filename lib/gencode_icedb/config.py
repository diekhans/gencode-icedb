"""
Configuration object.
"""

import os


class PathConfig(object):
    """
    Functions to provide access to data file paths
        root/data/data.db
        root/data/${organism}/genome/${assembly}/${assembly}.fa.gz
                             /rnaSeq/${rnaSeqSet}/
        root/analysis/${assembly}/${geneSet}/star/rl${realLength}
                                            /gencode/
                                            /rnaSeqEvid/${rnaSeqSet}/
    """
    def __init__(self, rootDir, organism, assembly=None, geneSet=None):
        """assembly and geneSet can be None if only using data part"""
        self.rootDir = rootDir
        self.organism = organism
        self.assembly = assembly
        self.geneSet = geneSet
        self.dataDir = os.path.join(rootDir, "data")
        self.analysisDir = os.path.join(rootDir, "analysis")

    def dataDatabase(self):
        return os.path.join(self.dataDir, "data.db")

    def dataOrganismDir(self):
        return os.path.join(self.dataDir, self.organism)

    def genomeDataDir(self):
        return os.path.join(self.dataOrganismDir(), "genome", self.assembly)

    def genomeFasta(self):
        return os.path.join(self.genomeDataDir(), '{}.fa.gz'.format(self.assembly))

    def rnaSeqSetDataDir(self, rnaSeqSet):
        return os.path.join(self.dataOrganismDir(), "rnaSeq", rnaSeqSet)

    def rnaSeqSetDataFile(self, rnaSeqSet, baseName):
        return os.path.join(self.rnaSeqSetDataDir(rnaSeqSet), baseName)

    def geneSetAnalysisDir(self):
        return os.path.join(self.analysisDir, self.assembly, self.geneSet)

    def starGenomeDir(self, readLength):
        return os.path.join(self.geneSetAnalysisDir(), "star", "rl{}".format(readLength))

    def starGenomeDoneFlag(self, readLength):
        return os.path.join(self.starGenomeDir(readLength), "done")

    def rnaSeqSetAnalysisDir(self, rnaSeqSet):
        return os.path.join(self.geneSetAnalysisDir(), "rnaSeqEvid", rnaSeqSet)

    def rnaSeqSetAnalysisSjOut(self, rnaSeqSet, rnaSeqRunName):
        return os.path.join(self.rnaSeqSetAnalysisDir(rnaSeqSet), "{}.sj.tab".format(rnaSeqRunName))

    def rnaSeqSetSupportDir(self):
        return os.path.join(self.geneSetAnalysisDir(), "rnaSeqSupport")

    def rnaSeqSetSupportTsv(self):
        return os.path.join(self.rnaSeqSetSupportDir(), "intron.evid.tsv")
