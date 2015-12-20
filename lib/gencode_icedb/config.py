"""
Configuration object.
"""

import os

class RnaSeqConfig(object):
    """
    Manage configuration and derived paths.
    """
    def __init__(self, dataRootDir, assembly):
        self.dataRootDir = dataRootDir
        self.assembly = assemlby

    @property
    def assemblyDir(self):
         return os.path.join(self.dataRootDir, self.assembly)

    @property
    def genomeSeqDir(self):
         return os.path.join(self.assemblyDir, "genome", "seq")

    def fastaSeqPath(self, fastaName):
        return os.path.join(self.genomeSeqDir, fastaName)

    def genomeStarDir(self, readLength, annotationSetName):
         return os.path.join(self.assemblyDir, "genome", "star", annotationSetName, str(readLength))

    def genomeStarDoneFlag(self, readLength, annotationSetName):
         return os.path.join(self.genomeStarDir, "done")


