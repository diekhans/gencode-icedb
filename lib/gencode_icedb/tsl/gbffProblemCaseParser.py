"""
Parse problem cases from set of genbank flat files.
"""
import re
from collections import namedtuple
from pycbio.sys import fileOps, symEnum
from gencode_icedb.tsl.genbankProblemCasesDb import GenbankProblemReason
import multiprocessing as mp


class Organism(symEnum.SymEnum):
    hs = 1
    mm = 2


class GenbankOrgProblemCase(namedtuple("GenbankOrgProblemCase",
                                       ("organism", "acc", "reason"))):
    """a problem case for an organism"""
    __slots__ = ()
    pass


class RecordInfoParser(object):
    "parse info from one record"
    __slots__ = ("accv", "organism", "isNedo", "isAthersysRage", "isOrestes")

    def __init__(self):
        self.accv = None
        self.organism = None
        self.isNedo = False
        self.isAthersysRage = False
        self.isOrestes = False

    def __str__(self):
        return "accv={accv} organism={organism} isNedo={isNedo} isAthersysRage={isAthersysRage} isOrestes=${isOrestes}".format(self)

    def __isCloneLibLine(self, line):
        return re.match("^ +/clone_lib=", line) is not None

    def __parseLineForVersion(self, line):
        m = re.match("^VERSION     (\S+)", line)
        if m is not None:
            assert self.accv is None
            self.accv = m.group(1)

    def __parseLineForOrganism(self, line):
        m = re.match("^  ORGANISM  (.*)$", line)
        if m is not None:
            assert self.organism is None
            if m.group(1) == "Homo sapiens":
                self.organism = Organism.hs
            elif m.group(1).startswith("Mus musculus"):
                self.organism = Organism.mm

    def parseLine(self, line):
        self.__parseLineForVersion(line)
        self.__parseLineForOrganism(line)
        if re.search("NEDO .*cDNA sequencing project", line):
            self.isNedo = True
        elif self.__isCloneLibLine(line):
            if line.find("Athersys RAGE Library") >= 0:
                self.isAthersysRage = True
            elif re.search("(\\s|^)ORESTES(\\s|$)", line) is not None:
                self.isOrestes = True

    def toProblemCase(self):
        "returns None if not a hit"
        if self.organism is None:
            return None
        elif self.isNedo:
            return GenbankOrgProblemCase(self.organism, self.accv, GenbankProblemReason.nedo)
        elif self.isAthersysRage:
            return GenbankOrgProblemCase(self.organism, self.accv, GenbankProblemReason.athRage)
        elif self.isOrestes:
            return GenbankOrgProblemCase(self.organism, self.accv, GenbankProblemReason.orestes)
        else:
            return None


class GbffParser(object):
    def __init__(self, gbffInFh):
        self.gbffInFh = gbffInFh

    def __skipToRecord(self):
        while True:
            line = self.gbffInFh.readline()
            if len(line) == 0:
                return None
            if line.startswith("LOCUS "):
                return line

    def __recordLineReader(self):
        while True:
            line = self.gbffInFh.readline()
            if ((len(line) == 0) or line.startswith("//")):
                break
            yield line

    def __scanRecord(self):
        """Scan the next GBFF record, this is crude pattern matching,
        not a true parser."""
        # note, LOCUS first line has been skipped
        recInfo = RecordInfoParser()
        for line in self.__recordLineReader():
            recInfo.parseLine(line)
        return recInfo.toProblemCase()

    def scanFile(self):
        problemCases = []
        while self.__skipToRecord():
            problemCase = self.__scanRecord()
            if problemCase is not None:
                problemCases.append(problemCase)
        return problemCases


class GbffProblemCaseParser(object):
    "parses problem cases and save"

    def __init__(self, gbffInputs, maxProcesses=2):
        self.problemEntries = []
        self.maxProcesses = maxProcesses
        for gbffIn in gbffInputs:
            self.__scanGbffFile(gbffIn)

    def __scanGbffFile(self, gbffIn):
        with fileOps.opengz(gbffIn) as gbffInFh:
            parser = GbffParser(gbffInFh)
            for problemEntry in parser.scanFile():
                self.problemEntries.append(problemEntry)


def parseGbff(gbffIn):
    with fileOps.opengz(gbffIn) as gbffInFh:
        return GbffParser(gbffInFh).scanFile()


def gbffProblemCaseParse(gbffInputs, maxProcesses):
    pool = mp.Pool(processes=maxProcesses)
    results = [pool.apply_async(parseGbff, args=(gbffIn,)) for gbffIn in gbffInputs]
    problemCases = []
    for result in results:
        problemCases.extend(result.get())
    return problemCases
