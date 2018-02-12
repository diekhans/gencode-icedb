"""
Parse problem cases from set of genbank flat files.
"""
import re
from pycbio.sys import fileOps
from gencode_icedb.tsl.supportDefs import Organism, GenbankProblemReason, EvidenceType
from gencode_icedb.tsl.genbankProblemCasesDb import GenbankProblemCase


class RecordInfoParser(object):
    "parse info from one record"
    __slots__ = ("organism", "etype", "accv", "isNedo", "isAthersysRage", "isOrestes")

    def __init__(self):
        self.etype = None
        self.accv = None
        self.organism = None
        self.isNedo = False
        self.isAthersysRage = False
        self.isOrestes = False

    def __str__(self):
        return "accv={accv} organism={organism} isNedo={isNedo} isAthersysRage={isAthersysRage} isOrestes=${isOrestes}".format(self)

    def _isCloneLibLine(self, line):
        return re.match("^ +/clone_lib=", line) is not None

    def _parseLineForEType(self, line):
        m = re.match("^LOCUS ", line)
        if m:
            self.etype = EvidenceType.EST if line.find(' EST ') > 0 else EvidenceType.RNA

    def _parseLineForVersion(self, line):
        m = re.match("^VERSION     (\S+)", line)
        if m is not None:
            assert self.accv is None
            self.accv = m.group(1)

    def _parseLineForOrganism(self, line):
        m = re.match("^  ORGANISM  (.*)$", line)
        if m is not None:
            assert self.organism is None
            if m.group(1) == "Homo sapiens":
                self.organism = Organism.hs
            elif m.group(1).startswith("Mus musculus"):
                self.organism = Organism.mm

    def parseLine(self, line):
        self._parseLineForEType(line)
        self._parseLineForVersion(line)
        self._parseLineForOrganism(line)
        if re.search("NEDO .*cDNA sequencing project", line):
            self.isNedo = True
        elif self._isCloneLibLine(line):
            if line.find("Athersys RAGE Library") >= 0:
                self.isAthersysRage = True
            elif re.search("(\\s|^)ORESTES(\\s|$)", line) is not None:
                self.isOrestes = True

    def _getReason(self):
        if self.isNedo:
            return GenbankProblemReason.nedo
        elif self.isAthersysRage:
            return GenbankProblemReason.athRage
        elif self.isOrestes:
            return GenbankProblemReason.orestes
        else:
            return None

    def toProblemCase(self):
        "returns None if not a hit"
        if (self.organism is None) or (self.organism not in Organism):
            return None
        else:
            reason = self._getReason()
            if reason is not None:
                return GenbankProblemCase(self.organism, self.etype, self.accv, self.accv, reason)
            else:
                return None


class GbffParser(object):
    def __init__(self, gbffIn, gbffInFh):
        self.gbffIn = gbffIn
        self.gbffInFh = gbffInFh

    def _skipToRecord(self):
        while True:
            line = self.gbffInFh.readline()
            if len(line) == 0:
                return None
            if line.startswith("LOCUS "):
                return line

    def _recordLineReader(self):
        while True:
            line = self.gbffInFh.readline()
            if ((len(line) == 0) or line.startswith("//")):
                break
            yield line

    def _scanRecord(self, locus_line):
        """Scan the next GBFF record, this is crude pattern matching,
        not a true parser."""
        # note, LOCUS first line has been skipped
        recInfo = RecordInfoParser()
        recInfo.parseLine(locus_line)
        for line in self._recordLineReader():
            recInfo.parseLine(line)

        # had bug where we got some bogus short ones once (AK1.1, AK941.1)
        if len(recInfo.accv) < 8:
            raise Exception("short acc parsed from {} ending at byte {}: {}".format(self.gbffIn, self.gbffInFh.tell(), recInfo))
        return recInfo.toProblemCase()

    def scanFile(self):
        problemCases = []
        while True:
            locus_line = self._skipToRecord()
            if locus_line is None:
                break
            problemCase = self._scanRecord(locus_line)
            if problemCase is not None:
                problemCases.append(problemCase)
        return problemCases


def parseGbff(gbffIn):
    with fileOps.opengz(gbffIn, encoding='latin-1', errors='replace') as gbffInFh:
        return GbffParser(gbffIn, gbffInFh).scanFile()
