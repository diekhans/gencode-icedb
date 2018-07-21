"""
GenBank problem cases management
"""
from intervaltree import IntervalTree
from gencode_icedb.tsl.genbankProblemCasesSqlite import GenbankProblemCaseSqliteTable


def accvToAcc(accv):
    """drop version"""
    return accv.split('.')[0]


class GenbankProblemCases(object):
    """Object that looks up GenBank problem cases stored in GenbankProblemCase
    table in database."""
    # Problem cases are stored by accession range.  Since the lists are
    # small, they are loaded into memory and an index by accession prefix

    def __init__(self, conn):
        self.itree = IntervalTree()
        dbTbl = GenbankProblemCaseSqliteTable(conn)
        for gpc in dbTbl.genAll():
            self._addCase(gpc)

    @staticmethod
    def _strAfter(s):
        """Return a accession string that is just after the current string for
        use a an half-open end coordinate"""
        # not general, doesn't handle overflow, but will work with ASCCI
        return s[:-1] + chr(ord(s[-1]) + 1)

    def _addCase(self, gpc):
        # must convert to half-open
        if len(gpc.startAcc) != len(gpc.endAcc):
            raise Exception("Accessions strings must be the same for indexing to work: {}".format(gpc))
        self.itree[gpc.startAcc:self._strAfter(gpc.endAcc)] = gpc

    def isProblem(self, accv):
        """does this accession have a problem """
        return len(self.itree.search(accvToAcc(accv))) > 0

    def getProblem(self, accv):
        """get the problem reason or None if no problem.  If multiple problems,
        a the lowest numbered one is returned."""
        probs = self.itree.search(accvToAcc(accv))
        if len(probs) == 0:
            return None
        else:
            return min([p.data.reason for p in probs])
