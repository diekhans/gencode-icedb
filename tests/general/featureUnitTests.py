from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.objDict import ObjDict
from pycbio.sys.testCaseBase import TestCaseBase
from gencode_icedb.general.genome import GenomeReaderFactory
from gencode_icedb.general.transFeatures import ExonFeature
from gencode_icedb.general.transFeatures import CdsRegionFeature, Utr3RegionFeature
from gencode_icedb.general.transFeatures import RnaInsertFeature, ChromInsertFeature
from gencode_icedb.general.evidFeatures import EvidencePslFactory
from gencode_icedb.general.annotFeatures import AnnotationGenePredFactory
from pycbio.hgdata.hgLite import sqliteConnect
from pycbio.hgdata.hgLite import PslDbTable, GenePredDbTable
from pycbio.hgdata.genePred import GenePredReader
from pycbio.sys.pprint2 import nswpprint


debugResults = False   # print out results for updated expected
noCheckResults = False  # don't check results

if debugResults or noCheckResults:
    print("WARNING: debug variables set", file=sys.stderr)


def getInputFile(base):
    "from input relative to test file"
    return os.path.join(os.path.dirname(__file__), "input", base)


class GenomeSeqSrc(object):
    "caching genome reader"

    srcs = {
        "hg19": "/hive/data/genomes/hg19/hg19.2bit",
        "mm10": "/hive/data/genomes/mm10/mm10.2bit",
    }
    readers = {}

    @classmethod
    def obtain(cls, db):
        reader = cls.readers.get(db)
        if reader is None:
            reader = cls.readers[db] = GenomeReaderFactory(cls.srcs[db]).obtain()
        return reader


class PslDbSrc(object):
    "caching psl sqlite database"
    srcs = {
        "set1": "set1.ucsc-mrna.psl",
        "hg38-mm10.transMap": "hg38-mm10.transMap.psl"
    }
    pslTbls = {}

    @classmethod
    def __loadDbTbl(cls, name):
        conn = sqliteConnect(None)
        dbTbl = PslDbTable(conn, "psls", create=True)
        dbTbl.loadPslFile(getInputFile(cls.srcs[name]))
        cls.pslTbls[name] = dbTbl
        return dbTbl

    @classmethod
    def obtain(cls, name):
        dbTbl = cls.pslTbls.get(name)
        if dbTbl is None:
            dbTbl = cls.__loadDbTbl(name)
        return dbTbl

    @classmethod
    def obtainPsl(cls, name, acc):
        psls = cls.obtain(name).getByQName(acc)
        if len(psls) == 0:
            raise Exception("psl not found: {}".format(acc))
        return psls[0]


class GenePredDbSrc(object):
    "caching genpred sqlite database"
    srcs = {
        "set1": "set1.gencodeCompV19.gp"
    }
    genePredTbls = {}

    @classmethod
    def __loadDbTbl(cls, name):
        conn = sqliteConnect(None)
        dbTbl = GenePredDbTable(conn, name, create=True)
        dbTbl.loadGenePredFile(getInputFile(cls.srcs[name]))
        cls.genePredTbls[name] = dbTbl
        return dbTbl

    @classmethod
    def obtain(cls, name):
        dbTbl = cls.genePredTbls.get(name)
        if dbTbl is None:
            dbTbl = cls.__loadDbTbl(name)
        return dbTbl

    @classmethod
    def obtainGenePred(cls, name, acc):
        gps = list(cls.obtain(name).getByName(acc))
        if len(gps) == 0:
            raise Exception("genePred not found: {}".format(acc))
        return gps[0]


class FeatureTestBase(TestCaseBase):
    def _assertFeatures(self, trans, expect):
        if debugResults:
            print("==== {} ==== ".format(self.id()))
            nswpprint(trans.toStrTree(), indent=2)
        if not noCheckResults:
            self.assertEqual(expect, trans.toStrTree())


class EvidenceTests(FeatureTestBase):
    def __getSet1Psl(self, acc):
        return PslDbSrc.obtainPsl("set1", acc)

    def __pslToEvidTranscript(self, psl):
        return EvidencePslFactory(GenomeSeqSrc.obtain("hg19")).fromPsl(psl)

    def __checkRnaAln(self, trans):
        "validate that full RNA is covered by alignment"
        prevChromEnd = trans.chrom.start
        prevRnaEnd = trans.rna.start
        for feat in trans.features:
            if feat.alignFeatures is not None:
                for alnFeat in feat.alignFeatures:
                    if alnFeat.chrom is not None:
                        self.assertEqual(alnFeat.chrom.start, prevChromEnd,
                                         msg="trans: {} feat: {} aln: {} chrom.start {} != prevChromEnd {}".format(trans.rna.name, feat, alnFeat, alnFeat.chrom.start, prevChromEnd))
                        prevChromEnd = alnFeat.chrom.end
                    if alnFeat.rna is not None:
                        self.assertEqual(alnFeat.rna.start, prevRnaEnd,
                                         msg="trans: {} feat: {} aln: {} rna.start {} != prevRnaEnd {}".format(trans.rna.name, feat, alnFeat, alnFeat.rna.start, prevRnaEnd))
                        prevRnaEnd = alnFeat.rna.end
        self.assertEqual(prevChromEnd, trans.chrom.end,
                         msg="trans: {} alignments don't cover chrom range".format(trans.rna.name))
        self.assertEqual(prevRnaEnd, trans.rna.end,
                         msg="trans: {} alignments don't cover RNA range".format(trans.rna.name))

    def testAF010310(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18900294-18905926/+, rna=AF010310.1:13-901/- 901',
                              (('exon 18900294-18900875 rna=13-618',
                                (('aln 18900294-18900452 rna=13-171',),
                                 ('rins None-None rna=171-174',),
                                 ('cins 18900452-18900453 rna=None-None',),
                                 ('aln 18900453-18900535 rna=174-256',),
                                 ('rins None-None rna=256-257',),
                                 ('aln 18900535-18900565 rna=257-287',),
                                 ('rins None-None rna=287-288',),
                                 ('aln 18900565-18900567 rna=288-290',),
                                 ('rins None-None rna=290-291',),
                                 ('aln 18900567-18900577 rna=291-301',),
                                 ('rins None-None rna=301-302',),
                                 ('aln 18900577-18900579 rna=302-304',),
                                 ('rins None-None rna=304-305',),
                                 ('aln 18900579-18900606 rna=305-332',),
                                 ('rins None-None rna=332-335',),
                                 ('cins 18900606-18900607 rna=None-None',),
                                 ('aln 18900607-18900611 rna=335-339',),
                                 ('rins None-None rna=339-340',),
                                 ('aln 18900611-18900616 rna=340-345',),
                                 ('rins None-None rna=345-346',),
                                 ('aln 18900616-18900621 rna=346-351',),
                                 ('rins None-None rna=351-352',),
                                 ('aln 18900621-18900628 rna=352-359',),
                                 ('rins None-None rna=359-369',),
                                 ('cins 18900628-18900634 rna=None-None',),
                                 ('aln 18900634-18900643 rna=369-378',),
                                 ('rins None-None rna=378-379',),
                                 ('aln 18900643-18900645 rna=379-381',),
                                 ('rins None-None rna=381-382',),
                                 ('aln 18900645-18900664 rna=382-401',),
                                 ('rins None-None rna=401-402',),
                                 ('aln 18900664-18900669 rna=402-407',),
                                 ('rins None-None rna=407-408',),
                                 ('aln 18900669-18900762 rna=408-501',),
                                 ('rins None-None rna=501-502',),
                                 ('aln 18900762-18900790 rna=502-530',),
                                 ('rins None-None rna=530-531',),
                                 ('aln 18900790-18900803 rna=531-544',),
                                 ('rins None-None rna=544-545',),
                                 ('aln 18900803-18900871 rna=545-613',),
                                 ('rins None-None rna=613-614',),
                                 ('aln 18900871-18900875 rna=614-618',))),
                               ('intron 18900875-18900950 rna=618-618 sjBases=GT...AG (GT_AG)',
                                (('cins 18900875-18900950 rna=None-None',),)),
                               ('exon 18900950-18901039 rna=618-707',
                                (('aln 18900950-18901039 rna=618-707',),)),
                               ('intron 18901039-18904402 rna=707-707 sjBases=GT...AG (GT_AG)',
                                (('cins 18901039-18904402 rna=None-None',),)),
                               ('exon 18904402-18904501 rna=707-806',
                                (('aln 18904402-18904501 rna=707-806',),)),
                               ('intron 18904501-18905828 rna=806-806 sjBases=GT...AG (GT_AG)',
                                (('cins 18904501-18905828 rna=None-None',),)),
                               ('exon 18905828-18905926 rna=806-901',
                                (('aln 18905828-18905889 rna=806-867',),
                                 ('rins None-None rna=867-868',),
                                 ('cins 18905889-18905893 rna=None-None',),
                                 ('aln 18905893-18905926 rna=868-901',))))))
        self.__checkRnaAln(trans)

    def testX96484(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("X96484.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18893922-18899592/+, rna=X96484.1:48-1067/+ 1080',
                              (('exon 18893922-18893997 rna=48-123',
                                (('aln 18893922-18893997 rna=48-123',),)),
                               ('intron 18893997-18894077 rna=123-123 sjBases=GT...AG (GT_AG)',
                                (('cins 18893997-18894077 rna=None-None',),)),
                               ('exon 18894077-18894238 rna=123-285',
                                (('aln 18894077-18894173 rna=123-219',),
                                 ('rins None-None rna=219-220',),
                                 ('aln 18894173-18894238 rna=220-285',))),
                               ('intron 18894238-18897684 rna=285-285 sjBases=GT...AG (GT_AG)',
                                (('cins 18894238-18897684 rna=None-None',),)),
                               ('exon 18897684-18897785 rna=285-386',
                                (('aln 18897684-18897785 rna=285-386',),)),
                               ('intron 18897785-18898400 rna=386-386 sjBases=GT...AG (GT_AG)',
                                (('cins 18897785-18898400 rna=None-None',),)),
                               ('exon 18898400-18898541 rna=386-527',
                                (('aln 18898400-18898541 rna=386-527',),)),
                               ('intron 18898541-18899052 rna=527-527 sjBases=GT...AG (GT_AG)',
                                (('cins 18898541-18899052 rna=None-None',),)),
                               ('exon 18899052-18899592 rna=527-1067',
                                (('aln 18899052-18899592 rna=527-1067',),)))))
        self.__checkRnaAln(trans)

    def testX96484NoSJ(self):
        trans = EvidencePslFactory(None).fromPsl(self.__getSet1Psl("X96484.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18893922-18899592/+, rna=X96484.1:48-1067/+ 1080',
                              (('exon 18893922-18893997 rna=48-123',
                                (('aln 18893922-18893997 rna=48-123',),)),
                               ('intron 18893997-18894077 rna=123-123 sjBases=None',
                                (('cins 18893997-18894077 rna=None-None',),)),
                               ('exon 18894077-18894238 rna=123-285',
                                (('aln 18894077-18894173 rna=123-219',),
                                 ('rins None-None rna=219-220',),
                                 ('aln 18894173-18894238 rna=220-285',))),
                               ('intron 18894238-18897684 rna=285-285 sjBases=None',
                                (('cins 18894238-18897684 rna=None-None',),)),
                               ('exon 18897684-18897785 rna=285-386',
                                (('aln 18897684-18897785 rna=285-386',),)),
                               ('intron 18897785-18898400 rna=386-386 sjBases=None',
                                (('cins 18897785-18898400 rna=None-None',),)),
                               ('exon 18898400-18898541 rna=386-527',
                                (('aln 18898400-18898541 rna=386-527',),)),
                               ('intron 18898541-18899052 rna=527-527 sjBases=None',
                                (('cins 18898541-18899052 rna=None-None',),)),
                               ('exon 18899052-18899592 rna=527-1067',
                                (('aln 18899052-18899592 rna=527-1067',),)))))
        self.__checkRnaAln(trans)

    def testTransMapDropExon(self):
        # internal exon was not mapped, causing an intron to contain unaligned
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        self._assertFeatures(trans,
                             ('t=chr4:148039043-148056154/+, rna=ENST00000641446:0-2820/+ 2820',
                              (('exon 148039043-148039141 rna=0-106',
                                (('aln 148039043-148039050 rna=0-7',),
                                 ('cins 148039050-148039051 rna=None-None',),
                                 ('aln 148039051-148039059 rna=7-15',),
                                 ('rins None-None rna=15-16',),
                                 ('aln 148039059-148039079 rna=16-36',),
                                 ('rins None-None rna=36-44',),
                                 ('aln 148039079-148039141 rna=44-106',))),
                               ('intron 148039141-148041583 rna=106-115 sjBases=aa...ag (unknown)',
                                (('rins None-None rna=106-115',),
                                 ('cins 148039141-148041583 rna=None-None',))),
                               ('exon 148041583-148041829 rna=115-364',
                                (('aln 148041583-148041649 rna=115-181',),
                                 ('rins None-None rna=181-184',),
                                 ('aln 148041649-148041829 rna=184-364',))),
                               ('intron 148041829-148043424 rna=364-364 sjBases=GT...AG (GT_AG)',
                                (('cins 148041829-148043424 rna=None-None',),)),
                               ('exon 148043424-148043663 rna=364-603',
                                (('aln 148043424-148043663 rna=364-603',),)),
                               ('intron 148043663-148044443 rna=603-603 sjBases=GT...AG (GT_AG)',
                                (('cins 148043663-148044443 rna=None-None',),)),
                               ('exon 148044443-148044554 rna=603-714',
                                (('aln 148044443-148044554 rna=603-714',),)),
                               ('intron 148044554-148048072 rna=714-714 sjBases=GT...AG (GT_AG)',
                                (('cins 148044554-148048072 rna=None-None',),)),
                               ('exon 148048072-148048266 rna=714-908',
                                (('aln 148048072-148048266 rna=714-908',),)),
                               ('intron 148048266-148051371 rna=908-908 sjBases=GT...AG (GT_AG)',
                                (('cins 148048266-148051371 rna=None-None',),)),
                               ('exon 148051371-148051622 rna=908-1159',
                                (('aln 148051371-148051622 rna=908-1159',),)),
                               ('intron 148051622-148051864 rna=1159-1159 sjBases=GT...AG (GT_AG)',
                                (('cins 148051622-148051864 rna=None-None',),)),
                               ('exon 148051864-148051999 rna=1159-1294',
                                (('aln 148051864-148051999 rna=1159-1294',),)),
                               ('intron 148051999-148052098 rna=1294-1294 sjBases=GT...AG (GT_AG)',
                                (('cins 148051999-148052098 rna=None-None',),)),
                               ('exon 148052098-148052279 rna=1294-1475',
                                (('aln 148052098-148052279 rna=1294-1475',),)),
                               ('intron 148052279-148052513 rna=1475-1475 sjBases=GT...AG (GT_AG)',
                                (('cins 148052279-148052513 rna=None-None',),)),
                               ('exon 148052513-148052696 rna=1475-1658',
                                (('aln 148052513-148052696 rna=1475-1658',),)),
                               ('intron 148052696-148053579 rna=1658-1658 sjBases=GT...AG (GT_AG)',
                                (('cins 148052696-148053579 rna=None-None',),)),
                               ('exon 148053579-148053681 rna=1658-1760',
                                (('aln 148053579-148053681 rna=1658-1760',),)),
                               ('intron 148053681-148054984 rna=1760-1869 sjBases=GT...AG (GT_AG)',
                                (('rins None-None rna=1760-1869',),
                                 ('cins 148053681-148054984 rna=None-None',))),
                               ('exon 148054984-148055104 rna=1869-1989',
                                (('aln 148054984-148055104 rna=1869-1989',),)),
                               ('intron 148055104-148055372 rna=1989-1989 sjBases=GT...AG (GT_AG)',
                                (('cins 148055104-148055372 rna=None-None',),)),
                               ('exon 148055372-148055597 rna=1989-2217',
                                (('aln 148055372-148055566 rna=1989-2183',),
                                 ('rins None-None rna=2183-2186',),
                                 ('aln 148055566-148055597 rna=2186-2217',))),
                               ('intron 148055597-148055703 rna=2217-2297 sjBases=tc...cc (unknown)',
                                (('rins None-None rna=2217-2297',),
                                 ('cins 148055597-148055703 rna=None-None',))),
                               ('exon 148055703-148056110 rna=2297-2750',
                                (('aln 148055703-148055739 rna=2297-2333',),
                                 ('cins 148055739-148055749 rna=None-None',),
                                 ('aln 148055749-148055803 rna=2333-2387',),
                                 ('rins None-None rna=2387-2388',),
                                 ('aln 148055803-148055878 rna=2388-2463',),
                                 ('rins None-None rna=2463-2464',),
                                 ('aln 148055878-148055884 rna=2464-2470',),
                                 ('rins None-None rna=2470-2500',),
                                 ('aln 148055884-148055903 rna=2500-2519',),
                                 ('rins None-None rna=2519-2523',),
                                 ('aln 148055903-148055932 rna=2523-2552',),
                                 ('rins None-None rna=2552-2568',),
                                 ('aln 148055932-148055976 rna=2568-2612',),
                                 ('rins None-None rna=2612-2617',),
                                 ('aln 148055976-148055994 rna=2617-2635',),
                                 ('rins None-None rna=2635-2637',),
                                 ('aln 148055994-148055997 rna=2637-2640',),
                                 ('cins 148055997-148055998 rna=None-None',),
                                 ('aln 148055998-148056008 rna=2640-2650',),
                                 ('rins None-None rna=2650-2651',),
                                 ('aln 148056008-148056013 rna=2651-2656',),
                                 ('cins 148056013-148056018 rna=None-None',),
                                 ('aln 148056018-148056040 rna=2656-2678',),
                                 ('cins 148056040-148056041 rna=None-None',),
                                 ('aln 148056041-148056046 rna=2678-2683',),
                                 ('cins 148056046-148056048 rna=None-None',),
                                 ('aln 148056048-148056065 rna=2683-2700',),
                                 ('rins None-None rna=2700-2703',),
                                 ('aln 148056065-148056089 rna=2703-2727',),
                                 ('rins None-None rna=2727-2729',),
                                 ('aln 148056089-148056110 rna=2729-2750',))),
                               ('intron 148056110-148056150 rna=2750-2816 sjBases=ca...tg (unknown)',
                                (('rins None-None rna=2750-2816',),
                                 ('cins 148056110-148056150 rna=None-None',))),
                               ('exon 148056150-148056154 rna=2816-2820',
                                (('aln 148056150-148056154 rna=2816-2820',),)))))
        self.__checkRnaAln(trans)

    def testTransMapDropExonRc(self):
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        transRc = trans.reverseComplement()
        self._assertFeatures(transRc,
                             ('t=chr4:8451962-8469073/-, rna=ENST00000641446:0-2820/- 2820',
                              (('exon 8451962-8451966 rna=0-4', (('aln 8451962-8451966 rna=0-4',),)),
                               ('intron 8451966-8452006 rna=4-70 sjBases=ca...tg (unknown)',
                                (('cins 8451966-8452006 rna=None-None',), ('rins None-None rna=4-70',))),
                               ('exon 8452006-8452413 rna=70-523',
                                (('aln 8452006-8452027 rna=70-91',),
                                 ('rins None-None rna=91-93',),
                                 ('aln 8452027-8452051 rna=93-117',),
                                 ('rins None-None rna=117-120',),
                                 ('aln 8452051-8452068 rna=120-137',),
                                 ('cins 8452068-8452070 rna=None-None',),
                                 ('aln 8452070-8452075 rna=137-142',),
                                 ('cins 8452075-8452076 rna=None-None',),
                                 ('aln 8452076-8452098 rna=142-164',),
                                 ('cins 8452098-8452103 rna=None-None',),
                                 ('aln 8452103-8452108 rna=164-169',),
                                 ('rins None-None rna=169-170',),
                                 ('aln 8452108-8452118 rna=170-180',),
                                 ('cins 8452118-8452119 rna=None-None',),
                                 ('aln 8452119-8452122 rna=180-183',),
                                 ('rins None-None rna=183-185',),
                                 ('aln 8452122-8452140 rna=185-203',),
                                 ('rins None-None rna=203-208',),
                                 ('aln 8452140-8452184 rna=208-252',),
                                 ('rins None-None rna=252-268',),
                                 ('aln 8452184-8452213 rna=268-297',),
                                 ('rins None-None rna=297-301',),
                                 ('aln 8452213-8452232 rna=301-320',),
                                 ('rins None-None rna=320-350',),
                                 ('aln 8452232-8452238 rna=350-356',),
                                 ('rins None-None rna=356-357',),
                                 ('aln 8452238-8452313 rna=357-432',),
                                 ('rins None-None rna=432-433',),
                                 ('aln 8452313-8452367 rna=433-487',),
                                 ('cins 8452367-8452377 rna=None-None',),
                                 ('aln 8452377-8452413 rna=487-523',))),
                               ('intron 8452413-8452519 rna=523-603 sjBases=tc...cc (unknown)',
                                (('cins 8452413-8452519 rna=None-None',),
                                 ('rins None-None rna=523-603',))),
                               ('exon 8452519-8452744 rna=603-831',
                                (('aln 8452519-8452550 rna=603-634',),
                                 ('rins None-None rna=634-637',),
                                 ('aln 8452550-8452744 rna=637-831',))),
                               ('intron 8452744-8453012 rna=831-831 sjBases=GT...AG (GT_AG)',
                                (('cins 8452744-8453012 rna=None-None',),)),
                               ('exon 8453012-8453132 rna=831-951',
                                (('aln 8453012-8453132 rna=831-951',),)),
                               ('intron 8453132-8454435 rna=951-1060 sjBases=GT...AG (GT_AG)',
                                (('cins 8453132-8454435 rna=None-None',),
                                 ('rins None-None rna=951-1060',))),
                               ('exon 8454435-8454537 rna=1060-1162',
                                (('aln 8454435-8454537 rna=1060-1162',),)),
                               ('intron 8454537-8455420 rna=1162-1162 sjBases=GT...AG (GT_AG)',
                                (('cins 8454537-8455420 rna=None-None',),)),
                               ('exon 8455420-8455603 rna=1162-1345',
                                (('aln 8455420-8455603 rna=1162-1345',),)),
                               ('intron 8455603-8455837 rna=1345-1345 sjBases=GT...AG (GT_AG)',
                                (('cins 8455603-8455837 rna=None-None',),)),
                               ('exon 8455837-8456018 rna=1345-1526',
                                (('aln 8455837-8456018 rna=1345-1526',),)),
                               ('intron 8456018-8456117 rna=1526-1526 sjBases=GT...AG (GT_AG)',
                                (('cins 8456018-8456117 rna=None-None',),)),
                               ('exon 8456117-8456252 rna=1526-1661',
                                (('aln 8456117-8456252 rna=1526-1661',),)),
                               ('intron 8456252-8456494 rna=1661-1661 sjBases=GT...AG (GT_AG)',
                                (('cins 8456252-8456494 rna=None-None',),)),
                               ('exon 8456494-8456745 rna=1661-1912',
                                (('aln 8456494-8456745 rna=1661-1912',),)),
                               ('intron 8456745-8459850 rna=1912-1912 sjBases=GT...AG (GT_AG)',
                                (('cins 8456745-8459850 rna=None-None',),)),
                               ('exon 8459850-8460044 rna=1912-2106',
                                (('aln 8459850-8460044 rna=1912-2106',),)),
                               ('intron 8460044-8463562 rna=2106-2106 sjBases=GT...AG (GT_AG)',
                                (('cins 8460044-8463562 rna=None-None',),)),
                               ('exon 8463562-8463673 rna=2106-2217',
                                (('aln 8463562-8463673 rna=2106-2217',),)),
                               ('intron 8463673-8464453 rna=2217-2217 sjBases=GT...AG (GT_AG)',
                                (('cins 8463673-8464453 rna=None-None',),)),
                               ('exon 8464453-8464692 rna=2217-2456',
                                (('aln 8464453-8464692 rna=2217-2456',),)),
                               ('intron 8464692-8466287 rna=2456-2456 sjBases=GT...AG (GT_AG)',
                                (('cins 8464692-8466287 rna=None-None',),)),
                               ('exon 8466287-8466533 rna=2456-2705',
                                (('aln 8466287-8466467 rna=2456-2636',),
                                 ('rins None-None rna=2636-2639',),
                                 ('aln 8466467-8466533 rna=2639-2705',))),
                               ('intron 8466533-8468975 rna=2705-2714 sjBases=aa...ag (unknown)',
                                (('cins 8466533-8468975 rna=None-None',),
                                 ('rins None-None rna=2705-2714',))),
                               ('exon 8468975-8469073 rna=2714-2820',
                                (('aln 8468975-8469037 rna=2714-2776',),
                                 ('rins None-None rna=2776-2784',),
                                 ('aln 8469037-8469057 rna=2784-2804',),
                                 ('rins None-None rna=2804-2805',),
                                 ('aln 8469057-8469065 rna=2805-2813',),
                                 ('cins 8469065-8469066 rna=None-None',),
                                 ('aln 8469066-8469073 rna=2813-2820',))))))
        self.__checkRnaAln(transRc)

    def testGetAlignmentFeaturesOfType(self):
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        feats = trans.getAlignmentFeaturesOfType((RnaInsertFeature, ChromInsertFeature))
        featStrs = tuple([str(f) for f in feats])
        self.assertEquals(featStrs,
                          ('cins 148039050-148039051 rna=None-None',
                           'rins None-None rna=15-16',
                           'rins None-None rna=36-44',
                           'rins None-None rna=106-115',
                           'cins 148039141-148041583 rna=None-None',
                           'rins None-None rna=181-184',
                           'cins 148041829-148043424 rna=None-None',
                           'cins 148043663-148044443 rna=None-None',
                           'cins 148044554-148048072 rna=None-None',
                           'cins 148048266-148051371 rna=None-None',
                           'cins 148051622-148051864 rna=None-None',
                           'cins 148051999-148052098 rna=None-None',
                           'cins 148052279-148052513 rna=None-None',
                           'cins 148052696-148053579 rna=None-None',
                           'rins None-None rna=1760-1869',
                           'cins 148053681-148054984 rna=None-None',
                           'cins 148055104-148055372 rna=None-None',
                           'rins None-None rna=2183-2186',
                           'rins None-None rna=2217-2297',
                           'cins 148055597-148055703 rna=None-None',
                           'cins 148055739-148055749 rna=None-None',
                           'rins None-None rna=2387-2388',
                           'rins None-None rna=2463-2464',
                           'rins None-None rna=2470-2500',
                           'rins None-None rna=2519-2523',
                           'rins None-None rna=2552-2568',
                           'rins None-None rna=2612-2617',
                           'rins None-None rna=2635-2637',
                           'cins 148055997-148055998 rna=None-None',
                           'rins None-None rna=2650-2651',
                           'cins 148056013-148056018 rna=None-None',
                           'cins 148056040-148056041 rna=None-None',
                           'cins 148056046-148056048 rna=None-None',
                           'rins None-None rna=2700-2703',
                           'rins None-None rna=2727-2729',
                           'rins None-None rna=2750-2816',
                           'cins 148056110-148056150 rna=None-None'))

    def testExonRnaOverlap(self):
        aln1 = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        aln2 = self.__pslToEvidTranscript(self.__getSet1Psl("AF120278.1"))
        exon1 = aln1.features[0]
        exon2 = aln2.features[0]
        self.assertTrue(exon1.rnaOverlaps(exon2))
        exon2b = aln2.features[3]
        self.assertFalse(exon1.rnaOverlaps(exon2b))
        self.__checkRnaAln(aln1)
        self.__checkRnaAln(aln2)


class AnnotationTests(FeatureTestBase):
    def __getSet1Gp(self, acc):
        return GenePredDbSrc.obtainGenePred("set1", acc)

    def __gpToAnnotTranscript(self, gp):
        factory = AnnotationGenePredFactory(GenomeSeqSrc.obtain("hg19"))
        return factory.fromGenePred(gp)

    def testENST00000215794(self):
        # + strand
        trans = self.__gpToAnnotTranscript(self.__getSet1Gp("ENST00000215794.7"))
        self._assertFeatures(trans,
                             ('t=chr22:18632665-18660164/+, rna=ENST00000215794.7:0-2129/+ 2129',
                              (('exon 18632665-18632989 rna=0-324',
                                (("5'UTR 18632665-18632989 rna=0-324",),)),
                               ('intron 18632989-18640324 rna=324-324 sjBases=GC...AG (GC_AG)',),
                               ('exon 18640324-18640587 rna=324-587',
                                (("5'UTR 18640324-18640430 rna=324-430",),
                                 ('CDS 18640430-18640587 rna=430-587 0',))),
                               ('intron 18640587-18642938 rna=587-587 sjBases=GT...AG (GT_AG)',),
                               ('exon 18642938-18643035 rna=587-684',
                                (('CDS 18642938-18643035 rna=587-684 1',),)),
                               ('intron 18643035-18644556 rna=684-684 sjBases=GT...AG (GT_AG)',),
                               ('exon 18644556-18644702 rna=684-830',
                                (('CDS 18644556-18644702 rna=684-830 2',),)),
                               ('intron 18644702-18650021 rna=830-830 sjBases=GT...AG (GT_AG)',),
                               ('exon 18650021-18650101 rna=830-910',
                                (('CDS 18650021-18650101 rna=830-910 1',),)),
                               ('intron 18650101-18650656 rna=910-910 sjBases=GT...AG (GT_AG)',),
                               ('exon 18650656-18650803 rna=910-1057',
                                (('CDS 18650656-18650803 rna=910-1057 0',),)),
                               ('intron 18650803-18652610 rna=1057-1057 sjBases=GT...AG (GT_AG)',),
                               ('exon 18652610-18652706 rna=1057-1153',
                                (('CDS 18652610-18652706 rna=1057-1153 0',),)),
                               ('intron 18652706-18653519 rna=1153-1153 sjBases=GT...AG (GT_AG)',),
                               ('exon 18653519-18653687 rna=1153-1321',
                                (('CDS 18653519-18653687 rna=1153-1321 0',),)),
                               ('intron 18653687-18655916 rna=1321-1321 sjBases=GT...AG (GT_AG)',),
                               ('exon 18655916-18656048 rna=1321-1453',
                                (('CDS 18655916-18656048 rna=1321-1453 0',),)),
                               ('intron 18656048-18656559 rna=1453-1453 sjBases=GT...AG (GT_AG)',),
                               ('exon 18656559-18656609 rna=1453-1503',
                                (('CDS 18656559-18656609 rna=1453-1503 0',),)),
                               ('intron 18656609-18659538 rna=1503-1503 sjBases=GT...AG (GT_AG)',),
                               ('exon 18659538-18660164 rna=1503-2129',
                                (('CDS 18659538-18659584 rna=1503-1549 2',),
                                 ("3'UTR 18659584-18660164 rna=1549-2129",))))))

        self.assertEqual(str(trans.toBed("0,1,2")),
                         "chr22\t18632665\t18660164\tENST00000215794.7\t0\t+\t18640430\t18659584\t0,1,2\t11\t324,263,97,146,80,147,96,168,132,50,626,\t0,7659,10273,11891,17356,17991,19945,20854,23251,23894,26873,")

    def testENST00000334029(self):
        # - strand
        trans = self.__gpToAnnotTranscript(self.__getSet1Gp("ENST00000334029.2"))
        self._assertFeatures(trans,
                             ('t=chr22:18900294-18923964/+, rna=ENST00000334029.2:0-1985/- 1985',
                              (('exon 18900294-18900875 rna=0-581',
                                (("3'UTR 18900294-18900687 rna=0-393",),
                                 ('CDS 18900687-18900875 rna=393-581 2',))),
                               ('intron 18900875-18900950 rna=581-581 sjBases=GT...AG (GT_AG)',),
                               ('exon 18900950-18901039 rna=581-670',
                                (('CDS 18900950-18901039 rna=581-670 0',),)),
                               ('intron 18901039-18904402 rna=670-670 sjBases=GT...AG (GT_AG)',),
                               ('exon 18904402-18904501 rna=670-769',
                                (('CDS 18904402-18904501 rna=670-769 2',),)),
                               ('intron 18904501-18905828 rna=769-769 sjBases=GT...AG (GT_AG)',),
                               ('exon 18905828-18906004 rna=769-945',
                                (('CDS 18905828-18906004 rna=769-945 1',),)),
                               ('intron 18906004-18906963 rna=945-945 sjBases=GC...AG (GC_AG)',),
                               ('exon 18906963-18907110 rna=945-1092',
                                (('CDS 18906963-18907110 rna=945-1092 0',),)),
                               ('intron 18907110-18907218 rna=1092-1092 sjBases=GT...AG (GT_AG)',),
                               ('exon 18907218-18907311 rna=1092-1185',
                                (('CDS 18907218-18907311 rna=1092-1185 0',),)),
                               ('intron 18907311-18908854 rna=1185-1185 sjBases=GT...AG (GT_AG)',),
                               ('exon 18908854-18908936 rna=1185-1267',
                                (('CDS 18908854-18908936 rna=1185-1267 1',),)),
                               ('intron 18908936-18909837 rna=1267-1267 sjBases=GT...AG (GT_AG)',),
                               ('exon 18909837-18909917 rna=1267-1347',
                                (('CDS 18909837-18909917 rna=1267-1347 1',),)),
                               ('intron 18909917-18910329 rna=1347-1347 sjBases=GT...AG (GT_AG)',),
                               ('exon 18910329-18910446 rna=1347-1464',
                                (('CDS 18910329-18910446 rna=1347-1464 0',),)),
                               ('intron 18910446-18910627 rna=1464-1464 sjBases=GT...AG (GT_AG)',),
                               ('exon 18910627-18910692 rna=1464-1529',
                                (('CDS 18910627-18910692 rna=1464-1529 2',),)),
                               ('intron 18910692-18912563 rna=1529-1529 sjBases=GT...AG (GT_AG)',),
                               ('exon 18912563-18912713 rna=1529-1679',
                                (('CDS 18912563-18912713 rna=1529-1679 1',),)),
                               ('intron 18912713-18913200 rna=1679-1679 sjBases=GT...AG (GT_AG)',),
                               ('exon 18913200-18913235 rna=1679-1714',
                                (('CDS 18913200-18913235 rna=1679-1714 0',),)),
                               ('intron 18913235-18918502 rna=1714-1714 sjBases=GT...AG (GT_AG)',),
                               ('exon 18918502-18918711 rna=1714-1923',
                                (('CDS 18918502-18918660 rna=1714-1872 1',),
                                 ("5'UTR 18918660-18918711 rna=1872-1923",))),
                               ('intron 18918711-18923902 rna=1923-1923 sjBases=GT...AG (GT_AG)',),
                               ('exon 18923902-18923964 rna=1923-1985',
                                (("5'UTR 18923902-18923964 rna=1923-1985",),)))))

    def testENST00000334029NoSJ(self):
        factory = AnnotationGenePredFactory(None)
        trans = factory.fromGenePred(self.__getSet1Gp("ENST00000334029.2"))
        self._assertFeatures(trans,
                             ('t=chr22:18900294-18923964/+, rna=ENST00000334029.2:0-1985/- 1985',
                              (('exon 18900294-18900875 rna=0-581',
                                (("3'UTR 18900294-18900687 rna=0-393",),
                                 ('CDS 18900687-18900875 rna=393-581 2',))),
                               ('intron 18900875-18900950 rna=581-581 sjBases=None',),
                               ('exon 18900950-18901039 rna=581-670',
                                (('CDS 18900950-18901039 rna=581-670 0',),)),
                               ('intron 18901039-18904402 rna=670-670 sjBases=None',),
                               ('exon 18904402-18904501 rna=670-769',
                                (('CDS 18904402-18904501 rna=670-769 2',),)),
                               ('intron 18904501-18905828 rna=769-769 sjBases=None',),
                               ('exon 18905828-18906004 rna=769-945',
                                (('CDS 18905828-18906004 rna=769-945 1',),)),
                               ('intron 18906004-18906963 rna=945-945 sjBases=None',),
                               ('exon 18906963-18907110 rna=945-1092',
                                (('CDS 18906963-18907110 rna=945-1092 0',),)),
                               ('intron 18907110-18907218 rna=1092-1092 sjBases=None',),
                               ('exon 18907218-18907311 rna=1092-1185',
                                (('CDS 18907218-18907311 rna=1092-1185 0',),)),
                               ('intron 18907311-18908854 rna=1185-1185 sjBases=None',),
                               ('exon 18908854-18908936 rna=1185-1267',
                                (('CDS 18908854-18908936 rna=1185-1267 1',),)),
                               ('intron 18908936-18909837 rna=1267-1267 sjBases=None',),
                               ('exon 18909837-18909917 rna=1267-1347',
                                (('CDS 18909837-18909917 rna=1267-1347 1',),)),
                               ('intron 18909917-18910329 rna=1347-1347 sjBases=None',),
                               ('exon 18910329-18910446 rna=1347-1464',
                                (('CDS 18910329-18910446 rna=1347-1464 0',),)),
                               ('intron 18910446-18910627 rna=1464-1464 sjBases=None',),
                               ('exon 18910627-18910692 rna=1464-1529',
                                (('CDS 18910627-18910692 rna=1464-1529 2',),)),
                               ('intron 18910692-18912563 rna=1529-1529 sjBases=None',),
                               ('exon 18912563-18912713 rna=1529-1679',
                                (('CDS 18912563-18912713 rna=1529-1679 1',),)),
                               ('intron 18912713-18913200 rna=1679-1679 sjBases=None',),
                               ('exon 18913200-18913235 rna=1679-1714',
                                (('CDS 18913200-18913235 rna=1679-1714 0',),)),
                               ('intron 18913235-18918502 rna=1714-1714 sjBases=None',),
                               ('exon 18918502-18918711 rna=1714-1923',
                                (('CDS 18918502-18918660 rna=1714-1872 1',),
                                 ("5'UTR 18918660-18918711 rna=1872-1923",))),
                               ('intron 18918711-18923902 rna=1923-1923 sjBases=None',),
                               ('exon 18923902-18923964 rna=1923-1985',
                                (("5'UTR 18923902-18923964 rna=1923-1985",),)))))

    def testENST00000334029Rc(self):
        # - strand
        trans = self.__gpToAnnotTranscript(self.__getSet1Gp("ENST00000334029.2"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans,
                             ('t=chr22:32380602-32404272/-, rna=ENST00000334029.2:0-1985/+ 1985',
                              (('exon 32380602-32380664 rna=0-62',
                                (("5'UTR 32380602-32380664 rna=0-62",),)),
                               ('intron 32380664-32385855 rna=62-62 sjBases=GT...AG (GT_AG)',),
                               ('exon 32385855-32386064 rna=62-271',
                                (("5'UTR 32385855-32385906 rna=62-113",),
                                 ('CDS 32385906-32386064 rna=113-271 0',))),
                               ('intron 32386064-32391331 rna=271-271 sjBases=GT...AG (GT_AG)',),
                               ('exon 32391331-32391366 rna=271-306',
                                (('CDS 32391331-32391366 rna=271-306 2',),)),
                               ('intron 32391366-32391853 rna=306-306 sjBases=GT...AG (GT_AG)',),
                               ('exon 32391853-32392003 rna=306-456',
                                (('CDS 32391853-32392003 rna=306-456 1',),)),
                               ('intron 32392003-32393874 rna=456-456 sjBases=GT...AG (GT_AG)',),
                               ('exon 32393874-32393939 rna=456-521',
                                (('CDS 32393874-32393939 rna=456-521 1',),)),
                               ('intron 32393939-32394120 rna=521-521 sjBases=GT...AG (GT_AG)',),
                               ('exon 32394120-32394237 rna=521-638',
                                (('CDS 32394120-32394237 rna=521-638 0',),)),
                               ('intron 32394237-32394649 rna=638-638 sjBases=GT...AG (GT_AG)',),
                               ('exon 32394649-32394729 rna=638-718',
                                (('CDS 32394649-32394729 rna=638-718 0',),)),
                               ('intron 32394729-32395630 rna=718-718 sjBases=GT...AG (GT_AG)',),
                               ('exon 32395630-32395712 rna=718-800',
                                (('CDS 32395630-32395712 rna=718-800 2',),)),
                               ('intron 32395712-32397255 rna=800-800 sjBases=GT...AG (GT_AG)',),
                               ('exon 32397255-32397348 rna=800-893',
                                (('CDS 32397255-32397348 rna=800-893 0',),)),
                               ('intron 32397348-32397456 rna=893-893 sjBases=GT...AG (GT_AG)',),
                               ('exon 32397456-32397603 rna=893-1040',
                                (('CDS 32397456-32397603 rna=893-1040 0',),)),
                               ('intron 32397603-32398562 rna=1040-1040 sjBases=GC...AG (GC_AG)',),
                               ('exon 32398562-32398738 rna=1040-1216',
                                (('CDS 32398562-32398738 rna=1040-1216 0',),)),
                               ('intron 32398738-32400065 rna=1216-1216 sjBases=GT...AG (GT_AG)',),
                               ('exon 32400065-32400164 rna=1216-1315',
                                (('CDS 32400065-32400164 rna=1216-1315 2',),)),
                               ('intron 32400164-32403527 rna=1315-1315 sjBases=GT...AG (GT_AG)',),
                               ('exon 32403527-32403616 rna=1315-1404',
                                (('CDS 32403527-32403616 rna=1315-1404 2',),)),
                               ('intron 32403616-32403691 rna=1404-1404 sjBases=GT...AG (GT_AG)',),
                               ('exon 32403691-32404272 rna=1404-1985',
                                (('CDS 32403691-32403879 rna=1404-1592 1',),
                                 ("3'UTR 32403879-32404272 rna=1592-1985",))))))

    def testENST00000334029NoSJRc(self):
        # no splice sites, just sizes
        factory = AnnotationGenePredFactory(chromSizeFunc=GenomeSeqSrc.obtain("hg19").getChromSize)
        trans = factory.fromGenePred(self.__getSet1Gp("ENST00000334029.2"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans,
                             ('t=chr22:32380602-32404272/-, rna=ENST00000334029.2:0-1985/+ 1985',
                              (('exon 32380602-32380664 rna=0-62',
                                (("5'UTR 32380602-32380664 rna=0-62",),)),
                               ('intron 32380664-32385855 rna=62-62 sjBases=None',),
                               ('exon 32385855-32386064 rna=62-271',
                                (("5'UTR 32385855-32385906 rna=62-113",),
                                 ('CDS 32385906-32386064 rna=113-271 0',))),
                               ('intron 32386064-32391331 rna=271-271 sjBases=None',),
                               ('exon 32391331-32391366 rna=271-306',
                                (('CDS 32391331-32391366 rna=271-306 2',),)),
                               ('intron 32391366-32391853 rna=306-306 sjBases=None',),
                               ('exon 32391853-32392003 rna=306-456',
                                (('CDS 32391853-32392003 rna=306-456 1',),)),
                               ('intron 32392003-32393874 rna=456-456 sjBases=None',),
                               ('exon 32393874-32393939 rna=456-521',
                                (('CDS 32393874-32393939 rna=456-521 1',),)),
                               ('intron 32393939-32394120 rna=521-521 sjBases=None',),
                               ('exon 32394120-32394237 rna=521-638',
                                (('CDS 32394120-32394237 rna=521-638 0',),)),
                               ('intron 32394237-32394649 rna=638-638 sjBases=None',),
                               ('exon 32394649-32394729 rna=638-718',
                                (('CDS 32394649-32394729 rna=638-718 0',),)),
                               ('intron 32394729-32395630 rna=718-718 sjBases=None',),
                               ('exon 32395630-32395712 rna=718-800',
                                (('CDS 32395630-32395712 rna=718-800 2',),)),
                               ('intron 32395712-32397255 rna=800-800 sjBases=None',),
                               ('exon 32397255-32397348 rna=800-893',
                                (('CDS 32397255-32397348 rna=800-893 0',),)),
                               ('intron 32397348-32397456 rna=893-893 sjBases=None',),
                               ('exon 32397456-32397603 rna=893-1040',
                                (('CDS 32397456-32397603 rna=893-1040 0',),)),
                               ('intron 32397603-32398562 rna=1040-1040 sjBases=None',),
                               ('exon 32398562-32398738 rna=1040-1216',
                                (('CDS 32398562-32398738 rna=1040-1216 0',),)),
                               ('intron 32398738-32400065 rna=1216-1216 sjBases=None',),
                               ('exon 32400065-32400164 rna=1216-1315',
                                (('CDS 32400065-32400164 rna=1216-1315 2',),)),
                               ('intron 32400164-32403527 rna=1315-1315 sjBases=None',),
                               ('exon 32403527-32403616 rna=1315-1404',
                                (('CDS 32403527-32403616 rna=1315-1404 2',),)),
                               ('intron 32403616-32403691 rna=1404-1404 sjBases=None',),
                               ('exon 32403691-32404272 rna=1404-1985',
                                (('CDS 32403691-32403879 rna=1404-1592 1',),
                                 ("3'UTR 32403879-32404272 rna=1592-1985",))))))

    def testGetStructureFeaturesOfType(self):
        trans = self.__gpToAnnotTranscript(self.__getSet1Gp("ENST00000334029.2"))
        feats = trans.getStructureFeaturesOfType(ExonFeature)
        featStrs = tuple([str(f) for f in feats])
        self.assertEquals(featStrs,
                          ('exon 18900294-18900875 rna=0-581',
                           'exon 18900950-18901039 rna=581-670',
                           'exon 18904402-18904501 rna=670-769',
                           'exon 18905828-18906004 rna=769-945',
                           'exon 18906963-18907110 rna=945-1092',
                           'exon 18907218-18907311 rna=1092-1185',
                           'exon 18908854-18908936 rna=1185-1267',
                           'exon 18909837-18909917 rna=1267-1347',
                           'exon 18910329-18910446 rna=1347-1464',
                           'exon 18910627-18910692 rna=1464-1529',
                           'exon 18912563-18912713 rna=1529-1679',
                           'exon 18913200-18913235 rna=1679-1714',
                           'exon 18918502-18918711 rna=1714-1923',
                           'exon 18923902-18923964 rna=1923-1985'))

    def testGetAnnotationFeaturesOfType(self):
        trans = self.__gpToAnnotTranscript(self.__getSet1Gp("ENST00000334029.2"))
        feats = trans.getAnnotationFeaturesOfType(CdsRegionFeature)
        featStrs = tuple([str(f) for f in feats])
        self.assertEquals(featStrs,
                          ('CDS 18900687-18900875 rna=393-581 2',
                           'CDS 18900950-18901039 rna=581-670 0',
                           'CDS 18904402-18904501 rna=670-769 2',
                           'CDS 18905828-18906004 rna=769-945 1',
                           'CDS 18906963-18907110 rna=945-1092 0',
                           'CDS 18907218-18907311 rna=1092-1185 0',
                           'CDS 18908854-18908936 rna=1185-1267 1',
                           'CDS 18909837-18909917 rna=1267-1347 1',
                           'CDS 18910329-18910446 rna=1347-1464 0',
                           'CDS 18910627-18910692 rna=1464-1529 2',
                           'CDS 18912563-18912713 rna=1529-1679 1',
                           'CDS 18913200-18913235 rna=1679-1714 0',
                           'CDS 18918502-18918660 rna=1714-1872 1'))

    def testENST00000434390(self):
        # non-coding, - strand
        trans = self.__gpToAnnotTranscript(self.__getSet1Gp("ENST00000434390.1"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans,
                             (('t=chr22:32615884-32643762/-, rna=ENST00000434390.1:0-1859/+ 1859',
                               (('exon 32615884-32616025 rna=0-141',
                                 (('NC 32615884-32616025 rna=0-141',),)),
                                ('intron 32616025-32622058 rna=141-141 sjBases=GT...AG (GT_AG)',),
                                ('exon 32622058-32622117 rna=141-200',
                                 (('NC 32622058-32622117 rna=141-200',),)),
                                ('intron 32622117-32622585 rna=200-200 sjBases=GT...AG (GT_AG)',),
                                ('exon 32622585-32622610 rna=200-225',
                                 (('NC 32622585-32622610 rna=200-225',),)),
                                ('intron 32622610-32625556 rna=225-225 sjBases=GT...AG (GT_AG)',),
                                ('exon 32625556-32625643 rna=225-312',
                                 (('NC 32625556-32625643 rna=225-312',),)),
                                ('intron 32625643-32630551 rna=312-312 sjBases=GT...AG (GT_AG)',),
                                ('exon 32630551-32630617 rna=312-378',
                                 (('NC 32630551-32630617 rna=312-378',),)),
                                ('intron 32630617-32632616 rna=378-378 sjBases=GT...AG (GT_AG)',),
                                ('exon 32632616-32632678 rna=378-440',
                                 (('NC 32632616-32632678 rna=378-440',),)),
                                ('intron 32632678-32633143 rna=440-440 sjBases=GT...AG (GT_AG)',),
                                ('exon 32633143-32633177 rna=440-474',
                                 (('NC 32633143-32633177 rna=440-474',),)),
                                ('intron 32633177-32637733 rna=474-474 sjBases=GT...AG (GT_AG)',),
                                ('exon 32637733-32638690 rna=474-1431',
                                 (('NC 32637733-32638690 rna=474-1431',),)),
                                ('intron 32638690-32643334 rna=1431-1431 sjBases=GT...AG (GT_AG)',),
                                ('exon 32643334-32643762 rna=1431-1859',
                                 (('NC 32643334-32643762 rna=1431-1859',),))))))

    def testGencodeV26Regress(self):
        "regression test for gencodeV26"
        factory = AnnotationGenePredFactory()
        names = ("ENST00000610542.1",  # 4-base gap that caused error
                 )
        i = 0
        for gp in GenePredReader(self.getInputFile("gencodeV26.gp")):
            trans = factory.fromGenePred(gp)
            self.assertEqual(trans.rna.name, names[i])

    def testAnnotToBed(self):
        """test for conversion to BED"""
        # had a bug non-code genes had only basic columns because thickStart/thickEnd were None
        factory = AnnotationGenePredFactory()
        for gp in GenePredReader(self.getInputFile("set1.gencodeCompV19.gp")):
            trans = factory.fromGenePred(gp)
            bed = trans.toBed("100,0,0")
            self.assertEqual(len(bed.getRow()), 12)

    def testAttrs(self):
        def getUtr3Feature(trans):
            # test case has one 3'UTR feature
            return trans.getAnnotationFeaturesOfType(Utr3RegionFeature)[0]

        factory = AnnotationGenePredFactory(GenomeSeqSrc.obtain("hg19"))
        tattrs = ObjDict(name="Fred")
        uattrs = ObjDict(name="Barney")
        trans = factory.fromGenePred(self.__getSet1Gp("ENST00000334029.2"), tattrs)
        # has one 3'UTR feature, so use it for the attrs
        utr3 = getUtr3Feature(trans)
        utr3.attrs = uattrs

        rcTrans = trans.reverseComplement()
        rcUtr3 = getUtr3Feature(rcTrans)
        self.assertEqual("Fred", rcTrans.attrs.name)
        self.assertEqual("Barney", rcUtr3.attrs.name)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts


if __name__ == '__main__':
    unittest.main()
