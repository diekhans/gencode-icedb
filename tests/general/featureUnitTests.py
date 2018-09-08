import sys
import os
if __name__ == '__main__':
    rootDir = "../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.objDict import ObjDict
from pycbio.sys.testCaseBase import TestCaseBase
from gencode_icedb.general.genome import GenomeReader
from gencode_icedb.general.transFeatures import ExonFeature
from gencode_icedb.general.transFeatures import AnnotationFeature, CdsRegionFeature, Utr3RegionFeature
from gencode_icedb.general.transFeatures import RnaInsertFeature, ChromInsertFeature
from gencode_icedb.general.evidFeatures import EvidencePslFactory
from gencode_icedb.general.genePredAnnotFeatures import GenePredAnnotationFactory
from gencode_icedb.general.ensemblDbAnnotFeatures import EnsemblDbAnnotationFactory
from gencode_icedb.general.ensemblDb import ensemblGeneQuery
from pycbio.db import sqliteOps
from pycbio.db import mysqlOps
from pycbio.hgdata.genePredSqlite import GenePredSqliteTable
from pycbio.hgdata.pslSqlite import PslSqliteTable
from pycbio.hgdata.psl import Psl
from pycbio.hgdata.genePred import GenePredReader
from pycbio.sys.pprint2 import nswpprint

# FIXME: be more consisted on how test cases are obtained, have DbSrc caches get
# genome too.

debugResults = False   # print out results for updated expected
noCheckResults = False  # don't check results
#debugResults=True

if debugResults or noCheckResults:
    print("WARNING: debug variables set", file=sys.stderr)


def getInputFile(base):
    "from input relative to test file"
    return os.path.join(os.path.dirname(__file__), "input", base)


class GenomeSeqSrc(object):
    "caching genome reader"

    srcs = {
        "hg38": "/hive/data/genomes/hg38/hg38.2bit",
        "mm10": "/hive/data/genomes/mm10/mm10.2bit",
        "grch38": "../data/grch38/GRCh38.fa.gz"
    }
    readers = {}

    @classmethod
    def obtain(cls, db):
        reader = cls.readers.get(db)
        if reader is None:
            reader = cls.readers[db] = GenomeReader.getFromFileName(cls.srcs[db])
        return reader


class PslDbSrc(object):
    "caching psl sqlite database"
    srcs = {
        "V28": "ucsc-mrnaV28.psl",
        "hg38-mm10.transMap": "hg38-mm10.transMap.psl"
    }
    pslTbls = {}

    @classmethod
    def _loadDbTbl(cls, name):
        conn = sqliteOps.connect(None)
        dbTbl = PslSqliteTable(conn, "psls", create=True)
        dbTbl.loadPslFile(getInputFile(cls.srcs[name]))
        cls.pslTbls[name] = dbTbl
        return dbTbl

    @classmethod
    def obtain(cls, name):
        dbTbl = cls.pslTbls.get(name)
        if dbTbl is None:
            dbTbl = cls._loadDbTbl(name)
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
        "V28": "gencodeCompV28.gp",
    }
    genePredTbls = {}

    @classmethod
    def _loadDbTbl(cls, name):
        conn = sqliteOps.connect(None)
        dbTbl = GenePredSqliteTable(conn, name, create=True)
        dbTbl.loadGenePredFile(getInputFile(cls.srcs[name]))
        cls.genePredTbls[name] = dbTbl
        return dbTbl

    @classmethod
    def obtain(cls, name):
        dbTbl = cls.genePredTbls.get(name)
        if dbTbl is None:
            dbTbl = cls._loadDbTbl(name)
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
    def _getV28Psl(self, acc):
        return PslDbSrc.obtainPsl("V28", acc)

    def _pslToEvidTranscript(self, psl):
        return EvidencePslFactory(GenomeSeqSrc.obtain("hg38")).fromPsl(psl)

    def _checkRnaAln(self, trans):
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
        trans = self._pslToEvidTranscript(self._getV28Psl("AF010310.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18912781-18918413/+, rna=AF010310.1:13-901/- 901 <->',
                              (('exon 18912781-18913362 rna=13-618',
                                (('aln 18912781-18912939 rna=13-171',),
                                 ('rins None-None rna=171-174',),
                                 ('cins 18912939-18912940 rna=None-None',),
                                 ('aln 18912940-18913022 rna=174-256',),
                                 ('rins None-None rna=256-257',),
                                 ('aln 18913022-18913052 rna=257-287',),
                                 ('rins None-None rna=287-288',),
                                 ('aln 18913052-18913054 rna=288-290',),
                                 ('rins None-None rna=290-291',),
                                 ('aln 18913054-18913064 rna=291-301',),
                                 ('rins None-None rna=301-302',),
                                 ('aln 18913064-18913066 rna=302-304',),
                                 ('rins None-None rna=304-305',),
                                 ('aln 18913066-18913093 rna=305-332',),
                                 ('rins None-None rna=332-335',),
                                 ('cins 18913093-18913094 rna=None-None',),
                                 ('aln 18913094-18913098 rna=335-339',),
                                 ('rins None-None rna=339-340',),
                                 ('aln 18913098-18913103 rna=340-345',),
                                 ('rins None-None rna=345-346',),
                                 ('aln 18913103-18913108 rna=346-351',),
                                 ('rins None-None rna=351-352',),
                                 ('aln 18913108-18913115 rna=352-359',),
                                 ('rins None-None rna=359-369',),
                                 ('cins 18913115-18913121 rna=None-None',),
                                 ('aln 18913121-18913130 rna=369-378',),
                                 ('rins None-None rna=378-379',),
                                 ('aln 18913130-18913132 rna=379-381',),
                                 ('rins None-None rna=381-382',),
                                 ('aln 18913132-18913151 rna=382-401',),
                                 ('rins None-None rna=401-402',),
                                 ('aln 18913151-18913156 rna=402-407',),
                                 ('rins None-None rna=407-408',),
                                 ('aln 18913156-18913249 rna=408-501',),
                                 ('rins None-None rna=501-502',),
                                 ('aln 18913249-18913277 rna=502-530',),
                                 ('rins None-None rna=530-531',),
                                 ('aln 18913277-18913290 rna=531-544',),
                                 ('rins None-None rna=544-545',),
                                 ('aln 18913290-18913358 rna=545-613',),
                                 ('rins None-None rna=613-614',),
                                 ('aln 18913358-18913362 rna=614-618',))),
                               ('intron 18913362-18913437 rna=618-618 sjBases=GT...AG (GT_AG)',
                                (('cins 18913362-18913437 rna=None-None',),)),
                               ('exon 18913437-18913526 rna=618-707',
                                (('aln 18913437-18913526 rna=618-707',),)),
                               ('intron 18913526-18916889 rna=707-707 sjBases=GT...AG (GT_AG)',
                                (('cins 18913526-18916889 rna=None-None',),)),
                               ('exon 18916889-18916988 rna=707-806',
                                (('aln 18916889-18916988 rna=707-806',),)),
                               ('intron 18916988-18918315 rna=806-806 sjBases=GT...AG (GT_AG)',
                                (('cins 18916988-18918315 rna=None-None',),)),
                               ('exon 18918315-18918413 rna=806-901',
                                (('aln 18918315-18918376 rna=806-867',),
                                 ('rins None-None rna=867-868',),
                                 ('cins 18918376-18918380 rna=None-None',),
                                 ('aln 18918380-18918413 rna=868-901',))))))
        self._checkRnaAln(trans)

    def testX96484(self):
        trans = self._pslToEvidTranscript(self._getV28Psl("X96484.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18906409-18912079/+, rna=X96484.1:48-1067/+ 1080 <+>',
                              (('exon 18906409-18906484 rna=48-123',
                                (('aln 18906409-18906484 rna=48-123',),)),
                               ('intron 18906484-18906564 rna=123-123 sjBases=GT...AG (GT_AG)',
                                (('cins 18906484-18906564 rna=None-None',),)),
                               ('exon 18906564-18906725 rna=123-285',
                                (('aln 18906564-18906660 rna=123-219',),
                                 ('rins None-None rna=219-220',),
                                 ('aln 18906660-18906725 rna=220-285',))),
                               ('intron 18906725-18910171 rna=285-285 sjBases=GT...AG (GT_AG)',
                                (('cins 18906725-18910171 rna=None-None',),)),
                               ('exon 18910171-18910272 rna=285-386',
                                (('aln 18910171-18910272 rna=285-386',),)),
                               ('intron 18910272-18910887 rna=386-386 sjBases=GT...AG (GT_AG)',
                                (('cins 18910272-18910887 rna=None-None',),)),
                               ('exon 18910887-18911028 rna=386-527',
                                (('aln 18910887-18911028 rna=386-527',),)),
                               ('intron 18911028-18911539 rna=527-527 sjBases=GT...AG (GT_AG)',
                                (('cins 18911028-18911539 rna=None-None',),)),
                               ('exon 18911539-18912079 rna=527-1067',
                                (('aln 18911539-18912079 rna=527-1067',),)))))
        self._checkRnaAln(trans)

    def testX96484NoSJ(self):
        trans = EvidencePslFactory(None).fromPsl(self._getV28Psl("X96484.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18906409-18912079/+, rna=X96484.1:48-1067/+ 1080 <+>',
                              (('exon 18906409-18906484 rna=48-123',
                                (('aln 18906409-18906484 rna=48-123',),)),
                               ('intron 18906484-18906564 rna=123-123',
                                (('cins 18906484-18906564 rna=None-None',),)),
                               ('exon 18906564-18906725 rna=123-285',
                                (('aln 18906564-18906660 rna=123-219',),
                                 ('rins None-None rna=219-220',),
                                 ('aln 18906660-18906725 rna=220-285',))),
                               ('intron 18906725-18910171 rna=285-285',
                                (('cins 18906725-18910171 rna=None-None',),)),
                               ('exon 18910171-18910272 rna=285-386',
                                (('aln 18910171-18910272 rna=285-386',),)),
                               ('intron 18910272-18910887 rna=386-386',
                                (('cins 18910272-18910887 rna=None-None',),)),
                               ('exon 18910887-18911028 rna=386-527',
                                (('aln 18910887-18911028 rna=386-527',),)),
                               ('intron 18911028-18911539 rna=527-527',
                                (('cins 18911028-18911539 rna=None-None',),)),
                               ('exon 18911539-18912079 rna=527-1067',
                                (('aln 18911539-18912079 rna=527-1067',),)))))
        self._checkRnaAln(trans)

    def testTransMapDropExon(self):
        # internal exon was not mapped, causing an intron to contain unaligned
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        self._assertFeatures(trans,
                             ('t=chr4:148039043-148056154/+, rna=ENST00000641446:0-2820/+ 2820 <+>',
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
        self._checkRnaAln(trans)

    def testTransMapDropExonRc(self):
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        transRc = trans.reverseComplement()
        self._assertFeatures(transRc,
                             ('t=chr4:8451962-8469073/-, rna=ENST00000641446:0-2820/- 2820 <+>',
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
        self._checkRnaAln(transRc)

    def testGetAlignmentFeaturesOfType(self):
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        feats = trans.getFeaturesOfType((RnaInsertFeature, ChromInsertFeature))
        featStrs = tuple([str(f) for f in feats])
        self.assertEqual(featStrs,
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
        aln1 = self._pslToEvidTranscript(self._getV28Psl("AF010310.1"))
        aln2 = self._pslToEvidTranscript(self._getV28Psl("AF120278.1"))
        exon1 = aln1.features[0]
        exon2 = aln2.features[0]
        self.assertTrue(exon1.rnaOverlaps(exon2))
        exon2b = aln2.features[3]
        self.assertFalse(exon1.rnaOverlaps(exon2b))
        self._checkRnaAln(aln1)
        self._checkRnaAln(aln2)

    def testStructPrevNext(self):
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)
        feat = trans.firstFeature(ExonFeature)
        self.assertIs(feat, trans.features[0])
        self.assertIsInstance(feat, ExonFeature)
        exonCnt = 1
        lastFeat = None
        while True:
            feat = feat.nextFeature(ExonFeature)
            if feat is None:
                break
            exonCnt += 1
            lastFeat = feat
        self.assertEqual(exonCnt, 14)

        feat = trans.lastFeature(ExonFeature)
        self.assertIs(feat, lastFeat)
        self.assertIsInstance(feat, ExonFeature)
        exonCnt = 1
        while True:
            feat = feat.prevFeature(ExonFeature)
            if feat is None:
                break
            exonCnt += 1
        self.assertEqual(exonCnt, 14)

    def testAlignPrevNext(self):
        psl = PslDbSrc.obtainPsl("hg38-mm10.transMap", "ENST00000641446")
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("mm10")).fromPsl(psl)

        feat = trans.firstFeature(RnaInsertFeature)
        self.assertIsInstance(feat, RnaInsertFeature)
        rinsCnt = 1
        while True:
            feat = feat.nextFeature(RnaInsertFeature)
            if feat is None:
                break
            rinsCnt += 1
        self.assertEqual(rinsCnt, 18)

        feat = trans.lastFeature(RnaInsertFeature)
        rinsCnt = 1
        while True:
            feat = feat.prevFeature(RnaInsertFeature)
            if feat is None:
                break
            rinsCnt += 1
        self.assertEqual(rinsCnt, 18)

    def testEst3PslMinus(self):
        # test handing 3' EST for negative strand gene
        estPsl = ["620", "11", "0", "8", "0", "0", "3", "6544", "--", "BX371226.2", "645", "6", "645", "chr22", "50818468", "41938778", "41945961", "4", "49,80,92,418,", "0,49,129,221,", "8872507,8873164,8874767,8879272,"]
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("hg38")).fromPsl(Psl(estPsl), orientChrom=False)
        self._assertFeatures(trans,
                             ('t=chr22:8872507-8879690/-, rna=BX371226.2:0-639/- 645 <->',
                              (('exon 8872507-8872556 rna=0-49', (('aln 8872507-8872556 rna=0-49',),)),
                               ('intron 8872556-8873164 rna=49-49 sjBases=GT...AG (GT_AG)',
                                (('cins 8872556-8873164 rna=None-None',),)),
                               ('exon 8873164-8873244 rna=49-129', (('aln 8873164-8873244 rna=49-129',),)),
                               ('intron 8873244-8874767 rna=129-129 sjBases=GT...AG (GT_AG)',
                                (('cins 8873244-8874767 rna=None-None',),)),
                               ('exon 8874767-8874859 rna=129-221',
                                (('aln 8874767-8874859 rna=129-221',),)),
                               ('intron 8874859-8879272 rna=221-221 sjBases=GT...AG (GT_AG)',
                                (('cins 8874859-8879272 rna=None-None',),)),
                               ('exon 8879272-8879690 rna=221-639',
                                (('aln 8879272-8879690 rna=221-639',),)))))
        self.assertEqual('-', trans.transcriptionStrand)
        # check the same PSL can be oriented chrom +
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("hg38")).fromPsl(Psl(estPsl), orientChrom=True)
        self._assertFeatures(trans,
                             ('t=chr22:41938778-41945961/+, rna=BX371226.2:6-645/+ 645 <->',
                              (('exon 41938778-41939196 rna=6-424',
                                (('aln 41938778-41939196 rna=6-424',),)),
                               ('intron 41939196-41943609 rna=424-424 sjBases=GT...AG (GT_AG)',
                                (('cins 41939196-41943609 rna=None-None',),)),
                               ('exon 41943609-41943701 rna=424-516',
                                (('aln 41943609-41943701 rna=424-516',),)),
                               ('intron 41943701-41945224 rna=516-516 sjBases=GT...AG (GT_AG)',
                                (('cins 41943701-41945224 rna=None-None',),)),
                               ('exon 41945224-41945304 rna=516-596',
                                (('aln 41945224-41945304 rna=516-596',),)),
                               ('intron 41945304-41945912 rna=596-596 sjBases=GT...AG (GT_AG)',
                                (('cins 41945304-41945912 rna=None-None',),)),
                               ('exon 41945912-41945961 rna=596-645',
                                (('aln 41945912-41945961 rna=596-645',),)))))
        self.assertEqual('-', trans.transcriptionStrand)

    def testEst3PslPlus(self):
        # test handing 3' EST for positive strand gene
        estPsl = ["646", "3", "0", "1", "0", "0", "3", "9480", "+-", "BM969800.1", "675", "15", "665", "chr22", "50818468", "45422286", "45432416", "4", "133,167,228,122,", "15,148,315,543,", "5386052,5387402,5392293,5396060,"]
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("hg38")).fromPsl(Psl(estPsl), orientChrom=False)
        self._assertFeatures(trans,
                             ('t=chr22:5386052-5396182/-, rna=BM969800.1:15-665/+ 675 <+>',
                              (('exon 5386052-5386185 rna=15-148', (('aln 5386052-5386185 rna=15-148',),)),
                               ('intron 5386185-5387402 rna=148-148 sjBases=GT...AG (GT_AG)',
                                (('cins 5386185-5387402 rna=None-None',),)),
                               ('exon 5387402-5387569 rna=148-315',
                                (('aln 5387402-5387569 rna=148-315',),)),
                               ('intron 5387569-5392293 rna=315-315 sjBases=GT...AG (GT_AG)',
                                (('cins 5387569-5392293 rna=None-None',),)),
                               ('exon 5392293-5392521 rna=315-543',
                                (('aln 5392293-5392521 rna=315-543',),)),
                               ('intron 5392521-5396060 rna=543-543 sjBases=GT...AG (GT_AG)',
                                (('cins 5392521-5396060 rna=None-None',),)),
                               ('exon 5396060-5396182 rna=543-665',
                                (('aln 5396060-5396182 rna=543-665',),)))))
        self.assertEqual('+', trans.transcriptionStrand)
        # check the same PSL can be oriented chrom +
        trans = EvidencePslFactory(GenomeSeqSrc.obtain("hg38")).fromPsl(Psl(estPsl), orientChrom=True)
        self._assertFeatures(trans,
                             ('t=chr22:45422286-45432416/+, rna=BM969800.1:10-660/- 675 <+>',
                              (('exon 45422286-45422408 rna=10-132',
                                (('aln 45422286-45422408 rna=10-132',),)),
                               ('intron 45422408-45425947 rna=132-132 sjBases=GT...AG (GT_AG)',
                                (('cins 45422408-45425947 rna=None-None',),)),
                               ('exon 45425947-45426175 rna=132-360',
                                (('aln 45425947-45426175 rna=132-360',),)),
                               ('intron 45426175-45430899 rna=360-360 sjBases=GT...AG (GT_AG)',
                                (('cins 45426175-45430899 rna=None-None',),)),
                               ('exon 45430899-45431066 rna=360-527',
                                (('aln 45430899-45431066 rna=360-527',),)),
                               ('intron 45431066-45432283 rna=527-527 sjBases=GT...AG (GT_AG)',
                                (('cins 45431066-45432283 rna=None-None',),)),
                               ('exon 45432283-45432416 rna=527-660',
                                (('aln 45432283-45432416 rna=527-660',),)))))
        self.assertEqual('+', trans.transcriptionStrand)


class AnnotationCheckMixin(object):
    """Mixin of functions to validate results that a common between GenePred
    and EnsemblDb annotation tests"""

    def _toEnsemblExpect(self, expect):
        return (expect[0].replace('chr', ''), expect[1])

    def checkENST00000215794(self, trans, ensChroms=False):
        # coding, + strand
        expect = ('t=chr22:18149898-18177397/+, rna=ENST00000215794.7:0-2129/+ 2129 <+>',
                  (('exon 18149898-18150222 rna=0-324',
                    (("5'UTR 18149898-18150222 rna=0-324",),)),
                   ('intron 18150222-18157557 rna=324-324 sjBases=GC...AG (GC_AG)',),
                   ('exon 18157557-18157820 rna=324-587',
                    (("5'UTR 18157557-18157663 rna=324-430",),
                     ('CDS 18157663-18157820 rna=430-587 0',))),
                   ('intron 18157820-18160171 rna=587-587 sjBases=GT...AG (GT_AG)',),
                   ('exon 18160171-18160268 rna=587-684',
                    (('CDS 18160171-18160268 rna=587-684 1',),)),
                   ('intron 18160268-18161789 rna=684-684 sjBases=GT...AG (GT_AG)',),
                   ('exon 18161789-18161935 rna=684-830',
                    (('CDS 18161789-18161935 rna=684-830 2',),)),
                   ('intron 18161935-18167254 rna=830-830 sjBases=GT...AG (GT_AG)',),
                   ('exon 18167254-18167334 rna=830-910',
                    (('CDS 18167254-18167334 rna=830-910 1',),)),
                   ('intron 18167334-18167889 rna=910-910 sjBases=GT...AG (GT_AG)',),
                   ('exon 18167889-18168036 rna=910-1057',
                    (('CDS 18167889-18168036 rna=910-1057 0',),)),
                   ('intron 18168036-18169843 rna=1057-1057 sjBases=GT...AG (GT_AG)',),
                   ('exon 18169843-18169939 rna=1057-1153',
                    (('CDS 18169843-18169939 rna=1057-1153 0',),)),
                   ('intron 18169939-18170752 rna=1153-1153 sjBases=GT...AG (GT_AG)',),
                   ('exon 18170752-18170920 rna=1153-1321',
                    (('CDS 18170752-18170920 rna=1153-1321 0',),)),
                   ('intron 18170920-18173149 rna=1321-1321 sjBases=GT...AG (GT_AG)',),
                   ('exon 18173149-18173281 rna=1321-1453',
                    (('CDS 18173149-18173281 rna=1321-1453 0',),)),
                   ('intron 18173281-18173792 rna=1453-1453 sjBases=GT...AG (GT_AG)',),
                   ('exon 18173792-18173842 rna=1453-1503',
                    (('CDS 18173792-18173842 rna=1453-1503 0',),)),
                   ('intron 18173842-18176771 rna=1503-1503 sjBases=GT...AG (GT_AG)',),
                   ('exon 18176771-18177397 rna=1503-2129',
                    (('CDS 18176771-18176817 rna=1503-1549 2',),
                     ("3'UTR 18176817-18177397 rna=1549-2129",)))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000334029(self, trans, ensChroms=False):
        # coding, - strand
        expect = ('t=chr22:18912781-18936451/+, rna=ENST00000334029.6:0-1985/- 1985 <->',
                  (('exon 18912781-18913362 rna=0-581',
                    (("3'UTR 18912781-18913174 rna=0-393",),
                     ('CDS 18913174-18913362 rna=393-581 1',))),
                   ('intron 18913362-18913437 rna=581-581 sjBases=GT...AG (GT_AG)',),
                   ('exon 18913437-18913526 rna=581-670',
                    (('CDS 18913437-18913526 rna=581-670 2',),)),
                   ('intron 18913526-18916889 rna=670-670 sjBases=GT...AG (GT_AG)',),
                   ('exon 18916889-18916988 rna=670-769',
                    (('CDS 18916889-18916988 rna=670-769 2',),)),
                   ('intron 18916988-18918315 rna=769-769 sjBases=GT...AG (GT_AG)',),
                   ('exon 18918315-18918491 rna=769-945',
                    (('CDS 18918315-18918491 rna=769-945 0',),)),
                   ('intron 18918491-18919450 rna=945-945 sjBases=GC...AG (GC_AG)',),
                   ('exon 18919450-18919597 rna=945-1092',
                    (('CDS 18919450-18919597 rna=945-1092 0',),)),
                   ('intron 18919597-18919705 rna=1092-1092 sjBases=GT...AG (GT_AG)',),
                   ('exon 18919705-18919798 rna=1092-1185',
                    (('CDS 18919705-18919798 rna=1092-1185 0',),)),
                   ('intron 18919798-18921341 rna=1185-1185 sjBases=GT...AG (GT_AG)',),
                   ('exon 18921341-18921423 rna=1185-1267',
                    (('CDS 18921341-18921423 rna=1185-1267 2',),)),
                   ('intron 18921423-18922324 rna=1267-1267 sjBases=GT...AG (GT_AG)',),
                   ('exon 18922324-18922404 rna=1267-1347',
                    (('CDS 18922324-18922404 rna=1267-1347 0',),)),
                   ('intron 18922404-18922816 rna=1347-1347 sjBases=GT...AG (GT_AG)',),
                   ('exon 18922816-18922933 rna=1347-1464',
                    (('CDS 18922816-18922933 rna=1347-1464 0',),)),
                   ('intron 18922933-18923114 rna=1464-1464 sjBases=GT...AG (GT_AG)',),
                   ('exon 18923114-18923179 rna=1464-1529',
                    (('CDS 18923114-18923179 rna=1464-1529 1',),)),
                   ('intron 18923179-18925050 rna=1529-1529 sjBases=GT...AG (GT_AG)',),
                   ('exon 18925050-18925200 rna=1529-1679',
                    (('CDS 18925050-18925200 rna=1529-1679 1',),)),
                   ('intron 18925200-18925687 rna=1679-1679 sjBases=GT...AG (GT_AG)',),
                   ('exon 18925687-18925722 rna=1679-1714',
                    (('CDS 18925687-18925722 rna=1679-1714 2',),)),
                   ('intron 18925722-18930989 rna=1714-1714 sjBases=GT...AG (GT_AG)',),
                   ('exon 18930989-18931198 rna=1714-1923',
                    (('CDS 18930989-18931147 rna=1714-1872 0',),
                     ("5'UTR 18931147-18931198 rna=1872-1923",))),
                   ('intron 18931198-18936389 rna=1923-1923 sjBases=GT...AG (GT_AG)',),
                   ('exon 18936389-18936451 rna=1923-1985',
                    (("5'UTR 18936389-18936451 rna=1923-1985",),))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000334029NoSJ(self, trans, ensChroms=False):
        # coding, - strand, no splice junctions
        expect = ('t=chr22:18912781-18936451/+, rna=ENST00000334029.6:0-1985/- 1985 <->',
                  (('exon 18912781-18913362 rna=0-581',
                    (("3'UTR 18912781-18913174 rna=0-393",),
                     ('CDS 18913174-18913362 rna=393-581 1',))),
                   ('intron 18913362-18913437 rna=581-581',),
                   ('exon 18913437-18913526 rna=581-670',
                    (('CDS 18913437-18913526 rna=581-670 2',),)),
                   ('intron 18913526-18916889 rna=670-670',),
                   ('exon 18916889-18916988 rna=670-769',
                    (('CDS 18916889-18916988 rna=670-769 2',),)),
                   ('intron 18916988-18918315 rna=769-769',),
                   ('exon 18918315-18918491 rna=769-945',
                    (('CDS 18918315-18918491 rna=769-945 0',),)),
                   ('intron 18918491-18919450 rna=945-945',),
                   ('exon 18919450-18919597 rna=945-1092',
                    (('CDS 18919450-18919597 rna=945-1092 0',),)),
                   ('intron 18919597-18919705 rna=1092-1092',),
                   ('exon 18919705-18919798 rna=1092-1185',
                    (('CDS 18919705-18919798 rna=1092-1185 0',),)),
                   ('intron 18919798-18921341 rna=1185-1185',),
                   ('exon 18921341-18921423 rna=1185-1267',
                    (('CDS 18921341-18921423 rna=1185-1267 2',),)),
                   ('intron 18921423-18922324 rna=1267-1267',),
                   ('exon 18922324-18922404 rna=1267-1347',
                    (('CDS 18922324-18922404 rna=1267-1347 0',),)),
                   ('intron 18922404-18922816 rna=1347-1347',),
                   ('exon 18922816-18922933 rna=1347-1464',
                    (('CDS 18922816-18922933 rna=1347-1464 0',),)),
                   ('intron 18922933-18923114 rna=1464-1464',),
                   ('exon 18923114-18923179 rna=1464-1529',
                    (('CDS 18923114-18923179 rna=1464-1529 1',),)),
                   ('intron 18923179-18925050 rna=1529-1529',),
                   ('exon 18925050-18925200 rna=1529-1679',
                    (('CDS 18925050-18925200 rna=1529-1679 1',),)),
                   ('intron 18925200-18925687 rna=1679-1679',),
                   ('exon 18925687-18925722 rna=1679-1714',
                    (('CDS 18925687-18925722 rna=1679-1714 2',),)),
                   ('intron 18925722-18930989 rna=1714-1714',),
                   ('exon 18930989-18931198 rna=1714-1923',
                    (('CDS 18930989-18931147 rna=1714-1872 0',),
                     ("5'UTR 18931147-18931198 rna=1872-1923",))),
                   ('intron 18931198-18936389 rna=1923-1923',),
                   ('exon 18936389-18936451 rna=1923-1985',
                    (("5'UTR 18936389-18936451 rna=1923-1985",),))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000334029Rc(self, trans, ensChroms=False):
        # coding, - strand, reverse-complemented
        expect = ('t=chr22:31882017-31905687/-, rna=ENST00000334029.6:0-1985/+ 1985 <->',
                  (('exon 31882017-31882079 rna=0-62',
                    (("5'UTR 31882017-31882079 rna=0-62",),)),
                   ('intron 31882079-31887270 rna=62-62 sjBases=GT...AG (GT_AG)',),
                   ('exon 31887270-31887479 rna=62-271',
                    (("5'UTR 31887270-31887321 rna=62-113",),
                     ('CDS 31887321-31887479 rna=113-271 0',))),
                   ('intron 31887479-31892746 rna=271-271 sjBases=GT...AG (GT_AG)',),
                   ('exon 31892746-31892781 rna=271-306',
                    (('CDS 31892746-31892781 rna=271-306 2',),)),
                   ('intron 31892781-31893268 rna=306-306 sjBases=GT...AG (GT_AG)',),
                   ('exon 31893268-31893418 rna=306-456',
                    (('CDS 31893268-31893418 rna=306-456 1',),)),
                   ('intron 31893418-31895289 rna=456-456 sjBases=GT...AG (GT_AG)',),
                   ('exon 31895289-31895354 rna=456-521',
                    (('CDS 31895289-31895354 rna=456-521 1',),)),
                   ('intron 31895354-31895535 rna=521-521 sjBases=GT...AG (GT_AG)',),
                   ('exon 31895535-31895652 rna=521-638',
                    (('CDS 31895535-31895652 rna=521-638 0',),)),
                   ('intron 31895652-31896064 rna=638-638 sjBases=GT...AG (GT_AG)',),
                   ('exon 31896064-31896144 rna=638-718',
                    (('CDS 31896064-31896144 rna=638-718 0',),)),
                   ('intron 31896144-31897045 rna=718-718 sjBases=GT...AG (GT_AG)',),
                   ('exon 31897045-31897127 rna=718-800',
                    (('CDS 31897045-31897127 rna=718-800 2',),)),
                   ('intron 31897127-31898670 rna=800-800 sjBases=GT...AG (GT_AG)',),
                   ('exon 31898670-31898763 rna=800-893',
                    (('CDS 31898670-31898763 rna=800-893 0',),)),
                   ('intron 31898763-31898871 rna=893-893 sjBases=GT...AG (GT_AG)',),
                   ('exon 31898871-31899018 rna=893-1040',
                    (('CDS 31898871-31899018 rna=893-1040 0',),)),
                   ('intron 31899018-31899977 rna=1040-1040 sjBases=GC...AG (GC_AG)',),
                   ('exon 31899977-31900153 rna=1040-1216',
                    (('CDS 31899977-31900153 rna=1040-1216 0',),)),
                   ('intron 31900153-31901480 rna=1216-1216 sjBases=GT...AG (GT_AG)',),
                   ('exon 31901480-31901579 rna=1216-1315',
                    (('CDS 31901480-31901579 rna=1216-1315 2',),)),
                   ('intron 31901579-31904942 rna=1315-1315 sjBases=GT...AG (GT_AG)',),
                   ('exon 31904942-31905031 rna=1315-1404',
                    (('CDS 31904942-31905031 rna=1315-1404 2',),)),
                   ('intron 31905031-31905106 rna=1404-1404 sjBases=GT...AG (GT_AG)',),
                   ('exon 31905106-31905687 rna=1404-1985',
                    (('CDS 31905106-31905294 rna=1404-1592 1',),
                     ("3'UTR 31905294-31905687 rna=1592-1985",)))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000334029NoSJRc(self, trans, ensChroms=False):
        # coding, - strand, reverse-complemented, no splice junctions
        expect = ('t=chr22:31882017-31905687/-, rna=ENST00000334029.6:0-1985/+ 1985 <->',
                  (('exon 31882017-31882079 rna=0-62',
                    (("5'UTR 31882017-31882079 rna=0-62",),)),
                   ('intron 31882079-31887270 rna=62-62',),
                   ('exon 31887270-31887479 rna=62-271',
                    (("5'UTR 31887270-31887321 rna=62-113",),
                     ('CDS 31887321-31887479 rna=113-271 0',))),
                   ('intron 31887479-31892746 rna=271-271',),
                   ('exon 31892746-31892781 rna=271-306',
                    (('CDS 31892746-31892781 rna=271-306 2',),)),
                   ('intron 31892781-31893268 rna=306-306',),
                   ('exon 31893268-31893418 rna=306-456',
                    (('CDS 31893268-31893418 rna=306-456 1',),)),
                   ('intron 31893418-31895289 rna=456-456',),
                   ('exon 31895289-31895354 rna=456-521',
                    (('CDS 31895289-31895354 rna=456-521 1',),)),
                   ('intron 31895354-31895535 rna=521-521',),
                   ('exon 31895535-31895652 rna=521-638',
                    (('CDS 31895535-31895652 rna=521-638 0',),)),
                   ('intron 31895652-31896064 rna=638-638',),
                   ('exon 31896064-31896144 rna=638-718',
                    (('CDS 31896064-31896144 rna=638-718 0',),)),
                   ('intron 31896144-31897045 rna=718-718',),
                   ('exon 31897045-31897127 rna=718-800',
                    (('CDS 31897045-31897127 rna=718-800 2',),)),
                   ('intron 31897127-31898670 rna=800-800',),
                   ('exon 31898670-31898763 rna=800-893',
                    (('CDS 31898670-31898763 rna=800-893 0',),)),
                   ('intron 31898763-31898871 rna=893-893',),
                   ('exon 31898871-31899018 rna=893-1040',
                    (('CDS 31898871-31899018 rna=893-1040 0',),)),
                   ('intron 31899018-31899977 rna=1040-1040',),
                   ('exon 31899977-31900153 rna=1040-1216',
                    (('CDS 31899977-31900153 rna=1040-1216 0',),)),
                   ('intron 31900153-31901480 rna=1216-1216',),
                   ('exon 31901480-31901579 rna=1216-1315',
                    (('CDS 31901480-31901579 rna=1216-1315 2',),)),
                   ('intron 31901579-31904942 rna=1315-1315',),
                   ('exon 31904942-31905031 rna=1315-1404',
                    (('CDS 31904942-31905031 rna=1315-1404 2',),)),
                   ('intron 31905031-31905106 rna=1404-1404',),
                   ('exon 31905106-31905687 rna=1404-1985',
                    (('CDS 31905106-31905294 rna=1404-1592 1',),
                     ("3'UTR 31905294-31905687 rna=1592-1985",)))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000434390(self, trans, ensChroms=False):
        # non-coding, - strand
        expect = ('t=chr22:18178037-18205915/+, rna=ENST00000434390.1:0-1859/- 1859 <->',
                  (('exon 18178037-18178465 rna=0-428',
                    (('NC 18178037-18178465 rna=0-428',),)),
                   ('intron 18178465-18183109 rna=428-428 sjBases=GT...AG (GT_AG)',),
                   ('exon 18183109-18184066 rna=428-1385',
                    (('NC 18183109-18184066 rna=428-1385',),)),
                   ('intron 18184066-18188622 rna=1385-1385 sjBases=GT...AG (GT_AG)',),
                   ('exon 18188622-18188656 rna=1385-1419',
                    (('NC 18188622-18188656 rna=1385-1419',),)),
                   ('intron 18188656-18189121 rna=1419-1419 sjBases=GT...AG (GT_AG)',),
                   ('exon 18189121-18189183 rna=1419-1481',
                    (('NC 18189121-18189183 rna=1419-1481',),)),
                   ('intron 18189183-18191182 rna=1481-1481 sjBases=GT...AG (GT_AG)',),
                   ('exon 18191182-18191248 rna=1481-1547',
                    (('NC 18191182-18191248 rna=1481-1547',),)),
                   ('intron 18191248-18196156 rna=1547-1547 sjBases=GT...AG (GT_AG)',),
                   ('exon 18196156-18196243 rna=1547-1634',
                    (('NC 18196156-18196243 rna=1547-1634',),)),
                   ('intron 18196243-18199189 rna=1634-1634 sjBases=GT...AG (GT_AG)',),
                   ('exon 18199189-18199214 rna=1634-1659',
                    (('NC 18199189-18199214 rna=1634-1659',),)),
                   ('intron 18199214-18199682 rna=1659-1659 sjBases=GT...AG (GT_AG)',),
                   ('exon 18199682-18199741 rna=1659-1718',
                    (('NC 18199682-18199741 rna=1659-1718',),)),
                   ('intron 18199741-18205774 rna=1718-1718 sjBases=GT...AG (GT_AG)',),
                   ('exon 18205774-18205915 rna=1718-1859',
                    (('NC 18205774-18205915 rna=1718-1859',),))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000434390Rc(self, trans, ensChroms=False):
        # non-coding, - strand, reverse complemented
        expect = ('t=chr22:32612553-32640431/-, rna=ENST00000434390.1:0-1859/+ 1859 <->',
                  (('exon 32612553-32612694 rna=0-141',
                    (('NC 32612553-32612694 rna=0-141',),)),
                   ('intron 32612694-32618727 rna=141-141 sjBases=GT...AG (GT_AG)',),
                   ('exon 32618727-32618786 rna=141-200',
                    (('NC 32618727-32618786 rna=141-200',),)),
                   ('intron 32618786-32619254 rna=200-200 sjBases=GT...AG (GT_AG)',),
                   ('exon 32619254-32619279 rna=200-225',
                    (('NC 32619254-32619279 rna=200-225',),)),
                   ('intron 32619279-32622225 rna=225-225 sjBases=GT...AG (GT_AG)',),
                   ('exon 32622225-32622312 rna=225-312',
                    (('NC 32622225-32622312 rna=225-312',),)),
                   ('intron 32622312-32627220 rna=312-312 sjBases=GT...AG (GT_AG)',),
                   ('exon 32627220-32627286 rna=312-378',
                    (('NC 32627220-32627286 rna=312-378',),)),
                   ('intron 32627286-32629285 rna=378-378 sjBases=GT...AG (GT_AG)',),
                   ('exon 32629285-32629347 rna=378-440',
                    (('NC 32629285-32629347 rna=378-440',),)),
                   ('intron 32629347-32629812 rna=440-440 sjBases=GT...AG (GT_AG)',),
                   ('exon 32629812-32629846 rna=440-474',
                    (('NC 32629812-32629846 rna=440-474',),)),
                   ('intron 32629846-32634402 rna=474-474 sjBases=GT...AG (GT_AG)',),
                   ('exon 32634402-32635359 rna=474-1431',
                    (('NC 32634402-32635359 rna=474-1431',),)),
                   ('intron 32635359-32640003 rna=1431-1431 sjBases=GT...AG (GT_AG)',),
                   ('exon 32640003-32640431 rna=1431-1859',
                    (('NC 32640003-32640431 rna=1431-1859',),))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)

    def checkENST00000538324(self, trans, ensChroms=False):
        # ABO error in genome
        expect = ('t=chr9:133255601-133275214/+, rna=ENST00000538324.2:0-1153/- 1153 <->',
                  (('exon 133255601-133256356 rna=0-755',
                    (('CDS 133255601-133255670 rna=0-69 0',),
                     ('gap 133255670-133255674 rna=69-73',),
                     ('CDS 133255674-133256356 rna=73-755 2',))),
                   ('intron 133256356-133257408 rna=755-755 sjBases=GT...AG (GT_AG)',),
                   ('exon 133257408-133257542 rna=755-889',
                    (('CDS 133257408-133257521 rna=755-868 0',),
                     ('gap 133257521-133257523 rna=868-870',),
                     ('CDS 133257523-133257542 rna=870-889 2',))),
                   ('intron 133257542-133258096 rna=889-889 sjBases=GT...AG (GT_AG)',),
                   ('exon 133258096-133258132 rna=889-925',
                    (('CDS 133258096-133258132 rna=889-925 2',),)),
                   ('intron 133258132-133259818 rna=925-925 sjBases=GT...AG (GT_AG)',),
                   ('exon 133259818-133259866 rna=925-973',
                    (('CDS 133259818-133259866 rna=925-973 2',),)),
                   ('intron 133259866-133261317 rna=973-973 sjBases=GT...AG (GT_AG)',),
                   ('exon 133261317-133261374 rna=973-1030',
                    (('CDS 133261317-133261374 rna=973-1030 2',),)),
                   ('intron 133261374-133262098 rna=1030-1030 sjBases=GT...AG (GT_AG)',),
                   ('exon 133262098-133262168 rna=1030-1100',
                    (('CDS 133262098-133262168 rna=1030-1100 1',),)),
                   ('intron 133262168-133275161 rna=1100-1100 sjBases=GT...AG (GT_AG)',),
                   ('exon 133275161-133275214 rna=1100-1153',
                    (('CDS 133275161-133275189 rna=1100-1128 0',),
                     ("5'UTR 133275189-133275214 rna=1128-1153",)))))
        if ensChroms:
            expect = self._toEnsemblExpect(expect)
        self._assertFeatures(trans, expect)


class GenePredAnnotationTests(FeatureTestBase, AnnotationCheckMixin):
    def _getTransAnnot(self, acc):
        factory = GenePredAnnotationFactory(GenomeSeqSrc.obtain("hg38"))
        return factory.fromGenePred(GenePredDbSrc.obtainGenePred("V28", acc))

    def _getTransAnnotNoSJ(self, acc):
        factory = GenePredAnnotationFactory(chromSizeFunc=GenomeSeqSrc.obtain("hg38").getChromSize)
        return factory.fromGenePred(GenePredDbSrc.obtainGenePred("V28", acc))

    def testENST00000215794(self):
        # coding, + strand
        trans = self._getTransAnnot("ENST00000215794.7")
        self.checkENST00000215794(trans)
        self.assertEqual("chr22\t18149898\t18177397\tENST00000215794.7\t0\t+\t18157663\t18176817\t0,1,2\t11\t324,263,97,146,80,147,96,168,132,50,626,\t0,7659,10273,11891,17356,17991,19945,20854,23251,23894,26873,",
                         str(trans.toBed("0,1,2")),)

    def testENST00000334029(self):
        # coding, - strand
        trans = self._getTransAnnot("ENST00000334029.6")
        self.checkENST00000334029(trans)

    def testENST00000334029NoSJ(self):
        # coding, - strand, no splice junctions
        trans = self._getTransAnnotNoSJ("ENST00000334029.6")
        self.checkENST00000334029NoSJ(trans)

    def testENST00000334029Rc(self):
        # coding, - strand, reverse-complemented
        trans = self._getTransAnnot("ENST00000334029.6")
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self.checkENST00000334029Rc(rcTrans)

    def testENST00000334029NoSJRc(self):
        # coding, - strand, reverse-complemented, no splice junctions
        trans = self._getTransAnnotNoSJ("ENST00000334029.6")
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self.checkENST00000334029NoSJRc(rcTrans)

    def testAnnotPrevNext(self):
        trans = self._getTransAnnot("ENST00000334029.6")

        feat = trans.firstFeature(AnnotationFeature)
        self.assertIs(feat, trans.features[0].annotFeatures[0])
        cdsCnt = 1 if isinstance(trans, CdsRegionFeature) else 0
        while True:
            feat = feat.nextFeature(CdsRegionFeature)
            if feat is None:
                break
            cdsCnt += 1
        self.assertEqual(cdsCnt, 13)

        feat = trans.lastFeature(AnnotationFeature)
        cdsCnt = 1 if isinstance(trans, CdsRegionFeature) else 0
        while True:
            feat = feat.prevFeature(CdsRegionFeature)
            if feat is None:
                break
            cdsCnt += 1
        self.assertEqual(cdsCnt, 13)

    def testGetStructureFeaturesOfType(self):
        # coding, - strand, exon features
        trans = self._getTransAnnot("ENST00000334029.6")
        feats = trans.getFeaturesOfType(ExonFeature)
        featStrs = tuple([str(f) for f in feats])
        self.assertEqual(('exon 18912781-18913362 rna=0-581',
                          'exon 18913437-18913526 rna=581-670',
                          'exon 18916889-18916988 rna=670-769',
                          'exon 18918315-18918491 rna=769-945',
                          'exon 18919450-18919597 rna=945-1092',
                          'exon 18919705-18919798 rna=1092-1185',
                          'exon 18921341-18921423 rna=1185-1267',
                          'exon 18922324-18922404 rna=1267-1347',
                          'exon 18922816-18922933 rna=1347-1464',
                          'exon 18923114-18923179 rna=1464-1529',
                          'exon 18925050-18925200 rna=1529-1679',
                          'exon 18925687-18925722 rna=1679-1714',
                          'exon 18930989-18931198 rna=1714-1923',
                          'exon 18936389-18936451 rna=1923-1985'),
                         featStrs)

    def testGetAnnotationFeaturesOfType(self):
        # coding, - strand, CDS features
        trans = self._getTransAnnot("ENST00000334029.6")
        feats = trans.getFeaturesOfType(CdsRegionFeature)
        featStrs = tuple([str(f) for f in feats])
        self.assertEqual(('CDS 18913174-18913362 rna=393-581 1',
                          'CDS 18913437-18913526 rna=581-670 2',
                          'CDS 18916889-18916988 rna=670-769 2',
                          'CDS 18918315-18918491 rna=769-945 0',
                          'CDS 18919450-18919597 rna=945-1092 0',
                          'CDS 18919705-18919798 rna=1092-1185 0',
                          'CDS 18921341-18921423 rna=1185-1267 2',
                          'CDS 18922324-18922404 rna=1267-1347 0',
                          'CDS 18922816-18922933 rna=1347-1464 0',
                          'CDS 18923114-18923179 rna=1464-1529 1',
                          'CDS 18925050-18925200 rna=1529-1679 1',
                          'CDS 18925687-18925722 rna=1679-1714 2',
                          'CDS 18930989-18931147 rna=1714-1872 0'),
                         featStrs)

    def testENST00000434390(self):
        # non-coding, - strand
        trans = self._getTransAnnot("ENST00000434390.1")
        self.checkENST00000434390(trans)

    def testENST00000434390Rc(self):
        # non-coding, - strand, reverse complemented
        trans = self._getTransAnnot("ENST00000434390.1")
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self.checkENST00000434390Rc(rcTrans)

    def testENST00000538324(self):
        # ABO error in genome
        trans = self._getTransAnnot("ENST00000538324.2")
        self.checkENST00000538324(trans)

    def testGencodeV28Regress(self):
        "regression test for V26 and V28"
        # 4-base gap that caused error, just try to parse
        trans = self._getTransAnnot("ENST00000610542.1")
        self.assertEqual(trans.rna.name, "ENST00000610542.1")

    def testAnnotToBed(self):
        """test for conversion to BED"""
        # had a bug non-code genes had only basic columns because thickStart/thickEnd were None
        factory = GenePredAnnotationFactory()
        for gp in GenePredReader(self.getInputFile("gencodeCompV28.gp")):
            trans = factory.fromGenePred(gp)
            bed = trans.toBed("100,0,0")
            self.assertEqual(len(bed.getRow()), 12)

    def testAttrs(self):
        def getUtr3Feature(trans):
            # test case has one 3'UTR feature
            return trans.getFeaturesOfType(Utr3RegionFeature)[0]

        factory = GenePredAnnotationFactory(GenomeSeqSrc.obtain("hg38"))
        tattrs = ObjDict(name="Fred")
        uattrs = ObjDict(name="Barney")
        trans = factory.fromGenePred(GenePredDbSrc.obtainGenePred("V28", "ENST00000334029.6"), tattrs)
        # has one 3'UTR feature, so use it for the attrs
        utr3 = getUtr3Feature(trans)
        utr3.attrs = uattrs

        rcTrans = trans.reverseComplement()
        rcUtr3 = getUtr3Feature(rcTrans)
        self.assertEqual("Fred", rcTrans.attrs.name)
        self.assertEqual("Barney", rcUtr3.attrs.name)


class EnsemblDbAnnotationTests(FeatureTestBase, AnnotationCheckMixin):
    # HOST = "ensembldb.ensembl.org"
    HOST = "useastdb.ensembl.org"
    CORE_DB = "homo_sapiens_core_92_38"

    conn = None
    annotFactory = None

    @classmethod
    def _getConnect(cls):
        if cls.conn is None:
            cls.conn = mysqlOps.connect(host=cls.HOST, port=5306, user="anonymous", password="", db=cls.CORE_DB)
        return cls.conn

    @classmethod
    def _getAnnotFactory(cls):
        if cls.annotFactory is None:
            cls.annotFactory = EnsemblDbAnnotationFactory(GenomeSeqSrc.obtain("grch38"))
        return cls.annotFactory

    def _getTransAnnot(self, transId):
        ensGenes = ensemblGeneQuery(self._getConnect(), transId)
        self.assertEqual(len(ensGenes), 1)
        self.assertEqual(len(ensGenes[0].transcripts), 1)
        return self._getAnnotFactory().fromEnsemblDb(ensGenes[0].transcripts[0])

    def testENST00000215794(self):
        # + strand
        trans = self._getTransAnnot("ENST00000215794.7")
        self.checkENST00000215794(trans, ensChroms=True)

    def testENST00000334029(self):
        # coding, - strand
        trans = self._getTransAnnot("ENST00000334029.6")
        self.checkENST00000334029(trans, ensChroms=True)

    def testENST00000434390(self):
        # non-coding, - strand
        trans = self._getTransAnnot("ENST00000434390.1")
        self.checkENST00000434390(trans, ensChroms=True)

    def testENST00000434390Rc(self):
        # non-coding, - strand, reverse complemented
        trans = self._getTransAnnot("ENST00000434390.1")
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self.checkENST00000434390Rc(rcTrans, ensChroms=True)

    def testENST00000538324(self):
        # ABO error in genome
        trans = self._getTransAnnot("ENST00000538324.2")
        self.checkENST00000538324(trans, ensChroms=True)

if __name__ == '__main__':
    unittest.main()
