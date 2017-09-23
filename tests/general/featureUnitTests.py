from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from gencode_icedb.general.genome import GenomeReaderFactory
from gencode_icedb.general.evidFeatures import EvidenceFeatureMap, EvidencePslFactory
from gencode_icedb.general.annotFeatures import AnnotationGenePredFactory
from pycbio.hgdata.hgLite import PslDbTable, GenePredDbTable
from pycbio.hgdata.genePred import GenePredReader
import sqlite3
import pprint


twoBitHg19 = "/hive/data/genomes/hg19/hg19.2bit"
mockReaderHg19PslTsv = "mockReaderHg19Psl.tsv"
mockReaderHg19GpTsv = "mockReaderHg19Gp.tsv"

updateMockReader = False   # set this to update from a real run
forceMockReader = False  # set this to check mock data

debugResults = False   # print out results for updated expected
noCheckResults = False  # don't check results

if updateMockReader or forceMockReader or debugResults or noCheckResults:
    print("Warning: debug variables set", file=sys.stderr)
if updateMockReader and forceMockReader:
    raise Exception("makes no sense to have both updateMockReader and forceMockReader set")


class FeatureTestBase(TestCaseBase):
    def _assertFeatures(self, trans, expect):
        if debugResults:
            print("==== {} ==== ".format(self.id()))
            pprint.pprint(trans.toStrTree(), width=1)
            print()
        if not noCheckResults:
            self.assertEqual(expect, trans.toStrTree())


class EvidenceTests(FeatureTestBase):
    genomeReaderFactory = None
    set1PslDbTbl = None

    def __obtainGenomeReader(self):
        if EvidenceTests.genomeReaderFactory is None:
            EvidenceTests.genomeReaderFactory = GenomeReaderFactory(twoBitHg19, self.getInputFile(mockReaderHg19PslTsv),
                                                                    updateMockReader, forceMockReader)
        return EvidenceTests.genomeReaderFactory.obtain()

    def __obtainSetPslDbTbl(self):
        if EvidenceTests.set1PslDbTbl is None:
            conn = sqlite3.connect(":memory:")
            EvidenceTests.set1PslDbTbl = PslDbTable(conn, "set1", create=True)
            EvidenceTests.set1PslDbTbl.loadPslFile(self.getInputFile("set1.ucsc-mrna.psl"))
        return EvidenceTests.set1PslDbTbl

    def __getSet1Psl(self, acc):
        pslDbTbl = self.__obtainSetPslDbTbl()
        psls = pslDbTbl.getByQName(acc)
        if len(psls) == 0:
            raise Exception("psl not found: {}".format(acc))
        return psls[0]

    def __pslToEvidTranscript(self, psl):
        factory = EvidencePslFactory(self.__obtainGenomeReader())
        return factory.fromPsl(psl)

    def testAF010310(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18900294-18905926/+, rna=AF010310.1:0-888/- 901',
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
                               ('intron 18900875-18900950 rna=618-618 sjBases=GT...AG (GT_AG)',),
                               ('exon 18900950-18901039 rna=618-707',
                                (('aln 18900950-18901039 rna=618-707',),)),
                               ('intron 18901039-18904402 rna=707-707 sjBases=GT...AG (GT_AG)',),
                               ('exon 18904402-18904501 rna=707-806',
                                (('aln 18904402-18904501 rna=707-806',),)),
                               ('intron 18904501-18905828 rna=806-806 sjBases=GT...AG (GT_AG)',),
                               ('exon 18905828-18905926 rna=806-901',
                                (('aln 18905828-18905889 rna=806-867',),
                                 ('rins None-None rna=867-868',),
                                 ('cins 18905889-18905893 rna=None-None',),
                                 ('aln 18905893-18905926 rna=868-901',))))))

    def testX96484(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("X96484.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18893922-18899592/+, rna=X96484.1:48-1067/+ 1080',
                              (('exon 18893922-18893997 rna=48-123',
                                (('aln 18893922-18893997 rna=48-123',),)),
                               ('intron 18893997-18894077 rna=123-123 sjBases=GT...AG (GT_AG)',),
                               ('exon 18894077-18894238 rna=123-285',
                                (('aln 18894077-18894173 rna=123-219',),
                                 ('rins None-None rna=219-220',),
                                 ('aln 18894173-18894238 rna=220-285',))),
                               ('intron 18894238-18897684 rna=285-285 sjBases=GT...AG (GT_AG)',),
                               ('exon 18897684-18897785 rna=285-386',
                                (('aln 18897684-18897785 rna=285-386',),)),
                               ('intron 18897785-18898400 rna=386-386 sjBases=GT...AG (GT_AG)',),
                               ('exon 18898400-18898541 rna=386-527',
                                (('aln 18898400-18898541 rna=386-527',),)),
                               ('intron 18898541-18899052 rna=527-527 sjBases=GT...AG (GT_AG)',),
                               ('exon 18899052-18899592 rna=527-1067',
                                (('aln 18899052-18899592 rna=527-1067',),)))))

    def testX96484NoSJ(self):
        factory = EvidencePslFactory(None)
        trans = factory.fromPsl(self.__getSet1Psl("X96484.1"))
        self._assertFeatures(trans,
                             ('t=chr22:18893922-18899592/+, rna=X96484.1:48-1067/+ 1080',
                              (('exon 18893922-18893997 rna=48-123',
                                (('aln 18893922-18893997 rna=48-123',),)),
                               ('intron 18893997-18894077 rna=123-123 sjBases=None',),
                               ('exon 18894077-18894238 rna=123-285',
                                (('aln 18894077-18894173 rna=123-219',),
                                 ('rins None-None rna=219-220',),
                                 ('aln 18894173-18894238 rna=220-285',))),
                               ('intron 18894238-18897684 rna=285-285 sjBases=None',),
                               ('exon 18897684-18897785 rna=285-386',
                                (('aln 18897684-18897785 rna=285-386',),)),
                               ('intron 18897785-18898400 rna=386-386 sjBases=None',),
                               ('exon 18898400-18898541 rna=386-527',
                                (('aln 18898400-18898541 rna=386-527',),)),
                               ('intron 18898541-18899052 rna=527-527 sjBases=None',),
                               ('exon 18899052-18899592 rna=527-1067',
                                (('aln 18899052-18899592 rna=527-1067',),)))))

    def testAF010310Rc(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans,
                             ('t=chr22:32398640-32404272/-, rna=AF010310.1:13-901/+ 901',
                              (('exon 32398640-32398738 rna=0-95',
                                (('aln 32398640-32398673 rna=0-33',),
                                 ('cins 32398673-32398677 rna=None-None',),
                                 ('rins None-None rna=33-34',),
                                 ('aln 32398677-32398738 rna=34-95',))),
                               ('intron 32398738-32400065 rna=95-95 sjBases=GT...AG (GT_AG)',),
                               ('exon 32400065-32400164 rna=95-194',
                                (('aln 32400065-32400164 rna=95-194',),)),
                               ('intron 32400164-32403527 rna=194-194 sjBases=GT...AG (GT_AG)',),
                               ('exon 32403527-32403616 rna=194-283',
                                (('aln 32403527-32403616 rna=194-283',),)),
                               ('intron 32403616-32403691 rna=283-283 sjBases=GT...AG (GT_AG)',),
                               ('exon 32403691-32404272 rna=283-888',
                                (('aln 32403691-32403695 rna=283-287',),
                                 ('rins None-None rna=287-288',),
                                 ('aln 32403695-32403763 rna=288-356',),
                                 ('rins None-None rna=356-357',),
                                 ('aln 32403763-32403776 rna=357-370',),
                                 ('rins None-None rna=370-371',),
                                 ('aln 32403776-32403804 rna=371-399',),
                                 ('rins None-None rna=399-400',),
                                 ('aln 32403804-32403897 rna=400-493',),
                                 ('rins None-None rna=493-494',),
                                 ('aln 32403897-32403902 rna=494-499',),
                                 ('rins None-None rna=499-500',),
                                 ('aln 32403902-32403921 rna=500-519',),
                                 ('rins None-None rna=519-520',),
                                 ('aln 32403921-32403923 rna=520-522',),
                                 ('rins None-None rna=522-523',),
                                 ('aln 32403923-32403932 rna=523-532',),
                                 ('cins 32403932-32403938 rna=None-None',),
                                 ('rins None-None rna=532-542',),
                                 ('aln 32403938-32403945 rna=542-549',),
                                 ('rins None-None rna=549-550',),
                                 ('aln 32403945-32403950 rna=550-555',),
                                 ('rins None-None rna=555-556',),
                                 ('aln 32403950-32403955 rna=556-561',),
                                 ('rins None-None rna=561-562',),
                                 ('aln 32403955-32403959 rna=562-566',),
                                 ('cins 32403959-32403960 rna=None-None',),
                                 ('rins None-None rna=566-569',),
                                 ('aln 32403960-32403987 rna=569-596',),
                                 ('rins None-None rna=596-597',),
                                 ('aln 32403987-32403989 rna=597-599',),
                                 ('rins None-None rna=599-600',),
                                 ('aln 32403989-32403999 rna=600-610',),
                                 ('rins None-None rna=610-611',),
                                 ('aln 32403999-32404001 rna=611-613',),
                                 ('rins None-None rna=613-614',),
                                 ('aln 32404001-32404031 rna=614-644',),
                                 ('rins None-None rna=644-645',),
                                 ('aln 32404031-32404113 rna=645-727',),
                                 ('cins 32404113-32404114 rna=None-None',),
                                 ('rins None-None rna=727-730',),
                                 ('aln 32404114-32404272 rna=730-888',))))))

    def testRangeMap1(self):
        # range is set1: chr22:18632931-19279166
        pslDbTbl = self.__obtainSetPslDbTbl()
        evidFeatureMap = EvidenceFeatureMap.dbFactory(pslDbTbl.conn, pslDbTbl.table,
                                                      "chr22", 18958026, 19109719,
                                                      self.__obtainGenomeReader())
        self.assertEqual(len(evidFeatureMap.transcripts), 30)
        overFeats = list(evidFeatureMap.overlapping("chr22", 18958026, 18982141))
        self.assertEqual(len(overFeats), 12)


class AnnotationTests(FeatureTestBase):
    genomeReaderFactory = None
    set1GpDbTbl = None

    def __obtainGenomeReader(self):
        if AnnotationTests.genomeReaderFactory is None:
            AnnotationTests.genomeReaderFactory = GenomeReaderFactory(twoBitHg19, self.getInputFile(mockReaderHg19GpTsv),
                                                                      updateMockReader, forceMockReader)
        return AnnotationTests.genomeReaderFactory.obtain()

    def __obtainSetGpDbTbl(self):
        if AnnotationTests.set1GpDbTbl is None:
            conn = sqlite3.connect(":memory:")
            AnnotationTests.set1GpDbTbl = GenePredDbTable(conn, "set1", create=True)
            AnnotationTests.set1GpDbTbl.loadGenePredFile(self.getInputFile("set1.gencodeCompV19.gp"))
        return AnnotationTests.set1GpDbTbl

    def __getSet1Gp(self, acc):
        gpDbTbl = self.__obtainSetGpDbTbl()
        gps = gpDbTbl.getByName(acc)
        if len(gps) == 0:
            raise Exception("gp not found: {}".format(acc))
        return gps[0]

    def __gpToEvidTranscript(self, gp):
        factory = AnnotationGenePredFactory(self.__obtainGenomeReader())
        return factory.fromGenePred(gp)

    def testENST00000215794(self):
        # + strand
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000215794.7"))
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
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000334029.2"))
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
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000334029.2"))
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
        factory = AnnotationGenePredFactory(chromSizeFunc=self.__obtainGenomeReader().getChromSize)
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

    def testENST00000434390(self):
        # non-coding, - strand
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000434390.1"))
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
            self.assertEqual(trans.rnaName, names[i])

    def testAnnotToBed(self):
        """test for conversion to BED"""
        # had a bug non-code genes had only basic columns because thickStart/thickEnd were None
        factory = AnnotationGenePredFactory()
        for gp in GenePredReader(self.getInputFile("set1.gencodeCompV19.gp")):
            trans = factory.fromGenePred(gp)
            bed = trans.toBed("100,0,0")
            self.assertEqual(len(bed.getRow()), 12)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts


if __name__ == '__main__':
    unittest.main()
