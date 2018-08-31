"""
Functions around reading and classifying sequence
"""
import os
from pycbio.hgdata import dnaOps
from twobitreader import TwoBitFile
from pysam.libcfaidx import FastaFile


class GenomeReader(object):
    """Base class for genome sequence file readers. Also serves as a factory."""

    @classmethod
    def getFromFileName(cls, genomeFile):
        """get read for *.2bit or *.fa file"""
        if genomeFile.endswith(".2bit"):
            return GenomeReaderTwoBit(genomeFile)
        elif genomeFile.endswith(".fa") or genomeFile.endswith(".fa.gz"):
            return GenomeReaderFasta(genomeFile)
        else:
            raise Exception("unrecognized genome sequence file extension: {}".format(genomeFile))

    @classmethod
    def getFromUcscDbName(cls, ucscDb):
        """create GenomeReader from a twobit in standard UCSC development
        server locations, checking for cluster location"""
        tb = "/scratch/data/{ucscDb}/{ucscDb}.2bit".format(ucscDb=ucscDb)
        if not os.path.exists(tb):
            tb = "/hive/data/genomes/{ucscDb}/{ucscDb}.2bit".format(ucscDb=ucscDb)
        return cls.getFromFileName(tb)

    @classmethod
    def addCmdOptions(cls, parser):
        """Add command options used to create a factory to the parser"""
        parser.add_argument('--genomeSeqs',
                            help="""Genome sequence twobit or false file to obtain splice sites""")

    @classmethod
    def getFromCmdOptions(cls, opts):
        """create a GenomeReader given option parse from command line"""
        # this will check sanity of options
        if opts.genomeSeqs is None:
            raise Exception("must specify --genomeSeqs")
        return GenomeReader.getFromFileName(opts.genomeSeqs)


class GenomeReaderTwoBit(GenomeReader):
    """
    Reads sequences from TwoBit format files.
    """
    def __init__(self, twoBitFile):
        self.twoBitFile = twoBitFile
        self.twoBitReader = TwoBitFile(twoBitFile)
        self.chrom = self.size = None

    def close(self):
        self.twoBitReader.close()
        self.twoBitReader = None
        self.chrom = self.size = None

    def getOptionArgs(self):
        "create a vector of options and values for passing to another program"
        return ["--genomeSeqs={}".format(self.twoBitFile)]

    def _obtainChrom(self, chrom):
        if self.chrom != chrom:
            self.twoBitSeq = self.twoBitReader[chrom]
            self.chrom = chrom
            self.size = len(self.twoBitSeq)

    def get(self, chrom, start, end, strand=None):
        self._obtainChrom(chrom)
        if strand == '-':
            start, end = dnaOps.reverseCoords(start, end, self.size)
        seq = self.twoBitSeq[start:end]
        if strand == '-':
            seq = dnaOps.reverseComplement(seq)
        return seq

    def haveChrom(self, chrom):
        return chrom in self.twoBitReader

    def getChroms(self):
        return sorted(self.twoBitReader.keys())

    def getChromSize(self, chrom):
        self._obtainChrom(chrom)
        return self.size


class GenomeReaderFasta(GenomeReader):
    """
    Reads sequences from indexed fasta format files.
    """
    def __init__(self, faFile):
        self.faFile = faFile
        self.faReader = FastaFile(faFile)
        self.chrom = self.size = None

    def close(self):
        self.faReader.close()
        self.faReader = None

    def getOptionArgs(self):
        "create a vector of options and values for passing to another program"
        return ["--genomeSeqs={}".format(self.faFile)]

    def get(self, chrom, start, end, strand=None):
        if strand == '-':
            start, end = dnaOps.reverseCoords(start, end, self.size)
        seq = self.faReader.fetch(chrom, start, end)
        if strand == '-':
            seq = dnaOps.reverseComplement(seq)
        return seq

    def haveChrom(self, chrom):
        return chrom in self.faReader

    def getChroms(self):
        return sorted(self.faReader.references)

    def getChromSize(self, chrom):
        return self.faReader.get_reference_length(chrom)
