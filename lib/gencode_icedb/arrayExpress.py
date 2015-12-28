import sys, os, re
from pycbio.tsv import TsvReader
from collections import defaultdict
from gencode_icedb import rnaSeqData
from peewee import SqliteDatabase

class ArrayExpressSdrfMapper(object):
    """Object to TSV column names to field names.  Counting occurrences of
    non-unique names and make them unique. Convert to lower case and make
    consistent `Comment [ENA_RUN]' and `Comment[ENA_RUN]'
    """
    def __init__(self):
        self.nameCount = dict()
    def __call__(self, cname):
        # replace ' *[' or ' ' with `_' and remove ']'
        cname = re.sub('( *\\[)|(\\])|( )',
                       lambda m: '' if m.group(0)==']' else '_', cname).lower()
        if cname in self.nameCount:
            self.nameCount[cname] += 1
            cname += str(self.nameCount[cname])
        else:
            self.nameCount[cname] = 1
        return cname

def arrayExpressSdrfReader(sdrfTsv, typeMap=None):
    """create a reader for an sdrf TSV file.  The column names
    are not consistent on all files and some are duplicated. They need to be mapped
    to valid python names (e.g. 'Factor Value[organism part] and made unique')
    """
    return TsvReader(sdrfTsv, columnNameMapper=ArrayExpressSdrfMapper(), typeMap=typeMap)

class ArrayExpressRuns(defaultdict):
    """Rows from ArrayExpress SDRF files by run.  Rows will be paired for pair-end
    reads, or one length arrays for unpaired, it will only have one entry per
    row
    """
    def __init__(self, sdrfTsv, scientificName):
        super(ArrayExpressRuns, self).__init__(list)
        self.scientificName = scientificName
        self.__load(sdrfTsv)

    def __keepRow(self, row):
        return row.characteristics_organism == self.scientificName

    def __load(self, sdrfTsv):
        for row in arrayExpressSdrfReader(sdrfTsv):
            if self.__keepRow(row):
                self.__addRow(row)

    def __addRow(self, row):
        run = row.comment_ena_run
        self[run].append(row)
        if len(self[run]) > 2:
            raise Exception("too many entries for run: " + run)

class ArrayExpressRegister(object):
    """Register RNA-Seq data from and ArrayExpress SDRF TSV files.
    filter by organism and locating files following PathConfig rules.
    """
    def __init__(self, pathConfig, scientificName, rnaSeqSetName, description,
                 skipExisting, verbose):
        self.pathConfig = pathConfig
        self.scientificName = scientificName
        self.rnaSeqSetName = rnaSeqSetName
        self.description = description
        self.skipExisting = skipExisting
        self.verbose = verbose
        self.database = SqliteDatabase(pathConfig.dataDatabase())

    def __verbPr(self, msg, rnaSeqRun):
        if self.verbose:
            sys.stderr.write(msg + ": " + self.rnaSeqSetName + "/" + rnaSeqRun + "\n")

    def __registerRun(self, loader, run):
        fq1 = run[0]
        fq2 = run[1] if len(run) > 1 else None
        status = loader.registerRun(fq1.comment_ena_run, self.description,
                                    fq1.characteristics_organism_part,
                                    fq1.comment_fastq_uri,
                                    os.path.split(fq1.comment_fastq_uri)[1],
                                    fq2.comment_fastq_uri if fq2 is not None else None,
                                    os.path.split(fq2.comment_fastq_uri)[1] if fq2 is not None else None)
        self.__verbPr(str(status), fq1.comment_ena_run)

    def load(self, sdrfTsv):
        loader = rnaSeqData.RnaSeqDataLoader(self.database, self.pathConfig,
                                             self.rnaSeqSetName, self.pathConfig.organism,
                                             self.skipExisting)
        runs = ArrayExpressRuns(sdrfTsv, self.scientificName)
        for run in runs.itervalues():
            self.__registerRun(loader, run)
