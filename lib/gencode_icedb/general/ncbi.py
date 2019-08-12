"""
Class and functions for dealing with NCBI data.
"""
from datetime import datetime
import csv
from pycbio.tsv import TsvReader, strOrNoneType
from pycbio.sys import fileOps


# run-info columns
#    Run ReleaseDate LoadDate spots bases spots_with_mates avgLength size_MB
#    AssemblyName download_path Experiment LibraryName LibraryStrategy
#    LibrarySelection LibrarySource LibraryLayout InsertSize InsertDev Platform
#    Model SRAStudy BioProject Study_Pubmed_id ProjectID Sample BioSample
#    SampleType TaxID ScientificName SampleName g1k_pop_code source
#    g1k_analysis_group Subject_ID Sex Disease Tumor Affection_Status
#    Analyte_Type Histological_Type Body_Site CenterName Submission
#    dbgap_study_accession Consent RunHash ReadHash


def datetimeParse(s):
    "parse a string that might be date or date and time"
    try:
        return datetime.strptime(s, '%Y-%m-%d %H:%M:%S')
    except ValueError:
        return datetime.strptime(s, '%Y-%m-%d')


datetimeType = (datetimeParse,
                lambda d: datetime.strftime(d, '%Y-%m-%d %H:%M:%S'))
yesNoType = (lambda s: True if s == "yes" else False,
             lambda v: "yes" if v else "no")
runInfoTypeMap = {
    "ReleaseDate": datetimeType,
    "LoadDate": datetimeType,
    "spots": int,
    "bases": int,
    "spots_with_mates": int,
    "avgLength": int,
    "size_MB": int,
    "InsertSize": int,
    "InsertDev": int,
    "Tumor": yesNoType,
}


def runInfoReader(runInfoTsv):
    """Return a TsvReader for runinfo, Can be CVS or TSV, but not piped"""
    with fileOps.opengz(runInfoTsv) as fh:
        dialect = csv.Sniffer().sniff(fh.readline())
        fh.seek(0)
    return TsvReader(runInfoTsv, typeMap=runInfoTypeMap,
                     defaultColType=strOrNoneType,
                     dialect=dialect)


# assembly constants
GRCh38 = "GRCh38"
GRCm38 = "GRCm37"

assemblyNameMap = {
    "grch38": GRCh38,
    "grcm38": GRCm38,
    "GCA_000001635": GRCm38,
    "GCF_000001635": GRCm38,
}


def normalizeAssemblyName(assembly):
    """Try heuristics to get one of the supported assemblies, return None if None, error if unknown"""
    if assembly is None:
        return None
    assem0 = assembly.split(',')
    name = assemblyNameMap.get(assem0)
    if name is None:
        name = assemblyNameMap.get(assembly.lower())
    if name is None:
        raise Exception("don't know how normalize assembly: {}".format(assembly))
    return name
