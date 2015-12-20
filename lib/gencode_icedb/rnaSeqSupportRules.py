"""
Toil rules for building RNA-Seq support
"""

def starGenerateGenome(refGeneomeName, refGeneomeFa, readLength, geneAnnotationSetName
                       geneAnnotationSetGtf, numThreads):
    """Generate STAR genome given reference genome, read, and annotation set.
    Create a flag file to indicate it is complete.
    """
    cmd = ["starGenerateGenome", "--numThreads="+str(numThreads), "refGeneomeFa",
           geneAnnotationSetGtf, readLength, genomeDir]


"""
"""
