"""functions for running commands and pipelines"""
import sys
from pycbio.sys import pipeline

def runCmd(cmds, stdin="/dev/null", stdout=None, verbose=False):
    "run a command or pipeline, capturing stderr in exceptions"
    pl = pipeline.Procline(cmds, stdin=stdin, stdout=stdout, stderr=pipeline.DataReader)
    if verbose:
        sys.stderr.write(str(pl)+"\n")
    pl.wait()

def callCmd(cmds, stdin="/dev/null", verbose=False):
    "run a command or pipeline returning output and capturing stderr in exceptions"
    out = pipeline.DataReader()
    runCmd(cmds, stdin=stdin, stdout=out, verbose=verbose)
    return out.get()

def openPipeline(cmds, mode='r', otherEnd=None, verbose=False):
    """open up a pipeline"""
    pl = pipeline.Pipeline(cmds, mode, otherEnd)
    if verbose:
        sys.stderr.write(str(pl)+"\n")
    return pl

