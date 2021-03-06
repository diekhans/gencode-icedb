"""
Import this to do standard library and executable path setup for python programs.
This is not a program, it is in the bin directory so that is will be
automatically on the path.
"""
import sys
import os

# Add our library to the path.  Using __file__ allows programs to be symlinked.
rootDir = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
binDir = os.path.join(rootDir, "bin")


def _addExtern(module, relDir):
    "add packages in extern directory to path"
    modDir = os.path.join(rootDir, "extern", module, relDir)
    if not os.path.exists(modDir):
        raise Exception("can't find {} directory: {}".format(module, modDir))
    sys.path.insert(0, modDir)


_addExtern("pycbio", "lib")
sys.path.insert(0, os.path.join(rootDir, "lib"))

# executable PATH, make sure we can override kent commands and
# we get ~markd commands rather than installed ones.
os.environ["PATH"] = ":".join([binDir, os.environ["PATH"]])
