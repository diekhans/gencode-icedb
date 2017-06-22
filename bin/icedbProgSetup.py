"""
Import this to do standard library and executable path setup for python programs.
This is not a program, it is in the bin directory so that is will be
automatically on the path.
"""
import sys
import os
rootDir = os.path.dirname(os.path.abspath(os.path.dirname(sys.argv[0])))


def _addExtern(module, relDir):
    modDir = os.path.join(rootDir, "extern", module, relDir)
    if not os.path.exists(modDir):
        raise Exception("can't find {} directory: {}".format(module, modDir))
    sys.path.insert(0, modDir)


_addExtern("pycbio", "lib")
_addExtern("pipettor", "build/lib")
_addExtern("ccds2", "output/lib/py")
sys.path.insert(0, os.path.join(rootDir, "lib"))

# executable PATH, make sure we can override kent commands and
# we get ~markd commands rather than installed ones.
os.environ["PATH"] = ":".join([os.path.join(rootDir, "bin"),
                               os.path.expanduser("~markd/opt/current/x86_64/bin"),
                               os.path.expanduser("~/kent/bin/x86_64"),
                               "/cluster/bin/x86_64",
                               os.environ["PATH"]])
