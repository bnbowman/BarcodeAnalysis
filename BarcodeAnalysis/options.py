from __future__ import absolute_import

import argparse, h5py, os, sys
from BarcodeAnalysis._version import __version__

__author__ = "Brett Bowman"

options = argparse.Namespace()

def consensusCoreVersion():
    try:
        import ConsensusCore
        return ConsensusCore.Version.VersionString()
    except:
        return None

def parseOptions():
    """
    Parse the options and perform some due diligence on them
    """
    desc = "Compute the likelihood of sequences coming from a list of barcodes."
    parser = argparse.ArgumentParser(description=desc, add_help=False)

    def canonicalizedFilePath(path):
        return os.path.abspath(os.path.expanduser(path))

    def checkInputFile(path):
        if not os.path.isfile(path):
            parser.error("Input file %s not found." % (path,))

    def checkOutputFile(path):
        try:
            f = open(path, "a")
            f.close()
        except:
            parser.error("Output file %s cannot be written." % (path,))

    basics = parser.add_argument_group("Basic required options")
    basics.add_argument(
        "inputFilename",
        type=canonicalizedFilePath,
        metavar="BAS.H5",
        help="The filename of the input Bas.H5 file")
    basics.add_argument(
        "barcodeFilename",
        type=canonicalizedFilePath,
        metavar="BARCODE_FASTA",
        help="The filename of the barcode FASTA")

    parameter = parser.add_argument_group("Parameter settings")
    parameter.add_argument(
        "-P", "--parametersFile",
        dest="parametersFile",
        type=str,
        default=None,
        help="Parameter set filename (QuiverParameters.ini), or directory D " + \
             "such that either D/*/GenomicConsensus/QuiverParameters.ini, "   + \
             "or D/GenomicConsensus/QuiverParameters.ini, is found.  In the " + \
             "former case, the lexically largest path is chosen.")
    parameter.add_argument(
        "-p", "--parametersSpec",
        action="store",
        dest="parametersSpec",
        type=str,
        default="auto",
        help="Name of parameter set (chemistry.model) to select from the "   + \
             "parameters file, or just the name of the chemistry, in which " + \
             "case the best available model is chosen.  Default is 'auto', " + \
             "which selects the best parameter set from the cmp.h5")
    parameter.add_argument(
        "-o", "--outputFilename",
        type=str,
        default=os.path.join(os.getcwd(), "output.csv"),
        metavar="CSV",
        help="The filename of the CSV to output barcode scoring data to.")
    parameter.add_argument(
        "-a", "--appendOutput",
        action='store_true',
        help="The")

    debugging = parser.add_argument_group("Verbosity and debugging/profiling")
    debugging.add_argument("--help", "-h",
                           action="help")
    debugging.add_argument(
        "--nZmws",
        default=-1,
        type=int,
        help="Label only the first N ZMWs for testing purposes.")
    debugging.add_argument(
        "--verbose", "-v",
        dest="verbosity",
        action="count",
        help="Set the verbosity level.")
    debugging.add_argument(
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Turn off all logging, including warnings")
    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "  BarcodeAnalysis version: %s" % __version__
            print "  ConsensusCore version: %s" % consensusCoreVersion()
            print "  h5py version: %s" % h5py.version.version
            print "  hdf5 version: %s" % h5py.version.hdf5_version
            sys.exit(0)
    debugging.add_argument("--version",
                           nargs=0,
                           action=PrintVersionAction)

    parser.parse_args(namespace=options)

    for path in (options.inputFilename, options.barcodeFilename):
        if path is not None:
            checkInputFile(path)

    for path in (options.outputFilename,):
        if path is not None:
            checkOutputFile(path)

    options.shellCommand = " ".join(sys.argv)
