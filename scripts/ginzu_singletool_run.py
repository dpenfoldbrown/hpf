#!/usr/bin/python 

# Script to execute the DDB command for rerunning on single tools. Modify as needed.
# Known to work for rerunning DDB (/Ginzu) elements: SignalP, Disopred
# Add tool names and methods below if desired
# dpb 08/28/2012

import argparse
import subprocess


def signalp(sequence_key, ginzu_version):
    cmd = list()
    cmd.append("ddb.pl")
    cmd.append("-mode", "sequence")
    cmd.append("-submode", "signalp")
    cmd.append("-sequence_key", str(args.sequence_key))
    cmd.append("-organismtype", "euk")
    cmd.append("-ginzu_version", str(args.ginzu_version))
    
    print "Executing command: {0}".format(" ".join(cmd))
    subprocess.check_call(cmd, shell=False)


def disopred(sequence_key, ginzu_version):
    cmd = list()
    cmd.append("ddb.pl")
    cmd.append("-mode", "sequence")
    cmd.append("-submode", "disopred")
    cmd.append("-sequence_key", str(args.sequence_key))
    cmd.append("-ginzu_version", str(args.ginzu_version))

    print "Executing command: {0}".format(" ".join(cmd))
    subprocess.check_call(cmd, )


def test(*args):
    print "Test tool method\nArgs:"
    for a in args:
        print a


functions = {
    'signalp': signalp,
    'disopred': disopred,
    'test': test,
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A tool for running supported single-tool DDB/Ginzu pipeline pieces")
    parser.add_argument("-t", "--tool", action="store", dest="tool", required=True, choices=["signalp", "disopred", "test"],
            help="The tool to run from the Ginzu pipeline")
    parser.add_argument("-m", "--mode", action="store", dest="mode", choices=["single", "batch"], default="single", 
            help="Program mode. If single, -s must be a sequence key. If batch, -s is a line-delim'ed file of sequence keys")
    parser.add_argument("-g", "--ginzu_version", action="store", dest="ginzu_version", type=int, required=True, 
            help="The ginzu version of the tool being run")
    parser.add_argument("-s", "--sequence", action="store", dest="sequence", required=True, 
            help="See mode. Either a sequence key or a file containing sequence keys")
    args = parser.parse_args()
   
    if args.mode == "single": 
        # Run the method corresponding to the tool name
        functions[args.tool](int(args.sequence), args.ginzu_version)
    else:
        with open(args.sequence) as handle:
            for line in handle:
                seqkey = int(line.rstrip())
                functions[args.tool](seqkey, args.ginzu_version)
    print "Completed {0} {1} run".format(args.mode, args.tool)

