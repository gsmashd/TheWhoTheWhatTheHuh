'''
This file includes anything involved in finding new flow cells to process.

Note that this includes anything done after a flow cell has been processed,
such as marking it as having been processed and sending emails.
'''

import os
import sys
import smtplib
import glob
from email.mime.text import MIMEText
import syslog
import xml.etree.ElementTree as ET
import flowcell_manager.flowcell_manager as fm


CUSTOM_OPTS = ['Organism', 'Libprep', 'SingleCell','RemoveHumanReads','SensitiveData','ReverseComplementIndexP5','ReverseComplementIndexP7','TrimAdapter']


#Returns True on processed, False on unprocessed
def flowCellProcessed(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    flowcells = fm.list_flowcell_all(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID")))
    path = "%s/%s%s/fastq.made" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes)
    if os.access(path, os.F_OK):
        return True
    elif not flowcells.empty:
        return True
    return False


# Get the number of lanes in the run. This might not match the number of lanes in the sampleSheet
def getNumLanes(d):
    try:
        tree = ET.parse("{}/RunInfo.xml".format(d))
        root = tree.getroot()[0]
        numLanes = root.findall("FlowcellLayout")[0]
        return int(numLanes.get("LaneCount"))
    except:
        return 1



def revComp(s):
    """Reverse complement a primer sequence"""
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    s = s[::-1]
    x = [d[c] for c in s]
    return "".join(x)


def formatHeaderLine(cols, colLabs, indexCols, storeLanes):
    """Format a single header line in a sample sheet"""
    l = []
    if storeLanes is True:
        l.append("Lane")
    if colLabs[1] is not None:
        l.append("Sample_ID")
    if colLabs[2] is not None:
        l.append("Sample_Name")
    if indexCols[0] is not None and len(cols[indexCols[0]]) > 0:
        l.append("index")
    if indexCols[1] is not None and len(cols[indexCols[1]]) > 0:
        l.append("index2")
    if colLabs[3] is not None:
        l.append("Sample_Project")
    return ",".join(l)


def formatLine(cols, colLabs, indexCols, storeLanes, rcI5=False):
    """Format a single line in a sample sheet"""
    l = []
    if storeLanes is True and colLabs[0] is not None:
        l.append(cols[colLabs[0]])
    if colLabs[1] is not None:
        l.append(cols[colLabs[1]])
    if colLabs[2] is not None:
        l.append(cols[colLabs[2]])
    if indexCols[0] is not None and len(cols[indexCols[0]]) > 0:
        l.append(cols[indexCols[0]])
    if indexCols[1] is not None and len(cols[indexCols[1]]) > 0:
        if rcI5 is True:
            l.append(revComp(cols[indexCols[1]]))
        else:
            l.append(cols[indexCols[1]])
    if colLabs[3] is not None:
        l.append(cols[colLabs[3]])
    return ",".join(l)


def reformatSS(rv,opts):
    """This is used in parseSampleSheet to reformat the output for upstream use"""
    ss = []
    laneOut = []
    bcLens = []
    nLanes = 0

    for k, v in rv.items():
        ss.append("\n".join(v[0]))
        lanes = ""
        if len(v[1]) > 0:
            nLanes += len(v[1])
            lanes = "_".join(["{}".format(x) for x in sorted(list(v[1]))])
        laneOut.append(lanes)
        bcLens.append(v[2])

    if len(ss) < 2 and nLanes == 8:
        laneOut = None
    return ss, laneOut, bcLens, opts

def parseSampleSheet(ss):
    """
    Return a dictionary with keys: (Barcode length 1, Barcode length 2)

    return ss, laneOut, bcLens
    """
    rv = dict()
    # If this is a NextSeq or a HiSeq 2500 rapid run, then don't store the incorrect Lane column
    storeLanes = True
    if getNumLanes(os.path.dirname(ss)) < 8:
        storeLanes = False

    # For NextSeq/MiSeq, reverse complement I5 sequences
    #rcI5 = mustRevComp(os.path.dirname(ss))

    f = open(ss)
    inData = False
    lastLine = None
    colLabs = [None, None, None, None] # Lane, Sample_ID, Sample_Name, Sample_Project
    indexCols = [None, None] # index, index2
    #custom_opts = ['Organism', 'Libprep', 'SingleCell','RemoveHumanReads','SensitiveData']
    opt_d = dict.fromkeys(CUSTOM_OPTS,False)
    opts_data = False
    for line in f:
        bcLen = '0,0'
        if line.startswith("[Data]"):
            if opts_data:
                opts_data=False
                inData = True
                continue
            else:
                return reformatSS(rv,opt_d)
        elif line.startswith("[CustomOptions]"):
            opts_data = True
            continue
        elif opts_data:
            key = line.split(',')[0]
            value = line.split(',')[-1]
            opt_d[key] = value.rstrip() if key in ['Organism','Libprep'] else str2bool(value)
        elif inData:
            cols = line.strip().split(",")
            if lastLine is True:
                if indexCols[0] is not None:
                    bcLen = "{}".format(len(cols[indexCols[0]]))
                    if indexCols[1] is not None:
                        bcLen = "{},{}".format(bcLen, len(cols[indexCols[1]]))
                    else:
                        bcLen += ",0"

                # Append to rv, with header
                if bcLen not in rv:
                    rv[bcLen] = [["[Data]", formatHeaderLine(cols, colLabs, indexCols, storeLanes)], set(), bcLen]

                # Add the lane to the set, if relevant
                if colLabs[0] is not None and storeLanes is True:
                    rv[bcLen][1].add(int(cols[colLabs[0]]))

                rv[bcLen][0].append(formatLine(cols, colLabs, indexCols, storeLanes))

            # Set columns for barcodes, etc.
            if lastLine is None:
                lastLine = True
                if "index" in cols:
                    indexCols[0] = cols.index("index")
                    if "index2" in cols:
                        indexCols[1] = cols.index("index2")

                if "Lane" in cols:
                    colLabs[0] = cols.index("Lane")

                if "Sample_ID" in cols:
                    colLabs[1] = cols.index("Sample_ID")

                if "Sample_Name" in cols:
                    colLabs[2] = cols.index("Sample_Name")

                if "Sample_Project" in cols:
                    colLabs[3] = cols.index("Sample_Project")
                continue
        elif line.startswith("[CustomOptions]"):
            opts_data = True
            continue
        elif line.startswith("[Data]"):
            inData = True
            continue
    #rv.update(opts_d)
    return reformatSS(rv,opt_d)

def str2bool(s):
    return s.lower() in ['true','1']


def getSampleSheets(d):
    """
    Provide a list of output directories and sample sheets
    """
    ss = glob.glob("%s/SampleSheet*.csv" % d)

    if len(ss) == 0:
        return ([None], [None], [''],None)

    laneOut = []
    bcLens = []
    ssUse = []
    for sheet in ss:
        ss_, laneOut_, bcLens_, opts = parseSampleSheet(sheet)
        nSS = 0
        if ss_ is not None and len(ss_) > 0:
            ssUse.extend(ss_)
            nSS = len(ss_)
        if nSS > 0 and laneOut_ is not None and len(laneOut_) > 0:
            laneOut.extend(laneOut_)
        elif nSS > 0:
            laneOut.extend([None] * nSS)
        if nSS > 0 and bcLens_ is not None and len(bcLens_) > 0:
            bcLens.extend(bcLens_)
        elif nSS > 0:
            bcLens.extend([None] * nSS)

    return ssUse, laneOut, bcLens, opts


'''
Iterate over all folders in config.baseDir from machine SN7001180. For each,
see if it contains an RTAComplete.txt file, as this signifies the completion
of a sequencing run.

If a run has finished, we then check to see if it's already been processed.
Runs are marked as having been processed if they appear in config.finalDir
and have a file called casava.finished or fastq.made. casava.finished is
produced by the old pipeline, while this one creates fastq.made.

This function always returns its configuration. If there's a new flow cell to
process, then the runID is filled in. Otherwise, that's set to None.
'''
def newFlowCell(config) :
    #HiSeq2500
    dirs = glob.glob("%s/*/data/*_SN7001334_*/ImageAnalysis_Netcopy_complete.txt" % config.get("Paths","baseDir"))
    #NextSeq 500
    dirs.extend(glob.glob("%s/*/data/*_NB501038_*/RunCompletionStatus.xml" % config.get("Paths","baseDir")))
    #MiSeq NTNU
    dirs.extend(glob.glob("%s/*/data/*_M026575*_*/ImageAnalysis_Netcopy_complete.txt" % config.get("Paths","baseDir")))
    #MiSeq St. Olav 
    dirs.extend(glob.glob("%s/*/data/*_M03942*_*/ImageAnalysis_Netcopy_complete.txt" % config.get("Paths","baseDir")))
    #HiSeq4000
    dirs.extend(glob.glob("%s/*/data/*_K00251*_*/SequencingComplete.txt" % config.get("Paths","baseDir")))
    for d in dirs :
        #Get the flow cell ID (e.g., 150416_SN7001180_0196_BC605HACXX)
        config.set('Options','runID',d.split("/")[-2])
        config.set('Options', 'sequencer',d.split("/")[-4])
        if flowCellProcessed(config):
            continue

        sampleSheet, lanes, bcLens, opts = getSampleSheets(os.path.dirname(d))

        if not opts:
            continue

        for ss, lane, bcLen in zip(sampleSheet, lanes, bcLens):
            config.set('Options','runID',d.split("/")[-2])
            config.set('Options', 'sequencer',d.split("/")[-4])
            lanesUse = ""
            if lane is not None and lane != "":
                config.set("Options","lanes",lane)
                lanesUse = "_lanes{}".format(lane)
            else:
                config.set("Options","lanes","")
            if ss is None:
                ss = ''
            if bcLen is not None and bcLen is not '':
                config.set("Options","bcLen",bcLen)
            else:
                config.set("Options","bcLen","0,0")

            if flowCellProcessed(config) is False:
                syslog.syslog("Found a new flow cell: %s\n" % config.get("Options","runID"))
                odir = "{}/{}{}".format(config.get("Paths", "outputDir"), config.get("Options", "runID"), lanesUse)
                if not os.path.exists(odir):
                    os.makedirs(odir)
                if ss is not None and not os.path.exists("{}/SampleSheet.csv".format(odir)):
                    o = open("{}/SampleSheet.csv".format(odir), "w")
                    for k,v in opts.items():
                        ss = '{},{}\n{}'.format(k,v,ss)
                    ss = '[CustomOptions]\n{}'.format(ss)
                    o.write(ss)
                    o.close()
                    ss = "{}/SampleSheet.csv".format(odir)
                config.set("Options","sampleSheet",ss)
                config = setConfFromOpts(config,opts)
                return config
            else :
                config.set("Options","runID","")
                config.set("Options","sequencer","")
                config = setConfFromOpts(config,opts,use_dict_values=False)
    config.set("Options","runID","")
    config.set("Options","sequencer","")
    config = setConfFromOpts(config,opts,use_dict_values=False)
    return config

def bool2strint(b):
    return '1' if b else '0'

def setConfFromOpts(config,opts,use_dict_values=True):
    if not opts:
        return config
    for k,v in opts.items():
        if use_dict_values:
            config.set("Options",k,v if v in ['Organism','Libprep'] else bool2strint(v))
        else:
            config.set("Options",k,"")
    return config


def markFinished(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    open("%s/%s%s/fastq.made" % (config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()

'''
This function needs to be run after newFlowCell() returns with config.runID
filled in. It creates the output directories.
'''
def MakeTargetDirs(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    assert(config["Paths"]["runID"] != None)
    os.mkdirs("%s/%s%s" % (config["Paths"]["outputDir"], config["Options"]["runID"], lanes))
