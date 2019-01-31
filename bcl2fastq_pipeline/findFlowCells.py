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
from shutil import copyfile
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


def parseSampleSheet(ss):
    """
    Return a dictionary with keys: (Barcode length 1, Barcode length 2)

    return ss, laneOut, bcLens
    """

    f = open(ss)
    opt_d = None
    opts_data = False
    for line in f:
        if line.startswith("[Data]"):
            return ss,opt_d
        elif line.startswith("[CustomOptions]"):
            opt_d = dict.fromkeys(CUSTOM_OPTS,False)
            opts_data = True
            continue
        elif opts_data:
            key = line.split(',')[0]
            value = line.split(',')[-1]
            opt_d[key] = value.rstrip() if key in ['Organism','Libprep'] else str2bool(value)
    return ss,opt_d

def str2bool(s):
    return s.lower() in ['true','1']


def getSampleSheets(d):
    """
    Provide a list of output directories and sample sheets
    """
    ss = glob.glob("%s/SampleSheet*.csv" % d)

    if len(ss) == 0:
        return ([None],None)

    for sheet in ss:
        ss_, opts = parseSampleSheet(sheet)
        if ss and opts:
            return ss, opts
    return None, None


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

        ss, opts = getSampleSheets(os.path.dirname(d))

        if not opts:
            continue

        syslog.syslog("Found a new flow cell: %s\n" % config.get("Options","runID"))
        odir = "{}/{}{}".format(config.get("Paths", "outputDir"), config.get("Options", "runID"), lanesUse)
        if not os.path.exists(odir):
            os.makedirs(odir)
        if ss is not None and not os.path.exists("{}/SampleSheet.csv".format(odir)):
            copyfile(ss,"{}/SampleSheet.csv".format(odir))
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
        opts = dict.fromkeys(CUSTOM_OPTS,False)
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
