#!/usr/bin/env python3
import sys
import os
import datetime
import time
import syslog
import bcl2fastq_pipeline.getConfig
import bcl2fastq_pipeline.findFlowCells
import bcl2fastq_pipeline.makeFastq
import bcl2fastq_pipeline.afterFastq
import bcl2fastq_pipeline.misc
import importlib
import signal
from threading import Event
import urllib3

# Disable excess warning messages if we disable SSL checks
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
gotHUP = Event()

def breakSleep(signo, _frame):
    gotHUP.set()

def sleep(config) :
    gotHUP.wait(timeout=float(config['Options']['sleepTime'])*60*60)
    gotHUP.clear()

signal.signal(signal.SIGHUP, breakSleep)

while True:
    #Reimport to allow reloading a new version
    importlib.reload(bcl2fastq_pipeline.getConfig)
    importlib.reload(bcl2fastq_pipeline.findFlowCells)
    importlib.reload(bcl2fastq_pipeline.makeFastq)
    importlib.reload(bcl2fastq_pipeline.afterFastq)
    importlib.reload(bcl2fastq_pipeline.misc)

    #Read the config file
    config = bcl2fastq_pipeline.getConfig.getConfig()
    if(config is None) :
        #There's no recovering from this!
        sys.exit("Error: couldn't read the config file!")

    #Get the next flow cell to process, or sleep
    config = bcl2fastq_pipeline.findFlowCells.newFlowCell(config)
    if(config.get('Options','runID') == "") :
        sleep(config)
        continue

    #Ensure we have sufficient space
    if(bcl2fastq_pipeline.misc.enoughFreeSpace(config) == False) :
        syslog.syslog("Error: insufficient free space!\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Error: insufficient free space!")
        sleep(config)
        continue

    startTime=datetime.datetime.now()

    #Make the fastq files, if not already done
    lanes = config["Options"]["lanes"]
    if lanes != "":
        lanes = "_lanes{}".format(lanes)
    if not os.path.exists("{}/{}{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes)):
        try:
            bcl2fastq_pipeline.makeFastq.bcl2fq(config)
            open("{}/{}{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()
        except :
            syslog.syslog("Got an error in bcl2fq\n")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in bcl2fq")
            sleep(config)
            continue



    if not os.path.exists("{}/{}{}/files.renamed".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes)):
        try:
            bcl2fastq_pipeline.makeFastq.fixNames(config)
            open("{}/{}{}/files.renamed".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()
        except :
            syslog.syslog("Got an error in fixNames\n")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in fixNames")
            sleep(config)
            continue  

    
    #Run post-processing steps
    try :
        message = bcl2fastq_pipeline.afterFastq.postMakeSteps(config)
    except :
        syslog.syslog("Got an error during postMakeSteps\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during postMakeSteps")
        sleep(config)
        continue

    #Get more statistics and create PDFs
    try :
        message += "\n\n"+bcl2fastq_pipeline.misc.parseConversionStats(config)
    except :
        syslog.syslog("Got an error during parseConversionStats\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during parseConversionStats")
        sleep(config)
        continue

    runTime = datetime.datetime.now()-startTime

    #Email finished message
    try :
        bcl2fastq_pipeline.misc.finishedEmail(config, message, runTime)
    except :
        #Unrecoverable error
        syslog.syslog("Couldn't send the finished email! Quiting")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during finishedEmail()")
        sleep(config)

    #Zip project archives
    try:
        bcl2fastq_pipeline.afterFastq.archive_worker(config)
    except Exception as e:
        syslog.syslog("Got an error during zipping\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), str(e))
        sleep(config)
        continue

    #md5sum final archive
    try:
        bcl2fastq_pipeline.afterFastq.md5sum_archive_worker(config)
    except Exception as e:
        syslog.syslog("Got an error during md5sum of archive\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), str(e))
        sleep(config)
        continue

    #Mark the flow cell as having been processed
    bcl2fastq_pipeline.findFlowCells.markFinished(config)



















