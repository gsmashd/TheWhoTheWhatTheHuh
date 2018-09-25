'''
This file includes code that actually runs FastQC and any other tools after the fastq files have actually been made. This uses a pool of workers to process each request.
'''
import multiprocessing as mp
import glob
import sys
import subprocess
import os
import os.path
import shlex
import shutil
import xml.etree.ElementTree as ET
import syslog
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as img

'''
Do we really need the md5sum?
'''

localConfig = None

def bgzip_worker(fname) :
    global localConfig
    config = localConfig
    cmd = "%s -r %s" % (
        config.get("bgzip","bgzip_command"),
        fname)
    syslog.syslog("[bgzip_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)


def fastq_screen_worker(fname) :
    global localConfig
    config = localConfig

    os.chdir(os.path.dirname(fname))

    #Skip read 2 when not single cell, skip if read 1 when single cell
    if config.get("Options","singleCell") == '0' and ("R2.fastq" in fname or "R2_001.fastq" in fname):
        return
    elif config.get("Options","singleCell") == '1' and ("R1.fastq" in fname or "R1_001.fastq" in fname):
        return

    ofile="{}/fastq_screen/{}".format(
            os.path.dirname(fname),
            os.path.basename(fname)
            )

    if os.path.exists(ofile.replace(".fastq.gz","_screen.html")):
        return
    
    os.makedirs(os.path.dirname(ofile), exist_ok=True)
    
    #fastq_screen
    cmd = "%s %s --outdir '%s' '%s'" % (
        config.get("fastq_screen", "fastq_screen_command"),
        config.get("fastq_screen", "fastq_screen_options"),
        os.path.dirname(ofile),
        fname)
    syslog.syslog("[fastq_screen_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

    #Unlink/rename
    #os.unlink(ofile)
    #os.rename(ofile.replace(".fastq","_screen.png"), fname.replace("_R1_001.fastq.gz", "_R1_001_screen.png"))


def suprDUPr_worker(fname) :
    global localConfig
    config = localConfig

    if config.get("Options","singleCell") == "1":
        return

    ofname = "{}/filtered/{}".format(os.path.dirname(fname),os.path.basename(fname).replace(".fastq.gz","_filtered.fastq.gz"))

    # WHEN ./filterfq is working cmd will be 
    #  "{cmd} {opts} {infile} | filterfq {infile} | tee {ofile}"
    cmd = "{cmd} {opts} {infile} > {ofile}".format(
          cmd = config.get("suprDUPr","suprdupr_command"),
          opts = config.get("suprDUPr","suprdupr_options"),
          infile = fname,
          ofile = ofname
          )

    # Skip if the output exists
    if os.path.exists(ofname):
        return

    os.makedirs(os.path.dirname(ofname), exist_ok=True)

    syslog.syslog("[suprDUPr_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def FastQC_worker(fname) :
    global localConfig
    config = localConfig
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    projectName = fname.split("/")[-3] #It's the penultimate directory
    libName = fname.split("/")[-2] #The last directory
    cmd = "%s %s -o %s/%s%s/FASTQC_%s/%s %s" % (
          config.get("FastQC","fastqc_command"),
          config.get("FastQC","fastqc_options"),
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          libName,
          fname)

    # Skip if the output exists
    if os.path.exists("%s/%s%s/FASTQC_%s/%s/%s_fastqc.zip" % (config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          libName,
          os.path.basename(fname)[:-9])):
        return

    os.makedirs("%s/%s%s/FASTQC_%s/%s" % (config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          libName), exist_ok=True)

    syslog.syslog("[FastQC_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def toDirs(files) :
    s = set()
    for f in files :
        d = os.path.dirname(f)
        s.add(d[:d.rfind('/')]) #We just want projects, not individual libraries
    return s

def md5sum_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    if os.path.exists("mdsums.txt"):
        return
    cmd = "md5sum */*.fastq.gz > md5sums.txt"
    syslog.syslog("[md5sum_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def multiqc_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    dname = d.split("/")
    dname[-1] = "FASTQC_{}".format(dname[-1])
    dname = os.path.join(os.path.dirname(d),dname[-1])
    cmd = "{} {} {}/*/*.zip {}/*/*".format(config.get("MultiQC", "multiqc_command"), config.get("MultiQC", "multiqc_options"), dname, d)
    syslog.syslog("[multiqc_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def clumpify_worker(d):
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)

    read1s = glob.glob("*_R1_001.fastq.gz")
    PE = 1
    for r1 in read1s:
        # This takes a while, don't duplicate work
        if os.path.exists("{}.duplicate.txt".format(r1[:-12])):
            os.chdir(oldWd)
            return

        r2 = "{}_R2.fastq.gz".format(r1[:-12])
        if os.path.exists(r2):
            cmd = "{} in={} in2={} out=temp.fq.gz {} dupedist={} threads={}".format(config.get("bbmap", "clumpify_command"),
                                                                                    r1, r2,
                                                                                    config.get("bbmap", "clumpify_options"),
                                                                                    config.get("bbmap", "clumpify_HiSeq3000_dist"),
                                                                                    config.get("bbmap", "clumpify_threads"))
        else:
            PE = 0
            cmd = "{} in={} out=temp.fq.gz {} dupedist={} threads={}".format(config.get("bbmap", "clumpify_command"),
                                                                             r1,
                                                                             config.get("bbmap", "clumpify_options"),
                                                                             config.get("bbmap", "clumpify_HiSeq3000_dist"),
                                                                             config.get("bbmap", "clumpify_threads"))
        syslog.syslog("[clumpify_worker] Processing %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        cmd = " ".join(["splitFastq", "temp.fq.gz", "{}".format(PE), r1[:-12], "{}".format(config.get("bbmap", "pigzThreads"))])
        syslog.syslog("[clumpify_worker] Splitting %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove("temp.fq.gz")
    os.chdir(oldWd)

def clumpifyNextSeq_worker(d):
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)

    read1s = glob.glob("*_R1_001.fastq.gz")
    PE = 1
    for r1 in read1s:
        # This takes a while, don't duplicate work
        if os.path.exists("{}.duplicate.txt".format(r1[:-12])):
            os.chdir(oldWd)
            return

        r2 = "{}_R2.fastq.gz".format(r1[:-12])
        if os.path.exists(r2):
            cmd = "{} in={} in2={} out=temp.fq.gz {} {} dupedist={} threads={}".format(config.get("bbmap", "clumpify_command"),
                                                                                       r1, r2,
                                                                                       config.get("bbmap", "clumpify_options"),
                                                                                       config.get("bbmap", "clumpify_NextSeq_options"),
                                                                                       config.get("bbmap", "clumpify_NextSeq_dist"),
                                                                                       config.get("bbmap", "clumpify_threads"))
        else:
            PE = 0
            cmd = "{} in={} out=temp.fq.gz {} {} threads={}".format(config.get("bbmap", "clumpify_command"),
                                                                                r1,
                                                                                config.get("bbmap", "clumpify_options"),
                                                                                config.get("bbmap", "clumpify_NextSeq_options"),
                                                                                config.get("bbmap", "clumpify_threads"))
        syslog.syslog("[clumpify_worker] Processing %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        cmd = " ".join(["splitFastq", "temp.fq.gz", "{}".format(PE), r1[:-16], "{}".format(config.get("bbmap", "pigzThreads"))])
        syslog.syslog("[clumpify_worker] Splitting %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove("temp.fq.gz")
    os.chdir(oldWd)

def parserDemultiplexStats(config) :
    '''
    Parse DemultiplexingStats.xml under outputDir/Stats/ to get the
    number/percent of undetermined indices.

    In particular, we extract the BarcodeCount values from Project "default"
    Sample "all" and Project "all" Sample "all", as the former gives the total
    undetermined and the later simply the total clusters
    '''
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    totals = [0,0,0,0,0,0,0,0]
    undetermined = [0,0,0,0,0,0,0,0]
    tree = ET.parse("%s/%s%s/Stats/DemultiplexingStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
    root = tree.getroot()
    for child in root[0].findall("Project") :
        if(child.get("name") == "default") :
            break
    for sample in child.findall("Sample") :
        if(sample.get("name") == "all") :
            break
    child = sample[0] #Get inside Barcode
    for lane in child.findall("Lane") :
        lnum = int(lane.get("number"))
        undetermined[lnum-1] += int(lane[0].text)

    for child in root[0].findall("Project") :
        if(child.get("name") == "all") :
            break
    for sample in child.findall("Sample") :
        if(sample.get("name") == "all") :
            break
    child = sample[0] #Get Inside Barcode
    for lane in child.findall("Lane") :
        lnum = int(lane.get("number"))
        totals[lnum-1] += int(lane[0].text)

    out = ""
    for i in range(8) :
        if(totals[i] == 0) :
            continue
        out += "\nLane %i: %i of %i reads/pairs had undetermined indices (%5.2f%%)" % (
            i+1,undetermined[i],totals[i],100*undetermined[i]/totals[i])
    return out

#All steps that should be run after `make` go here
def postMakeSteps(config) :
    '''
    Current steps are:
      1) Run FastQC on each fastq.gz file
      2) Run md5sum on the files in each project directory
    Other steps could easily be added to follow those. Note that this function
    will try to use a pool of threads. The size of the pool is set by config.postMakeThreads
    '''
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    projectDirs = glob.glob("%s/%s%s/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes))
    projectDirs.extend(glob.glob("%s/%s%s/*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes)))
    projectDirs = toDirs(projectDirs)
    sampleFiles = glob.glob("%s/%s%s/*/*R[12].fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
    sampleFiles.extend(glob.glob("%s/%s%s/*/*/*R[12].fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes)))
    sampleFiles.extend(glob.glob("%s/%s%s/*/*/*R[12]_001.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes)))

    global localConfig
    localConfig = config

    # skip clumpify
    """
    #Deduplicate if this is a HiSeq 3000 run
    if config.get("Options", "runID")[7] == "J":
        sampleDirs = glob.glob("%s/%s%s/Project_*/*_R1.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
        sampleDirs = [os.path.dirname(x) for x in sampleDirs]
        
        p = mp.Pool(int(config.get("Options", "deduplicateInstances")))
        p.map(clumpify_worker, sampleDirs)
        p.close()
        p.join()
        
    #Different deduplication for NextSeq samples
    elif config.get("Options", "runID")[7:9] == "NB":
        sampleDirs = glob.glob("%s/%s%s/Project_*/*_R1*.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
        #print("%s/%s%s/Project_*/*_R1*.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
        #sampleDirs = [os.path.dirname(x) for x in sampleDirs]
        sampleDirs = set([os.path.dirname(x) for x in sampleDirs])
        #print("SampleDirs")
        #print(sampleDirs)
        
        p = mp.Pool(int(config.get("Options", "deduplicateInstances")))
        p.map(clumpifyNextSeq_worker, sampleDirs)
        p.close()
        p.join()
        """

    # Avoid running post-processing (in case of a previous error) on optical duplicate files.
    #sampleFiles = [x for x in sampleFiles if "optical_duplicates" not in x]


    #suprDUPr

    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(suprDUPr_worker, sampleFiles)
    p.close()
    p.join()

    #FastQC

    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(FastQC_worker, sampleFiles)
    p.close()
    p.join()

    #md5sum
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(md5sum_worker, projectDirs)
    p.close()
    p.join()

    #fastq_screen
    p = mp.Pool(int(config.get("Options", "postMakeThreads")))
    p.map(fastq_screen_worker, sampleFiles)
    p.close()
    p.join()

    # multiqc
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(multiqc_worker, projectDirs)
    p.close()
    p.join()

    #disk usage
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024

    #Undetermined indices
    undeter = parserDemultiplexStats(config)

    message = "Current free space: %i of %i gigs (%5.2f%%)\n" % (
        free,tot,100*free/tot)
    message += undeter
    return(message)
