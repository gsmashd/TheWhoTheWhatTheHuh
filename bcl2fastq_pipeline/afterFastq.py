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
import yaml
import json
import re
'''
Do we really need the md5sum?
'''

localConfig = None

SEQUENCERS = {
        'NB501038' : 'NextSeq 500',
        'SN7001334' : 'HiSeq 2500',
        'K00251' : 'HiSeq 4000',
        'M02675' : 'MiSeq NTNU',
        'M03942' : 'MiSeq StOlav'
        }


def bgzip_worker(fname) :
    global localConfig
    config = localConfig
    cmd = "%s -r %s" % (
        config.get("bgzip","bgzip_command"),
        fname)
    syslog.syslog("[bgzip_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def clumpify_worker(fname):
    global localConfig
    config = localConfig

    if config.get("Options","singleCell") == "1":
        return

    #We use R1 files to construct R2 filenames for paired end
    if "R2.fastq.gz" in fname:
        return

    if os.path.exists(os.path.join(os.path.dirname(fname),"clumpify.done")):
        return

    r1 = fname
    r2 = fname.replace("R1.fastq.gz","R2.fastq.gz") if os.path.exists(fname.replace("R1.fastq.gz","R2.fastq.gz")) else None

    if r2:
        out_r1 = r1.replace(".fastq.gz","_clumped.fastq.gz") 
        out_r2 = r2.replace(".fastq.gz","_clumped.fastq.gz")
        cmd = "{clump_cmd} {clump_opts} in1={in1} in2={in2} out1={out1} out2={out2} rcomp=f rename=f overwrite=true".format(
                clump_cmd = config.get("clumpify","clumpify_cmd"),
                clump_opts = config.get("clumpify", "clumpify_opts"),
                in1 = r1,
                in2 = r2,
                out1 = out_r1, #yes, we want to overwrite the files
                out2 = out_r2
                )
    else:
        out_r1 = r1.replace(".fastq.gz","_clumped.fastq.gz") 
        cmd = "{clump_cmd} {clump_opts} in={in1} out={out1} rename=f overwrite=true".format(
                clump_cmd = config.get("clumpify","clumpify_cmd"),
                clump_opts = config.get("clumpify", "clumpify_opts"),
                in1 = r1,
                out1 = out_r1 
                )
    syslog.syslog("[clumpify_worker] Processing %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)
    os.rename(out_r1,r1)
    if r2:
        os.rename(out_r2,r2)

def decontaminate_worker(fname):
    global localConfig
    config = localConfig

    #We use R1 files to construct R2 filenames for paired end
    if "R2.fastq.gz" in fname:
        return

    if os.path.exists(fname.replace(os.path.basename(fname),"contaminated/{}".format(os.path.basename(fname)))):
        return

    r2 = fname.replace("R1.fastq.gz","R2.fastq.gz") if os.path.exists(fname.replace("R1.fastq.gz","R2.fastq.gz")) else None

    masked_path = config.get("MaskedGenomes","HGDir")
    if not masked_path:
        raise Exception("HGDir not set in config bcl2fastq.ini!\n")

    if not os.path.exists(masked_path):
        raise Exception("HGDir {} does not exist!\n".format(masked_path))
    
    os.makedirs("{}/contaminated/".format(os.path.dirname(fname)),exist_ok=True)

    if r2:
        cont_out = fname.replace(os.path.basename(fname),"contaminated/{}".format(os.path.basename(fname).replace("R1.fastq.gz","interleaved.fastq.gz")))
        cmd = "{bbmap_cmd} {bbmap_opts} path={masked_path} in={infile} in2={infile2} outu={clean_out} outm={contaminated_out}".format(
                bbmap_cmd = config.get("MaskedGenomes","bbmap_cmd"),
                bbmap_opts = config.get("MaskedGenomes","bbmap_opts"),
                masked_path = config.get("MaskedGenomes","HGDir"),
                infile = fname,
                infile2 = r2,
                clean_out = fname.replace("R1.fastq.gz","interleaved.fastq.gz"),
                contaminated_out = cont_out 
                )
    else:   
        cmd = "{bbmap_cmd} {bbmap_opts} path={masked_path} in={infile} outu={clean_out} outm={contaminated_out}".format(
                bbmap_cmd = config.get("MaskedGenomes","bbmap_cmd"),
                bbmap_opts = config.get("MaskedGenomes","bbmap_opts"),
                masked_path = config.get("MaskedGenomes","HGDir"),
                infile = fname,
                clean_out = fname,
                contaminated_out = fname.replace(os.path.basename(fname),"contaminated/{}".format(os.path.basename(fname)))
                )

    syslog.syslog("[decontaminate_worker] Processing %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)
    # clean up
    os.remove(fname)
    if r2:
        os.remove(r2)
        #If paired end, split interleaved files to R1 and R2
        cmd = "{rename_cmd} {rename_opts} in={interleaved} out1={out_r1} out2={out_r2}".format(
                rename_cmd = "rename.sh",
                rename_opts = "renamebymapping=t",
                interleaved = fname.replace("R1.fastq.gz","interleaved.fastq.gz"),
                out_r1 = fname,
                out_r2 = r2
                )
        syslog.syslog("[decontaminate_worker] De-interleaving %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove(fname.replace("R1.fastq.gz","interleaved.fastq.gz"))
        #De-interleave contaminated file
        cmd = "{rename_cmd} {rename_opts} in={interleaved} out1={out_r1} out2={out_r2}".format(
                rename_cmd = "rename.sh",
                rename_opts = "renamebymapping=t",
                interleaved = cont_out,
                out_r1 = cont_out.replace("interleaved.fastq.gz","R1.fastq.gz"),
                out_r2 = cont_out.replace("interleaved.fastq.gz","R2.fastq.gz")
                )
        syslog.syslog("[decontaminate_worker] De-interleaving %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove(cont_out)


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

def get_gcf_name(fname):
    for name in fname.split('/'):
        if re.match(r'GCF-[0-9]{4}-[0-9]{3,}',name):
            return name
    raise Exception('Unable to determine GCF project number for filename {}\n'.format(fname))



def FastQC_worker(fname) :
    global localConfig
    config = localConfig
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    #projectName = fname.split("/")[-3]
    projectName = get_gcf_name(fname)

    libName = fname.split("/")[-2] #The last directory
    cmd = "%s %s -o %s/%s%s/FASTQC_%s/ %s" % (
          config.get("FastQC","fastqc_command"),
          config.get("FastQC","fastqc_options"),
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          fname)

    # Skip if the output exists
    fastqc_fname = glob.glob("{}/{}{}/FASTQC_{}/*/{}".format(
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          os.path.basename(fname).replace(".fastq.gz","_fastqc.zip")
          ))

    fastqc_fname.extend(glob.glob("{}/{}{}/FASTQC_{}/{}".format(
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          os.path.basename(fname).replace(".fastq.gz","_fastqc.zip")
          )))

    if fastqc_fname:
        return
    
    """
    if os.path.exists("%s/%s%s/FASTQC_%s/%s/%s_fastqc.zip" % (config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          libName,
          os.path.basename(fname)[:-9])):
        return
    """

    os.makedirs("%s/%s%s/FASTQC_%s/" % (config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName), exist_ok=True)

    syslog.syslog("[FastQC_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def toDirs(files) :
    s = set()
    for f in files :
        d = os.path.dirname(f)
        if d.split('/')[-1].startswith("GCF-"):
            s.add(d)
        else:
            s.add(d[:d.rfind('/')]) 
    return s


def get_read_geometry(run_dir):
    stats_file = open('{}/Stats/Stats.json'.format(run_dir),'r')
    stats_json = json.load(stats_file)
    lane_info = stats_json['ReadInfosForLanes'][0].get('ReadInfos', None)
    if not lane_info:
        return 'Read geometry could not be automatically determined.'
    R1 = None
    R2 = None
    for read in lane_info:
        if read['IsIndexedRead'] == True:
            continue
        elif read['Number'] == 1:
            R1 = read['NumCycles']
        elif read['Number'] == 2:
            R2 = read['NumCycles']
    if R1 and R2:
        return 'Paired end - forward read length (R1): {}, reverse read length (R2): {}'.format(R1,R2)
    elif R1 and not R2:
        return 'Single end - read length (R1): {}'.format(R1)
    elif not R1 and not R2:
        return 'Read geometry could not be automatically determined.'

def get_sequencer(run_id):
	return SEQUENCERS.get(run_id.split('_')[1],'Sequencer could not be automatically determined.')

def md5sum_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    if os.path.exists("md5sums.txt"):
        return
    if glob.glob("*.fastq.gz"):
        cmd = "md5sum *.fastq.gz > md5sums.txt"
    else:
        cmd = "md5sum */*.fastq.gz > md5sums.txt"
    syslog.syslog("[md5sum_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def multiqc_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    #os.chdir("{}/{}".format(config.get('Paths','outputDir'), config.get('Options','runID')))
    dname = d.split("/")
    pname = dname[-1]

    conf_name = ".{}_multiqc_config.yaml".format(pname)
    in_conf = open("/root/multiqc_config.yaml","r")
    out_conf = open(conf_name,"w+")
    mqc_conf = yaml.load(in_conf)

    mqc_conf['title'] = pname

    read_geometry = get_read_geometry(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
    contact = config.get('MultiQC','report_contact')
    sequencer = get_sequencer(config.get('Options','runID'))

    report_header = [
    {'Contact E-mail': contact},
    {'Sequencing Platform': sequencer},
    {'Read Geometry': read_geometry}
    ]
    
    mqc_conf['report_header_info'] = report_header
    
    yaml.dump(mqc_conf,out_conf)
    in_conf.close()
    out_conf.close()

    

    cmd = "{multiqc_cmd} {multiqc_opts} --config {conf} {flow_dir}/FASTQC_{pname} {flow_dir}/{pname}".format(
            multiqc_cmd = config.get("MultiQC", "multiqc_command"), 
            multiqc_opts = config.get("MultiQC", "multiqc_options"), 
            conf = conf_name,
            flow_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
            pname=pname,  
            )
    syslog.syslog("[multiqc_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)
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


def clumpify_mark_done(project_dirs):
    for d in project_dirs:
        open(os.path.join(d,"clumpify.done"),"w+").close()


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
    sampleFiles.extend(glob.glob("%s/%s%s/*/*R[12]_001.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes)))
    sampleFiles.extend(glob.glob("%s/%s%s/*/*/*R[12]_001.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes)))

    sampleFiles = [sf for sf in sampleFiles if not 'contaminated' in sf]
    sampleFiles = [sf for sf in sampleFiles if not 'filtered' in sf]

    global localConfig
    localConfig = config

    #suprDUPr
    """
    SKIP SUPRDUPR FOR NOW - NEED ./filterfq to work
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(suprDUPr_worker, sampleFiles)
    p.close()
    p.join()
    """
    #Decontaminate with masked genome
    if config.get("Options","decontaminate") == "1":
        p = mp.Pool(int(1))
        p.map(decontaminate_worker, sampleFiles)
        p.close()
        p.join()
    
    #clumpify

    p = mp.Pool(int(config.get("Options","clumpifyWorkerThreads")))
    p.map(clumpify_worker, sampleFiles)
    p.close()
    p.join()
    clumpify_mark_done(projectDirs)

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
