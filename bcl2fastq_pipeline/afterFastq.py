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
import datetime as dt
import flowcell_manager.flowcell_manager as fm

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

    if config.get("Options","SingleCell") == "1":
        return

    #We use R1 files to construct R2 filenames for paired end
    if "R2.fastq.gz" in fname:
        return

    if os.path.exists(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"clumpify.done")):
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

def RemoveHumanReads_worker(fname):
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

    syslog.syslog("[RemoveHumanReads_worker] Processing %s\n" % cmd)
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
        syslog.syslog("[RemoveHumanReads_worker] De-interleaving %s\n" % cmd)
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
        syslog.syslog("[RemoveHumanReads_worker] De-interleaving %s\n" % cmd)
        subprocess.check_call(cmd, shell=True)
        os.remove(cont_out)


def fastq_screen_worker(fname) :
    global localConfig
    config = localConfig

    os.chdir(os.path.dirname(fname))

    project_nr = get_gcf_name(fname)

    #Skip read 2 when not single cell, skip if read 1 when single cell
    if config.get("Options","SingleCell") == '0' and ("R2.fastq" in fname or "R2_001.fastq" in fname):
        return
    elif config.get("Options","SingleCell") == '1' and ("R1.fastq" in fname or "R1_001.fastq" in fname):
        return

    ofile="{}/{}/QC_{}/fastq_screen/{}".format(
            os.path.join(config.get("Paths","outputDir")),
            config.get("Options","runID"),
            project_nr,
            os.path.basename(fname)
            )

    if os.path.exists(ofile.replace(".fastq.gz","_screen.html")) :
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

    if config.get("Options","SingleCell") == "1":
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
    cmd = "%s %s -o %s/%s%s/QC_%s/FASTQC %s" % (
          config.get("FastQC","fastqc_command"),
          config.get("FastQC","fastqc_options"),
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          fname)

    # Skip if the output exists
    fastqc_fname = glob.glob("{}/{}{}/QC_{}/FASTQC/*/{}".format(
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          os.path.basename(fname).replace(".fastq.gz","_fastqc.zip")
          ))

    fastqc_fname.extend(glob.glob("{}/{}{}/QC_{}/FASTQC/{}".format(
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          lanes,
          projectName,
          os.path.basename(fname).replace(".fastq.gz","_fastqc.zip")
          )))

    if fastqc_fname:
        return

    os.makedirs("%s/%s%s/QC_%s/FASTQC" % (config.get("Paths","outputDir"),
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
            R1 = int(read['NumCycles'])-1
        elif read['Number'] == 2:
            R2 = int(read['NumCycles'])-1
    if R1 and R2:
        return 'Paired end - forward read length (R1): {}, reverse read length (R2): {}'.format(R1,R2)
    elif R1 and not R2:
        return 'Single end - read length (R1): {}'.format(R1)
    elif not R1 and not R2:
        return 'Read geometry could not be automatically determined.'

def get_sequencer(run_id):
	return SEQUENCERS.get(run_id.split('_')[1],'Sequencer could not be automatically determined.')

def md5sum_worker(project_dirs) :
    global localConfig
    config = localConfig
    old_wd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
    pnames = get_project_names(project_dirs)
    for p in pnames:
        cmd = "find {} -type f -name '*.fastq.gz' | parallel md5sum > {}".format(
            p,
            'md5sums_{}.txt'.format(p)
        )
        syslog.syslog("[md5sum_worker] Processing %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),p))
        subprocess.check_call(cmd, shell=True)
    os.chdir(old_wd)

def set_mqc_conf_header(config, mqc_conf):

    read_geometry = get_read_geometry(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')))
    contact = config.get('MultiQC','report_contact')
    sequencer = get_sequencer(config.get('Options','runID'))
    prepkit = config.get('Options','Libprep')
    organism = config.get('Options','Organism')

    report_header = [
    {'Contact E-mail': contact},
    {'Sequencing Platform': sequencer},
    {'Read Geometry': read_geometry},
    ]
    if not organism in ['N/A', 'n/a', '', ' ']:
        report_header.append({'Organism': organism})
    if not prepkit in ['N/A', 'n/a', '', ' ']:
        report_header.append({'Lib prep kit': prepkit})

    mqc_conf['report_header_info'] = report_header

    return mqc_conf

def multiqc_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    #os.chdir("{}/{}".format(config.get('Paths','outputDir'), config.get('Options','runID')))
    dname = d.split("/")
    pname = dname[-1]

    conf_name = "{}/{}/QC_{}/.multiqc_config.yaml".format(config.get('Paths','outputDir'), config.get('Options','runID'),pname)
    in_conf = open("/root/multiqc_config.yaml","r")
    out_conf = open(conf_name,"w+")
    mqc_conf = yaml.load(in_conf)

    mqc_conf['title'] = pname
    mqc_conf = set_mqc_conf_header(config,mqc_conf)

    yaml.dump(mqc_conf,out_conf)
    in_conf.close()
    out_conf.close()

    cmd = "{multiqc_cmd} {multiqc_opts} --config {conf} {flow_dir}/QC_{pname} {flow_dir}/{pname} -o {flow_dir}/QC_{pname}".format(
            multiqc_cmd = config.get("MultiQC", "multiqc_command"), 
            multiqc_opts = config.get("MultiQC", "multiqc_options"), 
            conf = conf_name,
            flow_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
            pname=pname,  
            )
    syslog.syslog("[multiqc_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def multiqc_stats(project_dirs) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats'))

    #Illumina interop
    cmd = "interop_summary {} --csv=1 > {}".format(
            os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get("Options","runID")),
            os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats','interop_summary.csv'),
        )
    syslog.syslog("[multiqc_worker] Interop summary on %s\n" % os.path.join(config.get("Paths","baseDir"),config.get("Options","sequencer"),'data',config.get("Options","runID")))
    subprocess.check_call(cmd,shell=True)

    conf_name = "{}/{}/Stats/.multiqc_config.yaml".format(config.get('Paths','outputDir'), config.get('Options','runID'))
    in_conf = open("/root/multiqc_config.yaml","r")
    out_conf = open(conf_name,"w+")
    mqc_conf = yaml.load(in_conf)

    pnames = get_project_names(project_dirs)
    pnames = ' ,'.join(pnames)
    mqc_conf['title'] = pnames

    mqc_conf = set_mqc_conf_header(config,mqc_conf)

    yaml.dump(mqc_conf,out_conf)
    in_conf.close()
    out_conf.close()

    cmd = "{multiqc_cmd} {multiqc_opts} --config {conf} {flow_dir}/Stats -o {flow_dir}/Stats".format(
            multiqc_cmd = config.get("MultiQC", "multiqc_command"), 
            multiqc_opts = config.get("MultiQC", "multiqc_options"), 
            conf = conf_name,
            flow_dir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
            )
    syslog.syslog("[multiqc_worker] Processing %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'Stats'))
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def archive_worker(config):
    project_dirs = get_project_dirs(config)
    pnames = get_project_names(project_dirs)

    for p in pnames:
        if os.path.exists(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p))):
            os.remove(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p)))
        pw = None
        if config.get("Options","SensitiveData") == "1":
            pw = subprocess.check_output("xkcdpass -d '-'",shell=True).decode().strip('\n')
            with open(os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),"encryption.{}".format(p)),'w') as pwfile:
                pwfile.write('{}\n'.format(pw))
        opts = "-p{}".format(pw) if pw else ""
        cmd = "7za a {opts} {flowdir}/{pnr}.7za {flowdir}/{pnr}/ {flowdir}/QC_{pnr} {flowdir}/Stats {flowdir}/Undetermined*.fastq.gz {flowdir}/SampleSheet.csv".format(
                opts = opts,
                flowdir = os.path.join(config.get('Paths','outputDir'), config.get('Options','runID')),
                pnr = p
            )
        syslog.syslog("[archive_worker] Zipping %s\n" % os.path.join(config.get('Paths','outputDir'), config.get('Options','runID'),'{}.7za'.format(p)))
        subprocess.check_call(cmd, shell=True)

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


def clumpify_mark_done(config):
    global localConfig
    config = localConfig
    open(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),"clumpify.done"),"w+").close()

def get_project_names(dirs):
    gcf = set()
    for d in dirs:
        for catalog in d.split('/'):
            if catalog.startswith('GCF-'):
                gcf.add(catalog)
    return gcf

def get_project_dirs(config):
    projectDirs = glob.glob("%s/%s/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    projectDirs.extend(glob.glob("%s/%s/*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"))))
    return toDirs(projectDirs)


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

    projectDirs = get_project_dirs(config)
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
    if config.get("Options","RemoveHumanReads") == "1":
        p = mp.Pool(int(1))
        p.map(RemoveHumanReads_worker, sampleFiles)
        p.close()
        p.join()

    #clumpify

    p = mp.Pool(int(config.get("Options","clumpifyWorkerThreads")))
    p.map(clumpify_worker, sampleFiles)
    p.close()
    p.join()
    clumpify_mark_done(config)

    #FastQC

    p = mp.Pool(int(config.get("Options","fastqcThreads")))
    p.map(FastQC_worker, sampleFiles)
    p.close()
    p.join()

    #md5sum
    md5sum_worker(projectDirs)

    #fastq_screen
    p = mp.Pool(int(config.get("Options", "fastqScreenThreads")))
    p.map(fastq_screen_worker, sampleFiles)
    p.close()
    p.join()

    # multiqc
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(multiqc_worker, projectDirs)
    p.close()
    p.join()

    # multiqc_stats
    multiqc_stats(projectDirs)

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

    #save configfile to flowcell
    with open(os.path.join(config.get("Paths","outputDir"), config.get("Options","runID"),'bcl2fastq.ini'), 'w+') as configfile:
        config.write(configfile)


    return(message)
