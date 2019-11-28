#!/usr/bin/env python

import bcl2fastq_pipeline.getConfig
import sys
import pandas as pd
import os
import subprocess
import datetime

HELP_MESSAGE = """
flowcell_manager.py usage \n
\n
flowcell_manager.py add project-name flowcell-path timestamp --- adds project to inventory file\n
flowcell_manager.py archive-fowcell flowcell-path --- deletes the fastq.gz files and the .7za files for the flowcell\n\t\t use with --force to omit prompt\n
flowcell_manager.py rerun-fowcell flowcell-path --- deletes the flowcell and all containing projects\n\t\t use with --force to omit prompt\n
flowcell_manager.py list --- lists all processed projects in inventory file \n
flowcell_manager.py list-all --- lists all projects in inventory file (also unprocessed)\n
flowcell_manager.py list-project project-name --- lists all occurences of a specific project \n
flowcell_manager.py list-flowcell flowcell-path --- lists all occurences of a specific flowcell-path that is processed\n
flowcell_manager.py list-flowcell-all flowcell-path --- lists all occurences of a specific flowcell-path, even unprocessed flowcells \n
flowcell_manager.py help --- print this message \n
"""

pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_columns', 6)

def add_flowcell(project,path,timestamp):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    row_list = [
            {
            'project': project,
            'flowcell_path': path,
            'timestamp': timestamp,
            'archived': 0
            }
            ]
    
    df = pd.DataFrame(row_list)
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    flowcells_processed = flowcells_processed.append(df)
    flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp', 'archived'],
            )


def archive_flowcell(args,force=False):
    if '--force' in args:
        force = True
        args.remove('--force')
    flowcell = args[0]
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    fc_for_deletion = flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell]
    if fc_for_deletion.empty:
        print("No such flowcell in inventory!")
        return
    confirm = 'yes'
    if not force:
        print("Please confirm deletion of the following flowcell and the contained projects.\n")
        print(fc_for_deletion)
        confirm = input("Delete? (yes/no)")
    if confirm == 'yes':
        deletions = [os.path.join(flowcell,pid) for pid in fc_for_deletion['project']]
        deletions.append("{}/*.fastq.gz".format(flowcell))
        deletions.append("{}/*.7za".format(flowcell))
        cmd = "rm -rf {}".format(" ".join(deletions))
        print("DELETING: {}".format(cmd))
        subprocess.check_call(cmd, shell=True)
        flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell,'archived'] = datetime.datetime.now()
        flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp','archived'],
            )
 
def rerun_flowcell(args,force=False):
    if '--force' in args:
        force = True
        args.remove('--force')
    flowcell = args[0]
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    fc_for_deletion = flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell]
    if fc_for_deletion.empty:
        print("No such flowcell in inventory!")
        return
    confirm = 'yes'
    if not force:
        print("Please confirm deletion of the following flowcell and the contained projects.\n")
        print(fc_for_deletion)
        confirm = input("Delete? (yes/no)")
    if confirm == 'yes':
        cmd = "rm -rf {}".format(flowcell)
        print("DELETING FLOWCELL: {}".format(cmd))
        subprocess.check_call(cmd, shell=True)
        flowcells_processed = flowcells_processed.loc[flowcells_processed['flowcell_path'] != flowcell]
        flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp','archived'],
            )

def list_processed():
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['timestamp'] != '0') | (flowcells_processed['archived'] != '0')]

def list_project(project):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['project'] == project) & (flowcells_processed['timestamp'] != '0')]

def list_flowcell(flowcell):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['flowcell_path'] == flowcell) & (flowcells_processed['timestamp'] != '0')]

def list_flowcell_all(flowcell):
    #USED TO AVOID RUNNING OLD FLOWCELLS
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell]

def pretty_print(df):
    print("Project \t Flowcell path \t Timestamp \t Archived")
    for i, row in df.iterrows():
        print("{}\t{}\t{}\t{}".format(row['project'], row['flowcell_path'], row['timestamp'], row['archived']))

def main(argv):

    config = bcl2fastq_pipeline.getConfig.getConfig()
    if argv[0] == 'add':
        add_flowcell(*argv[1:])
    elif argv[0] == 'archive-flowcell':
        archive_flowcell(argv[1:])
    elif argv[0] == 'rerun-flowcell':
        rerun_flowcell(argv[1:])
    elif argv[0] == 'list':
        pretty_print(list_processed())
    elif argv[0] == 'list-all':
        df = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
        pretty_print(df)
    elif argv[0] == 'list-project':
        pretty_print(list_project(argv[1]))
    elif argv[0] == 'list-flowcell':
        pretty_print(list_flowcell(argv[1]))
    elif argv[0] == 'list-flowcell-all':
        pretty_print(list_flowcell_all(argv[1]))
    elif argv[0] == 'help':
        print(HELP_MESSAGE)
    else:
        print(HELP_MESSAGE)


if __name__=='__main__':
    #main(sys.argv[1:])
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultHelpFormatter)
    #parser.add_argument("--add", help="Add flowcell to inventory")
    subparsers = parse.add_subparsers()
    parser_add = subparsers.add_parser("add",help="adds project to inventory file")
    parser_add.set_defaults(func=add_flowcell)
    parser_add.add_argument("project",type=str,help="GCF project number")
    parser_add.add_argument("path",type=str,help="Flowcell path")
    parser_add.add_argument("timestamp",type=datetime.datetime.fromisoformat,help="A datetime isoformat string")

    args = parser.parse_args()
    args.func(args)


