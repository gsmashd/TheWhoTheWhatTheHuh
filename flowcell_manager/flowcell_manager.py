#!/usr/bin/env python

import bcl2fastq_pipeline.getConfig
import sys
import pandas as pd
import os
import subprocess

HELP_MESSAGE = """
flowcell_manager.py usage \n
\n
flowcell_manager.py add project-name flowcell-path timestamp --- adds project to inventory file\n
flowcell_manager.py delete-fowcell flowcell-path --- deletes the flowcell and all containing projects\n\t\t use with --force to omit prompt\n
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
            }
            ]
    
    df = pd.DataFrame(row_list)
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    flowcells_processed = flowcells_processed.append(df)
    flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp'],
            )


def delete_flowcell(args,force=False):
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
            columns = ['project','flowcell_path','timestamp'],
            )
 

def list_processed():
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['timestamp'] != '0')]

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
    print("Project \t Flowcell path \t Timestamp")
    for i, row in df.iterrows():
        print("{}\t{}\t{}".format(row['project'], row['flowcell_path'], row['timestamp']))

def main(argv):

    config = bcl2fastq_pipeline.getConfig.getConfig()
    if argv[0] == 'add':
        add_flowcell(*argv[1:])
    elif argv[0] == 'delete-flowcell':
        delete_flowcell(*argv[1:])
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
    main(sys.argv[1:])



