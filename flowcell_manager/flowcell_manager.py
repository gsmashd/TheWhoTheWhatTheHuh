#!/usr/bin/env python

import bcl2fastq_pipeline.getConfig
import sys
import pandas as pd
import os

HELP_MESSAGE = """
flowcell_manager.py usage \n
\n
flowcell_manager.py add project-name flowcell-path timestamp --- adds project to inventory file\n
flowcell_manager.py list --- lists all projects in inventory file \n
flowcell_manager.py list-project project-name --- lists all occurences of a specific project \n
flowcell_manager.py list-flowcell flowcell-paht --- lists all occurences of a specific flowcell-path \n
flowcell_manager.py help --- print this message \n
"""

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

def list_project(project):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    print(flowcells_processed.loc[flowcells_processed['project'] == project])

def list_flowcell(flowcell):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    print(flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell])

def main(argv):

    config = bcl2fastq_pipeline.getConfig.getConfig()
    if argv[0] == 'add':
        add_flowcell(*argv[1:])
    elif argv[0] == 'list':
        df = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
        print(df)
    elif argv[0] == 'list-project':
        list_project(argv[1])
    elif argv[0] == 'list-flowcell':
        list_flowcell(argv[1])
    elif argv[0] == 'help':
        print(HELP_MESSAGE)
    else:
        print(HELP_MESSAGE)


if __name__=='__main__':
    main(sys.argv[1:])



