#!/usr/bin/env python

import bcl2fastq_pipeline.getConfig
import sys
import pandas as pd
import os

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
    flowcell_processed = flowcells_processed.append(df)
    flowcells_processed.to_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))

def main(argv):

    print(argv)
    config = bcl2fastq_pipeline.getConfig.getConfig()

    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    print(flowcells_processed)




if __name__=='__main__':
    main(sys.argv[1:])

