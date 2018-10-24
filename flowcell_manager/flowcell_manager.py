#!/usr/bin/env python

import bcl2fastq_pipeline.getConfig
import sys
import pandas as pd
import os



def main(argv):

    print(argv)
    config = bcl2fastq_pipeline.getConfig.getConfig()
    manager_base = config.get("FlowCellManager","managerDir")
    flowcells_processed = pd.read_csv(os.path.join(manager_base,'flowcells.processed'), skiprows=1)
    print(flowcells_processed)




if __name__=='__main__':
    main(sys.argv[1:])

