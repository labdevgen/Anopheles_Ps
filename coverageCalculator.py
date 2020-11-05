import straw
from Expected_calculator import get_vaid_chrms_from_straw, get_contacts_using_juicer_dump
import pandas as pd

def getCoverage(file, resolution, juicer_tools_path, norm):
    # compute cis- and trans-coverage for each loci
    # returns dict chr-->coverage dataFrame
    # coverage dataFrame structure: loci start -- cis-coverage -- trans-coverage -- total coverage

    # file - hic-file
    # resolution - resolution to use
    strawObj = straw.straw(file)
    chrms = get_vaid_chrms_from_straw(strawObj, minchrsize=0, excludechrms=())
    contacts = []
    for ind,val in enumerate(chrms):
        for j in range(ind,len(chrms)):
            chrm1 = val
            chrm2 = chrms[j]
            contacts.append(get_contacts_using_juicer_dump(juicer_tools_path,file,chrm1,
                                                                    resolution,chrm2,datatype="observed",
                                                                    norm = norm))
    return pd.concat(contacts)