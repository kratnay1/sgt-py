""" Main decomposition algorithm """

import group_modules as gm
import requests
import sys

def main():
    # ITA number of grp to decompose
    gnum = sys.argv[1]
    # read file containing all Gamma_B < (Gamma)^Z
    subgroups = gm.get_subgroups('b_matrices/normal/b_{}_normal.dat'.format(gnum))
    # load B groups
    B = gm.get_b_groups(subgroups)
    # decompisition algorithm
    decomps = gm.decomp(gnum, subgroups, B)



if __name__ == "__main__":
    main()





