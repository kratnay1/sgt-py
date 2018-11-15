""" Main decomposition algorithm """

import group_modules_sym as gm
import requests
import sys

def main():
    # ITA number of grp to decompose
    gnum = sys.argv[1]
    # read file containing all Gamma_S < (Gamma)^Z
    fnames = gm.get_fnames(gnum)
    f_err = open('error_log', 'w')
    for fname in fnames:
        print('group: ' + sys.argv[1], f_err)
        print('file: ' + fname, f_err)
        subgroups = gm.get_subgroups('s_matrices/normal_split/s_{}_matrices/{}'.format(gnum, fname))
        # load S groups
        S = gm.get_s_groups(subgroups)
        # decompisition algorithm
        decomps = gm.decomp(gnum, subgroups, S, 'decomp/non_semi_split/decomp_{}/{}'.format(gnum, fname), f_err)
    f_err.close()



if __name__ == "__main__":
    main()






