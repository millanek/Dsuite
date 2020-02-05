#!/usr/bin/env python3

#General python libraries
import os, sys
import argparse
import pandas as pd
from matplotlib import pyplot as plt

#This library is deployed with Dsuite
import dtools




def main():
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description="Plot f-branch statistic as produced by Dsuite. "
                                                    "Produces .png and .svg files.",
                                                                                 add_help=True)
    argparser.add_argument("fbranch", type=argparse.FileType('r'),
                            help="Path to file containing f-branch matrix as produced by Dsuite Fbranch.")
    argparser.add_argument("tree", type = argparse.FileType('r'),
                            help="Path to .newick tree file as given to Dsuite Fbranch.")

    argparser.add_argument("-n", "--run-name", type=str,
                            help="Base file name for output plots.",default="fbranch")
    argparser.add_argument( "--outgroup", type=str,
                        help="Outgroup name in newick file.",default="Outgroup")

    args = argparser.parse_args()
    fb = pd.read_csv(args.fbranch, sep='\t',
                  index_col=[0,1])
    tree = dtools.HsTree(args.tree.read())
    fb1, tree_no_outgroup = dtools.align_fbranch_with_tree(fb, tree,
                                                     outgroup=args.outgroup)
    ax = dtools.plot_fbranch(fb1, tree_no_outgroup)
    plt.savefig(args.run_name+'.svg', bbox_inches='tight')
    plt.savefig(args.run_name+'.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    sys.exit(main())
