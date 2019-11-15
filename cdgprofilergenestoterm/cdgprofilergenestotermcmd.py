#!/usr/bin/env python

import os
import sys
import argparse
import json
from gprofiler import GProfiler
import cdgprofilergenestoterm

def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)
    parser.add_argument('input',
                        help='comma delimited list of genes in file')
    parser.add_argument('--maxpval', type=float, default=0.00001,
                        help='Max p value')
    parser.add_argument('--organism', default='hsapiens',
                        help='Organism to use')
    return parser.parse_args(args)


def read_inputfile(inputfile):
    """

    :param inputfile:
    :return:
    """
    with open(inputfile, 'r') as f:
        return f.read()


def run_gprofiler(inputfile, theargs,
                  gprofwrapper=GProfiler(user_agent='cdgprofilergenestoterm/' +
                                         cdgprofilergenestoterm.__version__,
                                         return_dataframe=True)):
    """
    todo

    :param inputfile:
    :param theargs:
    :param gprofwrapper:
    :return:
    """
    genes = read_inputfile(inputfile)
    genes = genes.strip(',').strip('\n').split(',')
    if genes is None or (len(genes) == 1 and len(genes[0].strip()) == 0):
        sys.stderr.write('No genes found in input')
        return None
    df_result = gprofwrapper.profile(query=genes, domain_scope="known",
                                     organism=theargs.organism,
                                     user_threshold=theargs.maxpval,
                                     no_evidences=False)
    if df_result.shape[0] == 0:
        return None

    df_result['Jaccard'] = 1.0 / (1.0 / df_result['precision'] +
                                  1.0 / df_result['recall'] - 1)
    df_result.sort_values(['Jaccard', 'p_value'],
                          ascending=[False, True], inplace=True)
    df_result.reset_index(drop=True, inplace=True)
    theres = {'name': df_result['name'][0],
              'source': df_result['source'][0],
              'p_value': df_result['p_value'][0],
              'description': df_result['description'][0],
              'intersections': df_result['intersections'][0]}

    return theres


def main(args):
    """
    Main entry point for program

    :param args: command line arguments usually :py:const:`sys.argv`
    :return: 0 for success otherwise failure
    :rtype: int
    """
    desc = """
        Running gprofiler-official 1.0.0, with python3!

        Takes file with comma delimited list of genes as input and
        outputs matching term if any
    """

    theargs = _parse_arguments(desc, args[1:])

    try:
        inputfile = os.path.abspath(theargs.input)
        theres = run_gprofiler(inputfile, theargs)
        if theres is None:
            sys.stderr.write('No terms found\n')
        else:
            json.dump(theres, sys.stdout)
        sys.stdout.flush()
        return 0
    except Exception as e:
        sys.stderr.write('Caught exception: ' + str(e))
        return 2


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
