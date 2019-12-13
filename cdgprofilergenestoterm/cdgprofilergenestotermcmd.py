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
    parser.add_argument('--maxpval', type=float, default=0.00000001,
                        help='Max p value')
    parser.add_argument('--omit_intersections', action='store_true',
                        help='If set, do NOT query for gene intersections')
    parser.add_argument('--maxgenelistsize', type=int,
                        default=500, help='Maximum number of genes that can'
                                          'be passed in via a query, '
                                          'exceeding this results in '
                                          'error')
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
    genelist_size = len(genes)
    if genes is None or (genelist_size == 1 and len(genes[0].strip()) == 0):
        sys.stderr.write('No genes found in input')
        return None
    if genelist_size > theargs.maxgenelistsize:
        sys.stderr.write('Gene list size of ' +
                         str(genelist_size) +
                         ' exceeds max gene list size of ' +
                         str(theargs.maxgenelistsize))
        return None
    df_result = gprofwrapper.profile(query=genes, domain_scope="known",
                                     organism=theargs.organism,
                                     user_threshold=theargs.maxpval,
                                     no_evidences=theargs.omit_intersections)
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
              'term_size': int(df_result['term_size'][0])}

    if theargs.omit_intersections is True:
        theres['intersections'] = []
    else:
        theres['intersections'] = df_result['intersections'][0]
    return theres


def main(args):
    """
    Main entry point for program

    :param args: command line arguments usually :py:const:`sys.argv`
    :return: 0 for success otherwise failure
    :rtype: int
    """
    desc = """
        Using gprofiler-official 1.0.0, Python module, this
        program takes a file with comma delimited list of genes 
        as input and outputs best matching term in JSON format to
        standard out. Any log messages and errors are output to
        standard error.
        
        Return 0 upon success otherwise error.
        
        Format of JSON output:
        
        {
         "name":"<TERM NAME>",
         "source":"<SOURCE, IF ANY, WHERE TERM NAME WAS OBTAINED>",
         "p_value":<PVALUE>,
         "description":"<DESCRIPTION, IF ANY, FOR TERM>",
         "term_size":<NUMBER OF GENES ASSOCIATED WITH TERM>,
         "intersections":["<LIST OF GENES USED TO GET TERM>"]
        }

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
