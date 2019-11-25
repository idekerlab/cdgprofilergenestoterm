#!/usr/bin/env python

import os
import sys
import argparse
from decimal import *
import tempfile
import json
import ijson
import shutil
import requests
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
    parser.add_argument('--maxgenelistsize', type=int,
                        default=5000, help='Maximum number of genes that can'
                                           'be passed in via a query, '
                                           'exceeding this results in '
                                           'error')
    parser.add_argument('--organism', default='hsapiens',
                        help='Organism to use')
    parser.add_argument('--tmpdir', default='/tmp',
                        help='Temp directory')
    parser.add_argument('--url', default='https://biit.cs.ut.ee/gprofiler/api/gost/profile/')
    return parser.parse_args(args)


def read_inputfile(inputfile):
    """

    :param inputfile:
    :return:
    """
    with open(inputfile, 'r') as f:
        return f.read()


def run_gprofiler(inputfile, theargs):
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

    user_agent = 'cdgprofilergenestoterm/' + cdgprofilergenestoterm.__version__

    result = requests.post(url=theargs.url, json={'organism': theargs.organism,
                                                  'user_threshold': theargs.maxpval,
                                                  'no_evidences': False,
                                                  'domain_scope': 'known',
                                                  'query': genes},
                           headers={'User_Agent': user_agent}, stream=True)
    sys.stdout.write(str(result.headers) + '\n')
    if result.status_code != 200:
        sys.stderr.write('Received non 200 status code: ' + str(result.status_code) + '\n')
        return None

    temp_dir = tempfile.mkdtemp(dir=theargs.tmpdir)
    try:

        tfile = os.path.join(temp_dir, 'raw.out')
        with open(tfile, 'wb') as fd:
            for chunk in result.iter_content(chunk_size=512):
                fd.write(chunk)
        sys.stderr.write('Received: ' + str(os.path.getsize(tfile)) + ' bytes from gprofiler\n')

        with open(tfile, 'r') as f:
            besthit = None
            for entry in ijson.items(f, 'result.item'):
                # print(entry)
                jaccard = Decimal(1.0) / (Decimal(1.0) / entry['precision'] + Decimal(1.0) / entry['recall'] -1)

                if besthit is None or jaccard > besthit['jaccard']:
                    besthit = entry
                    besthit['jaccard'] = jaccard
                    continue
                if jaccard == besthit['jaccard']:
                    if entry['p_value'] < besthit['p_value']:
                        besthit = entry
                        besthit['jaccard'] = jaccard

        theres = {'name': besthit['name'],
                  'source': besthit['source'],
                  'p_value': float(besthit['p_value']),
                  'description': besthit['description'],
                  'intersections': []}

        # besthit['intersections'] needs to be obtained from metadata specifically
        # genes = result['meta']['query'][result['result'][0]['query']]['ensgs']
        # besthit['intersections'] = ([gene for ev, gene in zip(result['result'][0]['intersections'], genes) if ev])
        return theres
    finally:
        pass
        # shutil.rmtree(temp_dir)


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
