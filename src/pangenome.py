# coding=utf-8
import argparse
import converter
import os
from re import fullmatch

def convert(args):
    if args.m and not fullmatch('\[\d+\:\d+\]', str(args.m)):
        parser.error('Wrong formatting for MERGE_BLOCKS')

    if (args.iter or args.hbmin or args.min_comp) and not args.c:
        parser.error("Options -iter, -hbmin, -min-comp requires -c.")

    try:
        if args.hbmin:
            hbmin = float(args.hbmin)
            if hbmin < 0 or hbmin > 1:
                parser.error('HBMIN must be a float in range [0,1].')
        if args.min_comp:
            min_comp = float(args.min_comp)
            if min_comp < 0 or hbmin > 1:
                parser.error('MIN_COMP must be in range [0,1].')
    except ValueError:
        parser.error('HBMIN and MIN_COMP must be a float in range [0,1].')

    if args.data != 'ebola' and args.data != 'mycoplasma':
        parser.error('METAVAR must be \'ebola\' or \'mycoplasma\'')

    file_abs_path = os.path.abspath(args.f)
    converter.convert_maf_to_po(file_abs_path, args.v, args.m, args.c, args.hbmin, args.fasta)


parser = argparse.ArgumentParser(description='PAN-GENOME tools')

subparsers = parser.add_subparsers(help='Tools available', dest='')
subparsers.required = True

parser_converter = subparsers.add_parser('mln', help='Processing multialignment')
parser_converter.add_argument('-f',
                              metavar='MAF_FILE',
                              type=str,
                              required=True,
                              help='path to the MAF (Multiple Alignment Format) file')
parser_converter.add_argument('-m',
                              metavar='MERGE_BLOCKS',
                              required=False,
                              help='''default behaviour is to merge all blocks, where possible; \nprovide MERGE_BLOCKS if special way of merging the blocks is required; \npass [idx1:idx2, idx3:idx4] to choose range of blocks to be merged; \nIDs begin from 1.''')
parser_converter.add_argument('-fasta',
                              action='store_true',
                              required=False,
                              help='generate FASTA files')
parser_converter.add_argument('-c',
                              action = 'store_true',
                              required = False,
                              help='generate consensus')
parser_converter.add_argument('-iter',
                              action='store_true',
                              required=False,
                              help='generate consensus iteratively')
parser_converter.add_argument('-hbmin',
                              metavar='HBMIN',
                              required=False,
                              help='HBMIN value for POA heaviest bundling alogrithm, min 0, max 1')
parser_converter.add_argument('-min_comp',
                              metavar='MINCOMP',
                              required=False,
                              help='minimum compatibility between source and consensus to match them')
parser_converter.add_argument('-v',
                              action='store_true',
                              required=False,
                              help='generate visualization')
parser_converter.add_argument('-data',
                              metavar="DATATYPE",
                              required=True,
                              help='ebola or mycoplasma')
parser_converter.set_defaults(func=convert)

args = parser.parse_args()
args.func(args)

