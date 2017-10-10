# coding=utf-8
import argparse
import converter
import os

def convert(args):
    file_abs_path = os.path.abspath(args.f)
    converter.convert_maf_to_po(file_abs_path, args.v, args.m, args.c, args.hbmin, args.fasta)


parser = argparse.ArgumentParser(description='PAN-GENOME tools')

subparsers = parser.add_subparsers(help='Tool', dest='')
subparsers.required = True

parser_converter = subparsers.add_parser('convert', help='Convert MAF files to PO')
parser_converter.add_argument('-f',
                              metavar='MAF_FILE',
                              type=str,
                              required=True,
                              help='Path to the MAF (Multiple Alignment Format) file.')
parser_converter.add_argument('-m',
                              metavar='MERGE_BLOCKS_OPTION',
                              required=False,
                              help='''Set if merging the blocks is required. \
                              Pass [idx1:idx2, idx3:idx4] to choose range of blocks to visualize or 
                              \'all\' to merge all blocks where possible. IDs begin from 1.''')
parser_converter.add_argument('-v',
                              action='store_true',
                              required=False,
                              help='Generate visualization')
parser_converter.add_argument('-c',
                              action = 'store_true',
                              required = False,
                              help='Calculate consensus')
parser_converter.add_argument('-fasta',
                              action='store_true',
                              required=False,
                              help='Generate FASTA files')
parser_converter.add_argument('-hbmin',
                              metavar='HBMIN',
                              required=False,
                              help='HBMIN value for POA heaviest bundling alogrithm, min 0, max 1')
parser_converter.set_defaults(func=convert)

args = parser.parse_args()
args.func(args)

