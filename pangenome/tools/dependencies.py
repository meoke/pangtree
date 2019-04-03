from argparse import Namespace

from datamodel.fasta_providers import FastaProvider
from datamodel.fasta_providers.ConstSymbolProvider import ConstSymbolProvider
from datamodel.fasta_providers.FromNCBI import FromNCBI
from datamodel.fasta_providers.FromFile import FromFile


def get_fasta_provider(args: Namespace) -> FastaProvider:
    if args.fasta_provider is None:
        return ConstSymbolProvider(args.missing_symbol)
    elif args.fasta_provider == 'ncbi':
        if args.email is None:
            raise Exception("Email address must be specified. It must be provided when fasta source is \'ncbi\'.")
        use_cache = args.cache if args.cache else False
        return FromNCBI(args.email, use_cache)
    elif args.fasta_provider == 'file':
        if args.fasta_file is None:
            raise Exception("Fasta file source must be specified. It must be provided when fasta source is \'local\'.")
        return FromFile(args.fasta_file)
    else:
        raise Exception("Not known fasta provider."
                        "Should be \'ncbi\' or \'file\' or None."
                        "Cannot build pangraph.")
