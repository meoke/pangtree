from datamodel.fasta_providers.FastaProvider import EmailAddress, FastaProvider


class FromNCBI(FastaProvider):
    def __init__(self, email_address: EmailAddress, use_cache: bool):
        pass