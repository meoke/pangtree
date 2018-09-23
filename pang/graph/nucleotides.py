_nucleotide_dictionary = {'-': 0,
                            'A': 1, 'C': 2, 'T': 3, 'G': 4,
                            'P': 5, 'K': 6, 'M': 7, 'I': 8, 'L': 9,
                            'R': 10, 'L': 11, 'Q': 12, 'H': 13, 'N': 14,
                            'E': 15, 'V': 16, 'S': 17, 'D': 18, 'F': 19,
                            'Y': 20}


_nucleotide_inverse_dictionary = {v: k for k, v in _nucleotide_dictionary.items()}


def code(nucleotide_base):
    return _nucleotide_dictionary[nucleotide_base]


def decode(nucleotide_code):
    return  _nucleotide_inverse_dictionary[nucleotide_code]
