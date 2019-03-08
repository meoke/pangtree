_nucleotide_dictionary = {'-': 0,
                          'A': 1, 'C': 2, 'T': 3, 'G': 4,
                          'B': 5,
                          'D': 6,
                          'E': 7,
                          'F': 8,
                          'H': 10,
                          'I': 11,
                          'J': 12,
                          'K': 13,
                          'L': 14,
                          'M': 15,
                          'N': 16,
                          'O': 17,
                          'P': 18,
                          'R': 19,
                          'S': 20,
                          'U': 22,
                          'W': 23,
                          'Y': 24,
                          'Z': 25,
                          'Q': 26,
                          'V': 27,
                          'X': 28}


_nucleotide_inverse_dictionary = {v: k for k, v in _nucleotide_dictionary.items()}


def code(nucleotide_base):
    return _nucleotide_dictionary[nucleotide_base]


def decode(nucleotide_code):
    return _nucleotide_inverse_dictionary[nucleotide_code]
