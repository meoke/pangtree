from subprocess import run
_mycoplasma_groups_dictionary = {0: '0', #brak grupy
                                 1: '1',
                                 2: '2'
                                }

_mycoplasma_dictionary = {
    'GCF_001509195': ('Mycoplasma pneumoniae', 'S355', 'NZ_CP013829.1', 2),
    'GCF_000143945': ('Mycoplasma pneumoniae', 'FH', 'NC_017504.1', 1),
    'GCF_002090295': ('Mycoplasma pneumoniae', 'S55-tet-R', 'NZ_CP020692.1', 0),
    'GCF_002127985': ('Mycoplasma pneumoniae', '685', 'NZ_CP017328.1',0),
    'GCF_002090235': ('Mycoplasma pneumoniae', 'FH-tet-R', 'NZ_CP020690.1', 0),
    'GCF_000027345': ('Mycoplasma pneumoniae', 'M129', 'NC_000912.1', 2),
    'GCF_001272815': ('Mycoplasma pneumoniae', '85138', 'NZ_CP010545.1', 2),
    'GCF_001272735': ('Mycoplasma pneumoniae', '51494', 'NZ_CP010541.1', 2),
    'GCF_001272795': ('Mycoplasma pneumoniae', '85084', 'NZ_CP010544.1', 2),
    'GCF_002090275': ('Mycoplasma pneumoniae', 'S68-tet-R', 'NZ_CP020691.1', 0),
    'GCF_002090315': ('Mycoplasma pneumoniae', 'S91-tet-R', 'NZ_CP020693.1',0),
    'GCF_002096035': ('Mycoplasma pneumoniae', 'S12-tet-R', 'NZ_CP020712.1',0),
    'GCF_002095995': ('Mycoplasma pneumoniae', 'S4-tet-R', 'NZ_CP020711.1',0),
    'GCF_002096015': ('Mycoplasma pneumoniae', 'S34-tet-R', 'NZ_CP020710.1',0),
    'GCF_002090215': ('Mycoplasma pneumoniae', 'S63-tet-R', 'NZ_CP020689.1',0),
    'GCF_000331085': ('Mycoplasma pneumoniae', 'M129-B7', 'NC_020076.2', 2),
    'GCF_002128065': ('Mycoplasma pneumoniae', 'E16', 'NZ_CP017332.1',0),
    'GCF_002128045': ('Mycoplasma pneumoniae', 'FL8', 'NZ_CP017331.1',0),
    'GCF_002128025': ('Mycoplasma pneumoniae', '549', 'NZ_CP017330.1',0),
    'GCF_000319675': ('Mycoplasma pneumoniae', 'PI1428', 'NZ_CP010538.1',0),
    'GCF_001558175': ('Mycoplasma pneumoniae', 'C267', 'NZ_CP014267.1', 2),
    'GCF_002147855': ('Mycoplasma pneumoniae', 'M129 2002', 'NZ_CP017343.1',0),
    'GCF_001272755': ('Mycoplasma pneumoniae', '54089', 'NZ_CP010542.1', 2),
    'GCF_001272775': ('Mycoplasma pneumoniae', '54524', 'NZ_CP010543.1', 2),
    'GCF_002128145': ('Mycoplasma pneumoniae', 'FL1', 'NZ_CP017333.1',0),
    'GCF_002128185': ('Mycoplasma pneumoniae', 'K27', 'NZ_CP017334.1',0),
    'GCF_001272855': ('Mycoplasma pneumoniae', 'M1139', 'NZ_CP010547.1', 1),
    'GCF_002355695': ('Mycoplasma pneumoniae', 'KCH-402', 'NZ_AP017318.1',0),
    'GCF_002355715': ('Mycoplasma pneumoniae', 'KCH-405', 'NZ_AP017319.1',0),
    'GCF_002128165': ('Mycoplasma pneumoniae', 'CO103', 'NZ_CP017335.1',0),
    'GCF_002128235': ('Mycoplasma pneumoniae', '1134', 'NZ_CP017338.1',0),
    'GCF_001272915': ('Mycoplasma pneumoniae', 'MAC Mac', 'NZ_CP010550.1', 1),
    'GCF_001272875': ('Mycoplasma pneumoniae', 'M2192', 'NZ_CP010548.1', 1),
    'GCF_002128105': ('Mycoplasma pneumoniae', '1801', 'NZ_CP017341.1',0),
    'GCF_000283755': ('Mycoplasma pneumoniae', '309', 'NC_016807.1', 1),
    'GCF_001272715': ('Mycoplasma pneumoniae', '39443', 'NZ_CP010540.1', 1),
    'GCF_001272895': ('Mycoplasma pneumoniae', 'M2592', 'NZ_CP010549.1', 1),
    'GCF_000319655': ('Mycoplasma pneumoniae', 'PO1', 'NZ_CP010551.1',0),
    'GCF_002128205': ('Mycoplasma pneumoniae', '1006', 'NZ_CP017337.1',0),
    'GCF_002128285': ('Mycoplasma pneumoniae', 'GA3', 'NZ_CP017336.1',0),
    'GCF_001901705': ('Mycoplasma pneumoniae', 'FH2009', 'NZ_CP017327.1',0),
    'GCF_000387745': ('Mycoplasma pneumoniae', '19294', 'NZ_CP010539.1',0),
    'GCF_002128085': ('Mycoplasma pneumoniae', 'RI3', 'NZ_CP017340.1',0),
    'GCF_002128005': ('Mycoplasma pneumoniae', 'E57', 'NZ_CP017329.1',0),
    'GCF_002128265': ('Mycoplasma pneumoniae', '519', 'NZ_CP017339.1',0),
    'GCF_002128125': ('Mycoplasma pneumoniae', 'CO3', 'NZ_CP017342.1',0),
    'GCF_000733995': ('Mycoplasma pneumoniae', 'M29', 'NZ_CP008895.1', 2)
}

def extract_ebola_code_part(original_ebola_src_name, index):
    return original_ebola_src_name.split('.')[index]


def get_strain_name(src_name):
    genbankID = src_name.split('.')[0]
    return _mycoplasma_dictionary[genbankID][1]

def get_group_name(src_name):
    genbankID = src_name.split('.')[0]
    groupID = _mycoplasma_dictionary[genbankID][3]
    return _mycoplasma_groups_dictionary[groupID]


def remove_shifted_or_reversed(mycoplasma_file_name):
    run(['sed', '\'/GCF_000733995\|GCF_000387745/d\|label=u0\|label=u1\|label=u2\'', mycoplasma_file_name])
