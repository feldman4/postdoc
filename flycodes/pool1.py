from ..sequence import get_genbank_features, reverse_translate_max, translate_dna
from ..sequence import reverse_complement as rc

import pandas as pd

# reverse are antisense to forward
features = dict(
    petcon_outer_forward = 'GGGTCGGCTTCGCATATG',
    petcon_outer_reverse = 'GGAACCTCCACCCTCGAG',
    pT14_outer_forward = 'aagaaggagagcagctATG',
    pT14_outer_reverse = 'gagactgccgctactgcg',
    petcon_inner_reverse_1 = 'CCATctgatacgggagcTT',
    petcon_inner_reverse_2 = 'CCATgacccttactgggTT',
    petcon_inner_reverse_3 = 'CCATgttgcccgtaagcTT',
    petcon_inner_reverse_4 = 'CCATtgagcacgacagcTT',
    petcon_inner_reverse_5 = 'CCATtggaactcgggcaTT',
    petcon_inner_reverse_6 = 'CCATtggacaaccagcgTT',
    petcon_inner_forward_1 = 'GCGAgagggtctacagaTT',
    petcon_inner_forward_2 = 'GCGAgcgaccttagagtTT',
    petcon_inner_forward_3 = 'GCGAcccactggcataaTT',
    petcon_inner_forward_4 = 'GCGAgtggcatggcaacTT',
    petcon_inner_forward_5 = 'GCGAccgcgtggttagaTT',
    petcon_inner_forward_6 = 'GCGAgccagttttaggcTT',
)


def fix_adapters(s):
    d = features
    return (s.replace(d['petcon_outer_forward'], d['pT14_outer_forward'])
             .replace(rc(d['petcon_outer_reverse']), rc(d['pT14_outer_reverse'])))


def load_pT14_parts():
    f = 'flycodes/wei_binders/pt14_wy.gb'
    features = get_genbank_features(f)
    assembled = features['assembled oligo']
    design = features['design']

    left, cterm = assembled.split(design)
    linker, right = [x for x in cterm.split('N') if x]
    return left, linker, right


def make_cterm_linker(length, base='GSK', repeat='GGS', left='G'):
    assert length >= len(base)
    backwards_repeat = ('GGS' * int((length + 1)/ 3))[::-1]
    n_to_add = length - len(base)
    assert n_to_add >= 0
    linker = backwards_repeat[:n_to_add][::-1] + base
    linker = 'G' + linker[1:]
    assert len(linker) == length
    return linker


def substitute_cterm_linker(df_agilent, linker='GSK', dna='GGATCC...'):
    """Would be much nicer to make a substitution pattern in mixed aa/dna characters.
    Or have defined components and a template to combine them.
    """
    assert len(linker)*3 == len(dna) # paranoid
    substitute = dna[:dna.index('.')] # up to first dot
    it = df_agilent[['assembly', 'second_insert', 'second', 'barcode', 'primer_4']].values
    fixed_assemblies = []
    fixed_second_inserts = []
    fixed_second = []
    for assembly, second_insert, second, barcode, primer_4 in it:
        # these count backwards from end of sequence
        start = (len(barcode) + len(linker)) * 3
        stop = start - len(substitute)
        # print(start, stop, len(primer_4), substitute)
        fixed_assemblies += [assembly[:-start] + substitute + assembly[-stop:]]
        fixed_second_inserts += [second_insert[:-start] + substitute + second_insert[-stop:]]

        start = (len(barcode) + len(linker)) * 3 + len(primer_4)
        stop = start - len(substitute)
        # print(start, stop)
        fixed_second += [second[:-start] + substitute + second[-stop:]]

        assert translate_dna(fixed_assemblies[-1]) == translate_dna(assembly)
        assert len(fixed_second_inserts[-1]) == len(second_insert)
        assert len(fixed_second[-1]) == len(second)

    return df_agilent.assign(assembly=fixed_assemblies, second_insert=fixed_second_inserts,
                             second=fixed_second)
        

    

def replace_barcodes(df_agilent, barcodes, num_repeats):
    """Expand number of barcodes after codon optimization/overlap design.
    """
    from collections import defaultdict
    barcodes_by_length = defaultdict(list)
    for s in barcodes:
        barcodes_by_length[len(s)] += [s]

    arr = []
    for _, row in df_agilent.iterrows():

        dna, seq, barcode = row[['assembly', 'aa_sequence', 'barcode']]
        barcode_dna = dna[-len(barcode)*3:]
        assert seq.count(barcode) == 1
        keep = 'name', 'overlap', 'overlap_length', 'source'
        for _ in range(num_repeats):
            info = {'old_assembly': dna, 'old_assembly_aa': seq, 'old_barcode': barcode,
                    'old_first': row['first'], 'old_second': row['second']}
            info.update({k: row[k] for k in keep})
            new_barcode = barcodes_by_length[len(barcode)].pop()
            new_barcode_dna = reverse_translate_max(new_barcode)
            info['new_assembly_aa'] = seq[:-len(barcode)] + new_barcode
            info['new_assembly'] = dna[:-len(barcode)*3] + new_barcode_dna
            info['new_barcode'] = new_barcode

            assert row['second'].count(barcode_dna) == 1

            info['new_first'] = row['first']
            info['new_second'] = row['second'].replace(barcode_dna, new_barcode_dna)

            arr += [info]

    df_split = pd.DataFrame(arr)

    assert (0 == (df_split['old_assembly_aa'].str.len() 
                  - df_split['new_assembly_aa'].str.len()).sum())
    
    cols = ['name', 'source', 'new_barcode', 'old_barcode', 'overlap_length', 'overlap',
            'new_first', 'new_second', 'old_first', 'old_second', 
            'new_assembly', 'new_assembly_aa', 'old_assembly', 'old_assembly_aa']

    return df_split[cols]