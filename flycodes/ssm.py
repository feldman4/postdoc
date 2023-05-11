from Levenshtein import editops

aa_order = list('AVLIMFWDERKHSTYNQGPC')

def add_mutations(df_ssm, design, col='aa_seq'):
    all_codes = []
    for aa_seq in df_ssm[col]:
        edits = editops(design, aa_seq)
        codes = []
        for op, i, j in edits:
            assert op == 'replace'
            codes += [f'{design[i]}{i+1}{aa_seq[i]}']
        all_codes += [codes]
    return (df_ssm
            .assign(mutation=[','.join(x) for x in all_codes])
            .assign(num_mutations=[len(x) for x in all_codes])
           )


def to_wide(df_ssm, col):
    return (df_ssm
            .pivot_table(columns='position', index='mutant', values=col, aggfunc='first')
            )


def add_mutations2(df_ssm, design, col='aa_seq'):
    all_codes = []
    for aa_seq in df_ssm[col]:
        assert len(aa_seq) == len(design)
        codes = []
        for i, (c1, c2) in enumerate(zip(design, aa_seq)):
            if c1 != c2:
                codes += [f'{c1}{i+1}{c2}']
        all_codes += [codes]
    return (df_ssm
            .assign(mutation=[','.join(x) for x in all_codes])
            .assign(num_mutations=[len(x) for x in all_codes])
            )


def make_ssm(parent_seq):
    arr = [parent_seq[:i] + aa + parent_seq[i+1:] for i in range(len(parent_seq)) for aa in aa_order]
    assert 20 * len(parent_seq) == len(arr)
    return [x for x in arr if x != parent_seq]

