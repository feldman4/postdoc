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
