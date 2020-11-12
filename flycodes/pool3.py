import pandas as pd
from postdoc.sequence import translate_dna, print_alignment
from postdoc.flycodes import assembly
from Levenshtein import distance

def patch_bamhi_linker(df_agilent):
    arr = []
    for _, row in df_agilent.iterrows():
        x = row['new_second']
        n = -(len('GSK' + row['barcode'] + 'GSGGSG') * 3)
        y = x[:n] + 'GGATCC' + x[n+6:]
        a = translate_dna(x[n:])
        b = translate_dna(y[n:])
        assert a == b
        x_, y_ = x[n:-20], y[n:-20]
        assert row['new_assembly'].count(x_) == 1
        row['new_assembly'] = row['new_assembly'].replace(x_, y_)
        row['new_second'] = row['new_second'].replace(x_, y_)
        arr += [row]
    return pd.concat(arr, axis=1).T


def check_overlap_quality(df_new, num_to_align=3, k=12):
    gb = df_new.drop_duplicates('name').groupby(['primer_2', 'primer_3'])

    for (p2, p3), df in gb:
        x = assembly.calculate_overlap_fraction(df['overlap'], k)
        y = assembly.calculate_overlap_fraction(df['new_assembly'], k)
        print(f'subpool {p3}: {len(df)} designs')
        print(f'  shared {k}-mers in "overlap" region: {x:.2f}')
        print(f'  shared {k}-mers in entire assembly: {y:.2f}')

        xs = df['overlap'].pipe(list)
        arr = []
        for x in xs:
            for y in xs:
                if x == y:
                    continue
                arr += [(distance(x, y), x, y)]
        arr = sorted(arr)

        for example in arr[:num_to_align]:
            print_alignment(example[1], example[2])

        print('-' * 50)


def add_start_length(df_new):
    arr = []
    starts = []
    lengths = []
    seqs = []
    for a, b, c, d in df_new[['new_assembly', 'overlap', 'new_first', 'new_second']].values:
        arr += [a.count(b)]

        frame = a.index(b) % 3
        start = a.index(b) - frame
        length = len(b) + frame
        length += 2
        length -= length % 3
        # this better work
        aa = translate_dna(a[start:start + length])
        starts += [start]
        lengths += [length]
        seq = a[start + 3:start + length - 3]
        assert c.count(seq) == 1
        assert d.count(seq) == 1
        seqs += [seq]

    assert set(arr) == {1}
    return df_new.assign(opt_start=starts, opt_lengths=lengths, opt_seqs=seqs)


def substitute_changes(df):
    cols = ['name', 'new_first', 'new_second', 'opt_seqs', 'opt_polished']
    arr = []
    for name, a, b, c, d in df[cols].values:
        assert a.count(c) == 1
        assert b.count(c) == 1
        first = a.replace(c, d)
        second = b.replace(c, d)
        arr += [[name, first], [name, second]]
    return pd.DataFrame(arr, columns=['name', 'oligo'])
