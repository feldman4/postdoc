from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqUtils.MeltingTemp import Tm_NN


tev_tag_0 = 'MGSHHHHHHENLYFQGWS'
tev_tag_1 = 'MGHHHHHHGWSENLYFQGS'

bsai_tail = 'CTTACGCACTTACTCATGGTCTCc'
gg_protein_fwd = bsai_tail + 'AAGAGC'
gg_protein_rev = bsai_tail + 'GTTA'

patterns = {'pTL12': 'AGTCGC(.*)AAGAGC',
            'pTL10': 'AGTCGC(.*)AAGACG'}

def translate(s):
    return str(Seq(s, generic_dna).translate())

def edge_primers(dna, melt_temp):    
    dna_rev = str(Seq(dna, generic_dna).reverse_complement())

    fwd = dna[:10]
    while Tm_NN(fwd) < melt_temp:
        fwd = dna[:len(fwd) + 1]
        
    rev = dna_rev[:10]
    while Tm_NN(rev) < melt_temp:
        rev = dna_rev[:len(rev) + 1]
        
    return fwd, rev

def clip_tags(cds, tags):
    protein = translate(cds)
    for tag in tags:  
        if protein.startswith(tag):
            dna = cds[len(tag)*3:]
            break
    else:
        raise ValueError
    
    return dna

def is_a_cds(dna):
    protein = translate(dna)
    return (protein.startswith('M') & 
            protein.endswith('*'))

def design_primers_from_cds(cds, melt_temp=54):
    assert is_a_cds(cds)
    dna = clip_tags(cds, (tev_tag_0, tev_tag_1))
    fwd, rev = edge_primers(dna, melt_temp)
    return gg_protein_fwd + fwd, gg_protein_rev + rev

