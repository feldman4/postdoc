from loess.loess_1d import loess_1d


def predict_iRT_loess(df_peptides, sample=1000):
    x_var, y_var = 'RTime', 'iRT'
    x, y = df_peptides.sample(sample, random_state=0)[[x_var, y_var]].values.T
    xnew = df_peptides['RTime'].drop_duplicates().values
    xout, yout, wout = loess_1d(x, y, xnew=xnew, degree=1, frac=0.5,
                                npoints=None, rotate=False, sigy=None)
    return df_peptides['RTime'].map(dict(zip(xout, yout)))