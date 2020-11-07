#!/usr/bin/env python

"""
trR_design_v4 applications powered by Google Fire. Lazy imports so command line 
tool loads fast.
"""

import os
import sys
sys.path.append('/home/dfeldman/packages/rtRosetta/')

import fire
import numpy as np

# reduce tensorflow garbage
import logging
logging.getLogger('tensorflow').setLevel(logging.ERROR)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class DSSPDesigner:
    """
    Minimal application to hallucinate sequence with desired DSSP.
    
    Example usage:
    python /home/dfeldman/s/app.py dssp_design \\
        initialize --num_models=1 - \\
        design --out_dir=test/ --dssp_target=dssp_targets_3.fa --N=2

    """
    def initialize(self, num_models=1):
        """Load trR and DSSP models, connect layers, calculate gradients.
        :param num_models: number of trR models to load (1-5)
        """
        from rtRosetta.apis import v4 as v4_api
        from rtRosetta.apis import v4_losses
        import tensorflow as tf

        # digs GPU issue
        v4_api.clear_session()
    
        f = 'wfc/20200723_DSSP/1DConv_model.weights'
        dssp_model = tf.keras.models.load_model(f)

        inputs = v4_api.define_inputs(bkg=True, dssp=True)
        I_feat, I_seq, I_pssm = v4_api.connect_inputs(inputs)
        O_feat, models = v4_api.load_resnet(I_feat, num_models)

        dssp_pred, dssp_loss = v4_losses.loss_dssp(O_feat, inputs['dssp'], dssp_model)
        losses = {
            'bkgr': v4_losses.loss_background(O_feat, inputs['bkg']),
            'dssp': dssp_loss
        }

        loss, loss_order = v4_api.combine_losses(losses, inputs['loss_weights'])
        grad = v4_api.calculate_gradient(inputs['I'], loss)

        self.loss_order = loss_order

        outputs = {'grad': grad, 'loss': loss, 'O_feat': O_feat,
                'dssp_pred': dssp_pred}
        self.model = tf.keras.Model(inputs.values(), outputs)

        return self # allows fire to chain commands


    def design(self, out_dir, dssp_target, loss_bkgr=1, loss_dssp=0.2, N=1, seed=0, 
        remove_aa='C'):
        """Hallucinate a batch of sequences against DSSP targets.
        :param out_dir: name of output directory ending in / (will be created if absent)
        :param dssp_target: either a DSSP string, or a fasta file of DSSP strings (HEL)
        :param loss_bkgr: weight for KL loss
        :param loss_dssp: weight for DSSP CCE loss
        :param N: number of sequences to hallucinate per DSSP target
        :param seed: random seed
        """

        if os.path.basename(out_dir) == '':
            os.makedirs(out_dir, exist_ok=True)
        else:
            raise ValueError(f'{out_dir} not a directory name')
        
        if os.path.exists(dssp_target):
            print(f'Reading dssp targets from fasta file {dssp_target}')
            dssp_targets = dict(read_fasta(dssp_target))
        elif set(dssp_target) - set('HEL'):
            raise ValueError(f'dssp_target must be HEL string or fasta file, not {dssp_target}')
        else:
            dssp_targets = {'noname': dssp_target}
            
            
        for i in range(N):
            for target_name, dssp_target in dssp_targets.items():
                prefix = f'{out_dir}{target_name}_{i:04d}'
                print(f'Designing {prefix}')
                print(f'seed={seed}')
                pred, history = self._design_one(
                    dssp_target, loss_bkgr, loss_dssp, seed, remove_aa)
                self._save(prefix, dssp_target, history, pred, seed)
                seed += 1


    def _design_one(self, dssp_target, loss_bkgr, loss_dssp, seed, remove_aa):
        from rtRosetta.apis import v4 as v4_api
        
        L = len(dssp_target) + 2
        
        loss_mods = {'bkgr': loss_bkgr, 'dssp': loss_dssp}

        input_data = v4_api.create_input(
            loss_mods, self.loss_order, L=L, seed=seed)

        remove = encode_oh(remove_aa).argmax(axis=-1)
        input_data['I'][..., remove] = -1e9

        dssp_oh = encode_oh(dssp_target, 'HEL')
        input_data['dssp'] = dssp_oh[None]

        minimizer = v4_api.OGMinimizer(opt_iter=50)
        pred, history = v4_api.design(self.model, input_data, minimizer)

        return pred, history


    def _save(self, prefix, dssp_target, history, pred, seed):
        from rtRosetta.apis import v4 as v4_api
        import seaborn as sns
        import matplotlib.pyplot as plt

        seq = decode_oh(history[-1]['I'][0, 0])
        dssp_pred = decode_oh(history[-1]['dssp_pred'], 'HEL')

        # save MDS backbone pdb, and scwrl4 sidechain pdb
        xyz, dm = v4_api.feat_to_xyz(pred['O_feat'])
        v4_api.save_PDB(f'{prefix}.mds.pdb', xyz, dm, seq)
        v4_api.do_scwrl(f'{prefix}.mds.pdb', f'{prefix}.scwrl4.pdb')

        # save npz format
        np.savez(f'{prefix}_pred.npz', **pred)

        # save CB-CB dist max heatmap
        dist_max = v4_api.split_feat(pred['O_feat'])['dist'].argmax(axis=-1)
        fig, ax = plt.subplots()
        sns.heatmap(dist_max, cmap='viridis', ax=ax)
        fig.savefig(f'{prefix}_dist.png')
        plt.close(fig)

        log = '\n'.join([
            '# seq and predicted DSSP',
            seq, f'-{dssp_pred}-', 
            '',
            '# target DSSP vs final DSSP', 
            print_alignment(dssp_target, dssp_pred, as_string=True),
            '',
            '# seed', str(seed),
            ])

        with open(f'{prefix}_result.txt', 'w') as fh:
            fh.write(log)
        print('-'*10, f'{prefix}_result.txt', '-'*10)
        print(log)

            


class SimplePredict:
    """
    Test application.
    """
    def predict(self, seq_or_fasta, num_models=1, diag=0.4):
        """
        Predict distance and angle map from sequence. 
        
        :param seq_or_fasta: protein sequence provided as string 
            (name defaults to "seq") or fasta file of protein sequences
        :param num_models: number of trRosetta models to use
        :param diag: diagonal value for 2D inverse covariance input
        """
        from rtRosetta.apis import v4 as v4_api
        from postdoc.sequence import read_fasta

        model = self._load_predict_model(num_models, diag)
        
        if os.path.exists(seq_or_fasta):
            sequences = read_fasta(seq_or_fasta)
        else:
            sequences = (('seq', seq_or_fasta),)

        loss_mods, loss_order = {}, []
        self.pred = {}
        for name, seq in sequences:
            input_data = v4_api.create_input(loss_mods, loss_order, seq=seq)
            self.pred[name] = v4_api.predict_and_format(model, input_data)

        return self

    def _load_predict_model(self, num_models, diag):
        from rtRosetta.apis import v4 as v4_api
        # clear namespace and fix digs GPU issue
        v4_api.clear_session()

        # create tensorflow Input layers
        inputs = v4_api.define_inputs()
        I_feat, _, _ = v4_api.connect_inputs(inputs, diag=diag)
        O_feat, _ = v4_api.load_resnet(I_feat, num_models)
        outputs = {'O_feat': O_feat}
        model = v4_api.tf.keras.Model(inputs.values(), outputs)
        return model

    def save(self, prefix, dataset=True, numpy=True):
        """
        Save results after calculation. 
        For example:
            app.py predict misc/app_test/seqs.fa - save misc/app_test/
        :param dataset: save predicted features as xarray.Dataset (Ofeat.nc)
        :param numpy: save model input and output (pred.npz)
        """
        from postdoc.wfc.dataset import feat_to_ds
        import numpy as np

        if os.path.basename(prefix) == '':
            os.makedirs(prefix, exist_ok=True)
        else:
            raise ValueError(f'{prefix} not a directory name')
        
        for name, pred in self.pred.items():
            if dataset:
                filename = f'{prefix}{name}_Ofeat.nc'
                feat_pred = feat_to_ds(pred['O_feat'])
                feat_pred.to_netcdf(filename)
            if numpy:
                filename = f'{prefix}{name}_pred.npz'
                np.savez(filename, **pred)


### HELPERS

aa_code_gap = 'ARNDCQEGHILKMFPSTWYV-'
aa_code = aa_code_gap[:-1]

def encode_oh(value, code=aa_code):
    if isinstance(value, str):
        value = np.array(list(value))
    n = len(code)
    index = np.zeros_like(value, dtype=int)
    index[:] = [code.index(v) for v in value.flat[:]]
    return np.eye(n)[index]


def decode_oh(data, code=aa_code):
    ix = data.argmax(-1)
    decoded = np.array(list(code))[ix]
    if data.ndim == 2:
        return ''.join(decoded)
    elif data.ndim == 3:
        return [''.join(x) for x in decoded]
    else:
        return decoded


def print_alignment(a, b, width=60, as_string=False):
    """Levenshtein alignment.
    """
    import edlib
    alignment = edlib.align(a, b, task='path')
    d = edlib.getNiceAlignment(alignment, a, b)
    lines = []
    for i in range(0, max(map(len, d.values())), width):
        lines += [str(i)]
        for x in d.values():
            lines += [x[i:i+width]]

    txt = '\n'.join(lines)
    if as_string:
        return txt
    else:
        print(txt)



def read_fasta(f):
    if f.endswith('gz'):
        fh = gzip.open(f)
        txt = fh.read().decode()
    else:
        fh = open(f, 'r')
        txt = fh.read()
    fh.close()
    return parse_fasta(txt)


def parse_fasta(txt):
    entries = []
    for raw in txt.split('>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries


if __name__ == '__main__':
    fire.Fire({'predict': SimplePredict, 'dssp_design': DSSPDesigner})
