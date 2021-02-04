import hashlib
import os
import pandas as pd
from tqdm.auto import tqdm

from ..pyrosetta import diy
from ..utils import nglob, timestamp, hash_set
from ..lab.rpxdock import generate_command, run_command

rocuronium_order_keys = 'ntf2/201203_rocuronium_designs/order_keys.txt'
rocuronium_designs = 'ntf2/rocuronium_designs.csv'
name_width = 10

num_designs = 10 # number to process?

rpx_beam_size = 50000
rpx_nout_top = 10
rpx_session = '20210202_test'

def rpxdock(name):
    output_prefix = f'ntf2/rpxdock/{rpx_session}/{name}'
    temp_name = f'{os.environ["TMPDIR"]}/{name}.pdb'
    df_inputs = pd.read_csv(rocuronium_designs)
    f = df_inputs.query('name == @name').iloc[0]['file']
    diy.read_pdb(f).pipe(diy.write_pdb, temp_name)
    cmd = generate_command(temp_name, output_prefix=output_prefix, 
            architecture='C2', beam_size=rpx_beam_size, nout_top=rpx_nout_top)
    result = run_command(cmd)
    os.remove(temp_name)
    log_file = timestamp(f'{output_prefix}_rpxlog.txt')
    with open(log_file, 'w') as fh:
        fh.write(result['stdout'])


def load_monomer_table():
    """Load pdb sequences and assign index for simple naming.
    """
    # search = 'ntf2/201203_rocuronium_designs/*pdb'
    # files = nglob(search)
    

    order_keys = pd.read_csv(rocuronium_order_keys, header=None)[0].str.strip()
    files = [f'ntf2/201203_rocuronium_designs/{f}' for f in order_keys]

    df = (pd.DataFrame({'file': files})
     .assign(design=[diy.read_pdb_sequences(f)['A'] for f in tqdm(files)])
    )
    cut = df['design'].duplicated()
    if cut.sum():
        print(f'Dropping {cut.sum()} duplicate design sequences...')
        print(df.loc[cut, 'file'].values)
        df = df.drop_duplicates(subset='design')
    df['name'] = hash_set(df['design'], name_width)
    df.to_csv(rocuronium_designs, index=False)


