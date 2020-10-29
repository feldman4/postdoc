import postdoc.sequence
import os

def update_sanger():
    from postdoc.drive import Drive
    drive = Drive()
    df_sanger = postdoc.sequence.sanger_database(drive)
    f = os.path.join(os.environ['HOME'], 'flycodes/sanger/database.csv')
    df_sanger.to_csv(f, index=None)

    print(f'Wrote {len(df_sanger)} entries to {f}')

if __name__ == '__main__':
    import fire
    fire.Fire()