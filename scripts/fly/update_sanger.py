#!/usr/bin/env python

if __name__ == '__main__':
    from postdoc.sequence import sanger_database
    from postdoc.drive import Drive
    import os

    drive = Drive()
    df_sanger = sanger_database(drive)
    f = os.path.join(os.environ['HOME'], 'flycodes/sanger/database.csv')
    df_sanger.to_csv(f, index=None)

    print(f'Wrote {len(df_sanger)} entries to {f}')

