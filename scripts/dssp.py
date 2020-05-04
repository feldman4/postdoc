#!/usr/bin/env python

from pyrosetta import init, pose_from_pdb

init(options='-mute all', silent=True)

import pyrosetta.rosetta.core.scoring.dssp
import sys

def pdb_to_dssp(filename):
    pose = pose_from_pdb(filename)
    dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
    return dssp.get_dssp_secstruct()

if __name__ == '__main__':
    for f in sys.stdin:
        print(pdb_to_dssp(f.rstrip()))



