%chk=305_trj10_frag1.chk
# hf/3-21g geom=connectivity

MODEL        1

0 2
 C(PDBName=C,ResName=UNL,ResNum=1)                    -1.75300000    2.16400000    1.41400000
 C(PDBName=C,ResName=UNL,ResNum=1)                    -0.96700000    4.95200000    0.44900000
 C(PDBName=C,ResName=UNL,ResNum=1)                     0.99800000    2.60900000   -0.52200000
 C(PDBName=C,ResName=UNL,ResNum=1)                     0.87600000    1.95300000    3.45300000
 O(PDBName=O,ResName=UNL,ResNum=1)                     1.36900000    1.48400000    4.30400000
 O(PDBName=O,ResName=UNL,ResNum=1)                     0.97800000    2.15500000    2.13000000
 Si(PDBName=Si,ResName=UNL,ResNum=1)                  -0.11800000    3.24200000    0.82100000
 H(PDBName=H,ResName=UNL,ResNum=1)                    -2.44900000    1.80200000    0.52500000
 H(PDBName=H,ResName=UNL,ResNum=1)                    -2.44500000    2.83700000    1.62400000
 H(PDBName=H,ResName=UNL,ResNum=1)                    -1.82800000    1.25900000    1.96700000
 H(PDBName=H,ResName=UNL,ResNum=1)                    -0.16300000    4.98400000   -0.38200000
 H(PDBName=H,ResName=UNL,ResNum=1)                    -1.54400000    5.03000000    1.35200000
 H(PDBName=H,ResName=UNL,ResNum=1)                    -1.87700000    4.10500000   -0.03500000
 H(PDBName=H,ResName=UNL,ResNum=1)                     1.91800000    2.10400000   -0.54000000
 H(PDBName=H,ResName=UNL,ResNum=1)                     1.31100000    3.42400000   -1.09100000
 H(PDBName=H,ResName=UNL,ResNum=1)                     0.27200000    2.02300000   -1.08700000

 1 8 1.0 7 1.0 9 1.0 10 1.0
 2 11 1.0 13 1.0 7 1.0 12 1.0
 3 15 1.0 16 1.0 14 1.0 7 1.0
 4 6 1.5 5 3.0
 5
 6 7 1.0
 7
 8
 9
 10
 11
 12
 13
 14
 15
 16


