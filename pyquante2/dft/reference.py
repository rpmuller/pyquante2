"""
This is reference data taken from the Density Functional Repository.
ftp://ftp.dl.ac.uk/qcg/dft_library/index.html

Need to think of a more comprehensive way to do this. Maybe keep data
files in the repository, and run the function against these.

These are parsed in to a table with 11 columns: the first 5 are inputs
and the next six are outputs. The field names are given by the array
'field_names'
"""

import numpy as np

field_names = ['rhoa','rhob','sigmaaa','sigmaab','sigmabb',
               'zk','vrhoa','vrhob','vsigmaaa','vsigmaab','vsigmabb']

x_lda_data = """\
 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
 zk            = -0.377592720836E+01
 vrhoa         = -0.148075576798E+01
 vrhob         = -0.148075576798E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.290344268232E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.290344268232E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
 zk            = -0.377592720836E+01
 vrhoa         = -0.148075576798E+01
 vrhob         = -0.148075576798E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.290344268232E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.290344268232E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
 zk            = -0.319555819038E+01
 vrhoa         = -0.142024808461E+01
 vrhob         = -0.142024808461E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.315610685470E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.315610685470E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
 zk            = -0.728453690414E-01
 vrhoa         = -0.551858856374E+00
 vrhob         = -0.551858856374E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.209037445596E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.209037445596E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
 zk            = -0.407494475320E+05
 vrhoa         = -0.150923879748E+02
 vrhob         = -0.150923879748E+02
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.279488666200E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.279488666200E-02

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
 zk            = -0.407494475320E+05
 vrhoa         = -0.150923879748E+02
 vrhob         = -0.150923879748E+02
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.279488666200E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.279488666200E-02

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
 zk            = -0.348271841145E+05
 vrhoa         = -0.145113267144E+02
 vrhob         = -0.145113267144E+02
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.302319306550E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.302319306550E-02

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
 zk            = -0.308832394647E+00
 vrhoa         = -0.791877934993E+00
 vrhob         = -0.791877934993E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.101522812179E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.101522812179E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
 zk            = -0.370503980143E+07
 vrhoa         = -0.466042742318E+02
 vrhob         = -0.466042742318E+02
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.293108642967E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.293108642967E-03

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
 zk            = -0.315661921284E+07
 vrhoa         = -0.447747406077E+02
 vrhob         = -0.447747406077E+02
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.317551351828E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.317551351828E-03

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
 zk            = -0.148324672136E+00
 vrhoa         = -0.659220765051E+00
 vrhob         = -0.659220765051E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.146493503345E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.146493503345E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.494484233083E+01
 vrhoa         = -0.188374945936E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.179404710416E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.494484233083E+01
 vrhoa         = -0.188374945936E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.179404710416E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.402615103023E+01
 vrhoa         = -0.178940045788E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.198822273098E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.208913119508E-01
 vrhoa         = -0.480260044845E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.276011520026E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
 zk            = -0.657615683804E+03
 vrhoa         = -0.539020244480E+01
 vrhob         = -0.536820137364E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.219113920520E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.220913636775E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
 zk            = -0.241948147838E+03
 vrhoa         = -0.420747936684E+01
 vrhob         = -0.417120618800E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.359613621097E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.365895279649E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
 zk            = -0.101616142698E+00
 vrhoa         = -0.628513933519E+00
 vrhob         = -0.566119777958E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.161157418851E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.198638518582E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
 zk            = -0.400731073431E-01
 vrhoa         = -0.530109182127E+00
 vrhob         = -0.389751405963E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.226542385525E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.419087533293E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
 zk            = -0.338253135027E+03
 vrhoa         = -0.457078149734E+01
 vrhob         = -0.454010418713E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.304718766489E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.308850624975E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
 zk            = -0.254588260307E+03
 vrhoa         = -0.424313767179E+01
 vrhob         = -0.424313767179E+01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.353594805982E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.353594805982E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
 zk            = -0.982681500273E-01
 vrhoa         = -0.611966348389E+00
 vrhob         = -0.575882382297E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.169990652330E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.191960794099E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
 zk            = -0.230346081831E-01
 vrhoa         = -0.450900660715E+00
 vrhob         = -0.362783167860E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.313125458830E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.483710890480E+01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00"""

c_vwn5_data = """\
 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
 zk            = -0.278978177367E+00
 vrhoa         = -0.907896301530E-01
 vrhob         = -0.907896301530E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.129443214985E-01
 v2rhoab       = -0.182559901422E-01
 v2rhob2       =  0.129443214985E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
 zk            = -0.278978177367E+00
 vrhoa         = -0.907896301530E-01
 vrhob         = -0.907896301530E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.129443214985E-01
 v2rhoab       = -0.182559901422E-01
 v2rhob2       =  0.129443214985E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
 zk            = -0.242883397986E+00
 vrhoa         = -0.896613951966E-01
 vrhob         = -0.896613951966E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.144649744464E-01
 v2rhoab       = -0.204638346911E-01
 v2rhob2       =  0.144649744464E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
 zk            = -0.101483720780E-01
 vrhoa         = -0.653289336535E-01
 vrhob         = -0.653289336535E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.171179303519E+00
 v2rhoab       = -0.263237442473E+00
 v2rhob2       =  0.171179303519E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
 zk            = -0.532741477023E+03
 vrhoa         = -0.157944671704E+00
 vrhob         = -0.157944671704E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.223669588268E-04
 v2rhoab       = -0.279504750065E-04
 v2rhob2       =  0.223669588268E-04

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
 zk            = -0.532741477023E+03
 vrhoa         = -0.157944671704E+00
 vrhob         = -0.157944671704E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.223669588268E-04
 v2rhoab       = -0.279504750065E-04
 v2rhob2       =  0.223669588268E-04

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
 zk            = -0.469795648279E+03
 vrhoa         = -0.156761418492E+00
 vrhob         = -0.156761418492E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.249624818738E-04
 v2rhoab       = -0.312385564000E-04
 v2rhob2       =  0.249624818738E-04

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
 zk            = -0.344300981310E-01
 vrhoa         = -0.743196778205E-01
 vrhob         = -0.743196778205E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.673517043262E-01
 v2rhoab       = -0.999958453856E-01
 v2rhob2       =  0.673517043262E-01

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
 zk            = -0.193016440590E+05
 vrhoa         = -0.192271893744E+00
 vrhob         = -0.192271893744E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.934975367062E-06
 v2rhoab       = -0.112791439004E-05
 v2rhob2       =  0.934975367062E-06

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
 zk            = -0.170016041806E+05
 vrhoa         = -0.191043581788E+00
 vrhob         = -0.191043581788E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.104726097228E-05
 v2rhoab       = -0.126473945017E-05
 v2rhob2       =  0.104726097228E-05

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
 zk            = -0.185432270230E-01
 vrhoa         = -0.697024933328E-01
 vrhob         = -0.697024933328E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.108356325631E+00
 v2rhoab       = -0.163679368941E+00
 v2rhob2       =  0.108356325631E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.149673920800E+00
 vrhoa         = -0.471789714951E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.130468710291E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.149673920800E+00
 vrhoa         = -0.471789714951E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.130468710291E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.126255646553E+00
 vrhoa         = -0.464766057364E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.151540949651E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.151940115037E-02
 vrhoa         = -0.297993118612E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.663179036789E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
 zk            = -0.191815088538E+02
 vrhoa         = -0.126816376070E+00
 vrhob         = -0.127719732599E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.386783435451E-03
 v2rhoab       = -0.511451913589E-03
 v2rhob2       =  0.397050991054E-03

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
 zk            = -0.851077910672E+01
 vrhoa         = -0.119099058995E+00
 vrhob         = -0.120906044904E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.756836181702E-03
 v2rhoab       = -0.102861281830E-02
 v2rhob2       =  0.800136175083E-03

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
 zk            = -0.132928938310E-01
 vrhoa         = -0.615913447654E-01
 vrhob         = -0.739125478228E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.963811005993E-01
 v2rhoab       = -0.210790984367E+00
 v2rhob2       =  0.193655360428E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
 zk            = -0.551636081757E-02
 vrhoa         = -0.482590171567E-01
 vrhob         = -0.811494055217E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.864399166736E-01
 v2rhoab       = -0.417632287561E+00
 v2rhob2       =  0.710228207528E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
 zk            = -0.111786110384E+02
 vrhoa         = -0.121711427226E+00
 vrhob         = -0.123144247973E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.605386382334E-03
 v2rhoab       = -0.814117737252E-03
 v2rhob2       =  0.632128050135E-03

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
 zk            = -0.887177858837E+01
 vrhoa         = -0.120365709995E+00
 vrhob         = -0.120365709995E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.751587360250E-03
 v2rhoab       = -0.992735036010E-03
 v2rhob2       =  0.751587360250E-03

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
 zk            = -0.130284832590E-01
 vrhoa         = -0.637104594188E-01
 vrhob         = -0.708674029459E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.114871866539E+00
 v2rhoab       = -0.215499991137E+00
 v2rhob2       =  0.172160006882E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
 zk            = -0.360439923610E-02
 vrhoa         = -0.488659901421E-01
 vrhob         = -0.708954619404E-01
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.164372114971E+00
 v2rhoab       = -0.574640089589E+00
 v2rhob2       =  0.724506769540E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00"""

def parsedata(dataname):
    import re
    from pyquante2.utils import parseline
    fields = []
    f = iter(dataname.splitlines())
    while True:
        try:
            line = f.next()
        except:
            break
        if line.lstrip().startswith('rhoa'):
            rhoa,rhob,sigaa,sigab,sigbb = parseline(line,'xfxfxfxfxf')
            field = [rhoa,rhob,sigaa,sigab,sigbb]
            fields.append(field)
            for i in range(6):
                line = f.next()
                term = parseline(line,'xxf')
                field.append(term)
        else:
            continue
    return np.array(fields)

if __name__ == '__main__':
    xlda = parsedata(x_lda_data)
    cvwn = parsedata(c_vwn5_data)
    print (xlda[:,:5] == cvwn[:,:5]).all()
