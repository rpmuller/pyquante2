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

x_b88_data="""\
 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
 zk            = -0.377592720836E+01
 vrhoa         = -0.148075576798E+01
 vrhob         = -0.148075576798E+01
 vsigmaaa      = -0.207006537839E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.207006537839E-02

 v2rhoa2       = -0.290344268232E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.290344268232E+00

 v2rhoasigmaaa =  0.162358068893E-02
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.162358068893E-02

 v2sigmaaa2    =  0.253445241319E-04
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.253445241319E-04

 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
 zk            = -0.378289713911E+01
 vrhoa         = -0.147807268065E+01
 vrhob         = -0.147807268065E+01
 vsigmaaa      = -0.203114756676E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.203114756676E-02

 v2rhoa2       = -0.293917869186E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.293917869186E+00

 v2rhoasigmaaa =  0.153738619102E-02
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.153738619102E-02

 v2sigmaaa2    =  0.208765215311E-04
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.208765215311E-04

 rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
 zk            = -0.334570325134E+01
 vrhoa         = -0.136817307919E+01
 vrhob         = -0.136817307919E+01
 vsigmaaa      = -0.185634599412E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.185634599412E-02

 v2rhoa2       = -0.371678164240E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.371678164240E+00

 v2rhoasigmaaa =  0.105687034683E-02
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.105687034683E-02

 v2sigmaaa2    =  0.926898408248E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.926898408248E-05

 rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
 zk            = -0.851611545044E-01
 vrhoa         = -0.501899165865E+00
 vrhob         = -0.501899165865E+00
 vsigmaaa      = -0.543404155466E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.543404155466E-01

 v2rhoa2       = -0.253321231850E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.253321231850E+01

 v2rhoasigmaaa =  0.239754146866E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.239754146866E+00

 v2sigmaaa2    =  0.221360010652E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.221360010652E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
 zk            = -0.407494475322E+05
 vrhoa         = -0.150923879747E+02
 vrhob         = -0.150923879747E+02
 vsigmaaa      = -0.191816494659E-06
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.191816494659E-06

 v2rhoa2       = -0.279488666210E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.279488666210E-02

 v2rhoasigmaaa =  0.142086292324E-09
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.142086292324E-09

 v2sigmaaa2    =  0.201646090402E-16
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.201646090402E-16

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
 zk            = -0.407494508312E+05
 vrhoa         = -0.150923867529E+02
 vrhob         = -0.150923867529E+02
 vsigmaaa      = -0.191816321255E-06
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.191816321255E-06

 v2rhoa2       = -0.279488824600E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.279488824600E-02

 v2rhoasigmaaa =  0.142085906983E-09
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.142085906983E-09

 v2sigmaaa2    =  0.201644008560E-16
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.201644008560E-16

 rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
 zk            = -0.362648637930E+05
 vrhoa         = -0.140333722784E+02
 vrhob         = -0.140333722784E+02
 vsigmaaa      = -0.174646643568E-06
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.174646643568E-06

 v2rhoa2       = -0.352244394438E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.352244394438E-02

 v2rhoasigmaaa =  0.971067113036E-10
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.971067113036E-10

 v2sigmaaa2    =  0.785386351399E-17
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.785386351399E-17

 rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
 zk            = -0.321148637763E+00
 vrhoa         = -0.766539815464E+00
 vrhob         = -0.766539815464E+00
 vsigmaaa      = -0.198197408319E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.198197408319E-01

 v2rhoa2       = -0.117949870779E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.117949870779E+01

 v2rhoasigmaaa =  0.685130252745E-01
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.685130252745E-01

 v2sigmaaa2    =  0.115351801846E-01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.115351801846E-01

 rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
 zk            = -0.370503980183E+07
 vrhoa         = -0.466042742267E+02
 vrhob         = -0.466042742267E+02
 vsigmaaa      = -0.210967131116E-08
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.210967131116E-08

 v2rhoa2       = -0.293108643192E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.293108643192E-03

 v2rhoasigmaaa =  0.530734919752E-13
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.530734919752E-13

 v2sigmaaa2    =  0.268272614874E-22
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.268272614874E-22

 rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
 zk            = -0.328152696735E+07
 vrhoa         = -0.433514250199E+02
 vrhob         = -0.433514250199E+02
 vsigmaaa      = -0.194182330561E-08
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.194182330561E-08

 v2rhoa2       = -0.368694289315E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.368694289315E-03

 v2rhoasigmaaa =  0.372175421272E-13
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.372175421272E-13

 v2sigmaaa2    =  0.108604300970E-22
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.108604300970E-22

 rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
 zk            = -0.161367392847E+00
 vrhoa         = -0.619947650806E+00
 vrhob         = -0.619947650806E+00
 vsigmaaa      = -0.341862053361E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.341862053361E-01

 v2rhoa2       = -0.179936512300E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.179936512300E+01

 v2rhoasigmaaa =  0.148255198864E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.148255198864E+00

 v2sigmaaa2    =  0.547109233247E-01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.547109233247E-01

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.494484233083E+01
 vrhoa         = -0.188374945936E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.790360507210E-03
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

 v2sigmaaa2    =  0.141061224494E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.494752158228E+01
 vrhoa         = -0.188273473796E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.785719874797E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.180074584533E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.132207800833E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.419401965265E+01
 vrhoa         = -0.172996985225E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.753968712720E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.231847776994E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.631038474677E-06
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.259998774808E-01
 vrhoa         = -0.428541964579E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.782798087404E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.312635380295E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.690680500539E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
 zk            = -0.740814331850E+03
 vrhoa         = -0.498247046111E+01
 vrhob         = -0.495532429311E+01
 vsigmaaa      = -0.678262141094E-05
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.682517937334E-05

 v2rhoa2       = -0.270432929879E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.272490060294E-01

 v2rhoasigmaaa =  0.426066228043E-07
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.425046982664E-07

 v2sigmaaa2    =  0.424725929436E-12
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.432961117719E-12

 rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
 zk            = -0.277987329958E+03
 vrhoa         = -0.385951846654E+01
 vrhob         = -0.381309494319E+01
 vsigmaaa      = -0.172434478018E-04
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.173712338362E-04

 v2rhoa2       = -0.441426807406E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.447245742260E-01

 v2rhoasigmaaa =  0.201415922856E-06
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.195961359539E-06

 v2sigmaaa2    =  0.700742719647E-11
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.718678968862E-11

 rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
 zk            = -0.120208982576E+00
 vrhoa         = -0.583637510333E+00
 vrhob         = -0.501672871724E+00
 vsigmaaa      = -0.379227871606E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.367802205908E-01

 v2rhoa2       = -0.199082326261E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.213690115573E+01

 v2rhoasigmaaa =  0.160652968405E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.609908850345E-01

 v2sigmaaa2    =  0.741970758035E-01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.737150455275E-01

 rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
 zk            = -0.416730804241E-01
 vrhoa         = -0.522700471817E+00
 vrhob         = -0.360524872778E+00
 vsigmaaa      = -0.111846105258E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.249411314940E+00

 v2rhoa2       = -0.245380935846E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.517371957290E+01

 v2rhoasigmaaa =  0.156984506289E+01
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.418857803825E+01

 v2sigmaaa2    =  0.244026452179E+01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.211148438264E+02

 rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
 zk            = -0.343037899309E+03
 vrhoa         = -0.451375531283E+01
 vrhob         = -0.448071957650E+01
 vsigmaaa      = -0.204625762085E-04
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.209266539812E-04

 v2rhoa2       = -0.327661718483E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.333090761894E-01

 v2rhoasigmaaa =  0.455875527750E-06
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.472402981076E-06

 v2sigmaaa2    =  0.153056541719E-10
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.162083837575E-10

 rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
 zk            = -0.260133861611E+03
 vrhoa         = -0.416254426190E+01
 vrhob         = -0.416322526434E+01
 vsigmaaa      = -0.262815715798E-04
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.263113502629E-04

 v2rhoa2       = -0.391759770307E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.391492164486E-01

 v2rhoasigmaaa =  0.680016391145E-06
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.681990700893E-06

 v2sigmaaa2    =  0.297024234619E-10
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.298552512045E-10

 rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
 zk            = -0.112603177863E+00
 vrhoa         = -0.568498809942E+00
 vrhob         = -0.520860206545E+00
 vsigmaaa      = -0.423139064219E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.436372569148E-01

 v2rhoa2       = -0.209994430549E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.229642938280E+01

 v2rhoasigmaaa =  0.195292854203E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.150061257405E+00

 v2sigmaaa2    =  0.103073123098E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.115652366462E+00

 rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
 zk            = -0.253983617946E-01
 vrhoa         = -0.431645138108E+00
 vrhob         = -0.325991019475E+00
 vsigmaaa      = -0.175455092294E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.260075409570E+00

 v2rhoa2       = -0.374315005133E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.566470285703E+01

 v2rhoasigmaaa =  0.291762144791E+01
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.301407606104E+01

 v2sigmaaa2    =  0.765442610536E+01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.248245711494E+02"""

x_pbe_data="""\
 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
 zk            = -0.377592720836E+01
 vrhoa         = -0.148075576798E+01
 vrhob         = -0.148075576798E+01
 vsigmaaa      = -0.165665974842E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.165665974842E-02

 v2rhoa2       = -0.290344268232E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.290344268232E+00

 v2rhoasigmaaa =  0.129934097915E-02
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.129934097915E-02

 v2sigmaaa2    =  0.361615443758E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.361615443758E-05

 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
 zk            = -0.378154942017E+01
 vrhoa         = -0.147855914532E+01
 vrhob         = -0.147855914532E+01
 vsigmaaa      = -0.165052935245E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.165052935245E-02

 v2rhoa2       = -0.293340073168E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.293340073168E+00

 v2rhoasigmaaa =  0.128494322308E-02
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.128494322308E-02

 v2sigmaaa2    =  0.359610088925E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.359610088925E-05

 rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
 zk            = -0.332917118617E+01
 vrhoa         = -0.136704102604E+01
 vrhob         = -0.136704102604E+01
 vsigmaaa      = -0.175922831660E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.175922831660E-02

 v2rhoa2       = -0.383048644703E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.383048644703E+00

 v2rhoasigmaaa =  0.123846484420E-02
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.123846484420E-02

 v2sigmaaa2    =  0.508271342885E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.508271342885E-05

 rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
 zk            = -0.847500738867E-01
 vrhoa         = -0.498335317577E+00
 vrhob         = -0.498335317577E+00
 vsigmaaa      = -0.545109539268E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.545109539268E-01

 v2rhoa2       = -0.229469146464E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.229469146464E+01

 v2rhoasigmaaa =  0.154401191219E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.154401191219E+00

 v2sigmaaa2    =  0.254715375324E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.254715375324E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
 zk            = -0.407494475322E+05
 vrhoa         = -0.150923879748E+02
 vrhob         = -0.150923879748E+02
 vsigmaaa      = -0.153509482897E-06
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.153509482897E-06

 v2rhoa2       = -0.279488666208E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.279488666208E-02

 v2rhoasigmaaa =  0.113710728070E-09
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.113710728070E-09

 v2sigmaaa2    =  0.287708461721E-17
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.287708461721E-17

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
 zk            = -0.407494501724E+05
 vrhoa         = -0.150923869969E+02
 vrhob         = -0.150923869969E+02
 vsigmaaa      = -0.153509458156E-06
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.153509458156E-06

 v2rhoa2       = -0.279488792967E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.279488792967E-02

 v2rhoasigmaaa =  0.113710673089E-09
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.113710673089E-09

 v2sigmaaa2    =  0.287708392165E-17
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.287708392165E-17

 rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
 zk            = -0.360960910073E+05
 vrhoa         = -0.140305340798E+02
 vrhob         = -0.140305340798E+02
 vsigmaaa      = -0.163703325386E-06
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.163703325386E-06

 v2rhoa2       = -0.361179088671E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.361179088671E-02

 v2rhoasigmaaa =  0.111691289742E-09
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.111691289742E-09

 v2sigmaaa2    =  0.400996995891E-17
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.400996995891E-17

 rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
 zk            = -0.319679717401E+00
 vrhoa         = -0.766494428708E+00
 vrhob         = -0.766494428708E+00
 vsigmaaa      = -0.185240091115E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.185240091115E-01

 v2rhoa2       = -0.120781973536E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.120781973536E+01

 v2rhoasigmaaa =  0.783950736036E-01
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.783950736036E-01

 v2sigmaaa2    =  0.578030314073E-02
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.578030314073E-02

 rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
 zk            = -0.370503980175E+07
 vrhoa         = -0.466042742277E+02
 vrhob         = -0.466042742277E+02
 vsigmaaa      = -0.168835611841E-08
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.168835611841E-08

 v2rhoa2       = -0.293108643147E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.293108643147E-03

 v2rhoasigmaaa =  0.424743677402E-13
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.424743677402E-13

 v2sigmaaa2    =  0.382771132308E-23
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.382771132308E-23

 rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
 zk            = -0.326657718168E+07
 vrhoa         = -0.433502031071E+02
 vrhob         = -0.433502031071E+02
 vsigmaaa      = -0.181368847096E-08
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.181368847096E-08

 v2rhoa2       = -0.377435398713E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.377435398713E-03

 v2rhoasigmaaa =  0.425352935840E-13
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.425352935840E-13

 v2sigmaaa2    =  0.541929951945E-23
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.541929951945E-23

 rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
 zk            = -0.160448304498E+00
 vrhoa         = -0.616293740752E+00
 vrhob         = -0.616293740752E+00
 vsigmaaa      = -0.340347075448E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.340347075448E-01

 v2rhoa2       = -0.188013619179E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.188013619179E+01

 v2rhoasigmaaa =  0.179505894965E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.179505894965E+00

 v2sigmaaa2    =  0.432509198788E-01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.432509198788E-01

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.494484233083E+01
 vrhoa         = -0.188374945936E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.632520331341E-03
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

 v2sigmaaa2    =  0.201266029008E-06
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.494699173727E+01
 vrhoa         = -0.188293152324E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.631836581688E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.179948990030E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.200939766657E-06
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.417440506676E+01
 vrhoa         = -0.172954529118E+01
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.707320404038E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.237820356120E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.323948969032E-06
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.258503716227E-01
 vrhoa         = -0.433574571158E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      = -0.743604018862E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.182332674956E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    =  0.934224897782E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
 zk            = -0.736865469611E+03
 vrhoa         = -0.493897932478E+01
 vrhob         = -0.491161504563E+01
 vsigmaaa      = -0.685185279408E-05
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.689688466662E-05

 v2rhoa2       = -0.273332629968E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.273922219983E-01

 v2rhoasigmaaa =  0.455358003173E-07
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.445075229580E-07

 v2sigmaaa2    =  0.413408272915E-12
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.427862514991E-12

 rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
 zk            = -0.276589791995E+03
 vrhoa         = -0.382556082420E+01
 vrhob         = -0.378108116179E+01
 vsigmaaa      = -0.174145337536E-04
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.175120610339E-04

 v2rhoa2       = -0.429564214817E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.424802511645E-01

 v2rhoasigmaaa =  0.185237729809E-06
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.161839553501E-06

 v2sigmaaa2    =  0.740514207206E-11
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.786563034093E-11

 rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
 zk            = -0.119421209325E+00
 vrhoa         = -0.578832621965E+00
 vrhob         = -0.517389121001E+00
 vsigmaaa      = -0.382340404759E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.332973674309E-01

 v2rhoa2       = -0.204167554105E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.722156009564E+00

 v2rhoasigmaaa =  0.181184032538E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb = -0.177031055838E+00

 v2sigmaaa2    =  0.685619910117E-01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.104342795817E+00

 rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
 zk            = -0.415413960050E-01
 vrhoa         = -0.523382772310E+00
 vrhob         = -0.357399049244E+00
 vsigmaaa      = -0.975929514687E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.251904684435E+00

 v2rhoa2       = -0.245525965200E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.524279284728E+01

 v2rhoasigmaaa =  0.155938982496E+01
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.452016269281E+01

 v2sigmaaa2    =  0.776664232752E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.203904030316E+02

 rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
 zk            = -0.342262372442E+03
 vrhoa         = -0.451953851009E+01
 vrhob         = -0.448651692914E+01
 vsigmaaa      = -0.177274055848E-04
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.181825088154E-04

 v2rhoa2       = -0.327443513253E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.333030095027E-01

 v2rhoasigmaaa =  0.445584305337E-06
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.464801667902E-06

 v2sigmaaa2    =  0.462724606680E-11
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.500476145449E-11

 rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
 zk            = -0.259316875633E+03
 vrhoa         = -0.416761407365E+01
 vrhob         = -0.416832346828E+01
 vsigmaaa      = -0.234431309174E-04
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.234541298631E-04

 v2rhoa2       = -0.394086605352E-01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.393741681091E-01

 v2rhoasigmaaa =  0.708870190119E-06
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.709919266338E-06

 v2sigmaaa2    =  0.109950768377E-10
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.110028156783E-10

 rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
 zk            = -0.112021068338E+00
 vrhoa         = -0.563849973546E+00
 vrhob         = -0.518540640509E+00
 vsigmaaa      = -0.426508442574E-01
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.433850716387E-01

 v2rhoa2       = -0.215582132758E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.193599570809E+01

 v2rhoasigmaaa =  0.221089275399E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.555875642118E-01

 v2sigmaaa2    =  0.948033727977E-01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.140057158296E+00

 rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
 zk            = -0.252665336863E-01
 vrhoa         = -0.430778167566E+00
 vrhob         = -0.326346959743E+00
 vsigmaaa      = -0.168954851377E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      = -0.254284088546E+00

 v2rhoa2       = -0.388578815717E+01
 v2rhoab       =  0.000000000000E+00
 v2rhob2       = -0.420071526691E+01

 v2rhoasigmaaa =  0.349932953984E+01
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb = -0.344308375753E+00

 v2sigmaaa2    =  0.467162912424E+01
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.317975451941E+02"""

c_lyp_data = """\
rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
 zk            = -0.179175399535E+00
 vrhoa         = -0.567254370239E-01
 vrhob         = -0.567254370239E-01
 vsigmaaa      =  0.603063052247E-04
 vsigmaab      =  0.562668577012E-04
 vsigmabb      =  0.603063052247E-04

 v2rhoa2       =  0.133393844920E-01
 v2rhoab       = -0.152396936370E-01
 v2rhob2       =  0.133393844920E-01

 v2rhoasigmaaa = -0.139567528486E-03
 v2rhoasigmaab = -0.291609267169E-04
 v2rhoasigmabb =  0.811155700026E-04
 v2rhobsigmaaa =  0.811155700026E-04
 v2rhobsigmaab = -0.291609267169E-04
 v2rhobsigmabb = -0.139567528486E-03

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
 zk            = -0.178874704439E+00
 vrhoa         = -0.568743789287E-01
 vrhob         = -0.568743789287E-01
 vsigmaaa      =  0.603063052247E-04
 vsigmaab      =  0.562668577012E-04
 vsigmabb      =  0.603063052247E-04

 v2rhoa2       =  0.138019299001E-01
 v2rhoab       = -0.154679607816E-01
 v2rhob2       =  0.138019299001E-01

 v2rhoasigmaaa = -0.139567528486E-03
 v2rhoasigmaab = -0.291609267169E-04
 v2rhoasigmabb =  0.811155700026E-04
 v2rhobsigmaaa =  0.811155700026E-04
 v2rhobsigmaab = -0.291609267169E-04
 v2rhobsigmabb = -0.139567528486E-03

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
 zk            = -0.148704448499E+00
 vrhoa         = -0.607218516385E-01
 vrhob         = -0.607218516385E-01
 vsigmaaa      =  0.741080708899E-04
 vsigmaab      =  0.701383136250E-04
 vsigmabb      =  0.741080708899E-04

 v2rhoa2       =  0.303235212020E-01
 v2rhoab       = -0.246964871909E-01
 v2rhob2       =  0.303235212020E-01

 v2rhoasigmaaa = -0.192341863313E-03
 v2rhoasigmaab = -0.411266458866E-04
 v2rhoasigmabb =  0.111052541885E-03
 v2rhobsigmaaa =  0.111052541885E-03
 v2rhobsigmaab = -0.411266458866E-04
 v2rhobsigmabb = -0.192341863313E-03

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
 zk            = -0.465024100803E-02
 vrhoa         = -0.610987646298E-01
 vrhob         = -0.610987646298E-01
 vsigmaaa      =  0.694062751783E-02
 vsigmaab      =  0.876388739167E-02
 vsigmabb      =  0.694062751783E-02

 v2rhoa2       =  0.939780887292E+00
 v2rhoab       = -0.498868805878E+00
 v2rhob2       =  0.939780887292E+00

 v2rhoasigmaaa = -0.228361811495E+00
 v2rhoasigmaab = -0.807965355076E-01
 v2rhoasigmabb =  0.107205512562E+00
 v2rhobsigmaaa =  0.107205512562E+00
 v2rhobsigmaab = -0.807965355076E-01
 v2rhobsigmabb = -0.228361811495E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
 zk            = -0.237638430952E+03
 vrhoa         = -0.665993271139E-01
 vrhob         = -0.665993271139E-01
 vsigmaaa      =  0.540555769454E-09
 vsigmaab      =  0.260773591409E-09
 vsigmabb      =  0.540555769454E-09

 v2rhoa2       =  0.144305866044E-04
 v2rhoab       = -0.146537621392E-04
 v2rhob2       =  0.144305866044E-04

 v2rhoasigmaaa = -0.158850198364E-11
 v2rhoasigmaab = -0.124882367057E-12
 v2rhoasigmabb =  0.108678004324E-11
 v2rhobsigmaaa =  0.108678004324E-11
 v2rhobsigmaab = -0.124882367057E-12
 v2rhobsigmabb = -0.158850198364E-11

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
 zk            = -0.237638419413E+03
 vrhoa         = -0.665993325023E-01
 vrhob         = -0.665993325023E-01
 vsigmaaa      =  0.540555769454E-09
 vsigmaab      =  0.260773591409E-09
 vsigmabb      =  0.540555769454E-09

 v2rhoa2       =  0.144306117017E-04
 v2rhoab       = -0.146537792037E-04
 v2rhob2       =  0.144306117017E-04

 v2rhoasigmaaa = -0.158850198364E-11
 v2rhoasigmaab = -0.124882367057E-12
 v2rhoasigmabb =  0.108678004324E-11
 v2rhobsigmaaa =  0.108678004324E-11
 v2rhobsigmaab = -0.124882367057E-12
 v2rhobsigmabb = -0.158850198364E-11

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
 zk            = -0.204955819468E+03
 vrhoa         = -0.697313972428E-01
 vrhob         = -0.697313972428E-01
 vsigmaaa      =  0.658117255270E-09
 vsigmaab      =  0.319515846955E-09
 vsigmabb      =  0.658117255270E-09

 v2rhoa2       =  0.328352356048E-04
 v2rhoab       = -0.277614126197E-04
 v2rhob2       =  0.328352356048E-04

 v2rhoasigmaaa = -0.217196289305E-11
 v2rhoasigmaab = -0.172301824814E-12
 v2rhoasigmabb =  0.148472344173E-11
 v2rhobsigmaaa =  0.148472344173E-11
 v2rhobsigmaab = -0.172301824814E-12
 v2rhobsigmabb = -0.217196289305E-11

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
 zk            = -0.216471443930E-01
 vrhoa         = -0.526995644881E-01
 vrhob         = -0.526995644881E-01
 vsigmaaa      =  0.127229458316E-02
 vsigmaab      =  0.145290124256E-02
 vsigmabb      =  0.127229458316E-02

 v2rhoa2       =  0.134806861660E+00
 v2rhoab       = -0.118221268828E+00
 v2rhob2       =  0.134806861660E+00

 v2rhoasigmaaa = -0.160681338920E-01
 v2rhoasigmaab = -0.472571671370E-02
 v2rhoasigmabb =  0.827564425842E-02
 v2rhobsigmaaa =  0.827564425842E-02
 v2rhobsigmaab = -0.472571671370E-02
 v2rhobsigmabb = -0.160681338920E-01

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
 zk            = -0.712575807209E+04
 vrhoa         = -0.674210020682E-01
 vrhob         = -0.674210020682E-01
 vsigmaaa      =  0.190814948140E-11
 vsigmaab      =  0.816481365419E-12
 vsigmabb      =  0.190814948140E-11

 v2rhoa2       =  0.495469678644E-06
 v2rhoab       = -0.497965254694E-06
 v2rhob2       =  0.495469678644E-06

 v2rhoasigmaaa = -0.196249114596E-15
 v2rhoasigmaab = -0.130135394800E-16
 v2rhoasigmabb =  0.136188177073E-15
 v2rhobsigmaaa =  0.136188177073E-15
 v2rhobsigmaab = -0.130135394800E-16
 v2rhobsigmabb = -0.196249114596E-15

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
 zk            = -0.615255831041E+04
 vrhoa         = -0.703265205307E-01
 vrhob         = -0.703265205307E-01
 vsigmaaa      =  0.233162530668E-11
 vsigmaab      =  0.100028081622E-11
 vsigmabb      =  0.233162530668E-11

 v2rhoa2       =  0.110943737062E-05
 v2rhoab       = -0.946187279564E-06
 v2rhob2       =  0.110943737062E-05

 v2rhoasigmaaa = -0.270252437366E-15
 v2rhoasigmaab = -0.179871585341E-16
 v2rhoasigmabb =  0.187490309515E-15
 v2rhobsigmaaa =  0.187490309515E-15
 v2rhobsigmaab = -0.179871585341E-16
 v2rhobsigmabb = -0.270252437366E-15

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
 zk            = -0.106605459023E-01
 vrhoa         = -0.547299644084E-01
 vrhob         = -0.547299644084E-01
 vsigmaaa      =  0.303386607339E-02
 vsigmaab      =  0.365209975632E-02
 vsigmabb      =  0.303386607339E-02

 v2rhoa2       =  0.346499910453E+00
 v2rhoab       = -0.239377906462E+00
 v2rhob2       =  0.346499910453E+00

 v2rhoasigmaaa = -0.624221561613E-01
 v2rhoasigmaab = -0.201987433700E-01
 v2rhoasigmabb =  0.307385926277E-01
 v2rhobsigmaaa =  0.307385926277E-01
 v2rhobsigmaab = -0.201987433700E-01
 v2rhobsigmabb = -0.624221561613E-01

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            =  0.000000000000E+00
 vrhoa         =  0.000000000000E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.000000000000E+00
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
 zk            =  0.000000000000E+00
 vrhoa         =  0.000000000000E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.000000000000E+00
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
 zk            =  0.000000000000E+00
 vrhoa         =  0.000000000000E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.000000000000E+00
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
 zk            =  0.000000000000E+00
 vrhoa         =  0.000000000000E+00
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.000000000000E+00
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       =  0.000000000000E+00
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
 zk            = -0.903966286601E+01
 vrhoa         = -0.759264713358E-01
 vrhob         = -0.784370481382E-01
 vsigmaaa      =  0.907992031805E-07
 vsigmaab      =  0.581118313258E-07
 vsigmabb      =  0.100370432217E-06

 v2rhoa2       =  0.141362773226E-02
 v2rhoab       = -0.105147646057E-02
 v2rhob2       =  0.150485268125E-02

 v2rhoasigmaaa = -0.565167905921E-08
 v2rhoasigmaab = -0.515378386013E-09
 v2rhoasigmabb =  0.376492504615E-08
 v2rhobsigmaaa =  0.384361243685E-08
 v2rhobsigmaab = -0.747058939790E-09
 v2rhobsigmabb = -0.588190678666E-08

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
 zk            = -0.402158795173E+01
 vrhoa         = -0.762734644914E-01
 vrhob         = -0.830226435821E-01
 vsigmaaa      =  0.301052145436E-06
 vsigmaab      =  0.220298633297E-06
 vsigmabb      =  0.369624286402E-06

 v2rhoa2       =  0.331769729999E-02
 v2rhoab       = -0.248438749270E-02
 v2rhob2       =  0.384280348438E-02

 v2rhoasigmaaa = -0.398359773843E-07
 v2rhoasigmaab = -0.335415277613E-08
 v2rhoasigmabb =  0.263970784129E-07
 v2rhobsigmaaa =  0.275886078235E-07
 v2rhobsigmaab = -0.685474898360E-08
 v2rhobsigmabb = -0.433118929134E-07

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
 zk            = -0.535812609335E-02
 vrhoa         = -0.428946532738E-01
 vrhob         = -0.100396349615E+00
 vsigmaaa      =  0.178546727213E-02
 vsigmaab      =  0.678262565697E-02
 vsigmabb      =  0.827650024238E-02

 v2rhoa2       =  0.309477814856E+00
 v2rhoab       = -0.329371433841E+00
 v2rhob2       =  0.162055387526E+01

 v2rhoasigmaaa = -0.781343180939E-01
 v2rhoasigmaab = -0.432850540037E-02
 v2rhoasigmabb =  0.363514849034E-01
 v2rhobsigmaaa =  0.732832606483E-01
 v2rhobsigmaab = -0.108581862510E+00
 v2rhobsigmabb = -0.179156926286E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
 zk            = -0.303154148293E-02
 vrhoa         = -0.244020499979E-01
 vrhob         = -0.783890652209E-01
 vsigmaaa      = -0.352365411226E-02
 vsigmaab      =  0.373879285446E-01
 vsigmabb      =  0.434388664083E-01

 v2rhoa2       =  0.202330890903E+00
 v2rhoab       = -0.536912334084E+00
 v2rhob2       =  0.182166278214E+01

 v2rhoasigmaaa = -0.992049673308E-01
 v2rhoasigmaab = -0.485780496011E-01
 v2rhoasigmabb = -0.158164241941E+00
 v2rhobsigmaaa =  0.340893479361E+00
 v2rhobsigmaab = -0.162192597781E+01
 v2rhobsigmabb = -0.154746697775E+01

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
 zk            = -0.608750998131E+01
 vrhoa         = -0.643435848234E-01
 vrhob         = -0.656377677365E-01
 vsigmaaa      =  0.202448437183E-06
 vsigmaab      =  0.140789565804E-06
 vsigmabb      =  0.237924154032E-06

 v2rhoa2       =  0.644965086508E-03
 v2rhoab       = -0.630065421347E-03
 v2rhob2       =  0.684494690431E-03

 v2rhoasigmaaa = -0.207590458741E-07
 v2rhoasigmaab = -0.182473154038E-08
 v2rhoasigmabb =  0.137746587331E-07
 v2rhobsigmaaa =  0.142543841876E-07
 v2rhobsigmaab = -0.323521094196E-08
 v2rhobsigmabb = -0.221599858116E-07

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
 zk            = -0.485826625474E+01
 vrhoa         = -0.653859880589E-01
 vrhob         = -0.653239498026E-01
 vsigmaaa      =  0.314305582522E-06
 vsigmaab      =  0.205017935750E-06
 vsigmabb      =  0.314305582522E-06

 v2rhoa2       =  0.932011907377E-03
 v2rhoab       = -0.845990806989E-03
 v2rhob2       =  0.927946564209E-03

 v2rhoasigmaaa = -0.375898792310E-07
 v2rhoasigmaab = -0.455264241776E-08
 v2rhoasigmabb =  0.244483769850E-07
 v2rhobsigmaaa =  0.244483769850E-07
 v2rhobsigmaab = -0.455264241776E-08
 v2rhobsigmabb = -0.375898792310E-07

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
 zk            = -0.634688010938E-02
 vrhoa         = -0.486068357646E-01
 vrhob         = -0.739949355207E-01
 vsigmaaa      =  0.305337986175E-02
 vsigmaab      =  0.641200968640E-02
 vsigmabb      =  0.698202739086E-02

 v2rhoa2       =  0.433585531118E+00
 v2rhoab       = -0.365108886112E+00
 v2rhob2       =  0.100038106732E+01

 v2rhoasigmaaa = -0.102735565156E+00
 v2rhoasigmaab = -0.181538750084E-01
 v2rhoasigmabb =  0.501849138748E-01
 v2rhobsigmaaa =  0.732361795124E-01
 v2rhobsigmaab = -0.824126803195E-01
 v2rhobsigmabb = -0.165270268130E+00

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
 zk            = -0.172301075975E-02
 vrhoa         = -0.294822900234E-01
 vrhob         = -0.766946918598E-01
 vsigmaaa      =  0.127140367254E-02
 vsigmaab      =  0.519446045760E-01
 vsigmabb      =  0.608118658277E-01

 v2rhoa2       =  0.427000673980E+00
 v2rhoab       = -0.671140471911E+00
 v2rhob2       =  0.300158916489E+01

 v2rhoasigmaaa = -0.631566175836E+00
 v2rhoasigmaab = -0.540623931547E-01
 v2rhoasigmabb = -0.252818541902E-02
 v2rhobsigmaaa =  0.988123041938E+00
 v2rhobsigmaab = -0.292607572225E+01
 v2rhobsigmabb = -0.333300175870E+01

 v2sigmaaa2    =  0.000000000000E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00"""

c_pbe_data = """\
 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
 zk            = -0.277343302026E+00
 vrhoa         = -0.902545684170E-01
 vrhob         = -0.902545684170E-01
 vsigmaaa      =  0.828329874208E-03
 vsigmaab      =  0.165665974842E-02
 vsigmabb      =  0.828329874208E-03

 v2rhoa2       =  0.129726391337E-01
 v2rhoab       = -0.182703394917E-01
 v2rhob2       =  0.129726391337E-01

 v2rhoasigmaaa = -0.324835244787E-03
 v2rhoasigmaab = -0.649670489575E-03
 v2rhoasigmabb = -0.324835244787E-03
 v2rhobsigmaaa = -0.324835244787E-03
 v2rhobsigmaab = -0.649670489575E-03
 v2rhobsigmabb = -0.324835244787E-03

 v2sigmaaa2    = -0.649078703127E-05
 v2sigmaaaab   = -0.129815740625E-04
 v2sigmaaabb   = -0.649078703127E-05
 v2sigmaab2    = -0.259631481251E-04
 v2sigmaabbb   = -0.129815740625E-04
 v2sigmabb2    = -0.649078703127E-05

 rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
 zk            = -0.271855691853E+00
 vrhoa         = -0.923103473041E-01
 vrhob         = -0.923103473041E-01
 vsigmaaa      =  0.786385334368E-03
 vsigmaab      =  0.157277066874E-02
 vsigmabb      =  0.786385334368E-03

 v2rhoa2       =  0.141731701968E-01
 v2rhoab       = -0.168480493572E-01
 v2rhob2       =  0.141731701968E-01

 v2rhoasigmaaa = -0.281023204875E-03
 v2rhoasigmaab = -0.562046409750E-03
 v2rhoasigmabb = -0.281023204875E-03
 v2rhobsigmaaa = -0.281023204875E-03
 v2rhobsigmaab = -0.562046409750E-03
 v2rhobsigmabb = -0.281023204875E-03

 v2sigmaaa2    = -0.586202682069E-05
 v2sigmaaaab   = -0.117240536414E-04
 v2sigmaaabb   = -0.586202682069E-05
 v2sigmaab2    = -0.234481072827E-04
 v2sigmaabbb   = -0.117240536414E-04
 v2sigmabb2    = -0.586202682069E-05

 rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
 zk            = -0.156330536080E+00
 vrhoa         = -0.102948996370E+00
 vrhob         = -0.102948996370E+00
 vsigmaaa      =  0.378004824716E-03
 vsigmaab      =  0.756009649433E-03
 vsigmabb      =  0.378004824716E-03

 v2rhoa2       =  0.767829363534E-02
 v2rhoab       = -0.200905483619E-01
 v2rhob2       =  0.767829363534E-02

 v2rhoasigmaaa =  0.267671090070E-04
 v2rhoasigmaab =  0.535342180140E-04
 v2rhoasigmabb =  0.267671090070E-04
 v2rhobsigmaaa =  0.267671090070E-04
 v2rhobsigmaab =  0.535342180140E-04
 v2rhobsigmabb =  0.267671090070E-04

 v2sigmaaa2    = -0.171987026880E-05
 v2sigmaaaab   = -0.343974053761E-05
 v2sigmaaabb   = -0.171987026880E-05
 v2sigmaab2    = -0.687948107522E-05
 v2sigmaabbb   = -0.343974053761E-05
 v2sigmabb2    = -0.171987026880E-05

 rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
 zk            = -0.353114293615E-02
 vrhoa         = -0.644123523321E-01
 vrhob         = -0.644123523321E-01
 vsigmaaa      =  0.830804942063E-02
 vsigmaab      =  0.166160988413E-01
 vsigmabb      =  0.830804942063E-02

 v2rhoa2       = -0.294409365638E+00
 v2rhoab       = -0.534985274786E+00
 v2rhob2       = -0.294409365638E+00

 v2rhoasigmaaa =  0.746427566937E-01
 v2rhoasigmaab =  0.149285513387E+00
 v2rhoasigmabb =  0.746427566937E-01
 v2rhobsigmaaa =  0.746427566937E-01
 v2rhobsigmaab =  0.149285513387E+00
 v2rhobsigmabb =  0.746427566937E-01

 v2sigmaaa2    = -0.283776069209E-01
 v2sigmaaaab   = -0.567552138418E-01
 v2sigmaaabb   = -0.283776069209E-01
 v2sigmaab2    = -0.113510427684E+00
 v2sigmaabbb   = -0.567552138418E-01
 v2sigmabb2    = -0.283776069209E-01

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
 zk            = -0.531558156901E+03
 vrhoa         = -0.157670707589E+00
 vrhob         = -0.157670707589E+00
 vsigmaaa      =  0.767547413336E-07
 vsigmaab      =  0.153509482667E-06
 vsigmabb      =  0.767547413336E-07

 v2rhoa2       =  0.222158431130E-04
 v2rhoab       = -0.278288704562E-04
 v2rhob2       =  0.222158431130E-04

 v2rhoasigmaaa = -0.284276819003E-10
 v2rhoasigmaab = -0.568553638007E-10
 v2rhoasigmabb = -0.284276819003E-10
 v2rhobsigmaaa = -0.284276819003E-10
 v2rhobsigmaab = -0.568553638007E-10
 v2rhobsigmabb = -0.284276819003E-10

 v2sigmaaa2    = -0.526353544028E-16
 v2sigmaaaab   = -0.105270708806E-15
 v2sigmaaabb   = -0.526353544028E-16
 v2sigmaab2    = -0.210541417611E-15
 v2sigmaabbb   = -0.105270708806E-15
 v2sigmabb2    = -0.526353544028E-16

 rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
 zk            = -0.531555516738E+03
 vrhoa         = -0.157671685407E+00
 vrhob         = -0.157671685407E+00
 vsigmaaa      =  0.767529308359E-07
 vsigmaab      =  0.153505861672E-06
 vsigmabb      =  0.767529308359E-07

 v2rhoa2       =  0.222164315920E-04
 v2rhoab       = -0.278281914345E-04
 v2rhob2       =  0.222164315920E-04

 v2rhoasigmaaa = -0.284258379026E-10
 v2rhoasigmaab = -0.568516758052E-10
 v2rhoasigmabb = -0.284258379026E-10
 v2rhobsigmaaa = -0.284258379026E-10
 v2rhobsigmaab = -0.568516758052E-10
 v2rhobsigmabb = -0.284258379026E-10

 v2sigmaaa2    = -0.526328718655E-16
 v2sigmaaaab   = -0.105265743731E-15
 v2sigmaaabb   = -0.526328718655E-16
 v2sigmaab2    = -0.210531487462E-15
 v2sigmaabbb   = -0.105265743731E-15
 v2sigmabb2    = -0.526328718655E-16

 rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
 zk            = -0.204845328577E+03
 vrhoa         = -0.139570172131E+00
 vrhob         = -0.139570172131E+00
 vsigmaaa      =  0.609549582413E-08
 vsigmaab      =  0.121909916483E-07
 vsigmabb      =  0.609549582413E-08

 v2rhoa2       = -0.845856478622E-05
 v2rhoab       = -0.379435495181E-04
 v2rhob2       = -0.845856478622E-05

 v2rhoasigmaaa =  0.183107776239E-11
 v2rhoasigmaab =  0.366215552479E-11
 v2rhoasigmabb =  0.183107776239E-11
 v2rhobsigmaaa =  0.183107776239E-11
 v2rhobsigmaab =  0.366215552479E-11
 v2rhobsigmabb =  0.183107776239E-11

 v2sigmaaa2    = -0.402310862302E-18
 v2sigmaaaab   = -0.804621724604E-18
 v2sigmaaabb   = -0.402310862302E-18
 v2sigmaab2    = -0.160924344921E-17
 v2sigmaabbb   = -0.804621724604E-18
 v2sigmabb2    = -0.402310862302E-18

 rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
 zk            = -0.257202574589E-01
 vrhoa         = -0.866997760784E-01
 vrhob         = -0.866997760784E-01
 vsigmaaa      =  0.582853357983E-02
 vsigmaab      =  0.116570671597E-01
 vsigmabb      =  0.582853357983E-02

 v2rhoa2       =  0.728321951048E-01
 v2rhoab       = -0.756102345474E-01
 v2rhob2       =  0.728321951048E-01

 v2rhoasigmaaa = -0.305698872548E-02
 v2rhoasigmaab = -0.611397745095E-02
 v2rhoasigmabb = -0.305698872548E-02
 v2rhobsigmaaa = -0.305698872548E-02
 v2rhobsigmaab = -0.611397745095E-02
 v2rhobsigmabb = -0.305698872548E-02

 v2sigmaaa2    = -0.233980425019E-02
 v2sigmaaaab   = -0.467960850038E-02
 v2sigmaaabb   = -0.233980425019E-02
 v2sigmaab2    = -0.935921700077E-02
 v2sigmaabbb   = -0.467960850038E-02
 v2sigmabb2    = -0.233980425019E-02

 rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
 zk            = -0.192839940466E+05
 vrhoa         = -0.192144028270E+00
 vrhob         = -0.192144028270E+00
 vsigmaaa      =  0.844177976352E-09
 vsigmaab      =  0.168835595270E-08
 vsigmabb      =  0.844177976352E-09

 v2rhoa2       =  0.930079636518E-06
 v2rhoab       = -0.112362987226E-05
 v2rhob2       =  0.930079636518E-06

 v2rhoasigmaaa = -0.106185890697E-13
 v2rhoasigmaab = -0.212371781394E-13
 v2rhoasigmabb = -0.106185890697E-13
 v2rhobsigmaaa = -0.106185890697E-13
 v2rhobsigmaab = -0.212371781394E-13
 v2rhobsigmabb = -0.106185890697E-13

 v2sigmaaa2    = -0.216237891268E-21
 v2sigmaaaab   = -0.432475782536E-21
 v2sigmaaabb   = -0.216237891268E-21
 v2sigmaab2    = -0.864951565071E-21
 v2sigmaabbb   = -0.432475782536E-21
 v2sigmabb2    = -0.216237891268E-21

 rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
 zk            = -0.621606173598E+04
 vrhoa         = -0.145186193191E+00
 vrhob         = -0.145186193191E+00
 vsigmaaa      =  0.239893933515E-10
 vsigmaab      =  0.479787867030E-10
 vsigmabb      =  0.239893933515E-10

 v2rhoa2       = -0.415550006172E-06
 v2rhoab       = -0.137798903649E-05
 v2rhob2       = -0.415550006172E-06

 v2rhoasigmaaa =  0.270195164309E-15
 v2rhoasigmaab =  0.540390328617E-15
 v2rhoasigmabb =  0.270195164309E-15
 v2rhobsigmaaa =  0.270195164309E-15
 v2rhobsigmaab =  0.540390328617E-15
 v2rhobsigmabb =  0.270195164309E-15

 v2sigmaaa2    = -0.210724100250E-24
 v2sigmaaaab   = -0.421448200501E-24
 v2sigmaaabb   = -0.210724100250E-24
 v2sigmaab2    = -0.842896401002E-24
 v2sigmaabbb   = -0.421448200501E-24
 v2sigmabb2    = -0.210724100250E-24

 rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
 zk            = -0.103495051751E-01
 vrhoa         = -0.808320886304E-01
 vrhob         = -0.808320886304E-01
 vsigmaaa      =  0.782485884453E-02
 vsigmaab      =  0.156497176891E-01
 vsigmabb      =  0.782485884453E-02

 v2rhoa2       =  0.974468235880E-02
 v2rhoab       = -0.195273591427E+00
 v2rhob2       =  0.974468235880E-02

 v2rhoasigmaaa =  0.120240976646E-01
 v2rhoasigmaab =  0.240481953292E-01
 v2rhoasigmabb =  0.120240976646E-01
 v2rhobsigmaaa =  0.120240976646E-01
 v2rhobsigmaab =  0.240481953292E-01
 v2rhobsigmabb =  0.120240976646E-01

 v2sigmaaa2    = -0.914758040599E-02
 v2sigmaaaab   = -0.182951608120E-01
 v2sigmaaabb   = -0.914758040599E-02
 v2sigmaab2    = -0.365903216240E-01
 v2sigmaabbb   = -0.182951608120E-01
 v2sigmabb2    = -0.914758040599E-02

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.149930438056E+00
 vrhoa         = -0.472559970726E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.632520331340E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.130590959396E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    = -0.735327076770E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.147821297915E+00
 vrhoa         = -0.480325520180E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.608456906274E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.815520844173E-03
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    = -0.681211100497E-05
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.596933203014E-01
 vrhoa         = -0.494521116430E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.163547485507E-03
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.918226192391E-02
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    = -0.752364152605E-06
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
 zk            = -0.134201626930E-03
 vrhoa         = -0.113336793038E-01
 vrhob         =  0.000000000000E+00
 vsigmaaa      =  0.418957314771E-02
 vsigmaab      =  0.000000000000E+00
 vsigmabb      =  0.000000000000E+00

 v2rhoa2       = -0.639691764043E+00
 v2rhoab       =  0.000000000000E+00
 v2rhob2       =  0.000000000000E+00

 v2rhoasigmaaa =  0.000000000000E+00
 v2rhoasigmaab =  0.000000000000E+00
 v2rhoasigmabb =  0.000000000000E+00
 v2rhobsigmaaa =  0.000000000000E+00
 v2rhobsigmaab =  0.000000000000E+00
 v2rhobsigmabb =  0.000000000000E+00

 v2sigmaaa2    = -0.192576380057E+00
 v2sigmaaaab   =  0.000000000000E+00
 v2sigmaaabb   =  0.000000000000E+00
 v2sigmaab2    =  0.000000000000E+00
 v2sigmaabbb   =  0.000000000000E+00
 v2sigmabb2    =  0.000000000000E+00

 rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
 zk            = -0.466827360026E+01
 vrhoa         = -0.911386831086E-01
 vrhob         = -0.914525751657E-01
 vsigmaaa      =  0.195526696493E-06
 vsigmaab      =  0.391053392985E-06
 vsigmabb      =  0.195526696493E-06

 v2rhoa2       = -0.535201657866E-03
 v2rhoab       = -0.850604423798E-03
 v2rhob2       = -0.538228029079E-03

 v2rhoasigmaaa =  0.214289755089E-08
 v2rhoasigmaab =  0.428579510177E-08
 v2rhoasigmabb =  0.214289755089E-08
 v2rhobsigmaaa =  0.215378864780E-08
 v2rhobsigmaab =  0.430757729559E-08
 v2rhobsigmabb =  0.215378864780E-08

 v2sigmaaa2    = -0.128838364018E-13
 v2sigmaaaab   = -0.257676728036E-13
 v2sigmaaabb   = -0.128838364018E-13
 v2sigmaab2    = -0.515353456071E-13
 v2sigmaabbb   = -0.257676728036E-13
 v2sigmabb2    = -0.128838364018E-13

 rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
 zk            = -0.184442072405E+01
 vrhoa         = -0.814334534280E-01
 vrhob         = -0.820182123795E-01
 vsigmaaa      =  0.510839298939E-06
 vsigmaab      =  0.102167859788E-05
 vsigmabb      =  0.510839298939E-06

 v2rhoa2       = -0.124297349784E-02
 v2rhoab       = -0.183505806584E-02
 v2rhob2       = -0.125767116982E-02

 v2rhoasigmaaa =  0.134850158624E-07
 v2rhoasigmaab =  0.269700317248E-07
 v2rhoasigmabb =  0.134850158624E-07
 v2rhobsigmaaa =  0.136189478240E-07
 v2rhobsigmaab =  0.272378956480E-07
 v2rhobsigmabb =  0.136189478240E-07

 v2sigmaaa2    = -0.216571369852E-12
 v2sigmaaaab   = -0.433142739704E-12
 v2sigmaaabb   = -0.216571369852E-12
 v2sigmaab2    = -0.866285479407E-12
 v2sigmaabbb   = -0.433142739704E-12
 v2sigmabb2    = -0.216571369852E-12

 rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
 zk            = -0.416748767197E-02
 vrhoa         = -0.595487649525E-01
 vrhob         = -0.658275108773E-01
 vsigmaaa      =  0.500446791291E-02
 vsigmaab      =  0.100089358258E-01
 vsigmabb      =  0.500446791291E-02

 v2rhoa2       = -0.254320619062E+00
 v2rhoab       = -0.454911463637E+00
 v2rhob2       = -0.296322781759E+00

 v2rhoasigmaaa =  0.376458611019E-01
 v2rhoasigmaab =  0.752917222037E-01
 v2rhoasigmabb =  0.376458611019E-01
 v2rhobsigmaaa =  0.432477759464E-01
 v2rhobsigmaab =  0.864955518927E-01
 v2rhobsigmabb =  0.432477759464E-01

 v2sigmaaa2    = -0.873874773980E-02
 v2sigmaaaab   = -0.174774954796E-01
 v2sigmaaabb   = -0.873874773980E-02
 v2sigmaab2    = -0.349549909592E-01
 v2sigmaabbb   = -0.174774954796E-01
 v2sigmabb2    = -0.873874773980E-02

 rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
 zk            = -0.446674124729E-02
 vrhoa         = -0.576230218023E-01
 vrhob         = -0.884299955306E-01
 vsigmaaa      =  0.565507026894E-01
 vsigmaab      =  0.113101405379E+00
 vsigmabb      =  0.565507026894E-01

 v2rhoa2       =  0.211475716365E+00
 v2rhoab       = -0.313442480821E+00
 v2rhob2       =  0.661715882717E+00

 v2rhoasigmaaa = -0.366634765948E+00
 v2rhoasigmaab = -0.733269531896E+00
 v2rhoasigmabb = -0.366634765948E+00
 v2rhobsigmaaa = -0.183738128575E+00
 v2rhobsigmaab = -0.367476257149E+00
 v2rhobsigmabb = -0.183738128575E+00

 v2sigmaaa2    = -0.113703722438E+01
 v2sigmaaaab   = -0.227407444876E+01
 v2sigmaaabb   = -0.113703722438E+01
 v2sigmaab2    = -0.454814889752E+01
 v2sigmaabbb   = -0.227407444876E+01
 v2sigmabb2    = -0.113703722438E+01

 rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
 zk            = -0.853822762402E+01
 vrhoa         = -0.136330468738E+00
 vrhob         = -0.137562342876E+00
 vsigmaaa      =  0.397018264645E-05
 vsigmaab      =  0.794036529289E-05
 vsigmabb      =  0.397018264645E-05

 v2rhoa2       =  0.505383535819E-03
 v2rhoab       = -0.719176570320E-03
 v2rhob2       =  0.520080723382E-03

 v2rhoasigmaaa =  0.220509385365E-09
 v2rhoasigmaab =  0.441018770730E-09
 v2rhoasigmabb =  0.220509385365E-09
 v2rhobsigmaaa =  0.606911398719E-09
 v2rhobsigmaab =  0.121382279744E-08
 v2rhobsigmabb =  0.606911398719E-09

 v2sigmaaa2    = -0.518845125296E-11
 v2sigmaaaab   = -0.103769025059E-10
 v2sigmaaabb   = -0.518845125296E-11
 v2sigmaab2    = -0.207538050118E-10
 v2sigmaabbb   = -0.103769025059E-10
 v2sigmabb2    = -0.518845125296E-11

 rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
 zk            = -0.615277551535E+01
 vrhoa         = -0.133907518541E+00
 vrhob         = -0.133907518541E+00
 vsigmaaa      =  0.414743637478E-05
 vsigmaab      =  0.829487274956E-05
 vsigmabb      =  0.414743637478E-05

 v2rhoa2       =  0.435534219905E-03
 v2rhoab       = -0.969294804880E-03
 v2rhob2       =  0.435534219905E-03

 v2rhoasigmaaa =  0.125079346911E-07
 v2rhoasigmaab =  0.250158693823E-07
 v2rhoasigmabb =  0.125079346911E-07
 v2rhobsigmaaa =  0.125079346911E-07
 v2rhobsigmaab =  0.250158693823E-07
 v2rhobsigmabb =  0.125079346911E-07

 v2sigmaaa2    = -0.710218297869E-11
 v2sigmaaaab   = -0.142043659574E-10
 v2sigmaaabb   = -0.710218297869E-11
 v2sigmaab2    = -0.284087319147E-10
 v2sigmaabbb   = -0.142043659574E-10
 v2sigmabb2    = -0.710218297869E-11

 rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
 zk            = -0.513978017858E-02
 vrhoa         = -0.677801008926E-01
 vrhob         = -0.720903652109E-01
 vsigmaaa      =  0.723862883206E-02
 vsigmaab      =  0.144772576641E-01
 vsigmabb      =  0.723862883206E-02

 v2rhoa2       = -0.172744794034E+00
 v2rhoab       = -0.396171978594E+00
 v2rhob2       = -0.188314565988E+00

 v2rhoasigmaaa =  0.403167927387E-01
 v2rhoasigmaab =  0.806335854774E-01
 v2rhoasigmabb =  0.403167927387E-01
 v2rhobsigmaaa =  0.445643432070E-01
 v2rhobsigmaab =  0.891286864140E-01
 v2rhobsigmabb =  0.445643432070E-01

 v2sigmaaa2    = -0.149411369506E-01
 v2sigmaaaab   = -0.298822739012E-01
 v2sigmaaabb   = -0.149411369506E-01
 v2sigmaab2    = -0.597645478025E-01
 v2sigmaabbb   = -0.298822739012E-01
 v2sigmabb2    = -0.149411369506E-01

 rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
 zk            = -0.203181802597E-02
 vrhoa         = -0.608626957301E-01
 vrhob         = -0.779575998152E-01
 vsigmaaa      =  0.578479124549E-01
 vsigmaab      =  0.115695824910E+00
 vsigmabb      =  0.578479124549E-01

 v2rhoa2       =  0.152335258351E-01
 v2rhoab       = -0.720403954936E+00
 v2rhob2       =  0.952448412124E-01

 v2rhoasigmaaa =  0.224592118966E+00
 v2rhoasigmaab =  0.449184237932E+00
 v2rhoasigmabb =  0.224592118966E+00
 v2rhobsigmaaa =  0.517340140054E+00
 v2rhobsigmaab =  0.103468028011E+01
 v2rhobsigmabb =  0.517340140054E+00

 v2sigmaaa2    = -0.238365551254E+01
 v2sigmaaaab   = -0.476731102509E+01
 v2sigmaaabb   = -0.238365551254E+01
 v2sigmaab2    = -0.953462205018E+01
 v2sigmaabbb   = -0.476731102509E+01
 v2sigmabb2    = -0.238365551254E+01
"""

def parsedata(dataname):
    import re
    from pyquante2.utils import parseline
    fields = []
    f = iter(dataname.splitlines())
    while True:
        try:
            line = next(f)
        except StopIteration:
            break
        if line.lstrip().startswith('rhoa'):
            rhoa,rhob,sigaa,sigab,sigbb = parseline(line,'xfxfxfxfxf')
            field = [rhoa,rhob,sigaa,sigab,sigbb]
            fields.append(field)
            for i in range(6):
                line = next(f)
                term = parseline(line,'xxf')
                field.append(term)
        else:
            continue
    return np.array(fields)

data = {
    'xlda' : parsedata(x_lda_data),
    'xb88' : parsedata(x_b88_data),
    'xpbe' : parsedata(x_pbe_data),
    'cvwn' : parsedata(c_vwn5_data),
    'clyp' : parsedata(c_lyp_data),
    'cpbe' : parsedata(c_pbe_data),
    }
# Redundant names
data['xs'] = data['xlda']
data['cvwn5'] = data['cvwn']

if __name__ == '__main__':
    xlda = parsedata(x_lda_data)
    cvwn = parsedata(c_vwn5_data)
    print((xlda[:,:5] == cvwn[:,:5]).all())
