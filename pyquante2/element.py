"""\
 element.py: Miscellaneous data about elements

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

symbol = [
    "X","H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl","Pb","Bi","Po","At","Rn"]


sym2no = {}
for i in range(len(symbol)):
    sym2no[symbol[i]] = i
    sym2no[symbol[i].lower()] = i
    
name = [
    "dummy",
    "hydrogen", "helium",
    "lithium","beryllium","boron","carbon","nitrogen",
    "oxygen","fluorine","neon","sodium","magnesium",
    "aluminum","silicon","phosphorus","sulfur","chlorine",
    "argon","potassium", "calcium", "scandium", "titanium",
    "vanadium", "chromium", "manganese", "iron",
    "cobalt", "nickel", "copper", "zinc",
    "gallium", "germanium", "arsenic", "selenium", "bromine",
    "krypton", "rubidium", "strontium", "yttrium", "zirconium",
    "niobium", "molybdenum", "technetium", "ruthenium","rhodium",
    "palladium", "silver", "cadmium","indium", "tin", "antimony",
    "tellerium", "iodine", "xenon","cesium", "barium",
    "lanthanum","cerium","praseodymium","neodymium","promethium",
    "samarium","europium","gadolinium","terbium","dysprosium",
    "holmium","erbium","thulium","ytterbium","lutetium",
    "halfnium","tantalum","tungsten","rhenium","osmium","iridium",
    "platinum","gold","mercury","thallium","lead","bismuth",
    "polonium","astatine","radon"]

name2no = {}
for i in range(len(name)):
    name2no[name[i]] = i
    name2no[name[i].upper()] = i

# British(?) spellings
name2no['PHOSPHOROUS'] = name2no['phosphorus']
name2no['PHOSPHOROUS'.lower()] = name2no['phosphorus']

mass = [
    0.00,
    1.0008, 4.0026,
    6.941,9.0122,
    10.811,12.011,14.007,15.999,18.998,20.179,
    22.990,24.305,
    26.982,28.086,30.974,32.066,35.453,39.948,
    39.098, 40.078,
    44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845,
    58.9332, 58.6934, 63.546,65.39,
    69.723, 72.61, 74.9216, 78.96, 79.904, 83.80,
    85.4678, 87.62,
    88.90686, 91.224, 92.90638, 95.94, 98, 101.07,
    102.90550, 106.42, 107.8682, 112.411,
    114.818, 118.710, 121.760, 127.60, 126.90447, 131.29,
    132.90545, 137.327, 138.9055, 140.11, 140.90765, 144.24,
    145.0, 150.36, 151.964,
    157.25, 158.92534, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967,
    178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59]


color = [(255, 20, 147),(250, 235, 215),(255, 192, 203),(178, 34, 34),
         (34, 139, 34),(0, 255, 0),(112, 128, 144),(0, 191, 255),
         (255, 0, 0),(218, 165, 32),(255, 105, 180),(0, 0, 255),
         (34, 139, 34),(190, 190, 190),(218, 165, 32),(255, 165, 0),
         (255, 255, 0),(0, 255, 0),(255, 192, 203),(255, 20, 147),
         (128, 128, 128),(190, 190, 190),(190, 190, 190),(190, 190, 190),
         (190, 190, 190),(190, 190, 190),(255, 165, 0),(165, 42, 42),
         (165, 42, 42),(165, 42, 42),(165, 42, 42),(165, 42, 42),
         (85, 107, 47),(253, 245, 230),(152, 251, 152),(165, 42, 42),
         (50, 205, 50),(165, 42, 42),(190, 190, 190),(190, 190, 190),
         (190, 190, 190),(190, 190, 190),(255, 127, 80),(190, 190, 190),
         (190, 190, 190),(190, 190, 190),(190, 190, 190),(190, 190, 190),
         (255, 140, 0),(190, 190, 190),(190, 190, 190),(190, 190, 190),
         (190, 190, 190),(160, 32, 240),(255, 105, 180),(165, 42, 42),
         (190, 190, 190),(190, 190, 190),
         None,None,None,None,None,None,None,None,None,None,None,None,None,
         (190, 190, 190),(190, 190, 190),(190, 190, 190),(64, 224, 208),
         (190, 190, 190),(190, 190, 190),(190, 190, 190),(190, 190, 190),
         (255, 215, 0),(190, 190, 190),(190, 190, 190),(190, 190, 190),
         (255, 181, 197)]

radius = [1.0000,1.2000,1.4000,1.8200,1.3725,0.7950,1.7000,
          1.5500,1.5200,1.4700,1.5400,2.2700,1.7300,1.7000,
          2.1000,1.8000,1.8000,1.7500,1.8800,2.7500,2.4500,
          1.3700,1.3700,1.3700,1.3700,1.3700,1.4560,0.8800, 
          0.6900,0.7200,0.7400,1.3700,1.9500,1.8500,1.9000, 
          1.8500,2.0200,1.5800,2.1510,1.8010,1.6020,1.4680, 
          1.5260,1.3600,1.3390,1.3450,1.3760,1.2700,1.4240, 
          1.6630,2.1000,2.0500,2.0600,1.9800,2.0000,1.8400, 
          2.2430,1.8770,None,None,None,None,None, 
          None,None,None,None,None,None,None, 
          None,2.1700,1.5800,1.4670,1.5340,1.3750,1.3530, 
          1.3570,1.7500,1.6600,1.5500,1.9600,2.0200,2.1500]

