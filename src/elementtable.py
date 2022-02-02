
element_table = {
    0: 0.00,  # Dummy
    1: 0.31,  # Hydrogen
    2: 0.28,  # Helium
    3: 1.28,  # Lithium
    4: 0.96,  # Beryllium
    5: 0.84,  # Boron
    6: 0.76,  # Carbon
    7: 0.71,  # Nitrogen
    8: 0.66,  # Oxygen
    9: 0.57,  # Fluorine
    10: 0.58,  # Neon
    11: 1.66,  # Sodium
    12: 1.41,  # Magnesium
    13: 1.21,  # Aluminium
    14: 1.11,  # Silicon
    15: 1.07,  # Phosphorus
    16: 1.05,  # Sulfur
    17: 1.02,  # Chlorine
    18: 1.06,  # Argon
    19: 2.03,  # Potassium
    20: 1.76,  # Calcium
    21: 1.70,  # Scandium
    22: 1.60,  # Titanium
    23: 1.53,  # Vanadium
    24: 1.39,  # Chromium
    25: 1.39,  # Manganese
    26: 1.32,  # Iron
    27: 1.26,  # Cobalt
    28: 1.24,  # Nickel
    29: 1.32,  # Copper
    30: 1.22,  # Zinc
    31: 1.22,  # Gallium
    32: 1.20,  # Germanium
    33: 1.19,  # Arsenic
    34: 1.20,  # Selenium
    35: 1.20,  # Bromine
    36: 1.16,  # Krypton
    37: 2.20,  # Rubidium
    38: 1.95,  # Strontium
    39: 1.90,  # Yttrium
    40: 1.75,  # Zirconium
    41: 1.64,  # Niobium
    42: 1.54,  # Molybdenum
    43: 1.47,  # Technetium
    44: 1.46,  # Ruthenium
    45: 1.42,  # Rhodium
    46: 1.39,  # Palladium
    47: 1.45,  # Silver
    48: 1.44,  # Cadmium
    49: 1.42,  # Indium
    50: 1.39,  # Tin
    51: 1.39,  # Antimony
    52: 1.38,  # Tellurium
    53: 1.39,  # Iodine
    54: 1.40,  # Xenon
    55: 2.44,  # Caesium
    56: 2.15,  # Barium
    57: 2.07,  # Lanthanum
    58: 2.04,  # Cerium
    59: 2.03,  # Praseodymium
    60: 2.01,  # Neodymium
    61: 1.99,  # Promethium
    62: 1.98,  # Samarium
    63: 1.98,  # Europium
    64: 1.96,  # Gadolinium
    65: 1.94,  # Terbium
    66: 1.92,  # Dysprosium
    67: 1.92,  # Holmium
    68: 1.89,  # Erbium
    69: 1.90,  # Thulium
    70: 1.87,  # Ytterbium
    71: 1.87,  # Lutetium
    72: 1.75,  # Hafnium
    73: 1.70,  # Tantalum
    74: 1.62,  # Tungsten
    75: 1.51,  # Rhenium
    76: 1.44,  # Osmium
    77: 1.41,  # Iridium
    78: 1.36,  # Platinum
    79: 1.36,  # Gold
    80: 1.32,  # Mercury
    81: 1.45,  # Thallium
    82: 1.46,  # Lead
    83: 1.48,  # Bismuth
    84: 1.40,  # Polonium
    85: 1.50,  # Astatine
    86: 1.50,  # Radon
    87: 2.60,  # Francium
    88: 2.21,  # Radium
    89: 2.15,  # Actinium
    90: 2.06,  # Thorium
    91: 2.00,  # Protactinium
    92: 1.96,  # Uranium
    93: 1.90,  # Neptunium
    94: 1.87,  # Plutonium
    95: 1.80,  # Americium
    96: 1.69,  # Curium
    97: 1.60,  # Berkelium
    98: 1.60,  # Californium
    99: 1.60,  # Einsteinium
    100: 1.60, # Fermium
    101: 1.60, # Mendelevium
    102: 1.60, # Nobelium
    103: 1.60, # Lawrencium
    104: 1.60, # Rutherfordium
    105: 1.60, # Dubnium
    106: 1.60, # Seaborgium
    107: 1.60, # Bohrium
    108: 1.60, # Hassium
    109: 1.60, # Meitnerium
    110: 1.60, # Darmstadtium
    111: 1.60, # Roentgenium
    112: 1.60, # Copernicium
    113: 1.60, # Nihonium
    114: 1.60, # Flerovium
    115: 1.60, # Moscovium
    116: 1.60, # Livermorium
    117: 1.60, # Tennessine
    118: 1.60, # Oganesson
}

'''
table generated from https://github.com/openbabel/openbabel/blob/master/src/elementtable.h using the code below to extract data.
'''
def readopenbabeltable():
    with open('text.txt') as f:
        lines = f.readlines()

    new_text = ""
    for line in lines:
        linev = line.split(",")
        id = linev[0][10:]
        covalent_radius = linev[3][1:]
        name = linev[-2][1:-1]

        new_line = f"    {id}: {covalent_radius},    #{name}\n"
        new_text += new_line
    print(new_text)
