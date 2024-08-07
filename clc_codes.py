clc_codes = {
    # 1: Artificial surfaces
    # 1.1: Urban fabric
    111: "Continuous urban fabric",
    112: "Discontinuous urban fabric",
    # 1.2: Industrial, commercial and transport units
    121: "Industrial or commercial units",
    122: "Road and rail networks and associated land",
    123: "Port areas",
    124: "Airports",
    # 1.3: Mine, dump and construction sites
    131: "Mineral extraction sites",
    132: "Dump sites",
    133: "Construction sites",
    # 1.4: Artificial, non-agricultural vegetated areas
    141: "Green urban areas",
    142: "Sport and leisure facilities",
    # 2: Agricultural areas
    # 2.1: Arable land
    211: "Non-irrigated arable land",
    212: "Permanently irrigated land",
    213: "Rice fields",
    # 2.2: Permanent crops
    221: "Vineyards",
    222: "Fruit trees and berry plantations",
    223: "Olive groves",
    # 2.3: Pastures
    231: "Pastures",
    # 2.4: Heterogeneous agricultural areas
    241: "Annual crops associated with permanent crops",
    242: "Complex cultivation patterns",
    243: "Land principally occupied by agriculture, with significant areas of natural vegetation",
    244: "Agro-forestry areas",
    # 3: Forest and seminatural areas
    # 3.1: Forest
    311: "Broad-leaved forest",
    312: "Coniferous forest",
    313: "Mixed forest",
    # 3.2: Shrub and/or herbaceous vegetation associations
    321: "Natural grassland",
    322: "Moors and heathland",
    323: "Sclerophyllous vegetation",
    324: "Transitional woodland/shrub",
    # 3.3: Open spaces with little or no vegetation
    331: "Beaches, dunes, sands",
    332: "Bare rock",
    333: "Sparsely vegetated areas",
    334: "Burnt areas",
    335: "Glaciers and perpetual snow",
    # 4: Wetlands
    # 4.1: Inland wetlands
    411: "Inland marshes",
    412: "Peatbogs",
    # 4.2: Coastal wetlands
    421: "Salt marshes",
    422: "Salines",
    423: "Intertidal flats",
    # 5: Water bodies
    # 5.1: Inland waters
    511: "Water courses",
    512: "Water bodies",
    # 5.2: Marine waters
    521: "Coastal lagoons",
    522: "Estuaries",
    523: "Sea and ocean",
}

clc_reverse = {
    'Continuous urban fabric': 111,
    'Discontinuous urban fabric': 112,
    'Industrial or commercial units': 121,
    'Road and rail networks and associated land': 122,
    'Port areas': 123,
    'Airports': 124,
    'Mineral extraction sites': 131,
    'Dump sites': 132,
    'Construction sites': 133,
    'Green urban areas': 141,
    'Sport and leisure facilities': 142,
    'Non-irrigated arable land': 211,
    'Permanently irrigated land': 212,
    'Rice fields': 213,
    'Vineyards': 221,
    'Fruit trees and berry plantations': 222,
    'Olive groves': 223,
    'Pastures': 231,
    'Annual crops associated with permanent crops': 241,
    'Complex cultivation patterns': 242,
    'Land principally occupied by agriculture, with significant areas of natural vegetation': 243,
    'Agro-forestry areas': 244,
    'Broad-leaved forest': 311,
    'Coniferous forest': 312,
    'Mixed forest': 313,
    'Natural grassland': 321,
    'Moors and heathland': 322,
    'Sclerophyllous vegetation': 323,
    'Transitional woodland/shrub': 324,
    'Beaches, dunes, sands': 331,
    'Bare rock': 332,
    'Sparsely vegetated areas': 333,
    'Burnt areas': 334,
    'Glaciers and perpetual snow': 335,
    'Inland marshes': 411,
    'Peatbogs': 412,
    'Salt marshes': 421,
    'Salines': 422,
    'Intertidal flats': 423,
    'Water courses': 511,
    'Water bodies': 512,
    'Coastal lagoons': 521,
    'Estuaries': 522,
    'Sea and ocean': 523,
}
