import numpy as np

channel_numbers_to_names = {
    0: "ACT0L",
    1: "ACT0R",
    2: "ACT1L",
    3: "ACT1R",
    4: "ACT2L",
    5: "ACT2R",
    6: "ACT3L",
    7: "ACT3R",
    8: "TOF00",
    9: "TOF01",
    10: "TOF02",
    11: "TOF03",
    12: "TOF10",
    13: "TOF11",
    14: "TOF12",
    15: "TOF13",
    16: "Hole0",
    17: "Hole1",
    18: "PbGlass",
}

channel_names_to_numbers = {v: k for k, v in channel_numbers_to_names.items()}
