"""
    This is an example code to save a python data array to a binary file
    and c++ program can read from it. The data is saved with single float
    precision (float32).
"""

import array

outputFile = open("O16_PGCM.bin.in", "wb")
for config_i in data:
    float_array = array.array('f', config_i)
    float_array.tofile(outputFile)
outputFile.close()
