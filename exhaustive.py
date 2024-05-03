"""
exhaustive.py
functions for dealing with exhaustive comparisons between Stylus genomes
"""

from itertools import permutations
from math import factorial

import numpy as np

import os

from xmlparse import loadRef, loadGeometryBases, getXmlScore, minXml

import stylusengine

stylusengine.setLogFile(b'errors.log')

stylusengine.setScope(
    b'file:///home/tulip/Documents/College/Stewart/stylusapp/hans',
    b'file:///home/tulip/Documents/College/Stewart/stylus/schemas'
)

def computeExhaustive(ref_char, f_read, data_dir, exhaust_dir = "Exhaustive", prog_interval = 100, save = True, xml_dir = "GenXml/Exhaustive", save_file = ""):
    ref_g, ref_l, output_size = loadRef(ref_char, "Reference")
    g_data, _, base_data, stroke_sets, _, f_names = loadGeometryBases(data_dir, output_size, f_read = f_read)
    n_strokes = len(ref_g)
    for i in range(len(g_data)):
        #print(f"Generating exhaustive scores for sample {f_read[i]}")
        bases = base_data[i]
        stroke_set = stroke_sets[i]
        exhaustive_alignments = permutations(range(1, n_strokes+1))
        exhaustive_scores = np.zeros(factorial(n_strokes))
        for j, p in enumerate(exhaustive_alignments):
            p_xml = minXml(ref_char, bases, stroke_set, p)
            exhaustive_scores[j] = getXmlScore(p_xml, f"{xml_dir}/{i}_{j}_{f_read[i]}", f"{xml_dir}/{i}_{j}_min_{f_read[i]}")
            #if j%prog_interval == 0:
            #    print(f"Scoring permutation {j} of {len(exhaustive_scores)}")
        if save:
            if save_file == "":
                f_name_cleaned = f_read[i].replace("/", "_")
                f"{exhaust_dir}/exhaust_{ref_char}_{f_name_cleaned}.npy"
            print(f"Wrote exhaustive scores to {save_file}")
            np.save(save_file, exhaustive_scores)
        yield exhaustive_scores

def computeExhaustiveAlign(ref_char, f_read, data_dir, exhaust_dir = "Exhaustive", prog_interval = 100, save = True, xml_dir = "GenXml/Exhaustive", save_file = ""):
    ref_g, ref_l, output_size = loadRef(ref_char, "Reference")
    g_data, _, base_data, stroke_sets, _, f_names = loadGeometryBases(data_dir, output_size, f_read = f_read)
    n_strokes = len(ref_g)
    for i in range(len(g_data)):
        #print(f"Generating exhaustive scores for sample {f_read[i]}")
        bases = base_data[i]
        stroke_set = stroke_sets[i]
        exhaustive_alignments = list(permutations(range(1, n_strokes+1)))
        exhaustive_scores = np.zeros(factorial(n_strokes))
        for j, p in enumerate(exhaustive_alignments):
            p_xml = minXml(ref_char, bases, stroke_set, p)
            exhaustive_scores[j] = getXmlScore(p_xml, f"{xml_dir}/{i}_{j}_{f_read[i]}", f"{xml_dir}/{i}_{j}_min_{f_read[i]}")
            #if j%prog_interval == 0:
            #    print(f"Scoring permutation {j} of {len(exhaustive_scores)}")
        if save:
            if save_file == "":
                f_name_cleaned = f_read[i].replace("/", "_")
                f"{exhaust_dir}/exhaust_{ref_char}_{f_name_cleaned}.npy"
            print(f"Wrote exhaustive scores to {save_file}")
            np.save(save_file, exhaustive_scores)
        yield exhaustive_scores, exhaustive_alignments

def readExhaustive(ref_char, f_name, exhaust_dir = "Exhaustive", save_file = ""):
    if save_file == "":
        f_name_cleaned = f_name.replace("/", "_")
        save_file = f"{exhaust_dir}/exhaust_{ref_char}_{f_name_cleaned}.npy"
    print(f"Read from {save_file}")
    data_read = np.load(save_file)
    return data_read

def exhaustScore(ref_char, f_name, data_dir, exhaust_dir = "Exhaustive", force_refresh = False, save = True, file_prefix = ""):
    f_name_cleaned = f_name.replace("/", "_")
    exhaust_name = f"{exhaust_dir}/exhaust_{file_prefix}{ref_char}_{f_name_cleaned}.npy"
    exhaust_maxes = []
    if not os.path.isfile(exhaust_name) or force_refresh:
        for e in computeExhaustive(ref_char, [f_name], data_dir, save = save, save_file = exhaust_name):
            exhaust_maxes.append(e.max())
    else:
        exhaust_maxes = readExhaustive(ref_char, f_name, exhaust_dir, exhaust_name)
    return np.max(exhaust_maxes)

def exhaustScoreAlignment(ref_char, f_name, data_dir, exhaust_dir = "Exhaustive", save = True):
    exhaust_maxes, exhaust_alignments = [], []
    for e in computeExhaustiveAlign(ref_char, [f_name], data_dir, prog_interval = prog_interval, save = save):
        exhaust_scores, alignments = list(e)
        exhaust_maxes = exhaust_scores.max()
        exhaust_alignments.append(alignments[np.argmax(exhaust_scores)])
    return exhaust_alignments[np.argmax(exhaust_maxes)]

