import os
import sys

import numpy as np
import pandas as pd

from itertools import permutations
from math import factorial, log
from time import time


from xmlparse import loadRef, loadGeometryBases, getXmlScore, minXml, loadScores
from score_strokes import alignStrokes, greedyAlign2
from exhaustive import computeExhaustive, exhaustScore, exhaustScoreAlignment


def matchError(heuristic_scores, standard_scores, timed=False):
    return log(standard_scores/heuristic_scores, 2)

def getScores(algorithm, ref_char, data_dir):
    """
    Calculates scores for a set of gene characters against a given reference with a heuristic algorithm
    Files from the directory specified with gene data are read in a deterministic order (same files = same order)
    Input:
    algorithm: A function with the below signature. Takes stroke geometry and archetype geometry and scores the stroke against the archetype
        Input:
        stroke_geometry: List of the strokes from the gene instance
        reference_geometry: List of the strokes from the archetype
        stroke_fractional_distance: List of fractional distances (a.k.a. progress percentages) for each stroke in the gene instance
        reference_fractional_distance: List of fractional distances (a.k.a. progress percentages) for each stroke in the reference
    ref_char: UTF-8 name of the archetype in question
    data_dir: Directory with the gene files for testing
    Output:
    heuristic_scores: The scores returned from the given algorithm for each of the genes in the directory
    heuristic_alignments: Alignments which the algorithm returned
    marks: Whether or not each respective score has a mark, which is an additional stroke which has no counterpart in the archetype
    """
    heuristic_alignments = []
    heuristic_scores = []
    marks = []
    ref_geometry, ref_progress_percentage, output_size = loadRef(ref_char, "Reference")
    g_data, han_chars, base_data, stroke_sets, stroke_orders, f_names = loadGeometryBases(data_dir, output_size)
    for (geometry_length, han_char, bases, stroke_set, stroke_order, f_name) in zip(g_data, han_chars, base_data, stroke_sets, stroke_orders, f_names):
        geometry, progress_percentage = geometry_length
        heuristic_alignment = np.array(algorithm(geometry, ref_geometry, progress_percentage, ref_progress_percentage))+1
        heuristic_alignments.append(heuristic_alignment)
        heuristic_xml = minXml(ref_char, bases, stroke_set, heuristic_alignment)
        heuristic_score = getXmlScore(heuristic_xml)
        heuristic_scores.append(heuristic_score)
        marks.append(len(geometry)!=len(ref_geometry))
    return heuristic_scores, heuristic_alignments, marks

def exhaustiveScores(ref_char, data_dir, timed = False, save = False, file_prefix=""):
    #ref_g, ref_l, output_size = loadRef(ref_char, "Reference")
    #g_data, han_chars, base_data, stroke_sets, stroke_orders, f_names = loadGeometryBases(data_dir, output_size)
    exhaustive_scores = []
    times_exhaustive = []
    file_names = []
    for _, _, files in os.walk(data_dir):
        for file in sorted(files):
            if file.endswith("gene"):
                file_names.append(file)
    for file in file_names:
        if timed:
            time_store = time()
            original_score = exhaustScore(ref_char, file, data_dir, force_refresh=True, save=save, file_pretfix=file_prefix)
            times_exhaustive.append(time()-time_store)
        else:
            original_score = exhaustScore(ref_char, file, data_dir, save=save, file_prefix=file_prefix)
        exhaustive_scores.append(original_score)
    if timed:
        print(times_exhaustive)
    return exhaustive_scores

