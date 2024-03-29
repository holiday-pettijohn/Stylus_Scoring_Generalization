#!/bin/python3
import numpy as np

"""
# old function with complex conflict resolution, archived for reference
def alignStrokes(strokes, ref, p_strokes, p_ref):
    stroke_map = np.full(len(strokes), -1)
    error_map = np.zeros((len(ref), len(strokes)), dtype=float)
    matches_tried = np.zeros(len(strokes), dtype=int)
    for i, ref_stroke, r_progresses in zip(range(len(ref)), ref, p_ref):
        for j, candidate_stroke, c_progresses in zip(range(len(strokes)), strokes, p_strokes):
            error_map[i, j] = strokeError(ref_stroke, candidate_stroke, r_progresses, c_progresses)
    for i, err in enumerate(error_map):
        candidate = np.argmin(err)
        stroke_map[candidate] = i
    # resolve conflicts until only one stroke is mapped to each reference stroke
    # the condition is based on the 'bad' candidates being set to -1, indicating no match
    while (np.unique(stroke_map).shape[0] != len(strokes) - (len(strokes)-len(ref)-1 if len(strokes)-len(ref) > 1 else 0)
           or (-1 in stroke_map) if len(strokes) == len(ref) else False):
        # conflict resolution is still rudimentery
        for i in range(len(stroke_map)):
            if stroke_map[i] == -1:
                prios = np.argsort(error_map[:, i])
                stroke_map[i] = prios[matches_tried[i]]
        for i, conflicted in enumerate(stroke_map):
            if conflicted != -1:
                candidates = np.argwhere(stroke_map == i).flatten()
                if candidates.shape[0] > 1:
                    best_candidate = error_map[conflicted, candidates].argmin()
                    stroke_map[candidates] = -1
                    matches_tried[candidates] += 1
                    stroke_map[candidates[best_candidate]] = conflicted
                    matches_tried[best_candidate] -= 1
        matches_tried = matches_tried%len(ref)
    return stroke_map
"""

def alignStrokes(strokes, ref, p_strokes, p_ref):

    error_maps = strokeErrorMatrix(strokes, ref, p_strokes, p_ref)

    # function to get stroke length given a stroke value (in this case, a stroke value is a 2d list
    # that contains an x coord in stroke[][0] and a y coord in stroke[][1])
    def getStrokeLen(stroke):
        length = 0 # adding all the lengths we get between two points to this variable
        # while it looks a little complicated, this is just the pythagorean theorem applied between two coordinates, just 
        # how one would calculate it on a graph: sqrt(a^2 + b^2)
        for i in range(len(stroke)-1):
            length += ((stroke[i][0] - stroke[i+1][0])**2 + (stroke[i][1] - stroke[i+1][1])**2)**0.5
        return length

    #get the lengths of each stroke for the order in the greedy algorithm
    ref_lengths = []
    for i in range(len(strokes)):
        ref_lengths.append(getStrokeLen(ref[i]))

    # -1 just means unmatched here since 0 (the other 'default' filler) is a meaningful number in this context
    stroke_map = np.full(len(strokes), -1)

    for each in ref_lengths:
        largestref = np.argmax(ref_lengths) # this is the index for the reference stroke that is largest
        smallerror = np.argmin(error_maps[largestref]) # access the error map from the largest stroke's index and see which error is smallest

        while(stroke_map[smallerror]!=-1):
            # change small error so that we do not repeat over indexes that are already taken
            # just keeps repeating until we land on an index that doesn't already have a value in its place
            error_maps[largestref][smallerror] = 10000
            smallerror = np.argmin(error_maps[largestref])

        stroke_map[smallerror] = largestref # set the index in the stroke_map to the reference stroke we designated

        ref_lengths[largestref] = 0 # set the length of the reference stroke to 0 so we never use it again

    return stroke_map

def alignStrokesResolve(strokes, ref, p_strokes, p_ref):
    
    # hack to ignore the last stroke due to dropout bug
    strokes, p_strokes = strokes[:len(ref)], p_strokes[:len(ref)]
    
    error_maps = strokeErrorMatrix(strokes, ref, p_strokes, p_ref)

    # function to get stroke length given a stroke value (in this case, a stroke value is a 2d list
    # that contains an x coord in stroke[][0] and a y coord in stroke[][1])
    def getStrokeLen(stroke):
        length = 0 # adding all the lengths we get between two points to this variable
        # while it looks a little complicated, this is just the pythagorean theorem applied between two coordinates, just 
        # how one would calculate it on a graph: sqrt(a^2 + b^2)
        for i in range(len(stroke)-1):
            length += ((stroke[i][0] - stroke[i+1][0])**2 + (stroke[i][1] - stroke[i+1][1])**2)**0.5
        return length

    #get the lengths of each stroke for the order in the greedy algorithm
    ref_lengths = []
    for i in range(len(ref)):
        ref_lengths.append(getStrokeLen(ref[i]))

    # -1 just means unmatched here since 0 (the other 'default' filler) is a meaningful number in this context
    stroke_map = np.full(len(ref), -1)

    for each in ref_lengths:
        largestref = np.argmax(ref_lengths) # this is the index for the reference stroke that is largest
        smallerror = np.argmin(error_maps[largestref]) # access the error map from the largest stroke's index and see which error is smallest

        while(stroke_map[smallerror]!=-1):
            # change small error so that we do not repeat over indexes that are already taken
            # just keeps repeating until we land on an index that doesn't already have a value in its place
            error_maps[largestref][smallerror] = 10000
            smallerror = np.argmin(error_maps[largestref])

        stroke_map[smallerror] = largestref # set the index in the stroke_map to the reference stroke we designated

        ref_lengths[largestref] = 0 # set the length of the reference stroke to 0 so we never use it again

    return stroke_map

def strokeErrorMatrix(strokes, ref, p_strokes, p_ref):
    error_map = np.zeros((len(ref), len(strokes)), dtype=float)
    matches_tried = np.zeros(len(strokes), dtype=int)
    for i, ref_stroke, r_progresses in zip(range(len(ref)), ref, p_ref):
        for j, candidate_stroke, c_progresses in zip(range(len(strokes)), strokes, p_strokes):
            error_map[i, j] = strokeError(ref_stroke, candidate_stroke, r_progresses, c_progresses)
    return error_map

def strokeError(stroke, ref_stroke, p_stroke, p_ref, mode="max"):
    forward_stroke_error, back_stroke_error = np.zeros(len(ref_stroke)), np.zeros(len(ref_stroke))
    forward_ref_error, back_ref_error = np.zeros(len(stroke)), np.zeros(len(stroke))
    for i, rpoint, rprogress in zip(range(len(ref_stroke)), ref_stroke, p_ref):
        forward_stroke_error[i] = np.linalg.norm((rpoint-strokeTrace(stroke, p_stroke, rprogress)))
    for i, rpoint, rprogress in zip(range(len(ref_stroke)), ref_stroke[::-1], p_ref[::-1]):
        back_stroke_error[i] = np.linalg.norm((rpoint-strokeTrace(stroke, p_stroke, 1-rprogress)))
    for i, point, progress in zip(range(len(stroke)), stroke, p_stroke):
        forward_ref_error[i] = np.linalg.norm((point-strokeTrace(ref_stroke, p_ref, progress)))
    for i, point, progress in zip(range(len(stroke)), stroke[::-1], p_stroke[::-1]):
        back_ref_error[i] = np.linalg.norm((point-strokeTrace(ref_stroke, p_ref, 1-progress)))
    final_error = min(max(forward_stroke_error.sum(), forward_ref_error.sum()), max(back_stroke_error.sum(), back_ref_error.sum()))
    return final_error

def strokeTrace(stroke, stroke_progresses, progress):
    if progress == 1:
        return stroke[-1]
    progress_line = len(stroke_progresses)-1
    for i in range(1, len(stroke_progresses)):
        if stroke_progresses[i] > progress:
            progress_line = i-1
            break
    startp, endp = stroke_progresses[progress_line], stroke_progresses[progress_line+1]
    norm_progress = (progress-startp)/endp
    if stroke[progress_line+1][0] == stroke[progress_line][0]:
        x = stroke[progress_line][0]
        y = norm_progress*(stroke[progress_line+1][1]-stroke[progress_line][1])+stroke[progress_line][1]
    else:
        slope = (stroke[progress_line+1][1]-stroke[progress_line][1])/(stroke[progress_line+1][0]-stroke[progress_line][0])
        intercept = stroke[progress_line][1]-slope*stroke[progress_line][0]
        x = norm_progress*(stroke[progress_line+1][0]-stroke[progress_line][0])+stroke[progress_line][0]
        y = slope*x + intercept
    return np.array((x, y))
