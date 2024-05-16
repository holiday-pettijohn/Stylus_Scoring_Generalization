#!/bin/python3
"""
Contains current implementations of heuristic matching algorithms
of the form which may be tested against exhaustive results
given the existing framework
"""

import numpy as np

def alignStrokes(strokes, ref, p_strokes, p_ref):
    """
    Initial greedy stroke alignment algorithm
    Input:
    strokes: list of points in the gene strokes
    ref: list of points in the reference strokes
    p_strokes: percentages of progress for each point in the gene strokes
    p_ref: percentages of progress for each point in the reference strokes
    Output:
    stroke_map: a mapping for the strokes from the gene to the archetype character
    """
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
    stroke_map = np.full(len(strokes), -1)

    for each in ref_lengths:
        largestref = np.argmax(ref_lengths) # this is the index for the reference stroke that is largest
        smallerror = np.argmin(error_maps[largestref]) # access the error map from the largest stroke's index and see which error is smallest
        #print(error_maps[largestref])
        while(stroke_map[smallerror]!=-1):
            # change small error so that we do not repeat over indexes that are already taken
            # just keeps repeating until we land on an index that doesn't already have a value in its place
            error_maps[largestref][smallerror] = 10000
            smallerror = np.argmin(error_maps[largestref])
        #print(largestref, smallerror)

        stroke_map[smallerror] = largestref # set the index in the stroke_map to the reference stroke we designated

        ref_lengths[largestref] = 0 # set the length of the reference stroke to 0 so we never use it again

    # get rid of invalid (-1) values
    stroke_map = np.array([n for n in stroke_map if n != -1])

    return stroke_map

def strokeError(stroke, ref_stroke, p_stroke, p_ref):
    """
    Calculates the error between two particular strokes
    Input:
    strokes: list of points in the gene stroke
    ref: list of points in the reference stroke
    p_strokes: percentage of progress for each point in the gene strokes
    p_ref: percentage of progress for each point in the reference strokes
    """
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
    #print("Stroke error options")
    #print(forward_stroke_error, forward_ref_error, back_stroke_error, back_ref_error)
    #TODO: This was the problem - we should have just taken the pure minimum of these, not the max between the stroke and ref
    #Apparently the geometry can lead to some messy scenarios, the original line is below
    #min(min(forward_stroke_error.max(), forward_ref_error.max()), min(back_stroke_error.max(), back_ref_error.max()))
    final_error = min(min(forward_stroke_error.max(), forward_ref_error.max()), min(back_stroke_error.max(), back_ref_error.max()))
    return final_error

def strokeErrorHybrid(stroke, ref_stroke, p_stroke, p_ref, w1=0.5, w2=0.5):
    """
    Calculates the error between two particular strokes
    Input:
    strokes: list of points in the gene stroke
    ref: list of points in the reference stroke
    p_strokes: percentage of progress for each point in the gene strokes
    p_ref: percentage of progress for each point in the reference strokes
    w1: weight to give to the first technique (max dev) in hybrid mode
    w2: weight to give to the second technique (mean dev) in hybrid mode
    """
    avg_error = strokeErrorScaled(stroke, ref_stroke, p_stroke, p_ref)
    max_error = strokeError(stroke, ref_stroke, p_stroke, p_ref)
    final_error = w1*avg_error+w2*max_error
    return final_error

def strokeErrorScaled(stroke, ref_stroke, p_stroke, p_ref):
    """
    Calculates the error between two particular strokes
    Scales the strokes on top of each other
    Input:
    strokes: list of points in the gene stroke
    ref: list of points in the reference stroke
    p_strokes: percentage of progress for each point in the gene strokes
    p_ref: percentage of progress for each point in the reference strokes
    """
    center_ref = np.mean(ref_stroke, axis=0)
    center_stroke = np.mean(stroke, axis=0)
    center_shift = center_ref-center_stroke
    stroke = stroke.copy()+center_shift
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

def strokeErrorMatrix(strokes, ref, p_strokes, p_ref, error_function=strokeError):
    """
    Calculates a matrix for errors between strokes to be used by alignment algorithms
    Input:
    strokes: list of points in the gene strokes
    ref: list of points in the reference strokes
    p_strokes: percentages of progress for each point in the gene strokes
    p_ref: percentages of progress for each point in the reference strokes
    Output:
    error_matrix: matrix containing the matching error
    between every stroke pair of form (archetype_stroke, gene_stroke)
    """
    error_map = np.zeros((len(ref), len(strokes)), dtype=float)
    for i, ref_stroke, r_progresses in zip(range(len(ref)), ref, p_ref):
        for j, candidate_stroke, c_progresses in zip(range(len(strokes)), strokes, p_strokes):
            error_map[i, j] = strokeError(ref_stroke, candidate_stroke, r_progresses, c_progresses)
    #print(error_map)
    return error_map

def strokeTrace(stroke, stroke_progresses, progress):
    """
    Gives the point along a stroke after tracing a percentage of progress
    Input:
    stroke: list of points in the stroke
    stroke_progress: percentages of progress for each point in the stroke
    progress: percentage of the way to trace the stroke to find the point
    """
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


"""
The following functions are not currently working with the current pipeline
They were meant to be different approaches to the heuristic matching algorithm
and may still be useful as references for such, however.
They are not fully functional or documented
though most have input/output structures defined.
"""

def alignStrokeOptions(strokes, ref, p_strokes, p_ref, remaining_options=2):
    """
    Greedy stroke alignment algorithm which gives options for alignments,
    filtering out low scores with heuristics
    TODO: Implementation not complete
    Input:
    strokes: list of points in the gene strokes
    ref: list of points in the reference strokes
    p_strokes: percentages of progress for each point in the gene strokes
    p_ref: percentages of progress for each point in the reference strokes
    mode: the technique to use for error calculation
    Output:
    stroke_maps: a set of potential mappings for the strokes from the gene to the archetype character
    """
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
    stroke_options = []

    for _ in ref_lengths:
        largestref = np.argmax(ref_lengths) # this is the index for the reference stroke that is largest
        largest_errors = np.argmax(error_maps[largestref]) # access the error map from the largest stroke's index and see which error is smallest

        for _ in range(len(ref)-remaining_options):
            # change small error so that we do not repeat over indexes that are already taken
            # just keeps repeating until we land on an index that doesn't already have a value in its place
            error_maps[largestref][largest_errors] = 10000
            largest_errors = np.argmax(error_maps[largestref])

    # get combinations of the remaining strokes
    stroke_map[smallerror] = largestref # set the index in the stroke_map to the reference stroke we designated

    ref_lengths[largestref] = 0 # set the length of the reference stroke to 0 so we never use it again

    return stroke_map

def alignOptions(stroke_errors, chosen = None, priority=None, noptions = 2):
    nstrokes = stroke_errors.shape[0]
    if -1 not in chosen:
        return chosen
    if chosen is None:
        chosen = np.empty(nstrokes)
        chosen.fill(-1)
    if priority is None:
        priority = range(stroke_errors.shape[0])
    for i in priority:
        for _ in range(noptions):
            small_error = np.argmin(stroke_errors[i])
            stroke_errors[i][small_error] = 10000
            chosen[smallerror] = i
            alignOptions(stroke_errors, chosen, priority, noptions)

def greedyAlign2(strokes, ref, p_strokes, p_ref, mode="max"):
    """
    Second greedy algorithm which uses a different priority for matching strokes
    Input:
    strokes: list of points in the gene strokes
    ref: list of points in the reference strokes
    p_strokes: percentages of progress for each point in the gene strokes
    p_ref: percentages of progress for each point in the reference strokes
    mode: the technique to use for error calculation
    Output:
    stroke_map: a mapping for the strokes from the gene to the archetype character
    """
    error_maps = strokeErrorMatrix(strokes, ref, p_strokes, p_ref)

    # -1 just means unmatched here since 0 (the other 'default' filler) is a meaningful number in this context
    stroke_map = np.full(len(strokes), -1)

    for _ in stroke_map:
        ref_i, smallerror = np.unravel_index(np.argmin(error_maps), error_maps.shape)
        # access the error map from the largest stroke's index and see which error is smallest, and to which ref stroke it belongs
        # make all of the selected reference index's errors large to prevent reuse
        error_maps[ref_i] = 10000

        stroke_map[smallerror] = ref_i # set the index in the stroke_map to the reference stroke we chose

    # get rid of invalid (-1) values
    stroke_map = np.array([n for n in stroke_map if n != -1])

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
