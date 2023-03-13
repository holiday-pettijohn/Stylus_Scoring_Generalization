#!/bin/python3

import argparse
import os
import xmltodict

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import Line2D
from matplotlib import collections as mc
from os import path
from PIL import Image, ImageDraw

def scanDir(xdir, out_dir="./HanBitmap", recursive=False, verbose=0, from_recursive=False):
    """
    Iterates over a directory and scans all of the gene files within, outputting them in the form of bitmap (X) and score (y)
    """
    if (int(verbose) >= 1 and not from_recursive) or (int(verbose) >= 2):
        print(f"Scanning directiory {xdir}...")
    dir_list = os.listdir(xdir)
    files_X = {}
    files_y = {}
    gene_count = 0
    if not path.isdir(out_dir):
        os.mkdir(out_dir)
    for f in dir_list:
        if not f.startswith("."):
            if path.isdir(os.path.join(xdir, f)) and recursive:
                if int(verbose) >= 1:
                    print(f"Recursively scanning child directory {f}...")
                files_X_r, files_y_r = scanDir(path.join(xdir, f), out_dir=out_dir, recursive=recursive, verbose=verbose, from_recursive=True)
                for han_char in files_X_r:
                    if han_char not in files_X:
                        files_X[han_char] = b""
                        files_y[han_char] = ""
                    files_X[han_char] += files_X_r[han_char]
                    files_y[han_char] += files_y_r[han_char]
            elif f.endswith(".gene"):
                if int(verbose) >= 2:
                    print(f"Scanning gene file {f}...")
                han_char, score, bitmap = xmlToBitmap(f"{xdir}/{f}")
                if han_char not in files_X or han_char not in files_y:
                    files_X[han_char] = b""
                    files_y[han_char] = ""
                files_X[han_char] += bitmap + b"\n"
                files_y[han_char] += str(score) + "\n"
                gene_count += 1
    if int(verbose) >= 1:
        print(f"Processed {gene_count} files in directory {xdir}")
    if not from_recursive:
        for han_char in files_X:
            assert len(files_X[han_char]) == len(files_X[han_char])
            fX = open(f"{out_dir}/{han_char}_X", "wb")
            fy = open(f"{out_dir}/{han_char}_y", "w")
            # don't do dupe checking since numpy concatenate crashes for unknown reasons when handling datasets of scale
            fX.write(files_X[han_char])
            fy.write(files_y[han_char])
    else:
        return files_X, files_y

def xmlToBitmap(xfile, output_size=(32, 32), border=4):
    """
    Convert an XML character genome file to a bitmap array representing the character
    Saves the fitness score of said genome along with the bitmap
    """
    xml_data = open(xfile, "r").read()
    root = xmltodict.parse(xml_data)
    score = root["genome"]["statistics"]["@score"]
    han_char = root["genome"]["genes"]["gene"]["hanReferences"]["hanReference"]["@unicode"]
    segments = root["genome"]["genes"]["gene"]["segments"]["segment"]
    drawn_char = drawXml(segments, output_size=output_size, border=border)
    return (han_char, score, drawn_char)

def drawXml(segments, output_size=(32, 32), border=4):
    """
    Process and render the coherent strokes in the xml data
    """
    img = Image.new(mode="1", size=output_size)
    draw = ImageDraw.Draw(img)
    minx, miny, maxx, maxy = None, None, None, None
    for segment in segments:
        if segment["@coherent"] == "true":
            for point in segment["point"]:
                if minx is None:
                    minx = float(point["@x"])
                if miny is None:
                    miny = float(point["@y"])
                if maxx is None:
                    maxx = float(point["@x"])
                if maxy is None:
                    maxy = float(point["@y"])
                if float(point["@x"]) > maxx:
                    maxx = float(point["@x"])
                if float(point["@x"]) < minx:
                    minx = float(point["@x"])
                if float(point["@y"]) > maxy:
                    maxy = float(point["@y"])
                if float(point["@y"]) < miny:
                    miny = float(point["@y"])
    for segment in segments:
        if segment["@coherent"] == "true":
            for i in range(len(segment["point"][:-1])):
                x1, y1, x2, y2 = ((float(segment["point"][i]["@x"])-minx)*((output_size[0]-border)/(maxx-minx))+border/2,
                                  (float(segment["point"][i]["@y"])-miny)*((output_size[1]-border)/(maxy-miny))+border/2,
                                  (float(segment["point"][i+1]["@x"])-minx)*((output_size[0]-border)/(maxx-minx))+border/2,
                                  (float(segment["point"][i+1]["@y"])-miny)*((output_size[1]-border)/(maxy-miny))+border/2)
                y1 = -(y1-output_size[1]/2)+output_size[1]/2
                y2 = -(y2-output_size[1]/2)+output_size[1]/2
                draw.line(((x1, y1), (x2, y2)), width=1, fill=1)
    return np.array(img).flatten().tobytes()

def extractBases(xfile):
    """
    Extract essential information including han charactere, bases, stroke base start/stop points, and stroke order from an XML gene file
    """
    xml_data = open(xfile, "r").read()
    root = xmltodict.parse(xml_data)
    bases = root["genome"]["bases"]
    strokes = root["genome"]["genes"]["gene"]["hanReferences"]["hanReference"]["strokes"]["stroke"]
    stroke_series = [(int(s["@baseFirst"]), int(s["@baseLast"])) for s in strokes]
    stroke_order = np.array([int(s["@correspondsTo"]) for s in strokes])
    han_char = root["genome"]["genes"]["gene"]["hanReferences"]["hanReference"]["@unicode"]
    return (han_char, bases, stroke_series, stroke_order)

def minXml(han_char, bases, stroke_bases, stroke_order):
    """
    Convert a set of minimal information for an xml Stylus input into an XML bytestring
    """
    xml = b"<?xml version='1.0' encoding='UTF-8' ?>\n<genome xmlns='http://biologicinstitute.org/schemas/stylus/1.5'>\n"
    xml += bytes(f"<bases>{bases}</bases>\n", "UTF-8")
    xml += bytes(f"<genes>\n<gene baseFirst='1' baseLast='{len(bases)}'>\n<hanReferences>\n<hanReference unicode='{han_char}'>\n<strokes>\n", "UTF-8")
    for i, base in enumerate(stroke_bases):
        first, last = base
        xml += bytes(f"<stroke baseFirst='{first}' baseLast='{last}' correspondsTo='{stroke_order[i]}' />\n", "UTF-8")
    xml += b"</strokes>\n</hanReference>\n</hanReferences>\n</gene>\n</genes>\n</genome>"
    return xml

def xmlToSegments(xfile, output_size=(32, 32), border=4):
    """
    Convert an XML character genome file to a bitmap array representing the character
    Saves the fitness score of said genome along with the bitmap
    """
    xml_data = open(xfile, "r").read()
    root = xmltodict.parse(xml_data)
    score = root["genome"]["statistics"]["@score"]
    han_char = root["genome"]["genes"]["gene"]["hanReferences"]["hanReference"]["@unicode"]
    segments = root["genome"]["genes"]["gene"]["segments"]["segment"]
    drawn_char = drawSegments(segments, output_size=output_size, border=border)
    return (han_char, score, drawn_char)

def drawSegments(segments, output_size=(32, 32), border=4):
    """
    Process and render the coherent strokes in the xml data
    Saves the segments as line seperated byte arrays
    """
    imgs = b""
    minx, miny, maxx, maxy = None, None, None, None
    for segment in segments:
        if segment["@coherent"] == "true":
            for point in segment["point"]:
                if minx is None:
                    minx = float(point["@x"])
                if miny is None:
                    miny = float(point["@y"])
                if maxx is None:
                    maxx = float(point["@x"])
                if maxy is None:
                    maxy = float(point["@y"])
                if float(point["@x"]) > maxx:
                    maxx = float(point["@x"])
                if float(point["@x"]) < minx:
                    minx = float(point["@x"])
                if float(point["@y"]) > maxy:
                    maxy = float(point["@y"])
                if float(point["@y"]) < miny:
                    miny = float(point["@y"])
    for segment in segments:
        if segment["@coherent"] == "true":
            img = Image.new(mode="1", size=output_size)
            draw = ImageDraw.Draw(img)
            for i in range(len(segment["point"][:-1])):
                x1, y1, x2, y2 = ((float(segment["point"][i]["@x"])-minx)*((output_size[0]-border)/(maxx-minx))+border/2,
                                  (float(segment["point"][i]["@y"])-miny)*((output_size[1]-border)/(maxy-miny))+border/2,
                                  (float(segment["point"][i+1]["@x"])-minx)*((output_size[0]-border)/(maxx-minx))+border/2,
                                  (float(segment["point"][i+1]["@y"])-miny)*((output_size[1]-border)/(maxy-miny))+border/2)
                y1 = -(y1-output_size[1]/2)+output_size[1]/2
                y2 = -(y2-output_size[1]/2)+output_size[1]/2
                draw.line(((x1, y1), (x2, y2)), width=1, fill=1)
            imgs += np.array(img).flatten().tobytes() + b"\n"
    return imgs

def xmlToGeometry(xfile, output_size=(32, 32), border=0):
    xml_data = open(xfile, "r").read()
    root = xmltodict.parse(xml_data)
    segments = root["genome"]["genes"]["gene"]["segments"]["segment"]
    seg_list = []
    frac_dists = []
    minx, miny, maxx, maxy = None, None, None, None
    for segment in segments:
        if segment["@coherent"] == "true":
            for point in segment["point"]:
                if minx is None:
                    minx = float(point["@x"])
                if miny is None:
                    miny = float(point["@y"])
                if maxx is None:
                    maxx = float(point["@x"])
                if maxy is None:
                    maxy = float(point["@y"])
                if float(point["@x"]) > maxx:
                    maxx = float(point["@x"])
                if float(point["@x"]) < minx:
                    minx = float(point["@x"])
                if float(point["@y"]) > maxy:
                    maxy = float(point["@y"])
                if float(point["@y"]) < miny:
                    miny = float(point["@y"])
    for segment in segments:
        if segment["@coherent"] == "true":
            points = []
            for p in segment["point"]:
                if (output_size[0] == -1 or output_size[1] == -1):
                    x, y = (float(p["@x"]),
                            float(p["@y"]))
                else:
                    x, y = ((float(p["@x"])-minx)*((output_size[0]-border)/(maxx-minx))+border/2,
                            (float(p["@y"])-miny)*((output_size[1]-border)/(maxy-miny))+border/2,)
                points.append((x, y))
            seg_list.append(np.array(points))
    for segment in seg_list:
        dists = [0]
        for i in range(len(segment)-1):
            dists.append(np.linalg.norm((segment[i]-segment[i+1]))+dists[-1])
        dists = np.array(dists)
        dists /= dists.max()
        frac_dists.append(dists)
    return seg_list, frac_dists

def loadRef(han_char, ref_dir = "Reference"):
    stroke_list = []
    frac_dists = []
    ref_path = f"{ref_dir}/{han_char[0]}000/{han_char}.han"
    ref_xml = open(ref_path, "r").read()
    root = xmltodict.parse(ref_xml)
    bounds = root["hanDefinition"]["bounds"]
    x_min, y_min, x_max, y_max = (float(bounds["@left"]), float(bounds["@bottom"]), float(bounds["@right"]), float(bounds["@top"]))
    scale = (int(x_max-x_min), int(y_max-y_min))
    strokes = root["hanDefinition"]["strokes"]["stroke"]
    for stroke in strokes:
        points = stroke["points"]["forward"]
        point_arr = []
        frac_arr = []
        for point in points["pointDistance"]:
            point_arr.append((float(point["@x"])-x_min,
                              float(point["@y"])-y_min))
            frac_arr.append(float(point["@fractionalDistance"]))
        stroke_list.append(np.array(point_arr))
        frac_dists.append(np.array(frac_arr))
    return stroke_list, frac_dists, scale

def loadGeometry(data_dir, han_char, output_size = (32, 32), f_read = None):
    """
    Loads geometric data about a gene character directly from the XML source
    """
    """
    if f_read is None:
        dir_list = os.listdir(f"{data_dir}/{han_char}")
        dir_list.sort()
    else:
        dir_list = f_read
    g_data = []
    f_names = []
    for f in dir_list:
        flines = open(f"{data_dir}/{han_char}/{f}", "rb").readlines()
        f_names.append(flines[0].decode()[:-1])
    """
    for f in f_read:
        g = xmlToGeometry(f, output_size)
        g_data.append(g)
    return g_data

def loadGeometryFNames(data_dir, han_char, output_size = (32, 32), f_read = None):
    """
    Loads geometric data about a gene characetr directly from the XML source
    """
    if f_read is None:
        dir_list = os.listdir(f"{data_dir}/{han_char}")
        dir_list.sort()
    else:
        dir_list = f_read
    g_data = []
    f_names = []
    for f in dir_list:
        flines = open(f"{data_dir}/{han_char}/{f}", "rb").readlines()
        f_names.append(flines[0].decode()[:-1])
    for f in f_names:
        g = xmlToGeometry(f, output_size)
        g_data.append(g)
    return g_data, f_names

def loadGeometryBases(data_dir, han_char, output_size = (32, 32), f_read = None):
    """
    Loads geometric data about a gene characetr directly from the XML source
    """
    if f_read is None:
        dir_list = os.listdir(f"{data_dir}/{han_char}")
        dir_list.sort()
    else:
        dir_list = f_read
    g_data, han_chars, base_data, stroke_sets, stroke_orders = [], [], [], [], []
    f_names = []
    for f in dir_list:
        flines = open(f"{data_dir}/{han_char}/{f}", "rb").readlines()
        f_names.append(flines[0].decode()[:-1])
    for f in f_names:
        g = xmlToGeometry(f, output_size)
        han_char, bases, strokes, stroke_order = extractBases(f)
        g_data.append(g)
        han_chars.append(han_char)
        base_data.append(bases)
        stroke_sets.append(strokes)
        stroke_orders.append(stroke_order)
    return g_data, han_chars, base_data, stroke_sets, stroke_orders, f_names

def loadScores(data_dir, han_char, f_read = None):
    """
    Loads gene score from XML source
    """
    if f_read is None:
        dir_list = os.listdir(f"{data_dir}/{han_char}")
        dir_list.sort()
    else:
        dir_list = f_read
    scores = []
    f_names = []
    for f in dir_list:
        flines = open(f"{data_dir}/{han_char}/{f}", "rb").readlines()
        f_names.append(flines[0].decode()[:-1])
    for f in f_names:
        xml_data = open(f, "r").read()
        root = xmltodict.parse(xml_data)
        score = root["genome"]["statistics"]["@score"]
        scores.append(score)
    return scores

def scanSegments(xdir, out_dir="./HanBitmap", recursive=False, verbose=0, from_recursive=False, han_i={}):
    """
    Iterates over a directory and scans all of the gene files within, outputting them in the form of bitmap (X) and score (y)
    """
    if (int(verbose) >= 1 and not from_recursive) or (int(verbose) >= 2):
        print(f"Scanning directiory {xdir}...")
    dir_list = os.listdir(xdir)
    han_i_new = {}
    if not path.isdir(out_dir):
        os.mkdir(out_dir)
    for f in dir_list:
        if not f.startswith("."):
            if path.isdir(os.path.join(xdir, f)) and recursive:
                if int(verbose) >= 1:
                    print(f"Recursively scanning child directory {f}...")
                han_i = scanSegments(path.join(xdir, f), out_dir=out_dir, recursive=recursive, verbose=verbose, from_recursive=True, han_i=han_i)
            elif f.endswith(".gene"):
                if int(verbose) >= 2:
                    print(f"Scanning gene file {f}...")
                han_char, score, bitmap = xmlToSegments(f"{xdir}/{f}")
                if not os.path.isdir(f"{out_dir}/{han_char}"):
                    os.mkdir(f"{out_dir}/{han_char}")
                if not han_char in han_i:
                    han_i[han_char] = 0
                if not han_char in han_i_new:
                    han_i_new[han_char] = 0
                dfile = open(f"{out_dir}/{han_char}/{han_i[han_char]}", "wb")
                dfile.write(bytes(f"{xdir}/{f}", "UTF-8")+b"\n"+bytes(score, "UTF-8")+b"\n"+bitmap)
                han_i[han_char] += 1
                han_i_new[han_char] += 1
    if int(verbose) >= 1:
        total_files = 0
        for han_char in han_i_new:
            total_files += han_i_new[han_char]
        print(f"Processed {total_files} files in directory {xdir}")
        if not from_recursive:
            final_total_files = 0
            for han_char in han_i:
                final_total_files += han_i[han_char]
            print(f"Processed a total of {final_total_files} files")
    return han_i

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse the XML gene files in a given directory and output the binary results in the output directory")
    parser.add_argument("--in_dir", nargs="?", const=1, default="Genes", help="The directory with the XML files you want to scan")
    parser.add_argument("--out_dir", nargs="?", const=1, default="HanBitmap", help="The directory where you want to output the binary data")
    parser.add_argument("--recursive", "-r", action=argparse.BooleanOptionalAction, help="Whether to recursively search through the input directory for files")
    parser.add_argument("--verbose", "-v", nargs="?", const=1, default=0, help="Whether to print additional information as the program runs")
    parser.add_argument("--whole", "-w", default=False, action=argparse.BooleanOptionalAction, help="Whether to render whole characters instead of segments")
    args = parser.parse_args()
    print(f"Processing gene directory {args.in_dir}, writing results to {args.out_dir}")
    if args.whole:
        scanDir(args.in_dir, args.out_dir, recursive=args.recursive, verbose=int(args.verbose))
    else:
        scanSegments(args.in_dir, args.out_dir, recursive=args.recursive, verbose=int(args.verbose))
    print("Done")

