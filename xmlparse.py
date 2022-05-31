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
    
def dupeCheck(out_dir="./HanBitmap", files=()):
    """
    Checks the given files in the output directory for duplicates
    If no files are specified, all files in the directory, both X and y, will be checked
    """
    if len(files) == 0:
        files = os.listdir(out_dir)
    for f in files:
        if f.endswith("_X"):
            han_char = f.split("_")[0]
            # check for duplicates
            fX = open(f"{out_dir}/{han_char}_X", "rb")
            fX_new = np.array(fX.readlines())
            uniques, u_filter = np.unique(fX_new, return_index=True, axis=0)
            if len(fX_new) == len(uniques):
                fX.close()
            else:
                fX_new = fX_new[u_filter]
                fX_str = b""
                for s in fX_new:
                    fX_str += s
                fX.close()
                fX = open(f"{out_dir}/{han_char}_X", "wb")
                fX.write(fX_str)
                fX.close()
                # only open fy if we made changes to fX
                fy = open(f"{out_dir}/{han_char}_y", "r")
                fy_new = np.array(fy.readlines())[u_filter]
                fy_str = ""
                for s in fy_new:
                    fy_str += s
                fy.close()
                fy = open(f"{out_dir}/{han_char}_y", "w")
                fy.write(fy_str)
                fy.close()
                
def dupeCheckNew(han_char, new_X, new_y, out_dir="./HanBitmap"):
    """
    Checks new data for duplicates and removes them while appending the new data to an existing file
    Seems to have some issues when presented with a large dataset
    """
    files = os.listdir(out_dir)
    if f"{han_char}_X" not in files:
        fX = open(f"{out_dir}/{han_char}_X", "wb")
        fX.write(new_X)
        fy = open(f"{out_dir}/{han_char}_y", "w")
        fy.write(new_y)
        fX.close()
        fy.close()
        return
    # check for duplicates
    fX = open(f"{out_dir}/{han_char}_X", "rb")
    fX_old = fX.readlines()
    fX_new = np.concatenate((np.array(fX_old), np.array(new_X.splitlines())))
    uniques, u_filter = np.unique(fX_new, return_index=True, axis=0)
    if len(fX_new) == len(uniques):
        fX = open(f"{out_dir}/{han_char}_X", "ab")
        fX.write(new_X)
        fy = open(f"{out_dir}/{han_char}_y", "a")
        fy.write(new_y)
        fX.close()
        fy.close()
    elif len(fX_old) == len(uniques):
        pass
    else:
        fX_new = fX_new[u_filter]
        fX_str = b""
        for s in fX_new:
            fX_str += s
        fX.close()
        fX = open(f"{out_dir}/{han_char}_X", "wb")
        fX.write(fX_str)
        fX.close()
        # only open fy if we made changes to fX
        fy = open(f"{out_dir}/{han_char}_y", "r")
        fy_new = np.concatenate((np.array(fy.readlines()), np.array(new_y.splitlines())))[u_filter]
        fy_str = ""
        for s in fy_new:
            fy_str += s
        fy = open(f"{out_dir}/{han_char}_y", "w")
        fy.write(fy_str)
        fX.close()
        fy.close()

def xmlToBitmap(xfile):
    """
    Convert an XML character genome file to a bitmap array representing the character
    Saves the fitness score of said genome along with the bitmap
    """
    xml_data = open(xfile, "r").read()
    root = xmltodict.parse(xml_data)
    score = root["genome"]["statistics"]["@score"]
    han_char = root["genome"]["genes"]["gene"]["hanReferences"]["hanReference"]["@unicode"]
    segments = root["genome"]["genes"]["gene"]["segments"]["segment"]
    drawn_char = drawXml(segments)
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
    img.save("test.png")
    return np.array(img).flatten().tobytes()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse the XML gene files in a given directory and output the binary results in the output directory")
    parser.add_argument("--in_dir", nargs="?", const=1, default="Genes", help="The directory with the XML files you want to scan")
    parser.add_argument("--out_dir", nargs="?", const=1, default="HanBitmap", help="The directory where you want to output the binary data")
    parser.add_argument("--recursive", "-r", action=argparse.BooleanOptionalAction, help="Whether to recursively search through the input directory for files")
    parser.add_argument("--verbose", "-v", nargs="?", const=1, default=0, help="Whether to print additional information as the program runs")
    args = parser.parse_args()
    print(f"Processing gene directory {args.in_dir}, writing results to {args.out_dir}")
    scanDir(args.in_dir, args.out_dir, recursive=args.recursive, verbose=args.verbose)
    print("Done")

