import os
import xmltodict

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import Line2D
from matplotlib import collections as mc
from PIL import Image, ImageDraw

def scanDir(xdir, out_dir="./HanBitmap"):
    """
    Iterates over a directory and scans all of the gene files within, outputting them in the form of bitmap (X) and score (y)
    """
    dir_list = os.listdir(xdir)
    files_X = {}
    files_y = {}
    for f in dir_list:
        if f.endswith(".gene"):
            han_char, score, bitmap = xmlToBitmap(f"{xdir}/{f}")
            if han_char not in files_X or han_char not in files_y:
                files_X[han_char] = b""
                files_y[han_char] = ""
            files_X[han_char] += bitmap + b"\n"
            files_y[han_char] += str(score) + "\n"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    for han_char in files_X:
        fX = open(f"{out_dir}/{han_char}_X", "wb")
        fX.write(files_X[han_char])
        fy = open(f"{out_dir}/{han_char}_y", "w")
        fy.write(files_y[han_char])
    
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
