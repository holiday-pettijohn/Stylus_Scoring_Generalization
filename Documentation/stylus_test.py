#!/bin/python3

import argparse
import os
import re
import xmltodict

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import Line2D
from matplotlib import collections as mc
from os import path
from PIL import Image, ImageDraw

import stylusengine

stylusengine.setLogFile(b'errors.log')

stylusengine.setScope(
    b'file:///home/tulip/Documents/College/Stewart/stylusapp/hans',
    b'file:///home/tulip/Documents/College/Stewart/stylus/schemas'
)

##stylusengine.setScope(
##    f'file:///home/{getlogin()}/Stylus_Scoring_Generalization/Reference'.encode("UTF8"),
##    f'file:///home/{getlogin()}/stylus/schemas'.encode("UTF8")
##)
