#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np

from os import getlogin

import stylusengine

stylusengine.setLogFile(b'errors.log')

stylusengine.setScope(
    f'file:///home/{getlogin()}/Stylus_Scoring_Generalization/Reference'.encode("UTF8"),
    f'file:///home/{getlogin()}/stylus/schemas'.encode("UTF8")
)
