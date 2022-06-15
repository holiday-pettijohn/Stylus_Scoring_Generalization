import os

import numpy as np
import tensorflow as tf

from PIL import Image
from sklearn.model_selection import train_test_split

def loadXy(data_dir):
    dir_list = os.listdir(data_dir)
    X_data = {}
    y_data = {}
    for f in dir_list:
        han_char = f.split("_")[0]
        if f.endswith("_X"):
            Xs = open(f"{data_dir}/{f}", "rb").readlines()
            X_data[han_char] = [np.frombuffer(X, dtype=bool)[:-1].reshape(32, 32).astype(int) for X in Xs]
        if f.endswith("_y"):
            ys = open(f"{data_dir}/{f}", "r").readlines()
            y_data[han_char] = np.array([float(a) for a in ys])
    return X_data, y_data

def loadXyCharSplits(data_dir, han_char, test_size = 0.2, random_state = 42):
    X_data, y_data = loadXy(data_dir)
    X_data_char, y_data_char = X_data[han_char], y_data[han_char]
    X_train, X_test, y_train, y_test = train_test_split(
        X_data_char, y_data_char, test_size=test_size, random_state=random_state
    )
    X_train, X_test, y_train, y_test = (tf.convert_to_tensor(X_train, tf.float16),
                                        tf.convert_to_tensor(X_test, tf.float16),
                                        tf.convert_to_tensor(y_train, tf.float16),
                                        tf.convert_to_tensor(y_test, tf.float16),
                                       )
    return X_train, X_test, y_train, y_test

def displayBinary(arr):
    img = Image.fromarray(arr)
    return img
