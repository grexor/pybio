import os
import glob
import pybio

def init():
    pybio.path.root_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))
    pybio.path.genomes_folder = os.path.join(pybio.path.root_folder, "genomes")
