import os
import glob
import pybio3

def init():
    pybio3.path.root_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))
    pybio3.path.genomes_folder = os.path.join(pybio3.path.root_folder, "genomes")
