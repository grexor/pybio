import os
import glob
import pybio

def init():
    if getattr(pybio.path, "root_folder", None)==None:
        pybio.path.root_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))
    if getattr(pybio.path, "genomes_folder", None)==None:
        pybio.path.genomes_folder = os.path.join(pybio.path.root_folder, "genomes")
