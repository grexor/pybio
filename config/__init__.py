import pybio
import os

version = 1.3

def init():
    config_file = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..", "config", "config.txt"))
    if os.path.exists(config_file):
        config_lines = open(config_file).readlines()        
        for cline in config_lines:
            exec(cline.replace("\r", "").replace("\n", ""))
