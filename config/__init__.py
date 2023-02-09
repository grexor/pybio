import pybio
import os
import sys

version = 1.3

def init():
    config_module = sys.modules[__name__]
    config_file = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..", "pybio.config"))
    if os.path.exists(config_file):
        config_lines = open(config_file).readlines()
        for cline in config_lines:
            k, v = cline.split("=")
            setattr(config_module, k, eval(v))
            #exec(cline.replace("\r", "").replace("\n", ""))
    else:
        print(f"[config] expected config file at {config_file} not found, using defaults")