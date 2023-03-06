import pybio
import os
import sys

version = 1.3

def init():
    config_module = sys.modules[__name__]
    pybio_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..")) 
    config_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config"))
    config_example_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config.example"))
    if not os.path.exists(config_file):
        os.system(f"cp {config_example_file} {config_file}")
    config_lines = open(config_file).readlines()
    for cline in config_lines:
        if cline.startswith("#"):
            continue
        k, v = cline.split("=")
        if k.find("folder")!=-1:
            if not v.startswith("/"):
                v = "\"" + os.path.join(pybio_folder, eval(v)) + "\""
        setattr(config_module, k, eval(v))
