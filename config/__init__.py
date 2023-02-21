import pybio
import os
import sys

version = 1.3

def init():
    config_module = sys.modules[__name__]
    pybio_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..")) 
    config_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config"))
    if os.path.exists(config_file):
        config_lines = open(config_file).readlines()
        for cline in config_lines:
            if cline.startswith("#"):
                continue
            k, v = cline.split("=")
            if k.find("folder")!="-1":
                if not v.startswith("/"):
                    v = "\"" + os.path.join(pybio_folder, eval(v)) + "\""
            setattr(config_module, k, eval(v))
            #exec(cline.replace("\r", "").replace("\n", ""))
    else:
        print(f"[config] expected config file at {config_file} not found, using defaults")