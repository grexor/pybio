import pybio
import os
import sys

version = 0.3

def init():
    config_module = sys.modules[__name__]
    pybio_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..")) 
    config_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config"))
    config_example_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config.example"))
    if not os.path.exists(config_file):
        print("[pybio config] Please first check configuration parameters.")
        change()
    config_lines = open(config_file).readlines()
    for cline in config_lines:
        if cline.startswith("#"):
            continue
        k, v = cline.split("=")
        if k.find("folder")!=-1:
            if not v.startswith("/"):
                v = "\"" + os.path.join(pybio_folder, eval(v)) + "\""
        setattr(config_module, k, eval(v))

def change():
    config_module = sys.modules[__name__]
    pybio_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..")) 
    config_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config"))
    config_example_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config.example"))
    if not os.path.exists(config_file):
        os.system(f"cp {config_example_file} {config_file}")
    config_lines = open(config_file).readlines()
    print(f"Note: all paths are relative to the pybio home folder {pybio_folder}.")
    print()
    new_config = []
    for cline in config_lines:
        if cline.startswith("#"):
            continue
        k, v = cline.split("=")
        v = v.replace("\n", "").replace("\r", "").replace("\"", "").replace("'", "")
        v = os.path.expanduser(v)
        v = input(f"{k} [{v}]: ") or v
        v = os.path.expanduser(v)
        new_config.append((k, str(v)))
        setattr(config_module, k, v)
    f = open(config_file, "wt")
    for (k, v) in new_config:
        f.write(f"{k}=\"{v}\"\n")
    f.close()
