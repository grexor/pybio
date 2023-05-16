import os
import sys

def init(genomes_folder=None):
    config_module = sys.modules[__name__]
    pybio_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", "..")) 
    config_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config"))
    config_example_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config.example"))
    if not os.path.exists(config_file):
        config_example_file = os.path.abspath(os.path.join(pybio_folder, "pybio.config.example"))
        os.system(f"cp {config_example_file} {config_file}")
    config_lines = open(config_file).readlines()
    new_config = []
    for cline in config_lines:
        if cline.startswith("#"):
            continue
        k, v = cline.split("=")
        if k=="genomes_folder" and genomes_folder!=None:
            v = genomes_folder
        v = v.replace("\n", "").replace("\r", "").replace("\"", "").replace("'", "")
        v = os.path.expanduser(v)
        if k.find("folder")!=-1:
            if not v.startswith("/"):
                v = "\"" + os.path.join(pybio_folder, eval(v)) + "\""
        setattr(config_module, k, v)
        new_config.append((k, str(v)))
    if genomes_folder!=None:
        f = open(config_file, "wt")
        for (k, v) in new_config:
            f.write(f"{k}=\"{v}\"\n")
        f.close()
        print("[pybio] Pybio genome folder changed to '{genomes_folder}'")

