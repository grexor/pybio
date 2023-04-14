from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_desc = fh.read()
    long_desc = "\n".join(long_desc.split("\n")[1:])

setup(
    name='pybio',
    version = "0.3.2",
    package_dir = {"":"src"},
    packages=find_packages("src"),
    description='pybio genomics',
    long_description = long_desc,
    long_description_content_type = "text/markdown",
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    zip_safe=False,
    author='Gregor Rot',
    scripts=["src/pybio/pybio"],
    author_email='gregor.rot@gmail.com',
    url='https://github.com/grexor/pybio',
    keywords=['pybio', 'bioinformatics'],
    include_package_data=True,
    package_data={
        'pybio': ['pybio.config.example'],
    },
    install_requires=["pysam", "numpy", "psutil", "bs4", "requests"],
)
