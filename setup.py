from setuptools import setup
from setuptools import find_packages

setup(
    name='pybio',
    version = "0.3",
    package_dir = {"":"src"},
    packages=find_packages("src"),
    description='pybio genomics',
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
)
