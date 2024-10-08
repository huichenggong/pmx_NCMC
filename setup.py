from setuptools import setup, find_packages
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name='pmxNCMC',
    version=get_version("pmxNCMC/__init__.py"),
    author='Your Name',
    author_email='chenggong.hui@mpinat.mpg.de',
    description='IO based implementation of NCMC in pmx',
    packages=find_packages(),
    install_requires=["python>=3.7",
                      "pymbar>=3, <5",
                      "numpy",
                      "pandas",
                      "matplotlib",
                      "scipy>=1.7.0"],
    # scripts=['script/pmx_mdrun.py',
    #          'script/analysis_bar.py',
    #          ],
    entry_points={
        'console_scripts': [
            'pmx_mdrun=pmxNCMC.pmx_mdrun:main',
            'analysis_bar=pmxNCMC.analysis_bar:main',
        ],
    },
    classifiers=['Programming Language :: Python :: 3',],
)
