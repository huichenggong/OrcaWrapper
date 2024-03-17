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
    name='orcawrapper',
    version=get_version("orcawrapper/__init__.py"),
    description='Tool sets for read/write orca input file',
    author='Chenggong Hui',
    author_email='chenggong.hui@mpinat.mpg.de',
    packages=find_packages(),
    scripts=['script/orcainp_2_pdb.py',
             'script/orca_update_xyz.py'
             ],
    install_requires=["numpy", "matplotlib", "MDAnalysis", "pandas"],
)
