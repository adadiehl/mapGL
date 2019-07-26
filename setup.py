import sys

if sys.version_info < (2, 6):
    sys.exit("ERROR: mapGL requires Python 2.6 or greater")
elif sys.version_info > (3, ) and sys.version_info < (3, 3):
    sys.exit("ERROR: mapGL requires Python 3.3 or greater")

try:
    from setuptools import setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages
from pathlib import Path
from glob import glob

def readme():
    with open(Path(__file__).parent.resolve() / 'README.md', encoding='utf-8') as md:
        return md.read()

def main():

    metadata = dict(
        name="mapGL",
        version="0.0.3",
        author="Adam Diehl",
        author_email="adadiehl@umich.edu",
        description="Prediction of lineage-specific gain and loss of sequence elements using phylogenetic maximum parsimony.",
        long_description=readme(),
        long_description_content_type="text/markdown",
        url="https://github.com/adadiehl/mapGL",
        packages = find_packages(),
        package_data = {
            'mapGL': [
                "LICENSE",
                "CODE_OF_CONDUCT.md"
            ]
        },
        setup_requires=[
            'numpy',
            'bx-python',
            'six',
            'cython'
        ],
        install_requires=[
            'numpy',
            'six',
            'bx-python',            
        ],
        entry_points={
            'console_scripts': [
                'mapGL.py = map_GL.mapGL:main'
            ]
        },
        classifiers=[
            "Development Status :: 4 - Beta",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Natural Language :: English"
        ],
        keywords = "phylogenetics, evolution, ancestral reconstruction",
        include_package_data=True,
        zip_safe=False,
    )
            
    setup(**metadata)

     
# ---- Monkey patches -------------------------------------------------------
def monkey_patch_numpy():
    # Numpy pushes its tests into every importers namespace, yeccch.
    try:
        import numpy
        numpy.test = None
    except:
        pass
        
if __name__ == "__main__":
    monkey_patch_numpy()
    main()
