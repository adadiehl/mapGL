import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mapGL",
    version="0.0.1",
    author="Adam Diehl",
    author_email="adadiehl@umich.edu",
    description="Prediction of lineage-specific gain and loss of sequence elements using phylogenetic maximum parsimony.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adadiehl/mapGL",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
