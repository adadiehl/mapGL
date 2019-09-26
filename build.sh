#!/bin/bash

# Create a new mapGL release and update all repositories. Takes one
# argument to indicate whether this is a major, minor, or patch release.

# Increment version numbers
NV=$(bump2version $1 | grep new_version | sed -r s,"^.*=",,)

# Update PyPI
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/mapGL-$NV*

# Update BioConda
