# setup.py
from setuptools import setup, find_packages

setup(
    name="prism",
    version="1.0.0",
    description="PCR primer design & optimization pipeline",
    author="Ao Wang",
    author_email="wang.ao@ufl.edu",
    url="https://github.com/William-A-Wang/PRISM",
    license="GPL-3.0",
    packages=find_packages(where="Main"),     
    package_dir={"": "Main"},                  
    install_requires=[
        "primer3-py",
        "numpy",
        "pandas",
        "tqdm",
        "numba",
        "joblib",
    ],
    entry_points={
        "console_scripts": [
            "prism=main:main",  
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
