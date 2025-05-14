from setuptools import setup, find_packages
from pathlib import Path

this_dir = Path(__file__).parent
long_desc = (this_dir / "README.md").read_text(encoding="utf‑8")

setup(
    name="prism-bio",
    version="1.0.5",                    
    description="PCR primer design & optimization pipeline",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    author="Ao Wang",
    author_email="wang.ao@ufl.edu",
    url="https://github.com/William-A-Wang/PRISM",
    license="GPL-3.0",

    
    package_dir={"": "Main"},
    packages=find_packages(where="Main"),  
    
    py_modules=[
        "main",              
        "badness_utils",
        "data_io",
        "design_primers",
        "iterative",
        "optimization",
    ],

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
