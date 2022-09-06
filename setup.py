import setuptools
import genet

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
    
VERSION = genet.__version__

setuptools.setup(
    name            = "genet",
    version         = VERSION,
    author          = "Goosang Yu",
    author_email    = "gsyu93@gmail.com",
    description     = "Genome Editing Toolkit",
    url             = "https://github.com/Goosang-Yu/genet",
    packages        = setuptools.find_packages(exclude = ['dev_ing', 'dev_ing.*']),
    python_requires = ">=3.6",
    
    install_requires = ['pandas', 'pyensembl', 'biopython'],
    long_description = long_description,
    long_description_content_type = "text/markdown",
    project_urls={"Bug Tracker": "https://github.com/Goosang-Yu/genet/issues"},
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
)