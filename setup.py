import setuptools
import genet
from glob import glob

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
    
VERSION = genet.__version__

setuptools.setup(
    name            = "genet",
    version         = VERSION,
    author          = "Goosang Yu",
    author_email    = "gsyu93@gmail.com",
    description     = "GenET: Genome Editing Toolkit",
    url             = "https://github.com/Goosang-Yu/genet",
    packages        = setuptools.find_packages(exclude = ['dev_ing', 'dev_ing.*']),
    
    python_requires = ">=3.7",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Operating System :: Unix",
        "Topic :: Software Development :: Libraries :: Application Frameworks",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],

    package_data={
        'genet':[
            # DeepSpCas9
            'genet/predict/models/DeepSpCas9/DeepSpCas9_model.data-00000-of-00001',
            'genet/predict/models/DeepSpCas9/DeepSpCas9_model.index',
            'genet/predict/models/DeepSpCas9/DeepSpCas9_model.meta',
        ]
    },
    
    install_requires = [
        'pandas',
        'regex',
        'biopython',
        'tensorflow==2.8.0',
        'protobuf==3.20.*',
        'silence_tensorflow',
        'pyarrow',
        'fastparquet',
        'tqdm',
        'editdistance',
        ],


    long_description = long_description,
    long_description_content_type = "text/markdown",
    project_urls={"Bug Tracker": "https://github.com/Goosang-Yu/genet/issues"},
    
    include_package_data=True,
)