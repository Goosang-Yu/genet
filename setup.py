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
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX"
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
        'torch==1.11.0+cu113',
        'torchvision==0.12.0+cu113',
        'torchaudio==0.11.0',
        'protobuf==3.20.*',
        'silence_tensorflow'
        ],

    dependency_links=[
        'https://download.pytorch.org/whl/cu113',
        ],

    long_description = long_description,
    long_description_content_type = "text/markdown",
    project_urls={"Bug Tracker": "https://github.com/Goosang-Yu/genet/issues"},
    
    include_package_data=True,
)