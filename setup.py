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
    ],

    package_data={
        'genet':[
            # DeepSpCas9
            'genet/predict/models/DeepSpCas9/DeepSpCas9_model.data-00000-of-00001',
            'genet/predict/models/DeepSpCas9/DeepSpCas9_model.index',
            'genet/predict/models/DeepSpCas9/DeepSpCas9_model.meta',

            # DeepPrime-base
            'genet/predict/models/DeepPrime/DeepPrime_base/model_0.pt',
            'genet/predict/models/DeepPrime/DeepPrime_base/model_1.pt',
            'genet/predict/models/DeepPrime/DeepPrime_base/model_2.pt',
            'genet/predict/models/DeepPrime/DeepPrime_base/model_3.pt',
            'genet/predict/models/DeepPrime/DeepPrime_base/model_4.pt',
            'genet/predict/models/DeepPrime/DeepPrime_base/PE2_mean.csv',
            'genet/predict/models/DeepPrime/DeepPrime_base/PE2_std.csv',

            # DeepPrime-FT: HEK293T-PE2max
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_0.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_1.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_2.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_3.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_4.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_5.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_6.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_7.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_8.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_9.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_10.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_11.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_12.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_13.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_14.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_15.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_16.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_17.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_18.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/final_model_19.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/PE2max_mean.csv',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_PE2max/PE2max_std.csv',

            # DeepPrime-FT: HEK293T-NRCH-PE2
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_0.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_1.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_2.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_3.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_4.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_5.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_6.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_7.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_8.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_9.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_10.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_11.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_12.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_13.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_14.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_15.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_16.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_17.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_18.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/final_model_19.pt',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/NRCH_PE2_mean.csv',
            'genet/predict/models/DeepPrime/DeepPrime_FT/DPFT_293T_NTCH_PE2/NRCH_PE2_std.csv',
        ]
    },

    
    install_requires = [
        'regex',
        'biopython',
        'ViennaRNA',
        'tensorflow==2.8.0',
        'torch==1.11.0+cu113',
        'torchvision==0.12.0+cu113',
        'torchaudio==0.11.0',
        ],


    dependency_links=[
        'https://download.pytorch.org/whl/cu113',
        ],

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