import setuptools
import genet

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
    
VERSION = genet.__version__

setuptools.setup(
    name="genet",
    version=VERSION,
    author="Goosang Yu",
    author_email="gsyu93@gmail.com",
    description="Genome Editing Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Goosang-Yu/genet",
    install_requires = [
        'pandas'
    ],
    project_urls={
        "Bug Tracker": "https://github.com/Goosang-Yu/genet/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    include_package_data=True,
)