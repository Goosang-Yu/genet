BINARY_NAME := genet

CONDA_EXISTS := $(shell if [ -x "$$(command -v conda)" ]; then echo "true"; else echo "false"; fi)

ifeq ($(OS),Windows_NT)
	CONDA_ACTIVATE = activate
	ifeq ($(shell where nvcc 2> nul),)
		USE_CPU = true
	endif
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CONDA_ACTIVATE = conda activate
		ifeq ($(shell which nvcc),)
			USE_CPU = true
		else ifeq ($(shell which nvcc | grep -c "cuda"),1)
			USE_CUDA = true
		endif
	else
		CONDA_ACTIVATE = conda activate
		ifeq ($(shell uname -p),arm)
			USE_ARM = true
		else
			ifeq ($(shell which nvcc | grep -c "cuda"),1)
				USE_CUDA = true
			else
				USE_CPU = true
			endif
		endif
	endif
endif

.PHONY: install

install:
ifeq ($(CONDA_EXISTS),false)
	$(error "conda is not installed. Please install conda and try again." \
		"See https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html")
endif
	conda create -n $(BINARY_NAME) python=3.11
ifeq ($(USE_CPU),true)
	$(CONDA_ACTIVATE) $(BINARY_NAME) && pip install torch==2.0.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu
else ifeq ($(USE_CUDA),true)
	$(CONDA_ACTIVATE) $(BINARY_NAME) && pip install torch==2.0.1+cu118 -f https://download.pytorch.org/whl/apple/cpu/torch_stable.html
else ifeq ($(USE_ARM),true)
	$(CONDA_ACTIVATE) $(BINARY_NAME) && pip install torch==2.0.1 -f https://download.pytorch.org/whl/torch_stable.html
endif
	$(CONDA_ACTIVATE) $(BINARY_NAME) && pip install viennarna && pip install .

troubleshoot:
	$(CONDA_ACTIVATE) $(BINARY_NAME) && conda install -c conda-forge mamba
	$(CONDA_ACTIVATE) $(BINARY_NAME) && mamba install libgcc