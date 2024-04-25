from genet.predict.functional import *
from genet.predict.Nuclease import *
from genet.predict.PrimeEditor import *

#### Under development ###
# from genet.predict.DeepSmallCas9 import *
# from genet.predict.functional_dev import *

#### Will be removed ###
# from genet.predict.DeepCas9Variants import *

from silence_tensorflow import silence_tensorflow
import warnings

silence_tensorflow()
warnings.filterwarnings(action='ignore')

print('''[Info] From GenET >= 0.15.0, the input format of DeepPrime has changed significantly. 
       Please refer to the GenET documentation (https://goosang-yu.github.io/genet/). 
       This message will be removed in future updates.''')