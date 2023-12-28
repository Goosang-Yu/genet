import genet
import genet.utils
from genet.predict.PredUtils import *


def pe_score(Ref_seq: str, 
            ED_seq: str, 
            sAlt: str,
            sID:str       = 'Sample',
            pe_system:str = 'PE2max',
            cell_type:str = 'HEK293T',
            pbs_min:int   = 7,
            pbs_max:int   = 15,
            rtt_max:int   = 40,
            show_features = False,
            silence:bool  = False,
            ):
    '''### genet.predict.pe_score is deprecated. 
    Please consider genet.predict.DeepPrime instead.'''

    raise DeprecationWarning('''genet.predict.pe_score is deprecated. Please consider genet.predict.DeepPrime instead.''')
    

def pecv_score(cv_record,
               sID:str       = 'Sample',
               pe_system:str = 'PE2max',
               cell_type:str = 'HEK293T',
               pbs_min:int   = 7,
               pbs_max:int   = 15,
               rtt_max:int   = 40
               ):

    '''### genet.predict.pecv_score is deprecated.'''

    raise DeprecationWarning('''genet.predict.pecv_score is deprecated.''')
    
