import pandas as pd


class MAGeCKAnalyzer:
    def __init__(self, summary_file:str, anotation:dict=None):
        """MAGeCK test file (summary file)을 기반으로 DataFrame 형태로 변환해주면서
        다양한 분석에 응용할 수 있는 객체

        Args:
            summary_file (str): MAGeCK test output (summary) file.
            anotation (dict, optional): _description_. Defaults to None.
        """        
        self.data   = pd.read_csv(summary_file, sep='\t').set_index('id')
        self.simple = self._make_df_simple(self.data)

    def _get_pvalue(self):
        '''MAGeCK test 결과 파일에서 p-value를 가져오는 method.
        rlqhswjrdmfh pos / neg score 중에서 작은 것을 가져온다.
        ToDo: p-value / FDR / 등 다른 것들에 대해서 가져오는 것도 찾아보기?
        '''
        pos  = self.data['pos|score']
        neg  = self.data['neg|score']
        pval = [min(pos.loc[i], neg.loc[i]) for i in self.data.index]

        return pval
    
    def _make_df_simple(self, data):
        '''MAGeCK test 결과 파일에서 핵심 값들만 가져온 것 (LFC / p-value).'''

        df_simple = data[['pos|lfc']].copy()
        df_simple['p-value'] = self._get_pvalue()
        df_simple.columns = ['LFC', 'p-value']

        return df_simple

    def volcano(self):
        '''Make volcano plot'''


    def kde(self):
        '''Make KDE joyplot'''


    def positive(self, lfc_cut:float, pvalue_cut:float):
        '''get positive hits depending LFC / p-values'''

    def negative(self, lfc_cut:float, pvalue_cut:float):
        '''get negative hits depending LFC / p-values'''




class MAGeCK:
    def __init__(self, ):
        print('MAGeCK object')




    def count(self,):
        '''MAGeCK count method'''


        