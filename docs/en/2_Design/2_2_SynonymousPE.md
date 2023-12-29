## Introduce Synonymous Mutation in pegRNA
#### Additioal edit for efficient prime editing

Prime editing은 RT-PBS 외에도 다양한 요인들에 의해 효율이 결정된다. 대표적으로 잘 알려져있는 것은 mismatch repair (MMR) system에 의한 저해효과이다. MMR에 의한 prime editing efficiency 저해를 방지하기 위해, pegRNA에 추가적인 mutation을 도입하는 전략을 사용할 수 있다 ([Chen et al., 2021, Cell](https://doi.org/10.1016/j.cell.2021.09.018)). 

하지만, 추가적인 mutation을 도입하는 것은 MMR의 영향을 줄이는 효과와 함께 pegRNA의 activity를 떨어뜨리는 효과도 있을 수 있다. 또한, 만약 genome editing을 하려는 위치가 protein coding sequence (CDS) 영역이라면 단백질의 기능에 영향을 주지 않는 synonymous (silent) mutation을 도입할 필요가 있다. 위와 같은 내용들을 고려했을 때, 어떤 mutation을 추가로 도입하는 pegRNA를 선택할지 디자인하고 선택하는 것은 꽤나 번거로운 작업이다.


#### SynonymousPE module in GenET

GenET의 `SynonymousPE`는 추가 mutation이 도입된 pegRNA를 손쉽게 디자인 할 수 있는 기능을 제공한다. 특히, CDS에 맞춰서 가능한 synonymous mutation을 찾아주므로, 생물학적 연구에 활용할 때 유용하게 사용될 수 있다. 

또한 `SynonymousPE`는 `DeepPrime`과 직접적으로 호환되기 때문에, `DeepPrime`에서 디자인한 pegRNA에 바로 synonymous mutation을 도입한 pegRNA를 만들 수 있다. 


```python
from genet.predict import DeepPrime
from genet.design import SynonymousPE

# DeepPrime pipeline
seq_wt   = 'CTTGCCTGTCTCTGTGGGCTGAAGGCTGTTCCCTGTTTCCTTCAGCTCTACGTCTCCTCCGAGAGCCGCTTCAACACCCTGGCCGAGTTGGTTCATCATCATTCAACGGTGGCCGACGGGC'
seq_ed   = 'CTTGCCTGTCTCTGTGGGCTGAAGGCTGTTCCCTGTTTCCTTCAGCTCTACGTCTCCTCCAAGAGCCGCTTCAACACCCTGGCCGAGTTGGTTCATCATCATTCAACGGTGGCCGACGGGC'

pegrna = DeepPrime('ABL1_ex4_pos21G_A', seq_wt, seq_ed, edit_type='sub', edit_len=1)

pe2max_output = pegrna.predict(pe_system='PE2max', cell_type='HEK293T')


# Select a pegRNA record that you want to add synonymous mutation 
dp_record = pe2max_output.iloc[9]

# Setup SynonymousPE input parameters 
synony_pegrna = design.SynonymousPE(dp_record, ref_seq=seq_wt,
                                    frame=0, cds_start=45, cds_end=121)

# print selected RTT containing synonymous mutation
print(synony_pegrna.extension)
```

DeepPrime의 사용법에 대해서는 `genet.predict` module의 [documentation](/docs/en/1_Predict/predict_pe.md)에서 더 자세한 내용을 볼 수 있다. 위 예시에서는 디자인 된 수 많은 pegRNA 중에서 한 개의 pegRNA를 선택해서 synonymous mutation을 도입하였다. 각각의 pegRNA마다 RTT의 영역과 길이가 다르므로, 도입할 수 있는 additional mutation의 종류도 달라질 수 있다. 우선 DeepPrime score를 기준으로 적절한 pegRNA를 선정한 후, `SynonymousPE`를 추가로 활용하여 optimization을 해서 사용할 것을 권장한다. 

#### SynonymousPE의 parameters







