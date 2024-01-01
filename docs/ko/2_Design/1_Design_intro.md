## Design CRISPR systems for research
#### Genomic sequences and Guide RNAs

Genome editing을 위해 CRISPR system을 디자인 한다면, 각 system마다의 규칙에 맞는 guideRNA (gRNA)를 설계하는 것이 중요하다. 예를 들어, 흔히 사용되는 SpCas9에 사용되는 gRNA는 NGG PAM을 찾아서 20nt의 spacer sequence를 찾아가도록 만들어줘야 한다. 반면, 다른 종류의 Cas protein (Cas12a, CjCas9 등)은 다른 종류의 PAM과 다른 길이의 spacer를 인식하는 gRNA를 디자인해줘야 한다. 

#### gRNA library design
단순히 하나의 gRNA만 디자인 하는 것은 그리 어렵지 않지만, genome-wide scale로 디자인해야 하는 경우에는 쉽지 않다. 인간의 유전자는 약 2만개 정도 존재하는데, 그 유전자들을 각각 knock-out (KO) 시키는 gRNA들을 선정하기 위한 작업을 해야 한다. 

#### Synonymous mutation 도입
Prime editing에 사용되는 pegRNA에는 추가적인 mutation을 도입하면 editing efficiency를 높일 수 있다. 이때, protein coding sequence에서 amino acid 서열에 영향을 주지 않도록 synonymous mutation을 도입하는 것이 필요하다면, amino acid 서열의 codon frame과 sequence 정보에 맞춰서 pegRNA를 디자인 해야 한다. 이러한 기능을 수행해주는 것이 `SynonymousPE`이다. 자세한 내용은 [synonymous PE documentation](2_2_SynonymousPE.md)에서 볼 수 있다. 
