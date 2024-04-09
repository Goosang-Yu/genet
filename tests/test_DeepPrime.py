from genet.predict import DeepPrime



def test_deepprime():
    """Test DeepPrime for substitution."""
    
    sID = ''
    ref_seq=''
    ed_seq =''
    
    pegrna = DeepPrime(sID=sID, Ref_seq=ref_seq, Ed_seq=ed_seq,
                       edit_type='sub', edit_len=1)
    
    results = pegrna.predict(pe_system='PE2max', cell_type='HEK293T')
    
    assert results == expected, "Einstein quote counted incorrectly!"