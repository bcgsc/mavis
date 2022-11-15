import pytest

from mavis.annotate.variant import IndelCall, call_protein_indel

from ..mock import Mock, MockFunction


class TestIndelCall:
    def test_deletion(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 4
        assert indel.cterm_aligned == len(indel.ref_seq) - 8 + 1
        assert indel.del_seq == 'ghj'
        assert indel.ins_seq == ''
        assert not indel.is_dup

    def test_insertion(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfmmmghjkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 4
        assert indel.cterm_aligned == len(indel.ref_seq) - 5 + 1
        assert indel.del_seq == ''
        assert indel.ins_seq == 'mmm'
        assert not indel.is_dup

    def test_dup(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfsdfghjkl'
        indel = IndelCall(refseq, mutseq)
        print(indel)
        assert indel.nterm_aligned == 4
        assert indel.cterm_aligned == len(indel.ref_seq) - 2 + 1
        assert indel.del_seq == ''
        assert indel.ins_seq == 'sdf'
        assert indel.is_dup

    def test_delins(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfmmmkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 4
        assert indel.cterm_aligned == len(indel.ref_seq) - 8 + 1
        assert indel.del_seq == 'ghj'
        assert indel.ins_seq == 'mmm'
        assert not indel.is_dup

    def test_delete_start(self):
        refseq = 'asdfghjkl'
        mutseq = 'fghjkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 0
        assert indel.cterm_aligned == 6
        assert indel.del_seq == 'asd'
        assert indel.ins_seq == ''
        assert not indel.is_dup

    def test_delete_start_repetition(self):
        refseq = 'asdafghjkl'
        mutseq = 'afghjkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 0
        assert indel.cterm_aligned == 7
        assert indel.del_seq == 'asd'
        assert indel.ins_seq == ''
        assert not indel.is_dup

    def test_delete_end(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfgh'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 6
        assert indel.cterm_aligned == 0
        assert indel.del_seq == 'jkl'
        assert indel.ins_seq == ''
        assert not indel.is_dup

    def test_ins_start(self):
        refseq = 'asdfghjkl'
        mutseq = 'mmasdfghjkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 0
        assert indel.cterm_aligned == 9
        assert indel.del_seq == ''
        assert indel.ins_seq == 'mm'
        assert not indel.is_dup

    def test_ins_end(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfghjklmmm'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 9
        assert indel.cterm_aligned == 0
        assert indel.del_seq == ''
        assert indel.ins_seq == 'mmm'
        assert not indel.is_dup

    def test_delins_start(self):
        refseq = 'asdfghjkl'
        mutseq = 'mmfghjkl'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 0
        assert indel.cterm_aligned == 6
        assert indel.del_seq == 'asd'
        assert indel.ins_seq == 'mm'
        assert not indel.is_dup

    def test_delins_end(self):
        refseq = 'asdfghjkl'
        mutseq = 'asdfghjmmm'
        indel = IndelCall(refseq, mutseq)
        assert indel.nterm_aligned == 7
        assert indel.cterm_aligned == 0
        assert indel.del_seq == 'kl'
        assert indel.ins_seq == 'mmm'
        assert not indel.is_dup


class TestHgvsProteinNotation:
    def test_homopolymer(self):
        indel = IndelCall('ASDFGHJKKLQWERTYUIOP', 'ASDFGHJKKKKLQWERTYUIOP').hgvs_protein_notation()
        assert indel == 'p.K8_K9dupKK'

    def test_dup(self):
        indel = IndelCall('ASDFGHJKL', 'ASDFSDFGHJKL').hgvs_protein_notation()
        assert indel == 'p.S2_F4dupSDF'


class TestCallProteinIndel:
    def test_large_start_deletion(self):
        ref_translation = Mock(
            get_aa_seq=MockFunction(
                'MGLKAAQKTLFPLRSIDDVVRLFAAELGREEPDLVLLSLVLGFVEHFLAVNRVIPTNVPE'
                'LTFQPSPAPDPPGGLTYFPVADLSIIAALYARFTAQIRGAVDLSLYPREGGVSSRELVKK'
                'VSDVIWNSLSRSYFKDRAHIQSLFSFITGTKLDSSGVAFAVVGACQALGLRDVHLALSED'
                'HAWVVFGPNGEQTAEVTWHGKGNEDRRGQTVNAGVAERSWLYLKGSYMRCDRKMEVAFMV'
                'CAINPSIDLHTDSLELLQLQQKLLWLLYDLGHLERYPMALGNLADLEELEPTPGRPDPLT'
                'LYHKGIASAKTYYRDEHIYPYMYLAGYHCRNRNVREALQAWADTATVIQDYNYCREDEEI'
                'YKEFFEVANDVIPNLLKEAASLLEAGEERPGEQSQGTQSQGSALQDPECFAHLLRFYDGI'
                'CKWEEGSPTPVLHVGWATFLVQSLGRFEGQVRQKVRIVSREAEAAEAEEPWGEEAREGRR'
                'RGPRRESKPEEPPPPKKPALDKGLGTGQGAVSGPPRKPPGTVAGTARGPEGGSTAQVPAP'
                'TASPPPEGPVLTFQSEKMKGMKELLVATKINSSAIKLQLTAQSQVQMKKQKVSTPSDYTL'
                'SFLKRQRKGL*'
            ),
            name='ref',
        )
        mut_translation = Mock(
            get_aa_seq=MockFunction(
                'MRCDRKMEVAFMV'
                'CAINPSIDLHTDSLELLQLQQKLLWLLYDLGHLERYPMALGNLADLEELEPTPGRPDPLT'
                'LYHKGIASAKTYYRDEHIYPYMYLAGYHCRNRNVREALQAWADTATVIQDYNYCREDEEI'
                'YKEFFEVANDVIPNLLKEAASLLEAGEERPGEQSQGTQSQGSALQDPECFAHLLRFYDGI'
                'CKWEEGSPTPVLHVGWATFLVQSLGRFEGQVRQKVRIVSREAEAAEAEEPWGEEAREGRR'
                'RGPRRESKPEEPPPPKKPALDKGLGTGQGAVSGPPRKPPGTVAGTARGPEGGSTAQVPAP'
                'TASPPPEGPVLTFQSEKMKGMKELLVATKINSSAIKLQLTAQSQVQMKKQKVSTPSDYTL'
                'SFLKRQRKGL*'
            )
        )
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == (
            'ref:p.M1_Y227del'
            'MGLKAAQKTLFPLRSIDDVVRLFAAELGREEPDLVLLSLVLGFVEHFLAVNRVIPTNVPE'
            'LTFQPSPAPDPPGGLTYFPVADLSIIAALYARFTAQIRGAVDLSLYPREGGVSSRELVKK'
            'VSDVIWNSLSRSYFKDRAHIQSLFSFITGTKLDSSGVAFAVVGACQALGLRDVHLALSED'
            'HAWVVFGPNGEQTAEVTWHGKGNEDRRGQTVNAGVAERSWLYLKGSY'
        )

    def test_deletion_rep_at_breaks(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ABCDEFKJFEDAGFLKJ'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ABCDE' 'AGFLKJ'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.F6_D11delFKJFED'

    def test_insertion(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKIIILQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8_L9insIII'

    def test_deletion(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8_L9delKL'

    def test_synonymous(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation is None

    def test_delins(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJIIIQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8_L9delKLinsIII'

    def test_transcript_name(self):
        ref_translation = Mock(
            get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'),
            name=None,
            reference_object=Mock(name='reft'),
        )
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJIIIQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'reft:p.K8_L9delKLinsIII'

    def test_delete_start(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('FGHJKLQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.A1_D3delASD'

    def test_delete_single_aa_start(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('SDFGHJKLQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.A1delA'

    def test_delete_end(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYU'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.I17_P19delIOP'

    def test_delete_single_aa_end(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIO'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.P19delP'

    def test_ins_start(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('IIASDFGHJKLQWERTYUIOP'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.A1ext-2'

    def test_ins_end(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOPII'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.P19ext2'

    def test_no_reference_obj(self):
        ref_translation = Mock(
            get_aa_seq=MockFunction('ASDFGHJKLQWERTYUIOP'), name=None, reference_object='thing'
        )
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJIIIQWERTYUIOP'))
        with pytest.raises(AttributeError):
            call_protein_indel(ref_translation, mut_translation)

    def test_fs(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKL'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJMMM'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8Mfs'

        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKL'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJCMMEF'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8Cfs'

    def test_fs_with_stops(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLT*'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJMMMHGFTTSBF*TUHG*'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8Mfs*12'

    def test_fs_immeadiate_stop(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDFGHJKLT*'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('ASDFGHJMMMHGFTTSBF*'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.K8Mfs*12'

    def test_delete_start_with_rep(self):
        ref_translation = Mock(get_aa_seq=MockFunction('ASDAFGHJKL'), name='ref')
        mut_translation = Mock(get_aa_seq=MockFunction('AFGHJKL'))
        notation = call_protein_indel(ref_translation, mut_translation)
        assert notation == 'ref:p.A1_D3delASD'
