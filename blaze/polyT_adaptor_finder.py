"""Searching for adaptor polyT or poly A in a Nanopore read.
"""
import numpy as np
from fast_edit_distance import sub_edit_distance

import blaze.helper as helper
from blaze.config import *

class Read(object):
    """class storing all information of read and locating the putative barcodes

    Initial Attributes:
        id: read id
        seq: read sequence
        phred_score: phred qscore
        _strand: strand (transcript strand: +)

    Methods:
        add()
            add attribute
        find_adaptor()
            find adaptor from a read
        get_strand_and_raw_bc()


    """
    def __init__(self, read_id, sequence, phred_score=None, strand = None, kit=None, **kwargs):
        self.id = read_id
        self.seq = sequence
        self.phred_score = phred_score
        self._strand = strand
        self._get_strand_and_raw_bc_flag = False # record whether the raw bc searching has been performed
        self.kit = kit
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
            
    def add(self, attr_name, attr_val, overwrite = True):
        if attr_name in __dict__.keys() and not overwrite:
            pass
        else:    
            self.__dict__[attr_name] = attr_val
    
    def find_adapter_5_prime(self, read=None, strand=None, adapter_seq=ADPT_SEQ,
                             num_nt=ADPT_WIN,max_ed=ADPT_MAC_MATCH_ED,
                             tso_seq=TSO_SEQ):
        
        strand = self._strand if not strand else strand
        read = self.seq if not read else read

        if strand == '-':
            seq = read[:num_nt]
            d1, d2 = SEQ_SUFFIX_AFT_ADPT
            _, tso_seq_end = sub_edit_distance(seq, tso_seq, max_ed)
            if tso_seq_end > 0:
                tso_seq_start = max(0, tso_seq_end-len(TSO_SEQ))

                _, sub_seq_end = sub_edit_distance(
                    seq[max(0,tso_seq_start-d2):max(0,tso_seq_start-d1)],
                    adapter_seq, max_ed)

                return {'-':[max(0,tso_seq_start-d2)+sub_seq_end+1]} if sub_seq_end > 0 else {}
            else:
                return {}
        
        if strand == '+':
            read = helper.reverse_complement(read[-num_nt:])
            adpt_ends = self.find_adapter_5_prime(
                read=read, strand='-', adapter_seq=adapter_seq, num_nt=num_nt,
                max_ed=max_ed, tso_seq=tso_seq).get('-', [])

            return {'+':adpt_ends} if len(adpt_ends) else {}
        else:
            
            fwd_strand = self.find_adapter_5_prime(
                strand='-', adapter_seq=adapter_seq, 
                num_nt=num_nt, tso_seq=tso_seq)
            rev_strand = self.find_adapter_5_prime(
                strand='+', adapter_seq=adapter_seq, 
                num_nt=num_nt, tso_seq=tso_seq)

            rst = {**{k:v for k,v in fwd_strand.items() if len(v)},
                   **{k:v for k,v in rev_strand.items() if len(v)}}
            return rst

    def find_adaptor_3_prime(self, read=None, strand=None, adaptor_seq=ADPT_SEQ, 
                        num_nt=ADPT_WIN, max_ed=ADPT_MAC_MATCH_ED):
        
        def find_poly_T(seq, poly_T_len=PLY_T_LEN, 
                              min_match_prop=SEQ_SUFFIX_MIN_MATCH_PROP):
            '''Find poly T in seq
            Parameters
            ----------
            seg : STR
            poly_T_len : INT
                How many Ts are used to match in the read. Default: 4.
            min_match_prop : float
                min proportion of T in the polyT. Default: 1
            Returns
                bool
                True if poly T find, otherwise False
            -------
            '''

            # check input
            try:
                assert 0<=min_match_prop <= 1
            except:
                raise AssertionError('value of min_match_prop must between 0-1.')

            # convert to np_array T -> 1 and ACG -> 0
            read_code = np.array([int(x == 'T') for x in seq])
            T_prop = helper.sliding_window_mean(read_code, poly_T_len)
            return np.where(T_prop >= min_match_prop)[0]  

        strand = self._strand if not strand else strand
        read = self.seq if not read else read

        if strand == '-':
            seq = read[:num_nt]
            # find polyT
            ply_T_idx = find_poly_T(seq)
            d1, d2 = SEQ_SUFFIX_AFT_ADPT
            ply_T_idx = ply_T_idx[ply_T_idx >= d1]
            adpt_ends = []
            searching_win = []
            if len(ply_T_idx):
                # get searching window
                win_start, win_end = ply_T_idx[0]-d2, ply_T_idx[0]-d1
                for idx in ply_T_idx:
                    if idx-d2 > win_end:
                        searching_win.append((win_start, win_end))
                        win_start, win_end = idx-d2, idx-d1
                    else: 
                        win_end = idx-d1
                searching_win.append((win_start, win_end))

                for start, end in searching_win:
                    # check adptor seqence d1~d2 upstream to the polyT
                    _, sub_seq_end = sub_edit_distance(seq[max(start,0):end], adaptor_seq, max_ed)
                    # adjust the adaptor end position in related to the read not the flanking sequencce
                    if sub_seq_end > 0:
                        adpt_ends.append(max(start,0)+sub_seq_end+1)

            return {'-':list(set(adpt_ends))} if len(adpt_ends) else {}    

        
        # take reverse complement if read is coming from transcript strand (with ployA instead ployT)
        if strand == '+':
            read = helper.reverse_complement(read[-num_nt:])
            adpt_ends = self.find_adaptor_3_prime(
                            read, strand='-', adaptor_seq=adaptor_seq, 
                            num_nt=num_nt).get('-',[])
            return {'+':adpt_ends} if len(adpt_ends) else {}
                
        else:
            T_strand = self.find_adaptor_3_prime(
                strand='-', adaptor_seq=adaptor_seq, 
                num_nt=num_nt)
            A_strand = self.find_adaptor_3_prime(
                strand='+', adaptor_seq=adaptor_seq, 
                num_nt=num_nt)
            rst = {**{k:v for k,v in T_strand.items() if len(v)},
                   **{k:v for k,v in A_strand.items() if len(v)}}
            return rst

    
    def find_adaptor(self,read=None, strand=None, adaptor_seq=ADPT_SEQ, 
                        num_nt=ADPT_WIN, max_ed=ADPT_MAC_MATCH_ED, tso_seq=TSO_SEQ):
        '''
        find adaptor from a read
    
        Parameters
        ----------
        read : STR
        num_nt : INT
            The adaptor will be searched within first num_nt bases in the reads.
        strand : '+' or '-'
            Transcript strand, this function will find adaptor in '-' poly T strand
            Will do a reverse complement if strand == '-'
        adaptor_seq : STR
            10X adaptor sequence, the full-length adaptor sequence is 
            'CTACACGACGCTCTTCCGATCT' (22nt). Last 12nt is used by default.
        max_ed : float
            max edit distance in the adaptor alignment. Default: 2
        Returns
        -------
        Dict
        '''

        if self.kit == '3v4' or self.kit == '3v3' or self.kit == '3v2':
            return self.find_adaptor_3_prime(read=read, strand=strand,
                adaptor_seq=adaptor_seq, num_nt=num_nt, max_ed=max_ed)
        else:
            return self.find_adapter_5_prime(read=read, strand=strand,
                adapter_seq=adaptor_seq, num_nt=num_nt, max_ed=max_ed,
                tso_seq=tso_seq)

    def get_strand_and_raw_bc(self):
        '''
        Found the strand of read, raw start (i.e. end of adaptor) and end (
            i.e. start of poly T) position of BC before any correction.

        Add following attribute:
            self.raw_bc_start
            self.raw_bc_end
            self.strand
            self.adaptor_polyT_pass:
                0: if an unambiguous adaptor and polyT found
                otherwise (output None for raw_bc_start, raw_bc_end and strand):
                1: poly T and adaptor not found in any strand
                2: poly T and adaptor found in both strand
                10: multiple adaptor found in 20~50 nt of poly T upstream
        
        Returns
        -------
        None.

        '''
        self._get_strand_and_raw_bc_flag = True
        adapt_dict=self.find_adaptor()

        self.adaptor_polyT_pass = 0
        num_of_strand_find = len(adapt_dict)
        num_of_adptor_find = 0

        # check ambiguous finding
        if num_of_strand_find == 0:
            self.adaptor_polyT_pass += 1
        elif num_of_strand_find == 1:
            # check single location 
            adptors = set(list(adapt_dict.values())[0])
            num_of_adptor_find = len(adptors)
            if num_of_adptor_find == 0:
                raise ValueError(
                    "Bug found: strand found but num_of_adptor_find==0.")
            elif num_of_adptor_find == 1:
                pass
            else:
                self.adaptor_polyT_pass += 10

        elif num_of_strand_find == 2:
            self.adaptor_polyT_pass += 2
        else:
            raise ValueError("Bug found: number of strands not in 0,1,2.")
        if self.adaptor_polyT_pass > 0:
            self.raw_bc_start = None
            self.raw_bc_end = None
            self.strand = None
            self.raw_bc = None
            self.raw_bc_min_q = None

        else:
            try:
                assert num_of_strand_find == 1
                assert num_of_adptor_find == 1
            except:
                raise AssertionError(
                    'Bug found: read passed but raw bc location is still ambiguous')

            self.raw_bc_start = list(adapt_dict.values())[0][0]
            self.strand = list(adapt_dict.keys())[0]
            
            if self.strand == '+':
                self.raw_bc = \
                    helper.reverse_complement(
                            self.seq)[self.raw_bc_start: self.raw_bc_start+DEFAULT_BC_SIZE]
            else: 
                self.raw_bc = self.seq[self.raw_bc_start: self.raw_bc_start+DEFAULT_BC_SIZE]

            if self.phred_score is not None:
                if self.strand == '+':
                    phred_score = self.phred_score[::-1]
                else:
                    phred_score = self.phred_score

                bc_phred_score = phred_score[self.raw_bc_start: self.raw_bc_start+DEFAULT_BC_SIZE]      

                self.raw_bc_min_q = min([ord(x) for x in bc_phred_score]) - 33

            else:
                self.raw_bc_min_q = None
        
        # get poly_T_start_raw and adaptor_end_raw
    
    @property
    def adator_trimming_idx(self):
        """sequencing after trimming the adaptor and UMI
        Returns:
            index of the end of the UMI, negative if it's polyA strand
        """
        if not self._get_strand_and_raw_bc_flag:
            self.get_strand_and_raw_bc()
        if not self._strand:
            return None
        elif self._strand == '+':
            return int(-self.raw_bc_start-DEFAULT_BC_SIZE-self.umi_len)
        else: 
            return int(self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len)

    @property
    def putative_UMI(self):
        if not self._get_strand_and_raw_bc_flag:
            self.get_strand_and_raw_bc()
        if not self._strand:
            return None
        elif self._strand == '+':
            return helper.reverse_complement(
                        self.seq)[self.raw_bc_start+DEFAULT_BC_SIZE: self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len]
        else: 
            return self.seq[
                self.raw_bc_start+DEFAULT_BC_SIZE: self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len]
    
    @property
    def pre_bc_flanking(self):
        if not self._get_strand_and_raw_bc_flag:
            self.get_strand_and_raw_bc()
        if not self._strand:
            return None
        elif self._strand == '+':
            return helper.reverse_complement(
                        self.seq)[self.raw_bc_start-DEFAULT_GRB_FLANKING_SIZE: self.raw_bc_start]
        else: 
            return self.seq[self.raw_bc_start-DEFAULT_GRB_FLANKING_SIZE: self.raw_bc_start]
    
    @property
    def post_umi_flanking(self):
        """sequencing after trimming the adaptor and UMI
        Returns:
            index of the end of the UMI, negative if it's polyA strand
        """
        if not self._get_strand_and_raw_bc_flag:
            self.get_strand_and_raw_bc()
        if not self._strand:
            return None
        elif self._strand == '+':
            return helper.reverse_complement(
                        self.seq)[self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len: \
                                  self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len+ DEFAULT_GRB_FLANKING_SIZE]
        else: 
            return self.seq[
                self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len: \
                    self.raw_bc_start+DEFAULT_BC_SIZE+self.umi_len+DEFAULT_GRB_FLANKING_SIZE]
    
    @property
    def polyT_trimming_idx(self):
        """sequencing after trimming the adaptor and UMI
        Returns:
            index of the end of the polyT, negative if it's polyA strand
        """
        
        
        umi_end_idx = self.adator_trimming_idx
        
        if umi_end_idx is None:
            return None
        
        # take reverse complement if read is coming from transcript strand (with ployA instead ployT)
        if umi_end_idx < 0:
            reversed = True
            seq = helper.reverse_complement(self.seq)
            umi_end_idx = abs(umi_end_idx)
        else: 
            reversed = False
            seq = self.seq

        # find the first appearance of TTTT in seq
        polyT_start = seq.find('TTTT', umi_end_idx)
        if polyT_start == -1:
            return umi_end_idx
        
        read_code = np.array([int(x == 'T') for x in seq])
        for idx, nt in enumerate(read_code[polyT_start:]):
            if nt == 1:
                polyT_start += 1
            elif sum(read_code[polyT_start:polyT_start+10]) >= 7: # hard code definition of polyT. plotT end when >3 of 10 consecutive letter are not T 
                polyT_start += 1
            else:
                break 
        return int(polyT_start) if not reversed else int(-polyT_start)


    @property
    def strand(self):
        if not self._strand:
            # set new strand
            helper.warning_msg(
                """
                No strand information is recorded, try 'get_strand_and_raw_bc' 
                first.
                """)
            return None
        return self._strand
    
    @strand.setter
    def strand(self, v):
        if v == None:
            pass
        elif not isinstance(v, str) or v not in '+-':
            helper.warning_msg("Error: strand must be '+' or '-'")
        else:
            self._strand = v

# test code
def main():
        # test reads
        read_id = 'id_1'
        read_seq1 = 'G'*100+'CTTCCGGTCT'+'C'*10+'TTTTT'+'G'*20 + 'TTTGTTTTGT' + 'G'*500
        read_seq2 = helper.reverse_complement(read_seq1)
        r1 = Read(read_id,read_seq1)
        r1.get_strand_and_raw_bc()

        r2 = Read(read_id,read_seq2)
        r2.get_strand_and_raw_bc()
        print(r1.__dict__)
        print(r2.__dict__)
        # [(x.start,x.end, x.score) for x in r1.find_adaptor()['-']]
        # [(x.start,x.end, x.score) for x in r1.find_adaptor().get('+', [])]
        # r1.find_poly_T(strand = '-')
        # r1.find_poly_T(strand = '+')
        # r1.find_poly_T()
        # r1.get_strand_and_raw_bc()
        # r1.raw_bc

if __name__ == '__main__':
    main()
    
    
