"""Searching for adaptor polyT or poly A in a Nanopore read.
"""
import numpy as np
import Bio
import Bio.pairwise2

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
    def __init__(self, read_id, sequence, phred_score=None, strand = None, **kwargs):
        self.id = read_id
        self.seq = sequence
        self.phred_score = phred_score
        self._strand = strand
        self._get_strand_and_raw_bc_flag = False # record whether the raw bc searching has been performed
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
            
    def add(self, attr_name, attr_val, overwrite = True):
        if attr_name in __dict__.keys() and not overwrite:
            pass
        else:    
            self.__dict__[attr_name] = attr_val
    
    
    def find_adaptor(self,read=None, strand=None, adaptor_seq=ADPT_SEQ, 
                        num_nt=ADPT_WIN, min_match_prop=ADPT_MIN_MATCH_PROP):
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
        min_match_prop : float
            min proportion of matches in the alignment. Default: 0.8
        Returns
        -------
        Bio.pairwise2.Alignment
        *(adaptor_start, adaptor_end): adaptor start and end position (0-based) in reads
        '''
        # check input
        try:
            assert 0<=min_match_prop <= 1
        except:
            raise AssertionError('value of min_match_prop must between 0-1.')

        strand = self._strand if not strand else strand
        read = self.seq if not read else read
        

        def min_score(len_adptor, min_match_prop):
            '''
            Calculation based on the following alignment metrics:
                score for a matched base: 2
                score for a mismatch : -1
                score for open a gap : -1
                score for extent a gap: -1
            The output is rearanged form of :
                2*len_adaptor - 3 * len_adaptor * (1-min_match_prop)                 
            '''
            return len_adptor* (3*min_match_prop - 1)

        def check_poly_T(seq, poly_T_len=PLY_T_LEN, 
                              min_match_prop=PLY_T_MIN_MATCH_PROP):
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
            return np.any(T_prop >= min_match_prop)
    
        if strand == '-':
            seq = read[:num_nt]
            
            # Align adaptor sequencing to seq 
                # no panelty for skipping the start and end of the seq
                # score for a matched base: 2
                # score for a mismatch : -1
                # score for open a gap : -1
                # score for extent a gap: -1
            align = Bio.pairwise2.align.localms(seq, adaptor_seq, 2,-1,-1,-1)
            
            # filter out candidate adaptor when no polyT find
            adp_cand = [a for a in align if a.score >= min_score(len(adaptor_seq), min_match_prop)]
            d1, d2 = PLY_T_NT_AFT_ADPT
            adp_cand = [a for a in adp_cand if check_poly_T(read[a.end+d1:a.end+d2])]
            return {'-':adp_cand} if len(adp_cand) else {}
        
        # take reverse complement if read is coming from transcript strand (with ployA instead ployT)
        if strand == '+':
            read = helper.reverse_complement(read)
            adp_cand = self.find_adaptor(
                            read, strand='-', adaptor_seq=adaptor_seq, 
                            num_nt=num_nt, 
                            min_match_prop=min_match_prop).get('-',[])
            return {'+':adp_cand} if len(adp_cand) else {}
                
        else:
            T_strand = self.find_adaptor(
                strand='-', adaptor_seq=adaptor_seq, 
                num_nt=num_nt, min_match_prop=min_match_prop)
            A_strand = self.find_adaptor(
                strand='+', adaptor_seq=adaptor_seq, 
                num_nt=num_nt, min_match_prop=min_match_prop)
            rst = {**{k:v for k,v in T_strand.items() if len(v)},
                   **{k:v for k,v in A_strand.items() if len(v)}}
            return rst


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
            adptor_align = list(adapt_dict.values())[0]
            num_of_adptor_find = len(set([a.end for a in adptor_align]))
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


            self.raw_bc_start = list(adapt_dict.values())[0][0].end
            self.strand = list(adapt_dict.keys())[0]
            
            if self.strand == '+':
                self.raw_bc = \
                    helper.reverse_complement(
                            self.seq)[self.raw_bc_start: self.raw_bc_start+16]
            else: 
                self.raw_bc = self.seq[self.raw_bc_start: self.raw_bc_start+16]

            if self.phred_score is not None:
                if self.strand == '+':
                    phred_score = self.phred_score[::-1]
                else:
                    phred_score = self.phred_score

                self.raw_bc_min_q = \
                    min(phred_score[self.raw_bc_start: self.raw_bc_start+16])
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
            return int(-self.raw_bc_start-16-DEFAULT_UMI_SIZE)
        else: 
            return int(self.raw_bc_start+16+DEFAULT_UMI_SIZE)

    @property
    def putative_UMI(self):
        if not self._get_strand_and_raw_bc_flag:
            self.get_strand_and_raw_bc()
        if not self._strand:
            return None
        elif self._strand == '+':
            return helper.reverse_complement(
                        self.seq)[self.raw_bc_start+16: self.raw_bc_start+16+DEFAULT_UMI_SIZE]
        else: 
            return self.seq[
                self.raw_bc_start+16: self.raw_bc_start+16+DEFAULT_UMI_SIZE]
    
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
                        self.seq)[self.raw_bc_start+16+DEFAULT_UMI_SIZE: \
                                  self.raw_bc_start+16+DEFAULT_UMI_SIZE+ DEFAULT_GRB_FLANKING_SIZE]
        else: 
            return self.seq[
                self.raw_bc_start+16+DEFAULT_UMI_SIZE: \
                    self.raw_bc_start+16+DEFAULT_UMI_SIZE+DEFAULT_GRB_FLANKING_SIZE]
    

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
    
    