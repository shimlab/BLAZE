"""Searching for adaptor polyT or poly A in a Nanopore read.
"""
import numpy as np
import Bio
import Bio.pairwise2

import helper

class Read(object):
    '''
    class storing all information about a read
    '''
    def __init__(self, sequence, phred_score=None, strand = None, **kwargs):
        self.seq = sequence
        self.phred_score = phred_score
        self._strand = strand
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
            
            
    def add(self, attr_name, attr_val, overwrite = True):
        if attr_name in __dict__.keys() and not overwrite:
            pass
        else:    
            self.__dict__[attr_name] = attr_val
    
    
    def find_adaptor(self,read=None, strand=None, adaptor_seq='TCTTCCGATCT', num_nt=150, min_match_prop=0.8):
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
    
        if strand == '-':
            seq = read[:num_nt]
            align = Bio.pairwise2.align.localms(seq, adaptor_seq,2,-1,-1,-1)
            return {'-':[a for a in align if a.score >= min_score(len(adaptor_seq), min_match_prop)]}
        
        if strand == '+':
            read = helper.reverse_complement(read)
            return {'+':self.find_adaptor(read, strand='-', adaptor_seq=adaptor_seq, 
                                num_nt=num_nt, min_match_prop=min_match_prop)['-']}
                
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
        


    # def find_poly_T(self, poly_T_len=15, num_nt=150, min_match_prop=0.8):
    #     '''
    #     find adaptor from a read
    
    #     Parameters
    #     ----------
    #     read : STR
    #     num_nt : INT
    #         The adaptor will be searched within first num_nt bases in the reads.
    #     strand : '+' or '-'
    #         Transcript strand, this function will find adaptor in '-' poly T strand
    #         Will do a reverse complement if strand == '-'
    #     poly_T_len : INT
    #         How many Ts are used to match in the read. Default: 15.
    #     min_match_prop : float
    #         min proportion of matches in the alignment. Default: 0.8
    #     Returns
    #     -------
    #     Bio.pairwise2.Alignment
    #     *(adaptor_start, adaptor_end): adaptor start and end position (0-based) in reads
    #     '''
    
    #     Ts = 'T' * poly_T_len
    #     # reuse the same method of finding adaptors to find polyT
    #     return  self.find_adaptor(adaptor_seq=Ts, num_nt=num_nt, min_match_prop=0.8)
        
    def find_poly_T(self, read=None, strand=None, 
                       poly_T_len=10, num_nt=200, min_match_prop=0.85):
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
        poly_T_len : INT
            How many Ts are used to match in the read. Default: 15.
        min_match_prop : float
            min proportion of matches in the alignment. Default: 0.8
        Returns
        -------
        Bio.pairwise2.Alignment
        *(adaptor_start, adaptor_end): adaptor start and end position (0-based) in reads
        '''
        strand = self._strand if not strand else strand
        read = self.seq if not read else read
        if strand == '-':
            seq = read[:num_nt]
            # convert to np_array T -> 1 and ACG -> 0
            read_code = np.array([int(x == 'T') for x in seq])
            T_prop = helper.sliding_window_mean(read_code, poly_T_len)
            wind_start_T = read_code[:-poly_T_len]
            first_index = np.argmax(T_prop*wind_start_T > min_match_prop)
            return {'-':first_index}
        
        if strand == '+':
            return {'+': self.find_poly_T(
                read=helper.reverse_complement(read), strand='-', 
                poly_T_len=poly_T_len, num_nt=num_nt, 
                min_match_prop=min_match_prop)['-']}
        else:
            A_strand = {'+': self.find_poly_T(
                read=helper.reverse_complement(read), strand='-', 
                poly_T_len=poly_T_len, num_nt=num_nt, 
                min_match_prop=min_match_prop)['-']}
            T_strand = self.find_poly_T(
                strand='-', poly_T_len=poly_T_len, num_nt=num_nt, 
                min_match_prop=min_match_prop)
            return {**{k:v for k,v in T_strand.items() if v != 0},
                    **{k:v for k,v in A_strand.items() if v != 0}}
    
    
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
        poly_T_dict = self.find_poly_T()
        adapt_dict=self.find_adaptor()
        self.adaptor_polyT_pass = 0
        # get strand
        all_strand = set(poly_T_dict.keys()) & set(adapt_dict.keys())
        rst_strand = []
        rst_t_pos = []
        rst_apt_end = []
        for strand in all_strand:
            t_pos = poly_T_dict.get(strand)
            apt_end_pos = list(set([x.end for x in adapt_dict.get(strand, [])]))
            apt_end_pos = np.array(apt_end_pos)
            try:
                adp_in_t_upstream = np.array(
                    [int(15< t_pos - x < 50) for x in apt_end_pos])
            except:
                print(t_pos, apt_end_pos, self.seq, strand)
                raise ValueError('...')
            if any(adp_in_t_upstream):
                rst_strand.append(strand)
            if sum(adp_in_t_upstream == 1):
                rst_t_pos.append(t_pos)
                rst_apt_end.append(apt_end_pos[np.argmax(adp_in_t_upstream)])

            if sum(adp_in_t_upstream) > 1:
                self.adaptor_polyT_pass += 10
        if len(rst_strand) == 0:
            self.adaptor_polyT_pass += 1
        if len(rst_strand) > 1:
            self.adaptor_polyT_pass += 2
        
        if self.adaptor_polyT_pass > 0:

            self.raw_bc_start = None
            self.raw_bc_end = None
            self.strand = None
            self.raw_bc = None
            # if self.adaptor_polyT_pass == 1:
            #     print(self.seq)
            #     #print(rst_strand,rst_apt_end, rst_t_pos,adp_in_t_upstream,self.adaptor_polyT_pass)
            #     print('\n'*3)
        else:
            try:
                assert len(rst_strand) == 1
                assert len(rst_apt_end) == 1
                assert len(rst_t_pos) == 1
            except:
                print(rst_strand,apt_end_pos,rst_apt_end, rst_t_pos,adp_in_t_upstream,self.adaptor_polyT_pass)
                raise AssertionError('..')

            self.raw_bc_start = apt_end_pos[0]
            self.raw_bc_end = rst_t_pos[0]
            self.strand = rst_strand[0]
            self.raw_bc = self.seq[self.raw_bc_start: self.raw_bc_start+16]
            
        
        # get poly_T_start_raw and adaptor_end_raw
    @property
    def poly_T_start_raw(self):
        pass
    
    @property
    def adaptor_end_raw(self):
        pass      
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
        read_seq = 'TACATGTATTGCTGCTCAAAGGCCATTACGGCCTACACGACGCCTTCCGATCTGTGGGAACACCTGTCTCAGAACCAGTAATTTTTTTTTTTTTTTTTTTTTTTTTTTTATCTAAAGTTTTTCAAATTCTTTTTAAGTGACAAAACTTGTACATGTGTATCGCTCAATATTTTGTAGTCGACAGTAATCTTGCTTTCGAGAATGTAAACCAGGCAACTTAGGAAAATGCAGACAGCACGCCTCTCTTCTGGGACCATGGCTCATACTTCGAAGTGCTCGGAGCCCTTCCTCCAGACCGCCCTCCCACACCCCGCTCCAGGGCCCTGGGAGTTACAAGCCTCGCTGCAGGCTCCTGGGAACCCAACGCGGTGTCAGAGTAGCTGGGTCCCCACGAGGGACCAGGAGCTCCGGGCGGGCAGCAGCCGCGGAAGAGTCATGCGAGGCTTTCCCAGAACCCGGCAGGGCGAAGGCAGGAGTGGGGAGGCGGAACCGGGACCTCCAGAGCCCGGTCCCTGCGCCCCACAAGCTTCCTGCTTCCTGCTAGGCCGGGCAAGGCCGGGTGGCAGGGCAAGGCCCAGGAGGAAGCCCGGGGCGAGCCTCCGCACGCTTCCCGGGCGGTCGGGGCTTCCAGCGGCGTTCAGTGGAGCTGGGCACGGGCAGCGGGCCGCGGAACACCAGCTGGCGCAGGCTTTCCTGGTCAGGAACGGTCCCGGGCCTCCTGCCCGCCTCCTCCAGCTCCTCCGGTCCCCTACTTCGCCCCCGCCAGGCCCCCACGACCCTACTTCCCGCGGCCCCGGACGCTCCTCGCCCAGTGAGCCGTCCCGGAAGCTCCTGCCGCCCCTATGTACTCTGCGCCGATACCACTGCTTGGCCATTACGGCCTGTAAAGCAATACGTAATGAACGAAGTACAA'
        
        r = Read(read_seq)
        [(x.start,x.end, x.score) for x in r.find_adaptor()['-']]
        [(x.start,x.end, x.score) for x in r.find_adaptor().get('+', [])]
        r.find_poly_T(strand = '-')
        r.find_poly_T(strand = '+')
        r.find_poly_T()
        r.get_strand_and_raw_bc()
        r.raw_bc

if __name__ == '__main__':
    main()
    
    