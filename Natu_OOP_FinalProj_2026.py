#OOP Project 


import re
import doctest

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}


class Seq:

    def __init__(self,sequence,gene,species,kmers=[]):
        self.sequence = sequence.upper().strip()     
        self.gene=gene
        self.species=species
        self.kmers=[]

    def __str__(self):
        #self.sequence=self.sequence.upper().strip()
        return self.sequence

    def print_record(self):
        return self.species + " " + self.gene + ": " + self.sequence
    
    def make_kmers(self, k=3):
        """
        returns kmers by default k=3 if no parameter given. 
        append kmers to list 
        >>> seq1 = Seq("ATTAGAA","GENE007","BC")
        >>> print(seq1.make_kmers())
        ['ATT', 'TTA', 'TAG', 'AGA', 'GAA']
        """
        self.kmers = []
        for i in range(len(self.sequence)-k+1):
            kmer = self.sequence[i:i+k]
            self.kmers.append(kmer)            
        return self.kmers
    
    def fasta(self):
        return(">" + self.species + " " + self.gene + "\n" + self.sequence)

class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence = re.sub('[^ATGCU]','N',self.sequence) 
        self.geneid=geneid
 
    #additional operator overload to get length of sequence
    def __len__(self):
        return len(self.sequence)
    
    def analysis(self):
        """
        Returns G and C count
        >>> DNA_1 = DNA("aGGCGggggfsdsgTAATT","gene13","baco","1307")
        >>> gc_count = DNA_1.analysis()
        >>> print(gc_count)
        9
        """
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def print_info(self):
        return " " + self.geneid + " " + self.species + " " + self.gene + ": " + self.sequence        

    def reverse_complement(self):
        """
        Returns the reverse complement of the entered sequence in string format. 
        Converts A to T, T to A, G to C, C to G, U to A. 
        If a nucleotide is missing or unknown (denoted by N), keeps it as N.
        >>> s2=DNA("ATGCN", "geneX", "Human", "2002")
        >>> s2.reverse_complement()
        'NGCAT'
        
        """
        reverse= self.sequence[::-1]
        replacement_dict = {"A":"T", "T":"A", "G":"C", "C": "G", "U":"A", "N":"N"}
        reverse_complement = ""
        for base in reverse:
            reverse_complement += replacement_dict[base]
        return reverse_complement

    def six_frames(self):
        frames=[]
        for i in range(3):
            frame = []
            for j in range(i, len(self.sequence)-2,3):
                frame.append(self.sequence[i:i+3])
            frames.append(frame)
        rc=self.reverse_complement()
        for i in range(3):
            frame = []
            for j in range(i, len(rc)-2,3):
                frame.append(rc[i:i+3])
                frames.append(frame)

        return frames
    #additional method to calculate gc percentage
    def gc_content(self):
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return "GC content of the sequence is " + str(gc/len(self.sequence)*100) + "%"
    #additional method to find repeats
    def repeat_finder(self, repeat_length):
        """
        Finds repeated nucleotides of the specified lenth or longer and returns a list.
        >>> s1=DNA("AAACCCTT", "geneA", "Mus musculus", "1001")
        >>> s1.repeat_finder(3)
        ['AAA', 'CCC']

        """
        repeated_sequences=[]
        current_repeat=self.sequence[0]
        
        for i in range(1,len(self.sequence)):
            if self.sequence[i] == self.sequence[i-1]:
                current_repeat += self.sequence[i]
                                
            else:
                if len(current_repeat) >= repeat_length:
                    repeated_sequences.append(current_repeat)
                current_repeat = self.sequence[i]
                            
        if len(current_repeat) >= repeat_length:
            repeated_sequences.append(current_repeat)
    
        return repeated_sequences
    
#s=DNA("AAAAACCC", "BRCA1", "Human", "12345")
#print(s.print_record())
#print(s.repeat_finder(3))
#d = DNA("AAATTTGCCCCA", "BRCA1", "Homo sapiens", "gene001")
#print(d.repeat_finder(3))  
#print(s.gc_content())
#print(len(s))
#print(s.six_frames())
#print(s.reverse_complement())
#print(s.make_kmers())

class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species, geneid)
        self.codons=[]
        self.rna_seq=self.sequence.replace("T", "U")
        self.make_codons()
        
    #def __str__(self):
     #   self.rna_seq=self.sequence.replace("T", "U")
      #  return self.rna_seq
    def __len__(self):
        """
        Length object overloader.
        Returns the length of the sequence 
        >>> TOY_RNA = RNA("AFTFGTT","gene1","XILO", "R12")
        >>> print(len(TOY_RNA))
        7
        """
        return len(self.sequence)

    def make_codons(self):

        for i in range(0,len(self.rna_seq),3):
            codon = self.rna_seq[i:(i+3)]
            if len(codon) == 3:
                self.codons.append(codon)
        return self.codons
       
    def translate(self):       
        prot_seq=""
        prot=[]
        for codon in self.codons:
            if "N" in codon:
                prot.append("X")
            elif standard_code[codon] == "*":
                break
            else:
                prot.append(standard_code[codon])            
        prot_seq = ''.join(prot)
        return prot_seq

#s1=RNA("ACTGTAUAGAGG", "BRCA1", "Homo sapiens", "1234")        
#s1.make_codons()
#print(s1.translate())
#print(s1.prot_seq)


class Protein(Seq):

    def __init__(self, sequence, gene, species,accession,avgscore=[],**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence = re.sub('[^a-zA-Z]','X',sequence).upper() 
        self.avgscore = avgscore


#s1=Protein("ATCGaN2", "BRCA1", "Homo sapiens","1234")    
#print(s1.sequence)            
    def total_hydro(self):
        """
        Calculates total hydrophobicity of the protein sequence.
        >>> s = Protein("ACWR", "gene1", "Homo sapiens", "1234")
        >>> s.total_hydro()
        -1.1
        
        """
        hydro=[]
        for amino_acid in self.sequence:
            hydro.append(kyte_doolittle[amino_acid])
        return sum(hydro)
    
    
#s2=Protein("AMNPQ", "BRCA1", "human", 1234)
#print(s2.total_hydro())
    
    def mol_weight(self):
        m_weight=[]
        for a_acid in self.sequence:
            m_weight.append(aa_mol_weights[a_acid])
        return sum(m_weight)
    
    def avg_hydrophobicity(self,r_frame):
        """
        Returns average hydrophobicity of entered sequence and numerical r_frame.
        Total frames is calculate by the difference of the length of sequence and reading frame(r_frame) and add 1.
        While i is less than total frames, an amino acid length of r_frame will be created.
        While iterate through each amino acid and use Kyte_doolittle dictionary to match scores.
        Append scores for each amino acid to score list. 
        The sum of the cores within the list, 'score', will then be divided by reading frame.
        This average will be added to list self.avgscore 
        >>> testp=Protein('VIKING','test','unknown',999)
        >>> x = testp.avg_hydrophobicity(3)
        >>> print(x)
        [1.6000000000000003, 1.7, -0.9666666666666667, 0.19999999999999998]
        """
        i=0
        tot_frames=((len(self.sequence) - r_frame)+1) #create the total frames that should have the same r_frame length
        while i < tot_frames:
            score = []
            aa_rf = self.sequence[i:r_frame+i]#amino acid reading frame
            for aa in aa_rf:
                score.append(kyte_doolittle[aa])#append to the score list
            self.avgscore.append((sum(score))/r_frame)
            i+=1
        return self.avgscore

#s2=Protein("AMNPQ", "BRCA1", "human", 1234)
#print(s2.mol_weight())   

#x=DNA("G","tmp","m",000)


doctest.testmod(verbose=True)

