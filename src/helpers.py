from math import log2, log, sqrt
from itertools import product

def getP_j(seqs,j,I,J,k,p,alphabet):
    #Function that calculates the k-mer frequency distributions for position j.
    #input parameters:
    #seqs = list of lists of all input sequences
    #j = position index
    #I = total number of sequences
    #J = length of sequences
    #k = length of k-mers
    #p = pseudocount mass
    #alphabet = alphabet used
    n_a = len(alphabet)
    P_j = {} #key = k-mer, value = frequency
    
    if p>0:
        #adding pseudocount mass p equally to each k-mer count
        kmers = [''.join(c) for c in product(alphabet,repeat=k)]
        for kmer in kmers: P_j[kmer] = p/pow(n_a,k)
        
        for i in range(0,I):
            kmer = seqs[i][j:j+k]
            P_j[kmer] += 1.0
            
        #normalize counts
        for kmer in P_j: P_j[kmer] /= (float(I)+p)
          
    else: 
        for i in range(0,I):
            kmer = seqs[i][j:j+k]
            if kmer not in P_j: P_j[kmer] = 1.0
            else: P_j[kmer] += 1.0
            
        #normalize counts
        for kmer in P_j: P_j[kmer] /= float(I)
       
    return P_j

def getHE_mn(seqs,m,n,I,J,P_j,k,p,alphabet):
    #Function that calculates the Hellinger distance between k-mer distributions starting at position j
    #and all other positions after j
    #input parameters:
    #seqs = list of lists of all input sequences
    #m = position index 1 for the pair
    #n = position index 2 for the pair
    #I = total number of sequences
    #J = length of sequences
    #P = singe site k-mer frequencies
    #k = k-mer length
    #p = pseudocount mass
    #alphabet = alphabet used
    n_a = len(alphabet)
    #possible pseudocount has been applied to position-specific k-mer distributions earlier
    #calculate Hellinger distance HE(P_i,P_j) = \sum_{a \in K} \sqrt{1-P_i(a)\times P_j(a)}
    #where K is the set of all possible k-mers
    HE = 0.0
    HE_mn = {} #dictionary containing all individual contributions to HE
    for kmer_m in P_j[m]:
        if kmer_m not in P_j[n]: continue
        HE_mn[kmer_m] = sqrt(P_j[m][kmer_m]*P_j[n][kmer_m])
        HE += HE_mn[kmer_m]

    HE = sqrt(1.0-HE)    
    return [(m,n),HE,HE_mn,None]

def getBC_mn(seqs,m,n,I,J,P_j,k,p,alphabet,inv=False):
    #Function that calculates the Battacharyya distance between k-mer distributions starting at position j
    #and all other positions after j
    #input parameters:
    #seqs = list of lists of all input sequences
    #m = position index 1 for the pair
    #n = position index 2 for the pair
    #I = total number of sequences
    #J = length of sequences
    #P = singe site k-mer frequencies
    #k = k-mer length
    #p = pseudocount mass
    #alphabet = alphabet used
    n_a = len(alphabet)
    #possible pseudocount has been applied to position-specific k-mer distributions earlier
    #calculate Bhattarachayya distance D_B(P_i,P_j) = -ln(\sum_{a \in K} \sqrt{P_i(a)\times P_j(a)})
    #where K is the set of all possible k-mers
    BC = 0.0
    BC_mn = {} #dictionary containing all individual contributions to BC
    for kmer_m in P_j[m]:
        if kmer_m not in P_j[n]: continue
        BC_mn[kmer_m] = sqrt(P_j[m][kmer_m]*P_j[n][kmer_m])
        BC += BC_mn[kmer_m]

    if not inv: BC = -1*log(BC)
    else: BC = -1*log(1.0-BC)
    
    return [(m,n),BC,BC_mn,None]

def getMI_mn(seqs,m,n,I,J,P_j,k,p,alphabet):
    #Function that calculates mutual information between k-mer distributions starting at position j
    #and all other positions after j
    #input parameters:
    #seqs = list of lists of all input sequences
    #m = position index 1 for the MI pair
    #n = position index 2 for the MI pair
    #I = total number of sequences
    #J = length of sequences
    #P = singe site k-mer frequencies
    #k = k-mer length
    #p = pseudocount mass
    #alphabet = alphabet used
    n_a = len(alphabet)
    #calculate pairwise frequencies for positions m and n
    P_mn = {} #key = k-mer pair, for example "AAA" and "ACG is coded as "AAAACG", value = frequency
    
    if p>0:
        #adding pseudocount p/n_a^2k to each k-mer count
        kmers = [''.join(c) for c in product(alphabet,repeat=2*k)]
        for kmer in kmers: P_mn[kmer] = p/pow(n_a,2*k)
        #P_m = {}
        #P_n = {}
        CAA_count = 0
        for i in range(0,I):
            kmer1 = seqs[i][m:m+k]
            kmer2 = seqs[i][n:n+k]
            if kmer2=='CAA': CAA_count += 1
            kmer = kmer1+kmer2
            P_mn[kmer] += 1.0
        
        #normalize counts
        for kmer in P_mn: P_mn[kmer] /= (float(I)+p)
        
        
    else:
        for i in range(0,I):
            kmer1 = seqs[i][m:m+k]
            kmer2 = seqs[i][n:n+k]
            if kmer1+kmer2 not in P_mn: P_mn[kmer1+kmer2] = 1.0
            else: P_mn[kmer1+kmer2] += 1.0
              
        #normalize counts
        for kmer in P_mn:
            P_mn[kmer] /= float(I)
    
    #calculate mutual information
    #MI(m,n) = \sum_{kmer_m \in K} \sum_{kmer_n in K} P_mn(kmer_m,kmer_n)*log_2(P_mn(kmer_m,kmer_n)/(P_m(kmer_m)*P_n(kmer_n)))
    #where K is the set of all observed k-mers at positions m and n
    MI = 0.0
    MI_mn = {} #dictionary containing all individual contributions to MI
    for kmer_m in P_j[m]:
        for kmer_n in P_j[n]:
            kmer_mn = kmer_m+kmer_n
            if kmer_mn not in P_mn: continue
            MI_mn[kmer_mn] = P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))
            MI += MI_mn[kmer_mn]#P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))
                   
    return [(m,n),MI,MI_mn,P_mn]
    
    
