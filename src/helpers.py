from math import log2
from itertools import product

def getP_j(seqs,j,I,J,k,p):
    #Function that calculates the k-mer frequency distributions for position j.
    #input parameters:
    #seqs = list of lists of all input sequences
    #j = position index
    #I = total number of sequences
    #J = length of sequences
    #k = length of k-mers
    #p = pseudocount mass
    P_j = {} #key = k-mer, value = frequency
    
    if p>0:
        #adding pseudocount p to each k-mer count
        alphabet = ['A','C','G','T']
        kmers = [''.join(c) for c in product(alphabet,repeat=k)]
        for kmer in kmers: P_j[kmer] = p/pow(4,k)#p
        
        for i in range(0,I):
            kmer = seqs[i][j:j+k]
            P_j[kmer] += 1.0
            
        #normalize counts
        summ = 0.0
        for kmer in P_j: summ += P_j[kmer]        
        for kmer in P_j: P_j[kmer] /= (float(I)+p)#(float(I)+p*pow(4,k))#(float(I)+p*4**k)
          
    else: 
        for i in range(0,I):
            kmer = seqs[i][j:j+k]
            if kmer not in P_j: P_j[kmer] = 1.0
            else: P_j[kmer] += 1.0
            
        #normalize counts
        for kmer in P_j:
            P_j[kmer] /= float(I)
      
    
    
    return P_j

def getMI_mn(seqs,m,n,I,J,P_j,k,p):
    #Function that calculates mutual information between k-mers starting at position j
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
    
    #calculate pairwise frequencies for positions m and n
    P_mn = {} #key = k-mer pair, for example "AAA" and "ACG is coded as "AAAACG", value = frequency
    
    if p>0:
        #adding pseudocount p to each k-mer count
        alphabet = ['A','C','G','T']
        kmers = [''.join(c) for c in product(alphabet,repeat=2*k)]
        for kmer in kmers: P_mn[kmer] = p/pow(4,2*k)#p
        #P_m = {}
        #P_n = {}
        CAA_count = 0
        for i in range(0,I):
            kmer1 = seqs[i][m:m+k]
            kmer2 = seqs[i][n:n+k]
            if kmer2=='CAA': CAA_count += 1
            kmer = kmer1+kmer2
            P_mn[kmer] += 1.0#/(float(I)+p)**2
            
            #if kmer1 not in P_m: P_m[kmer1] = 1.0
            #else: P_m[kmer1] += 1.0
            
            #if kmer2 not in P_n: P_n[kmer2] = 1.0
            #else: P_n[kmer2] += 1.0
            
        
        #normalize counts
        summ = 0.0
        for kmer in P_mn:
            #kmer1 = kmer[:k]
            #kmer2 = kmer[k:]
            #P_mn[kmer] += (P_mn[kmer]*(p/P_m[kmer1]+p/P_n[kmer2]))
            summ += P_mn[kmer]
        
        #print("summ="+str(summ)+" / "+str((float(I)+p*pow(4,2*k))))
        for kmer in P_mn: P_mn[kmer] /= (float(I)+p)#(float(I)+p*pow(4,2*k))#(float(I)+p*4**(2*k))#summ
        
        
    else:
        for i in range(0,I):
            kmer1 = seqs[i][m:m+k]
            kmer2 = seqs[i][n:n+k]
            if kmer1+kmer2 not in P_mn: P_mn[kmer1+kmer2] = 1.0
            else: P_mn[kmer1+kmer2] += 1.0
              
        #normalize counts
        summ = 0
        for kmer in P_mn:
            P_mn[kmer] /= float(I)
    #if m==0:
    #    print("m="+str(m)+", n="+str(n))
    #    print(P_mn)
    #calculate mutual information
    #MI(m,n) = \sum_{kmer_m \in K} \sum_{kmer_n in K} P_mn(kmer_m,kmer_n)*log_2(P_mn(kmer_m,kmer_n)/(P_m(kmer_m)*P_n(kmer_n)))
    #where K is the set of all observed k-mers at positions m and n
    MI_mn = 0.0
    MI_auxs = []
    MI_GTG_sum = 0.0
    MI_others_sum = 0.0
    for kmer_m in P_j[m]:
        for kmer_n in P_j[n]:
            kmer_mn = kmer_m+kmer_n
            if kmer_mn not in P_mn: continue
            if p>0: MI_auxs.append(P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n])))
            if p>0:
                if kmer_m=='GTG': MI_GTG_sum += MI_auxs[-1]
                else: MI_others_sum += MI_auxs[-1]
            MI_mn += P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))
            if m==5 and n==80 and kmer_m=='GTG' and kmer_n=='CAA' and p>0: print("MI_5,80^(GTG,CAA)="+str(P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))))
            if m==5 and n==80 and kmer_m=='CAA' and kmer_n=='CAA' and p>0: print("MI_5,80^(CAA,CAA)="+str(P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))))
            if m==60 and n==80 and kmer_m=='GTG' and kmer_n=='CAA' and p>0: print("MI_60,80^(GTG,CAA)="+str(P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))))
    
    if m==5 and n==80 and p>0:
            print("P_5,80(GTG,CAA)="+str(P_mn["GTGCAA"]))
            print("P_5(GTG)="+str(P_j[m]["GTG"]))
            print("P_80(CAA)="+str(P_j[n]["CAA"]))
            print("observed number of pairs = "+str(summ))
            print("CAA count="+str(CAA_count))
            print("MI_5,80="+str(MI_mn))
            print("MI_GTG_sum="+str(MI_GTG_sum))
            print("MI_others_sum"+str(MI_others_sum))
            #print(MI_auxs)
    if m==60 and n==80 and p>0:
            print("P_60,80(GTG,CAA)="+str(P_mn["GTGCAA"]))
            print("P_60(GTG)="+str(P_j[m]["GTG"]))
            print("P_80(CAA)="+str(P_j[n]["CAA"]))
            print("observed number of pairs = "+str(summ))
            print("CAA count="+str(CAA_count))
            print("MI_60,80="+str(MI_mn))
            print("sum(MI_auxs)="+str(sum(MI_auxs)))
    
    #for kmer_mn in P_mn:
    #    kmer_m = kmer_mn[:k]
    #    kmer_n = kmer_mn[k:]
    #    #if m==0 and n==3: print(kmer_mn+": "+kmer_m+", "+kmer_n+"| P_mn="+str(P_mn[kmer_mn]))
    #    #if kmer_mn not in P_mn: continue
    #    #if m==0 and n==3: print(P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n])))
    #    MI_mn += P_mn[kmer_mn]*log2(P_mn[kmer_mn]/(P_j[m][kmer_m]*P_j[n][kmer_n]))
            
    return [(m,n),MI_mn]
    
    
