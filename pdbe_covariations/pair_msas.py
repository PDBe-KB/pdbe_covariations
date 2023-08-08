import sys
import os
import numpy as np
import argparse
import logging


def parse_a3m(filename):
    '''parse an a3m files as a dictionary {label->sequence}'''
    seq = []
    lab = []
    for line in open(filename, "r"):
        if line[0] == '>':
            lab.append(line.split()[0][1:])
        else:
            seq.append(line.rstrip())
    return seq, lab

def uni2idx(ids):
    '''convert uniprot ids into integers according to the structure 
    of uniprot accession numbers'''
    ids2 = [i+'AAA0' if len(i)==6 else i for i in ids]
    
    arr = np.array([list(s) for s in ids2],dtype='|S1').view(np.uint8)
    
    for i in [1,5,9]:
        arr[:,i] -= ord('0')
    
    arr[arr>=ord('A')] -= ord('A')
    arr[arr>=ord('0')] -= ord('0')-26
    
    arr[:,0][arr[:,0]>ord('Q')-ord('A')] -= 3

    arr = arr.astype(np.int64)
    
    coef = np.array([23,10,26,36,36,10,26,36,36,1], dtype=np.int64)
    coef = np.tile(coef[None,:],[len(ids),1])
    
    c1 = [i for i,id_ in enumerate(ids) if id_[0] in 'OPQ' and len(id_)==6]
    c2 = [i for i,id_ in enumerate(ids) if id_[0] not in 'OPQ' and len(id_)==6]
    
    coef[c1] = np.array([3, 10,36,36,36,1,1,1,1,1])
    coef[c2] = np.array([23,10,26,36,36,1,1,1,1,1])
    
    for i in range(1,10):
        coef[:,-i-1] *= coef[:,-i]

    return np.sum(arr*coef,axis=-1)


def run_pair_alignment(input_a3m_msa1,input_a3m_msa2):
    msa1,lab1 = parse_a3m(input_a3m_msa1)
    msa2,lab2 = parse_a3m(input_a3m_msa2)
    ids1 = [l.split('|')[1] for l in lab1[1:]]
    ids2 = [l.split('|')[1] for l in lab2[1:]]
    
    hash1 = uni2idx(ids1)
    hash2 = uni2idx(ids2)

    # find pairs of uniprot ids which are separated by at most 10
    i,j = np.where(np.abs(hash1[:,None]-hash2[None,:]) < 10)

    print("pairs found: ", i.shape[0])

    # save paired alignment
    
    with open('paired.a3m','w') as f:

        # here we assume that the pair of first sequences from the two MSAs
        # is already paired (we know that these two proteins do interact)
        f.write('>query\n%s%s\n'%(msa1[0],msa2[0]))
        idx1 = []
        idx2 = []
        # save all other pairs
        for i,j in zip(idx1,idx2):
            f.write(">%s_%s\n%s%s\n"%(lab1[i],lab2[j],msa1[i+1],msa2[j+1]))

#def main():
#    parser = argparse.ArgumentParser()
#    parser.add_argument(
#        "-msa1",
#        "--input_a3m_msa1",
#        help="Input CIF file",
#        required=True,
#    )
#    parser.add_argument(
#        "-msa2",
#        "--input_a3m_msa2",
#	help="Input CIF file",
#        required=True,
#    )

#    args = parser.parse_args()
    
#    input_a3m_msa1=args.input_a3m_msa1
#    input_a3m_msa2=args.input_a3m_msa2

#    run_pair_alignment(input_a3m_msa1,input_a3m_msa2)
    
#if "__main__" in __name__:
#    main()
