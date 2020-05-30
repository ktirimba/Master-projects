# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import scipy.stats as sci
import numpy as np
from math import *

def trailing_zeroes(num):
  if num == 0:
    return 32
  p = 0
  while (num >> p) & 1 == 0:
    p += 1
  return p


def probabilistic_counting(values):
    np.random.seed(17)
    m=30
    hash_vals=np.matrix([[np.random.randint(0,2) for i in range(m)] for j in range(values)])
    R=0
    for value in np.arange(values):
        arr=np.array(hash_vals[value,:])[0]
        b=''.join(map(str,arr))
        h_value= int(b,2)
        R=max(R,trailing_zeroes(h_value))
    # print 'estimated value:',2**R
    return 2**R


def loglog(values, k):
  np.random.seed(17)
##from damn cool algorithms by nick johnson
  vals=[np.random.randint(0,2**32) for i in range(values)]  #values to be hashed
  num_buckets = 2 ** k
  max_zeroes = [0] * num_buckets
  for value in np.arange(values):
        a=vals[value]
        # arr = np.array(hash_vals[value, :])[0]
        # b = ''.join(map(str, arr))
        # h_value = int(b, 2)
        bucket = a & (num_buckets - 1)
        h_value= a>>k
        print bin(h_value)
        max_zeroes[bucket] = max(max_zeroes[bucket], trailing_zeroes(h_value))
  #print 'estimated value:',2 ** (float(sum(max_zeroes)) / num_buckets) * num_buckets * 0.79402
  return 2 ** (float(sum(max_zeroes)) / num_buckets) * num_buckets * 0.79402



def flajolet_martin(values):
    #np.random.seed(17)
    num_hash=2**7
    #num_hash gives the number of hash functions to be tested
    R_list=[0]*num_hash  #list of highest number of trailing zeroes
    for h in range(num_hash):
        k=30  #number of bits for each hash value
        hash_vals = np.matrix([[np.random.randint(0, 2) for i in range(k)] for j in range(values)])
        #randomy generated hash values
        R=0
        for value in np.arange(values):
            arr = np.array(hash_vals[value, :])[0] ##makes an array out of each matrix row
            b = ''.join(map(str, arr))  #turns array into string
            h_value = int(b, 2)  #turns string of binary to integer
            R_list[h] = max(R_list[h], trailing_zeroes(h_value))  #keeps track of highest R
    means=[]
    #here i split
    for j in np.arange(0,len(R_list),int(ceil(np.log2((values))))):
        part_list=R_list[j:j+int(ceil(np.log2((values))))]
        means.append(np.mean(part_list))
    #print 'estimated value',2**np.median(means)
    return 2**np.median(means)
#print ([num_unique*1.0 / flajolet_martin([np.random.randint(0, 2 ** 32) for i in range(num_unique)])])







if __name__ == '__main__':

    print 'Relative Approximation Errors for the 3 algorithms.'
    num_unique = np.random.randint(1000,2**17)
    print 'number of unique values',num_unique
    print 'probabilistic counting:'
    print '%.3f' % np.mean([np.abs(1 - (1/0.7751)*probabilistic_counting(num_unique)/num_unique) for j in range(1)])
    print 'loglog:'
    print '%.3f' % ((np.abs(1- np.mean(
      [loglog(num_unique, 10)/num_unique for j in range(1)]))))
    print 'Flajolet Martin:'
    print '%.3f'% np.abs(1-np.mean([(flajolet_martin(num_unique)/num_unique)
                          for j in range(1)]))
