import subprocess
from random import randint
from timeit import default_timer as timer
import time
import os
import sys
num=[0,1,2,3,4,5,6]
inf=1000
names=["hash3","hash3_omp","rabin_karp","rabin_karp_vect","hash3_seq","rabin_karp_seq","rabin_karp_seq_vect"]
os.environ["OMP_PLACES"]="threads"
os.environ["OMP_PROC_BIND"]="spread"
i=10;
while(i<inf):
    file="genome_samples/s"+str(i)+".txt"
    for j in num:
         subprocess.call(["build/"+names[j],file,"2"])
         tm=open("tmp.txt").read()
         sample=open("ans/s"+str(i)+".txt").read()
         if (sample!=tm): print(names[j],i,tm,sample)
    i+=i//10
