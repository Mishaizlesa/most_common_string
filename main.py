import subprocess
from random import randint
from timeit import default_timer as timer
import time
import os
import sys
num= list(map(int, sys.argv[1:]))
inf=num[0]
num=num[1:]
names=["hash3","hash3_omp","rabin_karp","rabin_karp_vect","hash3_seq","rabin_karp_seq","rabin_karp_seq_vect"]
os.environ["OMP_PLACES"]="threads"
os.environ["OMP_PROC_BIND"]="spread"
i=10;
sys.stdout.write("size ")
for j in num: sys.stdout.write(names[j]+" ")
sys.stdout.write("\n")
while(i<inf):
    sys.stdout.write(str(i)+" ")
    file="genome_samples/s"+str(i)+".txt"
    for j in num:
         subprocess.call(["build/"+names[j],file,"1"])
         tm=open("tmp.txt").read()
         sys.stdout.write(tm[:-1]+" ")
    i+=i//10
    sys.stdout.write("\n")
