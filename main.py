import subprocess
from random import randint
from timeit import default_timer as timer
import time
import os
import sys
inf=int(sys.argv[1])
num=int(sys.argv[2])#количество запускаемых алгоритмов
names=["hash3","rabin_karp","rabin_karp_vect","hash3_seq","rabin_karp_seq","rabin_karp_seq_vect"]#алгоритмы
os.environ["OMP_PLACES"]="threads"
os.environ["OMP_PROC_BIND"]="spread"
i=10;
out=open("time_result.txt",'w');
out.write("size ")
for j in range (num): out.write(names[j]+" ")
out.write("\n")
while(i<inf):
    out.write(str(i)+" ");
    file="genome_samples/s"+str(i)+".txt"
    for j in range (num):
         subprocess.run(["build/Release/"+names[j],file,"1"])
         tm=open("tmp.txt").read()
         out.write(tm[:-1]+" ")
    i+=i//10
    out.write('\n')
