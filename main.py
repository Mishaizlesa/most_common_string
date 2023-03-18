import subprocess
from random import randint
from timeit import default_timer as timer
import time
import os
import sys
inf=int(sys.argv[1])
num=int(sys.argv[2])#количество запускаемых алгоритмов
names=["hash3","rabin_karp","rabin_karp_vect","hash3_seq","rabin_karp_seq"]#алгоритмы
os.environ["OMP_PLACES"]="threads"
os.environ["OMP_PROC_BIND"]="spread"
i=10;
s=open("genome_samples/all_data.txt",'r').read()
for ch in s:
        if (ch!='A' and ch!='C' and ch!='G' and ch!='T'):
                s=s.replace(ch,"A")
                print(ch)
f=open("genome_samples/all_data.txt",'w')
f.write(s)
f.close()
out=open("time_result.txt",'w');
out.write("size ")
for j in range (num): out.write(names[j]+" ")
out.write("\n")
while(i<inf):
    #out.write(str(i)+" ")
    #print(i)
    file="genome_samples/s"+str(i)+".txt"
    f=open(file,'w');
    pos=randint(0,inf-i);
    tmp=s[pos:pos+i]
    f.write(tmp)
    f.close()
    i+=i//10
i=10
while(i<inf):
    out.write(str(i)+" ");
    file="genome_samples/s"+str(i)+".txt"
    for j in range (num):
        #print(name)
        #print(file)
         subprocess.run(["build/Release/"+names[j],file,"1"])
         tm=open("tmp.txt").read()
         out.write(tm[:-1]+" ")
        #print(time)
    i+=i//10
    out.write('\n')
print(i)
