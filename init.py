import subprocess
from random import randint
import os
import sys
inf=2000000
i=10;
os.environ["OMP_PLACES"]="threads"
os.environ["OMP_PROC_BIND"]="spread"
s=open("genome_samples/all_data.txt",'r').read()
for ch in s:
        if (ch!='A' and ch!='C' and ch!='G' and ch!='T'):
                s=s.replace(ch,"A")
                print(ch)
f=open("genome_samples/all_data.txt",'w')
f.write(s)
f.close()
while(i<inf):
    file="genome_samples/s"+str(i)+".txt"
    file_out="ans/s"+str(i)+".txt"
    f=open(file,'w');
    pos=randint(0,inf-i);
    tmp=s[pos:pos+i]
    f.write(tmp)
    f.close()
    if (i<1000): subprocess.call(["build/naive",file,"2",file_out])
    i+=i//10
