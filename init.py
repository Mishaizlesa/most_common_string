import subprocess
from random import randint
from timeit import default_timer as timer
import time
import os
import sys
inf=2000000
i=10;
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
    f=open(file,'w');
    pos=randint(0,inf-i);
    tmp=s[pos:pos+i]
    f.write(tmp)
    f.close()
    i+=i//10
