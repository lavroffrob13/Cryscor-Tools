#!/usr/bin/python3
import numpy as np
import math
import sys


def main(fcidump):
  with open(fcidump,'r+') as inf:

#fo = open("test1", "r+")
#first=fo.readline()
#norb=int(first.split()[2].split(",")[0])
#print(norb)



    for line in inf:
      if 'NORB=' in line:
        norb = int(line.split('NORB=')[1].split(',')[0])
      if 'NELEC=' in line:
        nelec = int(line.split('NELEC=')[1].split(',')[0])
        nocc = int(round(.5*nelec))
      if 'END' in line: break


    inf.seek(0,2)
    off=inf.tell()
    inf.seek(off-44)
    l=inf.readline()
    inf.seek(off-44*(norb+1))

    while(True):
      line1=inf.readline()
      nums=line1.split()[2:5]
      if (nums[0]=='0' and nums[1]=='0' and nums[2]=='0'):
        break

    off=inf.tell()
    inf.seek(off-44)
#off2=fo.tell()
#l2=fo.readline()
    inf.truncate()
    inf.write(l)


#    integrals = {}

#    for line in inf:
#      split = line.split()
#      val = np.float(split[0])
#      idx = tuple(int(i) for i in split[1:])
#      integrals[idx] = val
#    inf.close()
  write('posthf.inp', fcidump,nocc, norb, nelec, 'guess')
#  enuc = integrals[(0,0,0,0)]
#  print('NORB,NELEC,NOCC:',norb,nelec,nocc)
#  print('Enuc:',enuc)
  print('Input file for post-HF calculation written')

def write(molprf, fcidump, nocc, norb, nelec, orbfile):
    with open(molprf, 'w') as outf:
        outf.write('memory,500,m\naoint,c_final=0\nnosym\nnoextra\ngeometry={\n')
        nbasis = int(math.floor(norb/nocc))
        rest = int((norb - nbasis*nocc))
        for i in range(nocc):
            if i==0:
                outf.write('h0,,'+str(i)+',0,0\n')
            else:
                outf.write('h0,,'+str(i)+'0,0,0\n')
        if (rest!=0):
            outf.write('h1,,'+str(i+1)+'0,0,0\n')
        outf.write('}\n')
         
        outf.write('basis={\ns, h0,')
        for j in range(nbasis):
            outf.write(str(j+1))
            if (j+1!=nbasis):
                outf.write(', ')

        if (rest!=0):
            outf.write('\ns, h1, ')
            for l in range(rest):
                outf.write(str(l+1))
                if (l+1!=rest):
                    outf.write(', ')
        outf.write('\n}\n')
        outf.write('set,nelec='+str(nelec)+'\nset,spin=0\n')
        outf.write('GTHRESH,THROVL=-1\nint\n')
        outf.write('{HAMILTONIAN, FCIDUMP}\n')
        outf.write('{matrop\nread,ORB,ORB,CANONICAL,FILE='+str(orbfile)+'\nsave,ORB,2100.2,ORBITALS}\n')
        outf.write('{hf;noenest;\nstart,2100.2\n}\n')
        outf.write('ccsd(t)')
        outf.close()
    size = norb*norb
    guess = np.eye(norb)
    if (size%6 ==0):
        guess.resize(size//6, 6)
    else:
        guess.resize(size//6+1, 6)
    np.savetxt(str(orbfile),guess,fmt='%i')

if __name__=='__main__':
  for fcidump in sys.argv[1:]:
    main(fcidump)
