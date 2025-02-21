#!/usr/bin/python3 
 
import math 
import sys
import numpy as np

input=sys.argv[1]
ind_z=int(sys.argv[2])

output=input+'.r1'
with open(input,'r+') as inp, open(output, 'w') as outp:
    line1=(inp.readline()).strip()
    outp.write(line1+', converted from vasp\n')
    outp.write('SLAB\n')
    outp.write('1\n')

    factor=float(inp.readline())
    line1=(inp.readline()).split()

    a1=np.zeros((3,3))
    a1[0,:]=np.array(line1,dtype=float)
    line1=(inp.readline()).split()
    a1[1,:]=np.array(line1,dtype=float)
    line1=(inp.readline()).split()
    a1[2,:]=np.array(line1,dtype=float)
    a1*=factor

#    print(ind_z)
    if (ind_z==1):
       i1=1
       i2=2
       print(1,ind_z)
    elif (ind_z==2):
       i1=2
       i2=0
       print(2,ind_z)
    else:
       i1=0
       i2=1
       print(3,ind_z)
#    print(i1)
#    print(i2)
     
    v1=a1[i1,:]
    v2=a1[i2,:]
#    print(v1)
#    print(v2)
    v3=np.zeros((3))
    v3[0]=v1[1]*v2[2]-v1[2]*v2[1]      # vector product v3 = v1 x v2
    v3[1]=v1[2]*v2[0]-v1[0]*v2[2]
    v3[2]=v1[0]*v2[1]-v1[1]*v2[0]
    norm=math.sqrt(v3[0]**2+v3[1]**2+v3[2]**2)
    v3*=1/norm
 

    a2=v1
    a2=np.vstack([a2, v2])
    a2=np.vstack([a2, v3])
  
    a2_inv=np.linalg.inv(a2)

    transf=np.matmul(a1,a2_inv)
   
    b1=math.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
    b2=math.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)

#    print(transf)
    alpha=math.acos((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/(b1*b2))/math.pi*180
    outp.write(str(b1)+'   '+str(b2)+'   '+str(alpha)+'\n')
  
#    line1=inp.readline()
    line1=((inp.readline()).strip()).split()
    n_at_type=len(line1)
    line2=((inp.readline()).strip()).split()
    at_inds=[]
    at_ns=[]
    tot_n=0
    for ats,ns in zip(line1,line2):
        if (ats == 'In'):
            ind='249'
        elif (ats == 'O'):
            ind='8'
        elif (ats == 'H'):
            ind='1'
        elif (ats == 'C'):
            ind='6'
        elif (ats == 'Ca'):
            ind='30'
        else:
            print('Add atom '+ats+' in the script')
            quit()
        at_inds.append(ind)
        at_ns.append(ns)
        tot_n+=int(ns)

    outp.write(str(tot_n)+'\n')

    while(True):
        line1=inp.readline().strip()
        if line1 == 'Direct':
            break
    

    curr_n=0
    curr_ind=-1
    curr_change=0
    

    while(True):
        curr_n+=1
        if curr_n>tot_n:
            break
        if curr_n>curr_change:
            curr_ind+=1
            curr_change+=int(at_ns[curr_ind])
            at=at_inds[curr_ind]
        line1=((inp.readline()).strip()).split()
        c1=np.array(line1[0:3],dtype=float)
#        print(c1)
        c1_new=np.matmul(c1,transf)
        outp.write(str(at)+'   '+str(c1_new[0])+'   '+str(c1_new[1])+'   '+str(c1_new[2])+'\n')

    outp.write('END\n')
    

            


    

