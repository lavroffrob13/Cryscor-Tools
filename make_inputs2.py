import os
import sys
import time
import glob

file_geom="geom_mod"

file_crystal_input="INPUT"
file_new_crystal_input="hf1.inp"
file_new_crystal_output="hf1.out"
file_dual_crystal_input="hf2.inp"
file_dual_crystal_output="hf2.out"

file_prop_input="prop.inp"
file_prop_output="prop.out"

file_frag_input="frag.inp"
file_frag_output="frag.out"


launch_prop_scr_file="Prop.job"
launch_HF1_scr_file="HF1.job"
launch_HF2_scr_file="HF2.job"
launch_frag_scr_file="Frag.job"

modules="XE2016.0.3.210 ompi-202_xe2016.0.3.210"
crystal_path="/huge/usvyat/Crystal17/bin/Linux-ifort17_XE_emt64/std/"
n_proc="01"

prop_path_addon="/opt/intel/XE2016.U3/bin"
prop_mpirun_path="/opt/intel/XE2016.U3/compilers_and_libraries_2016.3.210/linux/mpi/bin64/"
prop_path="/huge/usvyat/Crystal14_new/bin/Linux-ifort-i64-new2/std/"

modules_cryscor="XE2018.0.2.199 ompi-213_xe2018.0.2.199"
cryscor_path="/backup/usvyat/Crystal_GIT/git_cryscor/bin/Linux-ifort64-mkl-test-part-t13/std/"

idle_int=10

def run_HF1(launch_scr_file,hf1_input,hf1_output,n_proc,modules,crystal_path):
        scr_file=open(launch_scr_file,'w')
        scr_file.write("#!/bin/bash\n")
        scr_file.write("#PBS -j eo\n")
        scr_file.write("#PBS -r n\n")
        scr_file.write("#PBS -l nodes="+n_proc+"\n")
        scr_file.write("\n")
        scr_file.write("module load "+modules+"\n")
        scr_file.write("\n")
        
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("TMPDIR=/scratch/${USER}/tmp$$\n")
        scr_file.write("mkdir -p $TMPDIR\n")
        scr_file.write("cp "+hf1_input+" $TMPDIR/INPUT\n")
        scr_file.write("cd $TMPDIR\n")
        scr_file.write("\n")
        
        scr_file.write("mpirun "+crystal_path+"Pcrystal >& ${PBS_O_WORKDIR}/"+hf1_output+"\n")
        scr_file.write("cp fort.9 ${PBS_O_WORKDIR}\n")
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("rm -r $TMPDIR\n")
        scr_file.close()

        os.system("mkdir HF1")
        os.chdir("HF1")
        print("cp ../"+hf1_input+" .")
        os.system("cp ../"+hf1_input+" .")

        print("cp ../"+launch_scr_file+" .")
        os.system("cp ../"+launch_scr_file+" .")

        if glob.glob(launch_scr_file+".e*"):
                os.system("rm "+launch_scr_file+".e*")
        
        print("qsub "+launch_scr_file)
        os.system("qsub "+launch_scr_file)
        while not glob.glob(launch_scr_file+".e*"):
                print("waiting")
                time.sleep(idle_int)
      

        converged=False
        hf1_energy=0
        geom=False
        coord_angs=[]
        with open(hf1_output) as file:
                for line in file:
#                        print(line)
#                        print(geom)
                        if not geom:
                                if 'X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)' in line:
                                        geom=True
                                if 'SCF ENDED' in line:
                                        converged=True
                                        words=line.split()
                                        hf1_energy=float(words[8])
                                        break
                        else:
                                words=line.strip()
                                if words:
                                        if '*' in line:
                                                continue
#                                        if 'ATOM' in line:
#                                                continue
                                        coord_angs.append(line)
                                else:
                                        geom=False
                                                        
        os.chdir("../")       
        os.system("exit")
        return converged,hf1_energy,coord_angs


def run_HF2(launch_scr_file,hf2_input,hf2_output,n_proc,modules,crystal_path):
        scr_file=open(launch_scr_file,'w')
        scr_file.write("#!/bin/bash\n")
        scr_file.write("#PBS -j eo\n")
        scr_file.write("#PBS -r n\n")
        scr_file.write("#PBS -l nodes="+n_proc+"\n")
        scr_file.write("\n")
        scr_file.write("module load "+modules+"\n")
        scr_file.write("\n")
        
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("TMPDIR=/scratch/${USER}/tmp$$\n")
        scr_file.write("mkdir -p $TMPDIR\n")
        scr_file.write("cp "+hf2_input+" $TMPDIR/INPUT\n")
        scr_file.write("cp ../HF1/fort.9 fort.20\n")
        scr_file.write("cp fort.20 $TMPDIR\n")

        scr_file.write("cd $TMPDIR\n")
        scr_file.write("\n")
        
        scr_file.write("mpirun "+crystal_path+"Pcrystal >& ${PBS_O_WORKDIR}/"+hf2_output+"\n")
        scr_file.write("cp fort.9 ${PBS_O_WORKDIR}\n")
        scr_file.write("cp fort.78 ${PBS_O_WORKDIR}\n")
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("rm -r $TMPDIR\n")
        scr_file.close()

        os.system("mkdir HF2")
        os.chdir("HF2")
        print("cp ../"+hf2_input+" .")
        os.system("cp ../"+hf2_input+" .")
        
        if glob.glob(launch_scr_file+".e*"):
                os.system("rm "+launch_scr_file+".e*")
 
        print("cp ../"+launch_scr_file+" .")
        os.system("cp ../"+launch_scr_file+" .")
        
        print("qsub "+launch_scr_file)
        os.system("qsub "+launch_scr_file)

        os.chdir("../")
        os.system("exit")
        return


def run_prop(launch_prop_file,prop_input,prop_output,path_addon,prop_path,mpirun_path):
        scr_file=open(launch_prop_file,'w')
        scr_file.write("#!/bin/bash\n")
        scr_file.write("#PBS -j eo\n")
        scr_file.write("#PBS -r n\n")
        scr_file.write("#PBS -l nodes=1\n")
        scr_file.write("\n")
        scr_file.write("export PATH="+path_addon+":${PATH}\n")
        scr_file.write("\n")
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("TMPDIR=/scratch/${USER}/tmp$$\n")
        scr_file.write("mkdir -p $TMPDIR\n")
        scr_file.write("cp "+prop_input+" $TMPDIR/INPUT\n")
        scr_file.write("cp ../HF1/fort.9 $TMPDIR\n")
        scr_file.write("cd $TMPDIR\n")
        scr_file.write("\n")
        
        scr_file.write(mpirun_path+"mpirun -n 1 "+prop_path+"Pproperties >& ${PBS_O_WORKDIR}/"+prop_output+"\n")
        scr_file.write("cp fort.80 ${PBS_O_WORKDIR}\n")
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("rm -r $TMPDIR\n")
        scr_file.close()

        os.system("mkdir Prop")
        os.chdir("Prop")
        print("cp ../"+prop_input+" .")
        os.system("cp ../"+prop_input+" .")
        
        if glob.glob(launch_prop_file+".e*"):
                os.system("rm "+launch_prop_file+".e*")

        print("cp ../"+launch_prop_file+" .")
        os.system("cp ../"+launch_prop_file+" .")
        
        print("qsub "+launch_prop_file)
        os.system("qsub "+launch_prop_file)

        os.chdir("../")
        os.system("exit")
        return

def run_frag(launch_scr_file,frag_input,frag_output,modules,cryscor_path):
        scr_file=open(launch_scr_file,'w')
        scr_file.write("#!/bin/bash\n")
        scr_file.write("#PBS -j eo\n")
        scr_file.write("#PBS -r n\n")
        scr_file.write("#PBS -l nodes=1\n")
        scr_file.write("\n")
        scr_file.write("module load "+modules+"\n")
        scr_file.write("\n")
        
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("TMPDIR=/scratch/${USER}/tmp$$\n")
        scr_file.write("mkdir -p $TMPDIR\n")

        scr_file.write("cp "+frag_input+" $TMPDIR/\n")
        scr_file.write("cp ../HF2/fort.9 $TMPDIR/\n")
        scr_file.write("cp ../HF2/fort.78 $TMPDIR/\n")
        scr_file.write("cp ../Prop/fort.80 $TMPDIR/\n")
        scr_file.write("cp "+cryscor_path+"cryscor $TMPDIR/\n")
        scr_file.write("cd $TMPDIR\n")
        
        scr_file.write("\n")        
        scr_file.write("export MKL_NUM_THREADS=1\n")  
      
        scr_file.write("./cryscor < "+frag_input+" > ${PBS_O_WORKDIR}/"+frag_output+" 2>"+frag_output+".err")
        scr_file.write("\n")
        
        scr_file.write("cp "+frag_output+".err ${PBS_O_WORKDIR}/\n")
        scr_file.write("cp FCIDUMP ${PBS_O_WORKDIR}/\n")
        scr_file.write("cd ${PBS_O_WORKDIR}\n")
        scr_file.write("rm -r $TMPDIR\n")
        scr_file.close()

        os.system("mkdir Frag")
        os.chdir("Frag")
        print("cp ../"+frag_input+" .")
        os.system("cp ../"+frag_input+" .")

        print("cp ../"+launch_scr_file+" .")
        os.system("cp ../"+launch_scr_file+" .")

        if glob.glob(launch_scr_file+".e*"):
                os.system("rm "+launch_scr_file+".e*")
        
        print("qsub "+launch_scr_file)
        os.system("qsub "+launch_scr_file)

#        while not glob.glob(launch_scr_file+".e*"):
#                print("waiting")
#                time.sleep(idle_int)
      


def get_prop(list_to_add,list_to_remove,n_atom_old,inp,n_k,coord_angs):
        angs2bohr=1.88972613392125187641
        inp.write("NEWK\n")
        inp.write(str(n_k)+" "+str(n_k)+" "+str(n_k)+"\n")
        inp.write("1 0\n")
        inp.write("POTC\n")
        n_pts=len(list_to_remove)+len(list_to_add)
        inp.write("0 "+str(n_pts)+" 0\n")
#        print("n_atom_old",n_atom_old)
        if list_to_remove:
                for atom_rm in list_to_remove:
                        for atoms in coord_angs:
                                words=atoms.split()
                                if int(words[0])==atom_rm:
                                        ind_st=4
                                        if words[1].isnumeric():
                                                ind_st=3       
#                                        print(ind_st)
#                                        print(words[1],words[1].isnumeric())
#                                        print(words)
                                        x=float(words[ind_st])*angs2bohr
                                        y=float(words[ind_st+1])*angs2bohr
                                        z=float(words[ind_st+2])*angs2bohr
                                        inp.write(str(x)+" "+str(y)+" "+str(z)+"\n")
        if list_to_add:
                ind=0
                for atoms in coord_angs:
                        ind+=1
                        if ind<=n_atom_old:
                                continue
                        words=atoms.split()
                        ind_st=4
                        if words[1].isnumeric():
                                ind_st=3
#                        print(ind_st)
#                        print(words[1],words[1].isnumeric())
#                        print(words)
                        x=float(words[ind_st])*angs2bohr
                        y=float(words[ind_st+1])*angs2bohr
                        z=float(words[ind_st+2])*angs2bohr
                        inp.write(str(x)+" "+str(y)+" "+str(z)+"\n")
        inp.write("LOCALWF\n")
        inp.write("FULLBOYS\n")
        inp.write("6\n")
        inp.write("WANDM\n")
        inp.write("-1 8\n")
        inp.write("END\n")
        inp.write("END\n")


def get_frag(list_to_add,list_to_remove,inp,n_k,n_init_at,epots,vacancies):
        inp.write("READC14\n")
        inp.write("DUALBAS\n")
        inp.write("KNET\n")
        inp.write(str(n_k)+"\n")
        inp.write("MEMORY\n")
        inp.write("60000\n")
        inp.write("NOSING\n")
        n_old=len(list_to_remove)
        n_new=len(list_to_add)
        inp.write("DEFECT\n")
        inp.write(str(n_old)+"\n")
        if n_old:
                for atom,vac in zip(list_to_remove,vacancies): #epots[0:n_old]):
#                        words=epot.split()
#                        inp.write(str(atom)+" "+words[4]+"\n")
                        inp.write(str(atom)+" "+str(vac)+"\n") 
        inp.write(str(n_new)+"\n")                
        if n_new:
                ind=n_init_at+1
                for atom,epot in zip(list_to_add,epots[n_old:n_old+n_new]):
                        words=epot.split()
                        inp.write(str(ind)+" "+words[4]+"\n")
                        ind+=1
        inp.write("RATOMCAS\n")                        
        inp.write("10\n")
        inp.write("MINPOP\n")
        inp.write("0.2\n")
        inp.write("DOMPUL\n")
        inp.write("0.8\n")
        inp.write("DFITTING\n")
        inp.write("DIRECT\n")
        inp.write("G-AVTZ\n")
        inp.write("ENDDF\n")
        inp.write("PRINPLOT\n")
        inp.write("2\n")
        inp.write("END\n")
        inp.write("END\n")



def analyze_geom(header,geometry,list_to_add,coords_to_add,inp1,inp2):
        geom_inp=header
        n_atoms_old=int(geometry[0])
        n_atoms_new=len(list_to_add)
        same_position=[]
        
        geom_inp.append(str(n_atoms_old+n_atoms_new))

        for lines in geometry[1:n_atoms_old+1]:
                geom_inp.append(lines)

        for atoms,coords in zip(list_to_add,coords_to_add):
                tol=0.0001
                x_new=float(coords[0])
                y_new=float(coords[1])
                z_new=float(coords[2])
                is_same=False
                for lines_old in geometry[1:n_atoms_old]:
                        words=lines_old.split()
                        x_old=float(words[1])
                        y_old=float(words[2])
                        z_old=float(words[3])
                        if abs(x_new-x_old)<tol and abs(y_new-y_old)<tol and abs(z_new-z_old)<tol:
                                is_same=True
                                break
                
                same_position.append(is_same)
                geom_inp.append(str(atoms)+" "+' '.join(coords))
                        
        ind=1
        for lines in geometry:
                if ind>1+n_atoms_old:
                       geom_inp.append(lines)
                ind+=1

        for lines in geom_inp:
                inp1.write(lines.strip()+"\n")
                inp2.write(lines.strip()+"\n")
        return same_position
                

def analyze_basis(basis,list_to_add,atom_shift,n_old_atoms,inp1,inp2,same_position):
        
        narrow_GTO1="0 0 1 0 1" 
        narrow_GTO2="10000000. 1."
#        narrow_GTO2_same="100000000. 1."
        narrow_GTO2_same="10000000. 1."
        basis_spec=True
        pseudo_spec=False
        shell_spec=False
        ending=False
        basis_list1=[]
        basis_list_new1=[]
        basis_list2=[]
        basis_list_new2=[]
        atom_inds_old=[]
        n_dual=0
        dual_at=[]
        dual_n_shells=[]
        old_at1=[]
        old_at2=[]
        new_at1=[]
        new_at2=[]
        for lines in basis:
                if basis_spec:
                        basis_spec=False
                        words=lines.split()
                        atom=int(words[0])
                        if atom==99:
                                ending=True
                                BS_end=[]
                                BS_end.append(lines)
                                continue
                        ind_shell=0
                        n_shells=int(words[1])

                        old_line=lines
                        old_line_plus=words[0]+" "+str(n_shells+1)
                        old_line_one=words[0]+" 1"
                        
                        atom_ind=atom%100
                        atom_inds_old.append(atom_ind)
                        atom_to_add=False
                        atom_fully_new=False
                        for atoms,same_pos in zip(list_to_add,same_position):
                                if atom_ind==atoms%100:
                                        atom_to_add=True
                                        n_dual+=1
                                        if same_pos:
                                                narrow_GTO2_new=narrow_GTO2_same
                                        else:
                                                narrow_GTO2_new=narrow_GTO2
                                        dual_at.append(atoms)
                                        dual_n_shells.append(str(n_shells))
                                        if atom_shift[atom_ind]==0:
                                                atom_fully_new=True
                                                break
                                        new_line_plus=str(atoms)+" "+str(n_shells+1)
                                        new_line_one=str(atoms)+" 1"
                                        new_at1.append(new_line_one)
                                        new_at2.append(new_line_plus)
                                        break
                                
                        if atom_to_add:
#                                old_at2.append(old_line_plus)
                                if atom_fully_new:
                                        old_at1.append(old_line_one)
                                        old_at2.append(old_line_plus)
                                else:
#                                        old_at1.append(old_line_plus)
                                        old_at1.append(old_line)
                                        old_at2.append(old_line)
                        else:
                                old_at1.append(old_line)
                                old_at2.append(old_line)
                                
                        if atom>200 and atom<1000:
                                pseudo_spec=True
                                shell_spec=False
                        else:
                                if atom_to_add:
#                                        old_at1.append(narrow_GTO1)
#                                        old_at2.append(narrow_GTO1)
                                        if not atom_fully_new:
#                                                old_at1.append(narrow_GTO2)
#                                                old_at2.append(narrow_GTO2)
                                                new_at1.append(narrow_GTO1)
                                                new_at1.append(narrow_GTO2_new)
                                                new_at2.append(narrow_GTO1)
                                                new_at2.append(narrow_GTO2_new)
                                        else:
                                                old_at1.append(narrow_GTO1)
                                                old_at2.append(narrow_GTO1)
                                                old_at1.append(narrow_GTO2_new)
                                                old_at2.append(narrow_GTO2_new)
                                shell_spec=True
                elif pseudo_spec:
                        words=lines.split()
                        if words[0]=="0" or words[0]=="1" or words[0]=="2":
                                pseudo_spec=False
                                shell_spec=True
                                ind_shell+=1
                                if ind_shell==n_shells:
                                        last_shell=True
                                        ind_prim=0
                                else:
                                        last_shell=False
                                if atom_to_add:
#                                        old_at1.append(narrow_GTO1)
#                                        old_at2.append(narrow_GTO1)
                                        if not atom_fully_new:
#                                                old_at1.append(narrow_GTO2)
#                                                old_at2.append(narrow_GTO2)
                                                new_at1.append(narrow_GTO1)
                                                new_at1.append(narrow_GTO2_new)
                                                new_at2.append(narrow_GTO1)
                                                new_at2.append(narrow_GTO2_new)
                                        else:
                                                old_at1.append(narrow_GTO1)
                                                old_at2.append(narrow_GTO1)
                                                old_at1.append(narrow_GTO2_new)
                                                old_at2.append(narrow_GTO2_new)
                        old_at2.append(lines)
                        if not shell_spec:
                                old_at1.append(lines)
                        if atom_to_add and not atom_fully_new:
                                new_at2.append(lines)
                                if not shell_spec:
                                        new_at1.append(lines)
                elif shell_spec:
                        words=lines.split()
                        if words[0]=="0" or words[0]=="1" or words[0]=="2":
                                ind_shell+=1
                                if ind_shell==n_shells:
                                        last_shell=True
                                        ind_prim=0
                                        n_prims=int(words[2])
                                else:
                                        last_shell=False
                        else:
                                if last_shell:
                                        ind_prim+=1
                                        if ind_prim==n_prims:
                                                shell_spec=False
                                                basis_spec=True
                        old_at2.append(lines)
                        if not atom_fully_new:
                                old_at1.append(lines)
                        if atom_to_add and not atom_fully_new:
                                new_at2.append(lines)
                else:
                        BS_end.append(lines)

        for atoms in list_to_add:
                atom_to_add_exists=False                
                for atoms_old in atom_inds_old:
                        if atoms_old==atoms%100:
                                atom_to_add_exists=True
                                break
                if not atom_to_add_exists:
                        sys.exit("Please specify the basis for all new atoms")

        new_ghosts=""
        for i in range(len(list_to_add)):
                new_ghosts=new_ghosts+" "+str(i+n_old_atoms+1)
                

        BS_ghosts_end=[]
        ghosts_in=False 
        ghost_line_ind=0
        for lines in BS_end:
                if 'GHOSTS' in lines: 
                        ghosts_in=True
                        ghost_line_ind=1
                        BS_ghosts_end.append(lines)
                elif ghosts_in and ghost_line_ind==1:
                        words=lines.split()
                        n_old_ghosts=int(words[0])
                        n_ghosts=len(list_to_add)+n_old_ghosts
                        BS_ghosts_end.append(str(n_ghosts))
                        ghost_line_ind+=1
                elif ghosts_in and ghost_line_ind==2:
                        line=lines.strip()+" "+new_ghosts
                        BS_ghosts_end.append(str(line))
                        ghost_line_ind+=1
                else:
                        BS_ghosts_end.append(lines)
                
        if not ghosts_in:
                BS_ghosts_end=[]
                for lines in BS_end:
                        if 'END' in lines:
                                BS_ghosts_end.append("GHOSTS")
                                BS_ghosts_end.append(str(len(list_to_add)))
                                BS_ghosts_end.append(new_ghosts)
                        BS_ghosts_end.append(lines)

#       print("basis")
#       print(basis_list1)
#       print("basis end")
        for lines in old_at1:
#               print(lines)
                inp1.write(lines.strip()+"\n")
        for lines in new_at1:
                inp1.write(lines.strip()+"\n")
        for lines in old_at2:
                inp2.write(lines.strip()+"\n")
        for lines in new_at2:
                inp2.write(lines.strip()+"\n")
        for lines in BS_ghosts_end:
                inp1.write(lines.strip()+"\n")
                inp2.write(lines.strip()+"\n")
 
        return n_dual,dual_at,dual_n_shells

                                
def analyze_rest(inp1,inp2,rest,n_dual,dual_at,dual_n_shells):
        rest_inp1=[]
        rest_inp2=[]
        rest_inp2.append("GUESDUAL")
        rest_inp2.append(str(n_dual)+" 0")
        for atom,n_shells in zip(dual_at,dual_n_shells):
                rest_inp2.append(str(atom)+" 1 "+str(n_shells))
        n_k=1
        for lines in rest:
                rest_inp1.append(lines)
                rest_inp2.append(lines)
                if n_k==0:
                        words=lines.split()
                        n_k=int(words[0])
                if 'SHRINK' in lines:
                        n_k=0
        for lines in rest_inp1:
                inp1.write(lines.strip()+"\n")
        for lines in rest_inp2:
                inp2.write(lines.strip()+"\n")

        return n_k
                
def is_HF2_prop_done(launch_prop_scr_file,launch_HF2_scr_file,file_prop_output,file_HF2_output,hf1_en):
        HF2_done=False
        prop_done=False
        HF2_correct=False
        prop_correct=False
        Epots=[]

        while not (HF2_done and prop_done):
                if not prop_done:
#                        print("waiting prop")
                        if glob.glob("Prop/"+launch_prop_scr_file+".e*"):
                                prop_done=True
                                ePotlines=False
                                got_epot=False
                                got_WFs=False

                                with open("Prop/"+file_prop_output) as file:
                                        for line in file:
                                                if not ePotlines:
                                                        if 'CONVERGENCE ACHIEVED' in line:
                                                                got_WFs=True
                                                        if 'TOTAL ELECTROSTATIC' in line:
                                                                ePotlines=True
                                                                n_point=0
                                                else:
                                                        words=line.strip()
                                                        if words:
                                                                if '*' in line:
                                                                        continue
                                                                if 'POINT' in line:
                                                                        continue
                                                                Epots.append(line)
                                                                n_point+=1
                                                        else:
                                                                if n_point>0:
                                                                        got_epot=True
                                                                        ePotlines=False
                                                if got_WFs and got_epot:
                                                        prop_correct=True
                                                        break
                if not HF2_done:
#                        print("waiting HF2")
                        if glob.glob("HF2/"+launch_HF2_scr_file+".e*"):
                                HF2_done=True
                                with open("HF2/"+file_HF2_output) as file:
                                        for line in file:
                                                if 'CYC   0' in line:
                                                        words=line.split()
                                                        hf2_en=float(words[3])
                                                        if abs(hf2_en-hf1_en)<0.00000001:
                                                                HF2_correct=True
                                                        break
                time.sleep(idle_int)

        return HF2_correct,prop_correct,Epots





        

header=[]
geometry=[]
basis=[]
rest=[]
n_header=5
ind=0
geometry_st=False
n_atom_old=100000
basis_sec=True
n_basis=0
atom_shift=[0]*100

with open(file_crystal_input) as file:
        for line in file:
#               print(line)
                ind+=1
#               print(ind)
                if 'CRYSTAL' in line:
                        n_header=6
                if 'MOLECULE' in line:
                        n_header=4
                if ind<n_header:
                        header.append(line)
#                       print("header"+line+" "+str(ind)+" "+str(n_header))
                elif ind==n_header:
                        words=line.split()
                        n_atom_old=int(words[0])
                        geometry.append(line)
                elif ind<n_header+n_atom_old+2:
                        geometry.append(line)
                        if 'END' not in line:
                                words=line.split()
                                atom_ind=int(words[0])
                                if atom_ind<100:
                                        atom_shift[atom_ind%100]=max(100,atom_shift[atom_ind%100])
                                elif atom_ind<200:
                                        atom_shift[atom_ind%100]=1000
                                elif atom_ind<300:
                                        atom_shift[atom_ind%100]=max(300,atom_shift[atom_ind%100])
                                elif atom_ind<400:
                                        atom_shift[atom_ind%100]=400
                                else:
                                        sys.exit("too many different basis sets for a given species")
                elif basis_sec:
                        basis.append(line)
                        n_basis+=1
                        if 'END' in line:
                                basis_sec=False
                else:
                        rest.append(line)
print(atom_shift)
file.close()
print("Header")
for stri in header:
        print(stri.strip())
print("Geom")
for stri in geometry:
    print(stri.strip())
print("Basis")
for stri in basis:
    print(stri.strip())
print("Rest")
for stri in rest:
    print(stri.strip())



list_to_remove=[]
vacancies=[]
list_to_add=[]
coord_to_add=[]
curr_mode=0
with open(file_geom) as file:
        for line in file:
                if 'add' in line:
                        curr_mode=1
                        continue
                if 'remove' in line:
                        curr_mode=2
                        continue
                if 'move' in line:
                        curr_mode=3
                        continue
                if 'substitute' in line:
                        curr_mode=4
                        continue
                if 'vacancy' in line:
                        curr_mode=5
                        continue
                if 'nothing' in line:
                        curr_mode=0
                        continue
                if curr_mode==0: 
                        continue
                words=line.split()
                if words:
                        if curr_mode==1:
                                atom=int(words[0])%100
                                atom_ad=atom_shift[atom%100]+atom
                                coords_ad=words[1:4]
                                atom_rm=0
                        elif curr_mode==2:
                                atom_ad=0
                                atom_rm=int(words[0])
                                vac=0
                        elif curr_mode==3:
                                atom_rm=int(words[0])
                                vac=0
                                curr_geom=geometry[atom_rm]
                                atom_data=curr_geom.split()
                                atom=int(atom_data[0])%100
                                atom_ad=atom_shift[atom]+atom
                                coords_ad=words[1:4]
                        elif curr_mode==4:
                                atom_rm=int(words[0])
                                vac=0
                                atom=int(words[1])%100
                                atom_ad=atom_shift[atom%100]+atom
                                if len(words)==2:
                                        curr_geom=geometry[atom_rm]
                                        atom_data=curr_geom.split()
                                        if len(words)==2:
                                                curr_geom=geometry[atom_rm]
                                                atom_data=curr_geom.split()
                                                coords_ad=atom_data[1:4]
                                        else:
                                                coords_ad=words[2:5]
                        elif curr_mode==5:
                                atom_ad=0
                                atom_rm=int(words[0])
                                vac=1
                                if len(words)>1:
                                        add_bas_vac=words[1]
                        if atom_ad:
                                list_to_add.append(atom_ad)
                                coord_to_add.append(coords_ad)
                        if atom_rm:
                                list_to_remove.append(atom_rm)
                                vacancies.append(vac)

file.close()

hf1_file=open(file_new_crystal_input,"w")
hf2_file=open(file_dual_crystal_input,"w")
#prop_file=open(file_properties_input,"w")
same_position=analyze_geom(header,geometry,list_to_add,coord_to_add,hf1_file,hf2_file)

n_dual,dual_at,dual_n_shells=analyze_basis(basis,list_to_add,atom_shift,n_atom_old,hf1_file,hf2_file,same_position)


#print(dual_at,dual_at[0],dual_at[1])
n_k=analyze_rest(hf1_file,hf2_file,rest,n_dual,dual_at,dual_n_shells)

hf1_file.close()
hf2_file.close()

hf1_converged,hf1_energy,coord_angs=run_HF1(launch_HF1_scr_file,file_new_crystal_input,file_new_crystal_output,n_proc,modules,crystal_path)

print("Initial HF calculation converged!")
print("Periodic HF energy: "+str(hf1_energy))
              
prop_file=open(file_prop_input,"w")
get_prop(list_to_add,list_to_remove,n_atom_old,prop_file,n_k,coord_angs)
prop_file.close()


run_HF2(launch_HF2_scr_file,file_dual_crystal_input,file_dual_crystal_output,n_proc,modules,crystal_path)
run_prop(launch_prop_scr_file,file_prop_input,file_prop_output,prop_path_addon,prop_path,prop_mpirun_path)
HF2_ok,prop_ok,epots=is_HF2_prop_done(launch_prop_scr_file,launch_HF2_scr_file,file_prop_output,file_dual_crystal_output,hf1_energy)
                       
print(HF2_ok,prop_ok)
frag_file=open(file_frag_input,"w")
get_frag(list_to_add,list_to_remove,frag_file,n_k,n_atom_old,epots,vacancies)
frag_file.close()
               
run_frag(launch_frag_scr_file,file_frag_input,file_frag_output,modules_cryscor,cryscor_path)
