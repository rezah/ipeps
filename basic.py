import pyUni10 as uni10
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import math
import copy
import time
import Move
import MoveCorboz
import MoveFull
import basicA
import basicB
import basicC
import numpy as np
from numpy import linalg as LA


def Short_TrotterSteps(N_iterF):
 List_delN=[]

# Delta_N=(1.0, N_iterF)
# List_delN.append(Delta_N)

# Delta_N=(0.08, N_iterF)
# List_delN.append(Delta_N)


# Delta_N=(0.07, N_iterF)
# List_delN.append(Delta_N)

# Delta_N=(0.06, N_iterF)
# List_delN.append(Delta_N)

# Delta_N=(0.05, N_iterF)
# List_delN.append(Delta_N)

 Delta_N=(0.04, N_iterF)
 List_delN.append(Delta_N)

 Delta_N=(0.03, N_iterF)
 List_delN.append(Delta_N)

 Delta_N=(0.02, N_iterF)
 List_delN.append(Delta_N)

 Delta_N=(0.01, N_iterF)
 List_delN.append(Delta_N)


 Delta_N=(0.009, N_iterF)
 List_delN.append(Delta_N)


 Delta_N=(0.008, N_iterF)
 List_delN.append(Delta_N)


 Delta_N=(0.007, N_iterF)
 List_delN.append(Delta_N)



 Delta_N=(0.006, N_iterF)
 List_delN.append(Delta_N)


 Delta_N=(0.005, N_iterF)
 List_delN.append(Delta_N)



 Delta_N=(0.004, N_iterF)
 List_delN.append(Delta_N)


 Delta_N=(0.003, N_iterF)
 List_delN.append(Delta_N)


 Delta_N=(0.002, N_iterF)
 List_delN.append(Delta_N)



# for i in xrange(5, 1, -1):
#  Delta_N=(i*(1.0/10),N_iterF)
#  List_delN.append(Delta_N)

# for i in xrange(10, 1, -1):
#  Delta_N=(i*(1.0/100),N_iterF)
#  List_delN.append(Delta_N)


# for i in xrange(5, 5, -1):
#  Delta_N=(i*(1.0/100),N_iterF)
#  List_delN.append(Delta_N)

# for i in xrange(10, 1, -1):
#  Delta_N=(i*(1.0/1000),N_iterF)
#  List_delN.append(Delta_N)

# for i in xrange(10, 0, -1):
#  Delta_N=(i*(1.0/10000),N_iterF)
#  List_delN.append(Delta_N)

 return List_delN


def full_make_bond(Model, D, chi, d_phys):
 ######################### No-symmetry #############################################
 if Model is "Heisenberg":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
  q_list=[q0_even]
  qchi_list=[q0_even]
  q_phys=[q0_even]*d_phys[0]

 ###################### Z(2) ######################################
 if Model is "Heisenberg_Z2":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
  q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
  q_list=[q0_even,q0_odd]
  qchi_list=[q0_even,q0_odd]
  #q_phys=[q0_even,q0_even,q0_even,q0_even,q0_odd,q0_odd,q0_odd,q0_odd]
  q_phys=[q0_even,q0_odd]

 ##########################  U(1)  ################################
 if Model is "Heisenberg_U1":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
  q1_even = uni10.Qnum(1,uni10.PRT_EVEN);
  q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
  q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
  q4_even = uni10.Qnum(4,uni10.PRT_EVEN);
  q5_even = uni10.Qnum(5,uni10.PRT_EVEN);

  q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
  q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
  q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
  q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);
  q_5_even = uni10.Qnum(-5,uni10.PRT_EVEN);

  #qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
  #qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
  qchi_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even]
  #qchi_list=[q_1_even,q0_even,q1_even]
  #qchi_list=[q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even]
  #qchi_list=[q_5_even,q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even,q5_even]
  #qchi_list=[q_1_even,q1_even]

  q_list=[q_1_even,q0_even,q1_even]
  #q_list=[q0_even,q1_even,q2_even]
  #q_list=[q_1_even,q0_even,q1_even,q2_even]
  #q_list=[q_1_even,q0_even]
  #q_list=[q_2_even,q_1_even,q0_even,q1_even, q2_even]
  #q_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even, q2_even,q3_even]
  #q_list=[q_1_even,q1_even]

  q_phys=[q_1_even,q1_even]
  #q_phys=[q_1_even,q0_even,q1_even]

 ##########################  Z2*U(1)  ################################
 if Model is "Heisenberg_U1Z2":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
  q1_even = uni10.Qnum(1,uni10.PRT_EVEN);
  q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
  q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
  q4_even = uni10.Qnum(4,uni10.PRT_EVEN);

  q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
  q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
  q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
  q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);


  q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
  q1_odd = uni10.Qnum(1,uni10.PRT_ODD);
  q2_odd = uni10.Qnum(2,uni10.PRT_ODD);
  q3_odd = uni10.Qnum(3,uni10.PRT_ODD);
  q4_odd = uni10.Qnum(4,uni10.PRT_ODD);

  q_1_odd = uni10.Qnum(-1,uni10.PRT_ODD);
  q_2_odd = uni10.Qnum(-2,uni10.PRT_ODD);
  q_3_odd = uni10.Qnum(-3,uni10.PRT_ODD);
  q_4_odd = uni10.Qnum(-4,uni10.PRT_ODD);

  qchi_list=[q_1_even,q_1_odd,q0_even,q0_odd,q1_even,q1_odd]
  q_list=[q_1_even,q_1_odd,q0_even,q0_odd,q1_even,q1_odd]
  q_phys=[q_1_odd,q1_even]
 #############################################################################

 q_D, q_chi=make_bond(D, q_list, chi, qchi_list)
 return q_D, q_chi, q_phys


def make_bond(D, q_list, chi, qchi_list):
 q_D=[]
 q_chi=[]
 for i in xrange(len(D)):
  for q in xrange(D[i]):
   q_D.append(q_list[i])

 for i in xrange(len(chi)):
  for q in xrange(chi[i]):
   q_chi.append(qchi_list[i])
 return  q_D, q_chi




def produce_GammaLanda(q_D, q_phys):
 bdi = uni10.Bond(uni10.BD_IN, q_D)
 bdo = uni10.Bond(uni10.BD_OUT, q_D)
 bdi_pys = uni10.Bond(uni10.BD_IN, q_phys)

 Gamma_a=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_a")
 Gamma_b=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_b")
 Gamma_c=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_c")
 Gamma_d=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_d")

 Landa_1=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa_2=uni10.UniTensor([bdi,bdo],"Landa_2")
 Landa_3=uni10.UniTensor([bdi,bdo],"Landa_3")
 Landa_4=uni10.UniTensor([bdi,bdo],"Landa_4")
 Landa_5=uni10.UniTensor([bdi,bdo],"Landa_5")
 Landa_6=uni10.UniTensor([bdi,bdo],"Landa_6")
 Landa_7=uni10.UniTensor([bdi,bdo],"Landa_7")
 Landa_8=uni10.UniTensor([bdi,bdo],"Landa_8")
 
 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8


def produce_GammaLanda_manual(q_D, q_phys):
 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q1_even = uni10.Qnum(1,uni10.PRT_EVEN);
 q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
 q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
 q4_even = uni10.Qnum(4,uni10.PRT_EVEN);
 q5_even = uni10.Qnum(5,uni10.PRT_EVEN);
 q6_even = uni10.Qnum(6,uni10.PRT_EVEN);
 q6_even = uni10.Qnum(7,uni10.PRT_EVEN);

 q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
 q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
 q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
 q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);
 q_5_even = uni10.Qnum(-5,uni10.PRT_EVEN);
 q_6_even = uni10.Qnum(-6,uni10.PRT_EVEN);
 q_7_even = uni10.Qnum(-7,uni10.PRT_EVEN);


# q_D=[q_2_even, q_2_even, q0_even, q0_even, q2_even, q2_even, q4_even]
# q_DD=[ q_2_even, q0_even, q0_even, q2_even, q2_even, q4_even,q4_even]

# q_D1=[q_3_even, q_1_even, q_1_even, q1_even, q1_even, q3_even, q3_even]
# q_DD1=[q_3_even,q_3_even, q_1_even, q_1_even, q1_even, q1_even, q3_even]


 q_D=[q_2_even, q_2_even, q0_even, q0_even, q2_even, q2_even, q4_even,q4_even]
 q_DD=[ q_2_even,q_2_even, q0_even, q0_even, q2_even, q2_even, q4_even,q4_even]

 q_D1=[q_3_even,q_3_even, q_1_even, q_1_even, q1_even, q1_even, q3_even, q3_even]
 q_DD1=[q_3_even,q_3_even, q_1_even, q_1_even, q1_even, q1_even, q3_even,q3_even]



 q_phys=[q_1_even, q1_even]

 bdi = uni10.Bond(uni10.BD_IN, q_D)
 bdo = uni10.Bond(uni10.BD_OUT, q_D)
 bdii = uni10.Bond(uni10.BD_IN, q_DD)
 bdoo = uni10.Bond(uni10.BD_OUT, q_DD)

 bdi1 = uni10.Bond(uni10.BD_IN, q_D1)
 bdo1 = uni10.Bond(uni10.BD_OUT, q_D1)
 bdii1 = uni10.Bond(uni10.BD_IN, q_DD1)
 bdoo1 = uni10.Bond(uni10.BD_OUT, q_DD1)

 bdi_pys = uni10.Bond(uni10.BD_IN, q_phys)

 Gamma_a=uni10.UniTensor([bdi_pys,bdi,bdi,bdo1,bdoo], "Gamma_a")
 Gamma_b=uni10.UniTensor([bdi_pys,bdi1,bdii,bdo,bdo], "Gamma_b")
 Gamma_c=uni10.UniTensor([bdi_pys,bdii,bdii,bdoo1,bdo], "Gamma_c")
 Gamma_d=uni10.UniTensor([bdi_pys,bdii1,bdi,bdoo,bdoo], "Gamma_d")

 Landa_1=uni10.UniTensor([bdi1,bdo1],"Landa_1")
 Landa_2=uni10.UniTensor([bdi,bdo],"Landa_2")
 Landa_3=uni10.UniTensor([bdi,bdo],"Landa_3")
 Landa_4=uni10.UniTensor([bdii,bdoo],"Landa_4")
 Landa_5=uni10.UniTensor([bdii,bdoo],"Landa_5")
 Landa_6=uni10.UniTensor([bdii1,bdoo1],"Landa_6")
 Landa_7=uni10.UniTensor([bdii,bdoo],"Landa_7")
 Landa_8=uni10.UniTensor([bdi,bdo],"Landa_8")
 
 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8



def produce_env_init(q_chi,q_D):
 c1, c2,c3,c4=makec1(q_chi,q_D)
 Ta1, Tb1=makeTab(q_chi,q_D)
 Ta2, Tb2=makeTab(q_chi,q_D)
 Ta3, Tb3=makeTab1(q_chi,q_D)
 Ta4, Tb4=makeTab1(q_chi,q_D)
 Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]
 
 return Env


def produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d):
 Landa1=[Landa[2],Landa[1],Landa[0],Landa[3]]
 a_u,a=makeab(Landa1,Gamma_a)
 Landa1=[Landa[0],Landa[6],Landa[2],Landa[7]]
 b_u,b=makeab(Landa1,Gamma_b)
 Landa1=[Landa[4],Landa[3],Landa[5],Landa[1]]
 c_u,c=makeab(Landa1,Gamma_c)
 Landa1=[Landa[5],Landa[7],Landa[4],Landa[6]]
 d_u,d=makeab(Landa1,Gamma_d)
 return a_u,b_u,c_u,d_u,a,b,c,d

def produce_gamma_abcd(a_u,b_u,c_u,d_u,Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8):
 Gamma_a=copy.copy(a_u)
 Gamma_b=copy.copy(b_u)
 Gamma_c=copy.copy(c_u)
 Gamma_d=copy.copy(d_u)


 Landa_1=uni10.UniTensor([b_u.bond(1),a_u.bond(3)])
 Landa_3=uni10.UniTensor([a_u.bond(1),b_u.bond(3)])

 Landa_2=uni10.UniTensor([a_u.bond(2),c_u.bond(4)])
 Landa_4=uni10.UniTensor([c_u.bond(2),a_u.bond(4)])


 Landa_5=uni10.UniTensor([c_u.bond(1),d_u.bond(3)])
 Landa_6=uni10.UniTensor([d_u.bond(1),c_u.bond(3)])


 Landa_7=uni10.UniTensor([b_u.bond(2),d_u.bond(4)])
 Landa_8=uni10.UniTensor([d_u.bond(2),b_u.bond(4)])
 
 Landa_1.identity();  Landa_2.identity();  Landa_3.identity();
 Landa_4.identity();  Landa_5.identity();  Landa_6.identity();
 Landa_7.identity();  Landa_8.identity();  


 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8


def Initialize_function(Gamma,Landa):
 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
 #Gamma[0].randomize()
 for i in xrange(len(Gamma)):
  Gamma[i].randomize()
  #Gamma[i].orthoRand()
  Gamma[i]=Gamma[i]*(1.00/MaxAbs(Gamma[i]))
  
 for i in xrange(len(Landa)):
  #Landa[i].identity()
  blk_qnums = Landa[i].blockQnum()
  for qnum in blk_qnums:
    M=Landa[i].getBlock(qnum)
    if qnum == q0_even:
     M[0]=1.00
     M.randomize()
     #M.orthoRand()
     #M=M*M
     #M.identity()
     #print "M0", M
    else: 
     M[0]=0.1
     M.randomize()
     #M.orthoRand()
     #M.identity()

    Landa[i].putBlock(qnum,M)
  Landa[i]=Landa[i]*(1.00/MaxAbs(Landa[i]))

 
def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(0.5)*uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0])
  return Mat 

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(0.5)*uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);
  return Mat 

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(0.5)*uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);
  return Mat 


def matIden():
    spin_t=0.5
    dimT = int(2*spin_t + 1)
    Mat=uni10.Matrix(dimT, dimT,[1,0,0,1])
    return Mat


#cx_mat X1(3,3);  X1.zeros();         X1(0,1)=1; X1(1,0)=1; X1(1,2)=1,X1(2,1)=1;
#cx_mat Z1(3,3);  Z1.zeros();         Z1(0,0)=1; Z1(1,1)=0; Z1(2,2)=-1;
#cx_mat Y1(3,3);  Y1.zeros();         Y1(0,1)=-1; Y1(1,0)=1; Y1(1,2)=-1,Y1(2,1)=1;


#def matSx():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=(1.0/(2.0**(0.5)))*uni10.Matrix(dimT, dimT,[0, 1.0, 0 ,1.0,0, 1.0,0,1.0,0])
#    return Mat 
#    
#def matSy():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=(1.0/(2.0**(0.5)))*uni10.Matrix(dimT, dimT,[0,-1.0,0,1.0,0,-1.0,0,1.0,0])
#    return Mat 
#   

#def matSz():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=uni10.Matrix(dimT, dimT,[1,0,0,0,0,0,0,0,-1])
#    return Mat

#def matIden():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=uni10.Matrix(dimT, dimT,[1,0,0,0,1,0,0,0,1])
#    return Mat



def Heisenberg0(h, d_phys):

    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden=matIden()
    iden = matIden()

    ham =(h[0])*(h[3]*uni10.otimes(sz,sz)+h[4]*uni10.otimes(sx,sx)+h[5]*(-1.0)*uni10.otimes(sy,sy))
    H.putBlock(ham)
    return H



def Heisenberg00(h, d_phys):

    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden=matIden()
    iden = matIden()

    ham =(h[1])*(h[3]*uni10.otimes(sz,sz)+h[4]*uni10.otimes(sx,sx)+h[5]*(-1.0)*uni10.otimes(sy,sy))
    H.putBlock(ham)
    return H




def Heisenberg1(h, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden=matIden()
    iden = matIden()

    ham =(h[2]*2)*(h[3]*uni10.otimes(sz,sz)+h[4]*uni10.otimes(sx,sx)+h[5]*(-1.0)*uni10.otimes(sy,sy))
    H.putBlock(ham)
    return H

def Heisenberg1_Z2(J2,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    blk_qnums=H.blockQnum()
    M=H.getBlock(blk_qnums[0])
    M[0]=J2
    M[1]=0.0
    M[2]=0.0
    M[3]=J2
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=(-J2)
    M[1]=2.0*J2
    M[2]=2.0*J2
    M[3]=(-J2)
    H.putBlock(blk_qnums[1],M)
##########################################################################
    return H



def threebody(h,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)

    H = uni10.UniTensor([bdi, bdi,bdi, bdo, bdo,bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    szt=uni10.otimes(sz,sz)
    sxt=uni10.otimes(sx,sx)
    syt=(-1.0)*uni10.otimes(sy,sy)
    #ident=uni10.otimes(iden,iden)
    iden = matIden()

    sztt=uni10.otimes(iden,sz)
    sxtt=uni10.otimes(iden,sx)
    sytt=(-1.0)*uni10.otimes(iden,sy)

    ham =(h[0])*(h[3]*uni10.otimes(szt,iden)+h[4]*uni10.otimes(sxt,iden)+h[5]*uni10.otimes(syt,iden))

    
    ham =ham + h[1]*(h[3]*uni10.otimes(iden,szt)+h[4]*uni10.otimes(iden,sxt)+h[5]*uni10.otimes(iden,syt))


    ham =ham + 2.0*h[2]*(h[3]*uni10.otimes(sz,sztt)+h[4]*uni10.otimes(sx,sxtt)+h[5]*uni10.otimes(sy,sytt))


    A=(h[3]*uni10.otimes(szt,iden)+h[4]*uni10.otimes(sxt,iden)+h[5]*uni10.otimes(syt,iden))*(h[3]*uni10.otimes(szt,iden)+h[4]*uni10.otimes(sxt,iden)+h[5]*uni10.otimes(syt,iden))
    

    ham= ham +  h[6]*A 

    B=(h[3]*uni10.otimes(iden,szt)+h[4]*uni10.otimes(iden,sxt)+h[5]*uni10.otimes(iden,syt))*(h[3]*uni10.otimes(iden,szt)+h[4]*uni10.otimes(iden,sxt)+h[5]*uni10.otimes(iden,syt))
    
    ham= ham +  h[6]*B 


    C=2.0*(h[3]*uni10.otimes(sz,sztt)+h[4]*uni10.otimes(sx,sxtt)+h[5]*uni10.otimes(sy,sytt))*(h[3]*uni10.otimes(sz,sztt)+h[4]*uni10.otimes(sx,sxtt)+h[5]*uni10.otimes(sy,sytt))
    
    ham= ham +  h[7]*C 

#    H.setRawElem(ham)
    H.putBlock(ham)
#    print ham, H
#    print H
    #print sx, sy, sz
    return H

def threebody_U1(h,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)

    H = uni10.UniTensor([bdi, bdi,bdi, bdo, bdo,bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    szt=uni10.otimes(sz,sz)
    sxt=uni10.otimes(sx,sx)
    syt=(-1.0)*uni10.otimes(sy,sy)
    #ident=uni10.otimes(iden,iden)
    iden = matIden()

    sztt=uni10.otimes(iden,sz)
    sxtt=uni10.otimes(iden,sx)
    sytt=(-1.0)*uni10.otimes(iden,sy)

    ham =(h[0])*(h[3]*uni10.otimes(szt,iden)+h[4]*uni10.otimes(sxt,iden)+h[5]*uni10.otimes(syt,iden))

    
    ham =ham + h[1]*(h[3]*uni10.otimes(iden,szt)+h[4]*uni10.otimes(iden,sxt)+h[5]*uni10.otimes(iden,syt))


    ham =ham + 2.0*h[2]*(h[3]*uni10.otimes(sz,sztt)+h[4]*uni10.otimes(sx,sxtt)+h[5]*uni10.otimes(sy,sytt))


    A=(h[3]*uni10.otimes(szt,iden)+h[4]*uni10.otimes(sxt,iden)+h[5]*uni10.otimes(syt,iden))*(h[3]*uni10.otimes(szt,iden)+h[4]*uni10.otimes(sxt,iden)+h[5]*uni10.otimes(syt,iden))
    

    ham= ham +  h[6]*A 

    B=(h[3]*uni10.otimes(iden,szt)+h[4]*uni10.otimes(iden,sxt)+h[5]*uni10.otimes(iden,syt))*(h[3]*uni10.otimes(iden,szt)+h[4]*uni10.otimes(iden,sxt)+h[5]*uni10.otimes(iden,syt))
    
    ham= ham +  h[6]*B 


    C=2.0*(h[3]*uni10.otimes(sz,sztt)+h[4]*uni10.otimes(sx,sxtt)+h[5]*uni10.otimes(sy,sytt))*(h[3]*uni10.otimes(sz,sztt)+h[4]*uni10.otimes(sx,sxtt)+h[5]*uni10.otimes(sy,sytt))
    
    ham= ham +  h[7]*C 

    H.setRawElem(ham)
#    print ham, H
#    print H
    #print sx, sy, sz
    return H

def Heisenberg0_U1(h, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden=matIden()
    iden = matIden()

    ham =(h[0])*(h[3]*uni10.otimes(sz,sz)+h[4]*uni10.otimes(sx,sx)+h[5]*(-1.0)*uni10.otimes(sy,sy))
    H.setRawElem(ham)
    #print H, ham

    return H

def Heisenberg00_U1(h, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden=matIden()

    ham =(h[1])*(h[3]*uni10.otimes(sz,sz)+h[4]*uni10.otimes(sx,sx)+h[5]*(-1.0)*uni10.otimes(sy,sy))
    H.setRawElem(ham)
    #print H, ham

    return H


#def threebody_U1_help(h,d_phys):
#    bdi = uni10.Bond(uni10.BD_IN, d_phys)
#    bdo = uni10.Bond(uni10.BD_OUT, d_phys)

#    H = uni10.UniTensor([bdi, bdi,bdi, bdo, bdo,bdo], "Heisenberg")
#    sx = matSx()
#    sy = matSy()
#    sz = matSz()
#    szt=uni10.otimes(sz,sz)
#    sxt=uni10.otimes(sx,sx)
#    syt=(-1.0)*uni10.otimes(sy,sy)
#    #ident=uni10.otimes(iden,iden)
#    iden = matIden()

#    sztt=uni10.otimes(iden,sz)
#    sxtt=uni10.otimes(iden,sx)
#    sytt=(-1.0)*uni10.otimes(iden,sy)


#    ham =(0.25*h[1])*(h[0]*uni10.otimes(szt,iden)+uni10.otimes(sxt,iden)+uni10.otimes(syt,iden))
#    ham =ham + 0.25*h[1]*(h[0]*uni10.otimes(iden,szt)+uni10.otimes(iden,sxt)+uni10.otimes(iden,syt))
#    ham =ham + 0.25*2*h*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

#    #ham = 0.25*2.00*h*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

#    H.setRawElem(ham)
#    #print ham, H
#    #print sx, sy, sz
#    return H



def threebody_Z2(h,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)

    H = uni10.UniTensor([bdi, bdi,bdi, bdo, bdo,bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden=matIden()
    szt=uni10.otimes(sz,sz)
    sxt=uni10.otimes(sx,sx)
    syt=(-1.0)*uni10.otimes(sy,sy)
    ident=uni10.otimes(iden,iden)

    sztt=uni10.otimes(iden,sz)
    sxtt=uni10.otimes(iden,sx)
    sytt=(-1.0)*uni10.otimes(iden,sy)


    ham =(0.25*h[1])*(h[0]*uni10.otimes(szt,iden)+uni10.otimes(sxt,iden)+uni10.otimes(syt,iden))
    ham =ham + 0.25*h[1]*(h[0]*uni10.otimes(iden,szt)+uni10.otimes(iden,sxt)+uni10.otimes(iden,syt))
    ham =ham + 0.25*2.00*h*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

    #ham = 0.25*2.00*h*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

    H.setRawElem(ham)
    #print ham, H

    return H


def Heisenberg0_Z2(h,J1,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    blk_qnums=H.blockQnum()
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    M[1]=0.0
    M[2]=0.0
    M[3]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=(-J1)*h
    M[1]=2.0*J1
    M[2]=2.0*J1
    M[3]=(-J1)*h
    H.putBlock(blk_qnums[1],M)
##########################################################################

    return H

def Heisenberg00_Z2(h,J1,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    blk_qnums=H.blockQnum()
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    M[1]=0.0
    M[2]=0.0
    M[3]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=(-J1)*h
    M[1]=2.0*J1
    M[2]=2.0*J1
    M[3]=(-J1)*h
    H.putBlock(blk_qnums[1],M)
##########################################################################

    return H








def Heisenberg0_U1Z2(h,J1, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    #H.randomize()
    #print transverseIsing(h).getBlock()
    #H.setRawElem(transverseIsing(h).getBlock().getElem());
    #H.setRawElem(Heisenberg().getBlock());
    blk_qnums=H.blockQnum()
    blk_qnums[0]
    
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=J1*(-h)
    M[1]=J1*2.0
    M[2]=J1*2.0
    M[3]=J1*(-h)
    H.putBlock(blk_qnums[1],M)


    M=H.getBlock(blk_qnums[2])
    M[0]=J1*h
    H.putBlock(blk_qnums[2],M)
    H=H*0.250
    #print H
    return H









def makeTab(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_IN, D)
 #bdi1.combine(copy.copy(bdi1))
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem0.randomize()
 Tem0*=(1.00/MaxAbs(Tem0))
 Tem1=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem1.randomize()
 Tem1*=(1.00/MaxAbs(Tem1))
 return Tem0, Tem1

def makeTab1(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_OUT, D)
 #bdi1.combine(copy.copy(bdi1))
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem0.randomize()
 Tem0*=(1.00/Tem0.norm())
 Tem1=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem1.randomize()
 Tem1*=(1.00/Tem1.norm())
 return Tem0, Tem1

def MaxAbs(c):
 blk_qnums = c.blockQnum()
 max_list=[]
 for qnum in blk_qnums:
    c_mat=c.getBlock(qnum)
    max_list.append(c_mat.absMax())
 #sv_mat = uni10.Matrix( len(max_list), len(max_list), max_list, True)
 #return sv_mat.absMax()
 max_list_f=[abs(x) for x in max_list]
 #print max_list_f, max(max_list_f)
 return max(max_list_f)

def makec1(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c1=uni10.UniTensor([bdi, bdo])
 c2=uni10.UniTensor([bdi, bdi])
 c3=uni10.UniTensor([bdi, bdo])
 c4=uni10.UniTensor([bdo, bdo])
 c1.randomize()
 c2.randomize()
 c3.randomize()
 c4.randomize()
 c1*=(1.00/MaxAbs(c1))
 c2*=(1.00/MaxAbs(c2))
 c3*=(1.00/MaxAbs(c3))
 c4*=(1.00/MaxAbs(c4))
 return c1,c2,c3,c4
 
def makeab(Landa,Gamma):
 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa)) ]
 Landa_sq=sqrt(Landa_cp)
 a_u=copy.copy(Gamma)
 Landa_sq[0].setLabel([-1,1])
 Landa_sq[1].setLabel([-2,2])
 Landa_sq[2].setLabel([3,-3])
 Landa_sq[3].setLabel([4,-4])


 a_u.setLabel([0,1,2,3,4])
 
 a_u=(((((a_u*Landa_sq[0])*Landa_sq[1])*Landa_sq[2])*Landa_sq[3]))
 a_u.permute([0,-1,-2,-3,-4],3)

 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-3,-4,0,-1,-2])

 a_u.setLabel([0,1,2,3,4])

 a=a_u*a_d
 a.permute([1,-1,2,-2,3,-3,4,-4],4)

 return a_u, a


 
def sqrt(Landa):
 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa))   ]
 for q in xrange(len(Landa_cp)): 
  blk_qnums=Landa_cp[q].blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp[q].getBlock(qnum).col())
   Landa_cpm=Landa_cp[q].getBlock(qnum)
   Landam=Landa[q].getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
   Landa_cp[q].putBlock(qnum,Landa_cpm)
 return Landa_cp

 
def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1)
  invLt = uni10.Matrix(D, D1)
  invLt=Landa2.getBlock(qnum)
  for i in xrange(D):
    for j in xrange(D1):
     invL2[i*D1+j] = 0 if ((invLt[i*D1+j].real) < 1.0e-12) else (1.00 / (invLt[i*D1+j].real))
  invLanda2.putBlock(qnum,invL2)
 return invLanda2


def make_ab(a_u):
 a_u.setLabel([0,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.transpose()
 a_uc.setLabel([-3,-4,0,-1,-2])
 result=a_uc*a_u
 result.permute([1,-1,2,-2,3,-3,4,-4], 4)
 return result




def corner_transfer_matrix_twosite(a,b,c,d,chi,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D):
 z1=copy.copy(a)
 z1.identity()
 z2=copy.copy(a)
 z2.randomize()
 z1=z2#+(1.0e-1)*z2
 
 Accuracy=1.00e-7
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 #print  '\n', '\n', 'CTM'
 
 while Loop_iter is 0: 

  c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4=Move.add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D)
  
  norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  norm1=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
  E0=E1
  if (abs(norm) > 1.00e-10):
   E1=abs(norm1)/abs(norm)
   if (abs((E0-E1)/E0) < Accuracy):Loop_iter=1;
  else:
   E1=abs(norm1)
   if (abs((E0-E1)) < Accuracy) : print 'Warning: norm~0', E1; Loop_iter=1;
  count+=1
  if (count > 12 ): print 'break! CTM'; break;
  print E1, abs((E0-E1)/E1),norm, count
  #print E1, Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()
 
 
 #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4)
 #print 'CTM', norm
 return c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4



def corner_transfer_matrix_twosite_CTMRG(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,direct_env,N_env):
 z1=copy.copy(a)
 z1.identity()
 z2=copy.copy(a)
 z2.randomize()
 z1=z2#+(1.0e-1)*z2
 
 Accuracy=N_env[1]
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 #print  '\n', '\n', 'CTM' 
 while Loop_iter is 0: 
  t0=time.time()
  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveCorboz.add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D)
  
  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveCorboz.add_left1(c1,c2,c3,c4,Tb1,Ta2,Tb3,Ta4,Ta1,Tb2,Ta3,Tb4,b,a,d,c,chi,D) 

  MoveCorboz.permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveCorboz.add_left1(c1,c4,c3,c2,Ta4,Ta3,Ta2,Ta1,Tb4,Tb3,Tb2,Tb1,a,c,b,d,chi,D)
  

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveCorboz.add_left1(c1,c4,c3,c2,Tb4,Ta3,Tb2,Ta1,Ta4,Tb3,Ta2,Tb1,c,a,d,b,chi,D)

  MoveCorboz.permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4=MoveCorboz.equall_norm(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)

  
  norm=MoveCorboz.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
#  norm1=MoveCorboz.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
  norm_print=norm
#  norm=norm[0]
#  norm1=norm1[0]
   
  norm=1.0
  if direct_env is 'h':
   norm1=MoveCorboz.Env_energy_h(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'v':
   norm1=MoveCorboz.Env_energy_v(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'D':
   norm1=MoveCorboz.Env_energy_D(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'D1':
   norm1=MoveCorboz.Env_energy_D1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'three':
   norm1=MoveCorboz.Env_energy_three(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,c_u,d_u,H0,D,d_phys)
  elif direct_env is 'three1':
   norm1=MoveCorboz.Env_energy_three1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,c_u,d_u,H0,D,d_phys)


  E0=E1
  if (abs(norm) > 1.00e-10):
   E1=abs(norm1)/abs(norm)
   if (abs((E0-E1)/E0) < Accuracy):Loop_iter=1;
  else:
   E1=abs(norm1)
   if (abs((E0-E1)) < Accuracy) : print 'Warning: norm~0', E1; Loop_iter=1;
  count+=1
  if (count > N_env[0] ): print 'break! CTMRG'; break;
  print E1, abs((E0-E1)/E1),norm_print, count,  time.time() - t0, "CTMRG"
  #print E1, Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()
  
 #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4)
 #print 'CTM', norm
 return c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4

#@profile 
def corner_transfer_matrix_twosite_CTMFull(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,direct_env,N_env):
 z1=copy.copy(a)
 z1.identity()
 z2=copy.copy(a)
 z2.randomize()
 z1=z2*100#+(1.0e-1)*z2
 #z1=200*z2
 
 Accuracy=N_env[1]
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 while Loop_iter is 0: 
  t0=time.time()
  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveFull.add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D)

  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveFull.add_left1(c1,c2,c3,c4,Tb1,Ta2,Tb3,Ta4,Ta1,Tb2,Ta3,Tb4,b,a,d,c,chi,D) 

  MoveFull.permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveFull.add_left1(c1,c4,c3,c2,Ta4,Ta3,Ta2,Ta1,Tb4,Tb3,Tb2,Tb1,a,c,b,d,chi,D)
  

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveFull.add_left1(c1,c4,c3,c2,Tb4,Ta3,Tb2,Ta1,Ta4,Tb3,Ta2,Tb1,c,a,d,b,chi,D)

  MoveFull.permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4=MoveFull.equall_norm(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  
  norm=MoveFull.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
#  norm1=MoveFull.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
  norm_print=norm[0]
#  norm=norm[0]
#  norm1=norm1[0]
   
  norm=1.0
  if direct_env is 'h':
   norm1=MoveFull.Env_energy_h(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'v':
   norm1=MoveFull.Env_energy_v(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'D':
   norm1=MoveFull.Env_energy_D(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'D1':
   norm1=MoveFull.Env_energy_D1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
  elif direct_env is 'three':
   norm1=MoveFull.Env_energy_three(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,c_u,d_u,H0,D,d_phys)
  elif direct_env is 'three1':
   norm1=MoveFull.Env_energy_three1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,c_u,d_u,H0,D,d_phys)



  E0=E1
  if (abs(norm) > 1.00e-10):
   E1=abs(norm1)/abs(norm)
   if (abs((E0-E1)/E0) < Accuracy):Loop_iter=1;
  else:
   E1=abs(norm1)
   if (abs((E0-E1)) < Accuracy) : print 'Warning: norm~0', E1; Loop_iter=1;
  count+=1
  if (count > N_env[0] ): print 'break! CTMFull'; break;
  #print E1, abs((E0-E1)/E1),norm, count, time.time() - t0,"CTMFull"
  print E1, abs((E0-E1)/E1),norm_print, count, time.time() - t0,"CTMFull"

  #print E1, Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()

 #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 #print 'CTM', norm
 return c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4
 



def E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model,N_env):

############################################################################
# E_val1=Energy_cb(c_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val1", E_val1

# E_val2=Energy_ad(a_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val2", E_val2

# E_val3=Energy_cb(d_u,a_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val3", E_val3

# E_val4=Energy_ad(b_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val4", E_val4
# 
# E_val5=Energy_cb(a_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val5", E_val5

# E_val6=Energy_ad(c_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val6", E_val6

# E_val7=Energy_cb(b_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val7", E_val7

# E_val8=Energy_ad(d_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val8", E_val8
############################################################################
 E_1=Energy_cab(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_2=Energy_cab(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_3=Energy_cab(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_4=Energy_cab(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_5=Energy_abd(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_6=Energy_abd(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_7=Energy_abd(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,N_env)
 E_8=Energy_abd(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)


 print '\n', E_1, '  ', E_2, '\n', E_3, '  ', E_4, '\n'

 print E_5, '  ', E_6, '\n', E_7, '  ', E_8, '\n'


 #return E_5
# return (E_5+E_6+E_7+E_8)/4.00

 #return (E_1+E_2+E_3+E_4)/4.00 
 return (E_1+E_2+E_3+E_4+E_5+E_6+E_7+E_8)/8.00
# return (E_ab+E_ba)/2.0
# return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc) / 4.00)#+((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc)/4.00)


def E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model,N_env):

############################################################################
# E_val1=Energy_cb(c_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val1", E_val1

# E_val2=Energy_ad(a_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val2", E_val2

# E_val3=Energy_cb(d_u,a_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val3", E_val3

# E_val4=Energy_ad(b_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val4", E_val4
# 
# E_val5=Energy_cb(a_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val5", E_val5

# E_val6=Energy_ad(c_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val6", E_val6

# E_val7=Energy_cb(b_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val7", E_val7

# E_val8=Energy_ad(d_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
# #print "E_val8", E_val8
############################################################################

 E_1=Energy_cab(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env)

 E_2=Energy_cab(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,N_env)

 E_3=Energy_cab(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,N_env)

 E_4=Energy_cab(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
 
 E_5=Energy_abd(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env)

 E_6=Energy_abd(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,N_env)

 E_7=Energy_abd(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,N_env)

 E_8=Energy_abd(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)


 print '\n', E_1, '  ', E_2, '\n', E_3, '  ', E_4, '\n'
 
 print E_5, '  ', E_6, '\n', E_7, '  ', E_8, '\n'

 #return E_5
# return (E_5+E_6+E_7+E_8)/4.00

 #return (E_1+E_2+E_3+E_4)/4.00 
 return (E_1+E_2+E_3+E_4+E_5+E_6+E_7+E_8)/8.00
# return (E_ab+E_ba)/2.0
# return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc) / 4.00)#+((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc)/4.00)



def Energy_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)
 
# Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

# if Corner_method is 'CTM':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
# if Corner_method is 'CTMRG':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
# if Corner_method is 'CTMFull':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)

# reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 
 if Model is "Heisenberg":
   H0=Heisenberg0(h,d_phys)
   H00=Heisenberg00(h,d_phys)
   H1=Heisenberg1(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h,d_phys)
   H00=Heisenberg00_Z2(h,d_phys)
   H1=Heisenberg1_Z2(h,d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h,d_phys)
   H00=Heisenberg0_U1(h,d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h,d_phys)
   H00=Heisenberg0_U1Z2(h,d_phys)




 E_ca=basicB.Energy_ca(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H00,c_u,a_u)

# E_ca=basicB.Energy_ca_positive(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H00,c_u,a_u)
# print "E_ca", E_ca
 return E_ca




def Energy_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env):

 if Model is "Heisenberg":
   H0=Heisenberg0(h,d_phys)
   H1=Heisenberg1(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h,d_phys)
   H1=Heisenberg1_Z2(h,d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h,d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h,d_phys)
   H00=Heisenberg0_U1Z2(h,d_phys)

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 #c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Reload_Env()

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D)
 if Corner_method is 'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h',N_env)
 if Corner_method is 'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h',N_env)
  #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 

 reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4, Tb4,D,d_phys)

 E_ab=basicB.Energy_ab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0, a_u, b_u)

# E_ab=basicB.Energy_ab_positive(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0, a_u, b_u)
# print "E_ab", E_ab
 return E_ab


def Energy_cab(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env):
 
 if Model is "Heisenberg":
   H0=Heisenberg0(h,d_phys)
   H00=Heisenberg00(h,d_phys)
   H1=Heisenberg1(h,d_phys)
   H2=threebody(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h,d_phys)
   H00=Heisenberg00_Z2(h,d_phys)
   H1=Heisenberg1_Z2(h,d_phys)
   H2=threebody_Z2(h,d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h,d_phys)
   H00=Heisenberg0_U1(h,d_phys)
   H2=threebody_U1(h,d_phys)
   #H2=threebody_U1_help(h,d_phys)
   
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h,d_phys)
   H00=Heisenberg0_U1Z2(h,d_phys)
 
 
 #H0.randomize()
 
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 #c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Reload_Env()


 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 #t0=time.time()
 if Corner_method is 'CTM':
#  c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1=basic.make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1)
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D,H0,d_phys)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H2,d_phys,'three',N_env)
  #basic.Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H2,d_phys,'three',N_env)
 #print time.time() - t0, "CTM-H, Left"
  
  
 reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)

 E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=basicC.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u)


 E=basicC.energy_cab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, H2)

 return E

def Energy_abd(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,N_env):
 
 if Model is "Heisenberg":
   H0=Heisenberg0(h,d_phys)
   H00=Heisenberg00(h,d_phys)
   H1=Heisenberg1(h,d_phys)
   H2=threebody(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h,d_phys)
   H00=Heisenberg00_Z2(h,d_phys)
   H1=Heisenberg1_Z2(h,d_phys)
   H2=threebody_Z2(h,d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h,d_phys)
   H00=Heisenberg0_U1(h,d_phys)
   H2=threebody_U1(h,d_phys)
   #H2=threebody_U1_help(h,d_phys)
   
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h,d_phys)
   H00=Heisenberg0_U1Z2(h,d_phys)
 
 
 #H0.randomize()
 
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 #c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Reload_Env()


 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

# #t0=time.time()
# if Corner_method is 'CTM':
##  c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1=basic.make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1)
#  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D,H0,d_phys)
# if Corner_method is'CTMRG':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H2,d_phys,'three1')
#  #basic.Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
# if Corner_method is'CTMFull':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H2,d_phys,'three1')
# #print time.time() - t0, "CTM-H, Left"

  
  
 reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)

 E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=basicC.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u)


 E=basicC.energy_abd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, H2)

 return E

 

def M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model):


 if Model is "Heisenberg":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden =matIden() 
  szt=uni10.otimes(sz,iden)
  szt1=uni10.otimes(iden,sz)
  H.putBlock(szt)
  H1.putBlock(szt1)
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
#  print E_1
  E_2=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H1)
#  print E_2
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
#  print E_3
  E_4=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H1)
  print "z"
  print E_1, "    ", E_2, "\n",  E_3, "    ", E_4 
  E_z=((abs(E_1)+abs(E_2)+abs(E_3)+abs(E_4)) / 4.00)
  szt=uni10.otimes(sx,iden)
  szt1=uni10.otimes(iden,sx)
  H.putBlock(szt)
  H1.putBlock(szt1)
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  #print E_1
  E_2=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H1)
  #print E_2
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  #print E_3
  E_4=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H1)
  print "x"
  print E_1, "    ", E_2, "\n",  E_3, "    ", E_4 
  E_x=((abs(E_1)+abs(E_2)+abs(E_3)+abs(E_4)) / 4.00)
  szt=uni10.otimes(sy,iden)
  szt1=uni10.otimes(iden,sy)
  H.putBlock(szt)
  H1.putBlock(szt1)
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
#  print E_1
  E_2=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H1)
#  print E_2
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
#  print E_3
  E_4=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H1)
  print "y"
  print E_1, "    ", E_2, "\n",  E_3, "    ", E_4 
  E_y=((abs(E_1)+abs(E_2)+abs(E_3)+abs(E_4)) / 4.00)
  return (E_x+E_y+E_z) 
 if Model is "Heisenberg_Z2":
  #print d_phys
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  iden = matIden()
  szt=uni10.otimes(sz,iden)
  H.setRawElem(szt)
  szt1=uni10.otimes(iden,sz)
  H1.setRawElem(szt1)
  #print szt,H,szt1,H1
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  print E_1
  E_2=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H1)
  print E_2
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  print E_3
  E_4=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H1)
  print E_4
  return ((abs(E_1)+abs(E_2)+abs(E_3)+abs(E_4)) / 4.00)
 if Model is "Heisenberg_U1":
  #print d_phys
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  iden=matIden()
  szt=uni10.otimes(sz,iden)
  H.setRawElem(szt)
  szt1=uni10.otimes(iden,sz)
  H1.setRawElem(szt1)
  #print szt,H,szt1,H1
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  #print 
  E_2=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H1)
  #print E_2
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  #print E_3
  E_4=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H1)
  print "z"
  print E_1, "    ", E_2, "\n",  E_3, "    ", E_4 
  return ((abs(E_1)+abs(E_2)+abs(E_3)+abs(E_4)) / 4.00)
 if Model is "Heisenberg_U1Z2":
  H0=Heisenberg0_U1Z2(h,d_phys)
  H00=Heisenberg0_U1Z2(h,d_phys)
  H1=Heisenberg1_U1(h,d_phys)



def Translational_sym(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model):


 if Model is "Heisenberg":
  H=Heisenberg0(h,d_phys)
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  E_2=M_h(b_u,a_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,H)
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  E_4=M_h(d_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,H)
  print 'bond-H'
  print E_1, "    ", E_2, "\n",  E_3, "    ", E_4 
  max_list=[E_1, E_2, E_3, E_4]
  print 'max', max(max_list) 
  print 'min', min(max_list) 
  H=Heisenberg00(h,d_phys)
  E_5=M_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  E_6=M_v(d_u,b_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,H)
  E_7=M_v(a_u,c_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  E_8=M_v(b_u,d_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,H)
  print 'bond-V'
  print E_5, "    ", E_6, "\n",  E_7, "    ", E_8 
  max_list1=[E_5, E_6, E_7, E_8]
  print 'max', max(max_list1) 
  print 'min', min(max_list1) 

  D_x=max(max_list) - min(max_list)
  D_y=max(max_list1) - min(max_list1)

  print   'D_x=', D_x
  print   'D_y=', D_y
  
  
 if Model is "Heisenberg_U1":
  H=Heisenberg0_U1(h,d_phys)
  E_1=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  E_2=M_h(b_u,a_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,H)
  E_3=M_h(c_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  E_4=M_h(d_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,H)
  print 'bond-H'
  print E_1, "    ", E_2, "\n",  E_3, "    ", E_4 
  max_list=[E_1, E_2, E_3, E_4]
  print 'max', max(max_list) 
  print 'min', min(max_list) 
  H=Heisenberg00_U1(h,d_phys)

  E_5=M_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H)
  E_6=M_v(d_u,b_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model,H)
  E_7=M_v(a_u,c_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model,H)
  E_8=M_v(b_u,d_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,H)
  print 'bond-V'
  print E_5, "    ", E_6, "\n",  E_7, "    ", E_8 
  max_list1=[E_5, E_6, E_7, E_8]
  print 'max', max(max_list1) 
  print 'min', min(max_list1) 

  D_x=max(max_list) - min(max_list)
  D_y=max(max_list1) - min(max_list1)

  print   'D_x=', D_x
  print   'D_y=', D_y
  
 if Model is "Heisenberg_U1Z2":
  H0=Heisenberg0_U1Z2(h,d_phys)
  H00=Heisenberg0_U1Z2(h,d_phys)
  H1=Heisenberg1_U1(h,d_phys)



def M_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H00):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)
 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 E_ca=basicB.Energy_ca(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H00,c_u,a_u)
 return E_ca


def M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H0):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)
 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 E_ab=basicB.Energy_ab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0,a_u,b_u)
 return E_ab












#@profile
def CorrelationH(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,distance_final,fileCorr,fileCorrLength):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 if Model is "Heisenberg":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)(-1.0)*uni10.otimes(sy,sy)
  H=uni10.otimes(sz,iden)+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1=uni10.otimes(iden,sz)+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 if Model is "Heisenberg_Z2":
  #print d_phys
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  iden = matIden()
  szt=uni10.otimes(sz,iden)
  H.setRawElem(szt)
  szt1=uni10.otimes(iden,sz)
  H1.setRawElem(szt1)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])
 if Model is "Heisenberg_U1":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH_tem=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy)
  H_tem=uni10.otimes(sz,iden)#+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1_tem=uni10.otimes(iden,sz)#+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setRawElem(HH_tem)
  H.setRawElem(H_tem)
  H1.setRawElem(H1_tem)
  Iden=copy.copy(H)
  Iden.identity()
  #print HH_tem, HH, H_tem, H, H1_tem, H1
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 vec_left=make_vleft(Tb4,Ta4,c1,c4)
 vec_right=make_vright(Ta2,Tb2,c2,c3)
 ap=make_ap_openindex(a_u)
 bp=make_ap_openindex(b_u)
 cp=make_ap_openindex(c_u)
 dp=make_ap_openindex(d_u)
 dis_val_list=[]
 Corr_val_list=[]
 dis_val_list1=[]
 Corr_val_list1=[]
 dis_val_list2=[]
 Corr_val_list2=[]
 dis_val_list3=[]
 Corr_val_list3=[]
 vec_right_copy=copy.copy(vec_right)
 Corr_length=ED_right(c2, Ta2, Tb2, c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3,vec_right_copy)
 print "Corr_length",Corr_length
 fileCorrLength.write(str(Corr_length)  + "\n")
 fileCorrLength.flush()



 print "\n"
#######################a-a#####################################################

 vec_left_1=Make_first_vecleft_a(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)
 vec_right_1=Make_first_vecright_a(vec_right, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list.append(2)
 Corr_val_list.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list.append(dis_val)
  Corr_val_list.append(Corr_val) 
###################################################################################

 print "\n"
#######################b-b#####################################################

 vec_left_1=Make_first_vecleft_b(vec_left, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)
 vec_right_1=Make_first_vecright_b(vec_right, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list1.append(2)
 Corr_val_list1.append(Corr_val)

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1.append(dis_val)
  Corr_val_list1.append(Corr_val) 
###################################################################################


 print "\n"
#######################c-c#####################################################

 vec_left_1=Make_first_vecleft_c(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)
 vec_right_1=Make_first_vecright_c(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list2.append(2)
 Corr_val_list2.append(Corr_val)

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2.append(dis_val)
  Corr_val_list2.append(Corr_val) 
###################################################################################

 print "\n"
#######################d-d#####################################################

 vec_left_1=Make_first_vecleft_d(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)
 vec_right_1=Make_first_vecright_d(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list3.append(2)
 Corr_val_list3.append(Corr_val)

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3.append(dis_val)
  Corr_val_list3.append(Corr_val) 
###################################################################################

############################################################################################################


 dis_val_listo=[]
 Corr_val_listo=[]
 dis_val_list1o=[]
 Corr_val_list1o=[]
 dis_val_list2o=[]
 Corr_val_list2o=[]
 dis_val_list3o=[]
 Corr_val_list3o=[]

 print "\n"

######################a-bo#####################################################

 vec_left_1=Make_first_vecleft_a(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)
 vec_right_1=Make_first_vecright_b(vec_right, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,3)
 dis_val_listo.append(3)
 Corr_val_listo.append(Corr_val) 

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_listo.append(dis_val)
  Corr_val_listo.append(Corr_val) 
##################################################################################


 print "\n"

######################b-ao#####################################################

 vec_left_1=Make_first_vecleft_b(vec_left, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)
 vec_right_1=Make_first_vecright_a(vec_right, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)

 #Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,1)
 #dis_val_list1.append(2)
 #Corr_val_list1.append(Corr_val)

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1o.append(dis_val)
  Corr_val_list1o.append(Corr_val) 
##################################################################################


 print "\n"

#######################c-do#####################################################

 vec_left_1=Make_first_vecleft_c(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)
 vec_right_1=Make_first_vecright_d(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,3)
 dis_val_list2o.append(3)
 Corr_val_list2o.append(Corr_val)

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2o.append(dis_val)
  Corr_val_list2o.append(Corr_val) 
###################################################################################

 print "\n"

#######################d-co#####################################################

 vec_left_1=Make_first_vecleft_d(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)
 vec_right_1=Make_first_vecright_c(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)

# Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
# dis_val_list3.append(2)
# Corr_val_list3.append(Corr_val)

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3o.append(dis_val)
  Corr_val_list3o.append(Corr_val) 
###################################################################################


 print dis_val_list,'\n,\n'
 print Corr_val_list,'\n,\n'
 print Corr_val_list1,'\n,\n'
 print Corr_val_list2,'\n,\n'
 print Corr_val_list3,'\n,\n'

 print dis_val_listo,'\n,\n'
 print Corr_val_listo,'\n,\n'
 print Corr_val_list1o,'\n,\n'
 print Corr_val_list2o,'\n,\n'
 print Corr_val_list3o,'\n,\n'


 Corr_val_list_ave=[ (sum(t)*(1.0/16.0)) for t in zip(Corr_val_list, Corr_val_list1, Corr_val_list2, Corr_val_list3)]
 print Corr_val_list_ave,'\n,\n'


 Corr_val_list_avo=[ (sum(t)*(1.0/16.0)) for t in zip(Corr_val_listo, Corr_val_list1o, Corr_val_list2o, Corr_val_list3o)]
 print Corr_val_list_avo,'\n,\n'


 Corr_val_list_final=Corr_val_list_ave+Corr_val_list_avo
 dis_val_list_final=dis_val_list+dis_val_listo

 print '\n,\n'
 print dis_val_list_final
 print Corr_val_list_final
 print '\n,\n'


 for i in xrange(len(dis_val_list_final)):
  fileCorr.write(str(dis_val_list_final[i])  + " " + str(Corr_val_list_final[i]) + "\n")
  fileCorr.flush()

def CorrelationV(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,distance_final,fileCorr,fileCorrLength):
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)


 if Model is "Heisenberg":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)(-1.0)*uni10.otimes(sy,sy)
  H=uni10.otimes(sz,iden)+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1=uni10.otimes(iden,sz)+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 if Model is "Heisenberg_Z2":
  #print d_phys
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  iden = matIden()
  szt=uni10.otimes(sz,iden)
  H.setRawElem(szt)
  szt1=uni10.otimes(iden,sz)
  H1.setRawElem(szt1)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])
 if Model is "Heisenberg_U1":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH_tem=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy)
  H_tem=uni10.otimes(sz,iden)#+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1_tem=uni10.otimes(iden,sz)#+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setRawElem(HH_tem)
  H.setRawElem(H_tem)
  H1.setRawElem(H1_tem)
  Iden=copy.copy(H)
  Iden.identity()
  #print HH_tem, HH, H_tem, H, H1_tem, H1
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 vec_down=make_down(c4,Ta3, Tb3,c3)
 vec_up=make_up(c1,Tb1, Ta1,c2)
 ap=make_ap_openindex(a_u)
 bp=make_ap_openindex(b_u)
 cp=make_ap_openindex(c_u)
 dp=make_ap_openindex(d_u)
 dis_val_list=[]
 Corr_val_list=[]
 dis_val_list1=[]
 Corr_val_list1=[]
 dis_val_list2=[]
 Corr_val_list2=[]
 dis_val_list3=[]
 Corr_val_list3=[]

# Corr_length=ED_right(c2, Ta2, Tb2, c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3,vec_right)
# print "Corr_length",Corr_length
 vec_up_copy=copy.copy(vec_up)
 Corr_length=ED_up(c1,Tb1, Ta1,c2, a, b, c, d, Tb2, Ta2, Ta4, Tb4,vec_up_copy)
 print "Corr_length",Corr_length
 fileCorrLength.write(str(Corr_length)  + "\n")
 fileCorrLength.flush()


 print "\n"
#######################a-a#####################################################
 vec_down_1=Make_first_vecdown_a(vec_down, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_a(vec_up, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
 dis_val_list.append(2)
 Corr_val_list.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list.append(dis_val)
  Corr_val_list.append(Corr_val) 
###################################################################################

 print "\n"
########################b-b#####################################################

 vec_down_1=Make_first_vecdown_b(vec_down, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_b(vec_up, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
 dis_val_list1.append(2)
 Corr_val_list1.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1.append(dis_val)
  Corr_val_list1.append(Corr_val) 
####################################################################################


 print "\n"
########################c-c#####################################################

 vec_down_1=Make_first_vecdown_c(vec_down, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_c(vec_up, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
 dis_val_list2.append(2)
 Corr_val_list2.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2.append(dis_val)
  Corr_val_list2.append(Corr_val) 
####################################################################################

 print "\n"
########################d-d#####################################################

 vec_down_1=Make_first_vecdown_d(vec_down, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_d(vec_up, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
 dis_val_list3.append(2)
 Corr_val_list3.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3.append(dis_val)
  Corr_val_list3.append(Corr_val) 
####################################################################################

#############################################################################################################


 dis_val_listo=[]
 Corr_val_listo=[]
 dis_val_list1o=[]
 Corr_val_list1o=[]
 dis_val_list2o=[]
 Corr_val_list2o=[]
 dis_val_list3o=[]
 Corr_val_list3o=[]

 print "\n"

#######################a-co#####################################################

 vec_down_1=Make_first_vecdown_a(vec_down, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_c(vec_up, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)

# Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
# dis_val_listo.append(2)
# Corr_val_listo.append(Corr_val) 

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_listo.append(dis_val)
  Corr_val_listo.append(Corr_val) 

###################################################################################


 print "\n"

#######################c-ao#####################################################
 vec_down_1=Make_first_vecdown_c(vec_down, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_a(vec_up, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,3)
 dis_val_list1o.append(3)
 Corr_val_list1o.append(Corr_val) 

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1o.append(dis_val)
  Corr_val_list1o.append(Corr_val) 
###################################################################################


 print "\n"

########################b-do#####################################################

 vec_down_1=Make_first_vecdown_b(vec_down, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_d(vec_up, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)

# Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
# dis_val_list2.append(2)
# Corr_val_list2.append(Corr_val) 

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2o.append(dis_val)
  Corr_val_list2o.append(Corr_val) 
####################################################################################

 print "\n"

########################d-bo#####################################################

 vec_down_1=Make_first_vecdown_d(vec_down, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_b(vec_up, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,3)
 dis_val_list3o.append(3)
 Corr_val_list3o.append(Corr_val) 

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3o.append(dis_val)
  Corr_val_list3o.append(Corr_val) 
####################################################################################


 print dis_val_list,'\n,\n'
 print Corr_val_list,'\n,\n'
 print Corr_val_list1,'\n,\n'
 print Corr_val_list2,'\n,\n'
 print Corr_val_list3,'\n,\n'

 print dis_val_listo,'\n,\n'
 print Corr_val_listo,'\n,\n'
 print Corr_val_list1o,'\n,\n'
 print Corr_val_list2o,'\n,\n'
 print Corr_val_list3o,'\n,\n'


 Corr_val_list_ave=[ (sum(t)*(1.0/16.0)) for t in zip(Corr_val_list, Corr_val_list1, Corr_val_list2, Corr_val_list3)]
 print Corr_val_list_ave,'\n,\n'


 Corr_val_list_avo=[ (sum(t)*(1.0/16.0)) for t in zip(Corr_val_listo, Corr_val_list1o, Corr_val_list2o, Corr_val_list3o)]
 print Corr_val_list_avo,'\n,\n'


 Corr_val_list_final=Corr_val_list_ave+Corr_val_list_avo
 dis_val_list_final=dis_val_list+dis_val_listo

 print '\n,\n'
 print dis_val_list_final
 print Corr_val_list_final
 print '\n,\n'


 for i in xrange(len(dis_val_list_final)):
  fileCorr.write(str(dis_val_list_final[i])  + " " + str(Corr_val_list_final[i]) + "\n")
  fileCorr.flush()


def Mat_np_to_Uni(Mat_np):
 d0=np.size(Mat_np,0)
 d1=np.size(Mat_np,1)
 Mat_uni=uni10.Matrix(d0,d1)
 for i in xrange(d0):
  for j in xrange(d1):
   Mat_uni[i*d1+j]=Mat_np[i,j]
 return  Mat_uni


def Mat_nptoUni(Mat_np):
 d0=np.size(Mat_np,0)
 Mat_uni=uni10.Matrix(d0,d0, True)
 for i in xrange(d0):
   Mat_uni[i]=Mat_np[i]
 return  Mat_uni

 
def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np

def eig_np(A):
 D_eig=[A]*2
 A_np=Mat_uni_to_np(A)
 w, v = LA.eig(A_np)
 #print w,"\n,\n"#, v
 D_eig[0]=Mat_nptoUni(w)
 D_eig[1]=Mat_np_to_Uni(v)
 return D_eig

def  make_Q(q_vec): 
 D=int(q_vec[0].row())
 m=len(q_vec)
 Q=uni10.Matrix(D, m)
 for i in xrange(m):
  for j in xrange(D):
    Q[j*m+i]=q_vec[i][j]
 return Q

def return_vec(A, index ):
 D=int(A.row())
 vec_tem=uni10.Matrix(D,1)
 for i in xrange(D): 
  vec_tem[i]=A[i*D+index].real
 return vec_tem


def find_maxindex(A):
 D=int(A.row())
 max_val=0
 index=0
 #print A
 for i in xrange(D):
  if (i == 0) or ( max_val < abs(A[i]) ):
   max_val = abs(A[i])
   index=i
 return max_val, index, A[index] 


##############################################################################################
def ED_right(c2, Ta2, Tb2, c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3,Vec_uni):

 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D=",  D

 m=10
 W=2
 num=0
 E1=0
 p=0

 while p  <  (W+1):
  #print "norm", p, Vec_F.norm(),Vec_F[0], Vec_F[1], Vec_F[2], Vec_F[3] 
  r=copy.copy(Vec_F)
  #r = r* (1.00/r.norm()) 
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()
  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3)
   for i in xrange(j+1):
    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 
  D_eig=eig_np(h)
  Lambda, index, Lambda_comp =find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'r0', Lambda, Lambda_comp
  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  Vec_FL=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m+=5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m+=5
     E1=Lambda
  p+=1


 E1L=copy.copy(E1)
 Vec_FL=Vec_FL*(1.00/Vec_FL.norm())
 m=10
 W=2
 num=0
 E1=0
 p=0
 Vec_F.randomize()
 Vec_F=Vec_F*(1.00/Vec_F.norm())
 while p  <  (W+1):
  #print "norm", p, Vec_F.norm(),Vec_F[0], Vec_F[1], Vec_F[2], Vec_F[3] 
  r=copy.copy(Vec_F)
 
  Vec_FL_trans=copy.copy(Vec_FL)
  Vec_FL_trans.transpose()
  dot_vec=Vec_FL_trans*r
  dot_val=dot_vec.trace()
  r=r+(-1.00*dot_val*Vec_FL)
  
  
  #r = r* (1.00/r.norm()) 
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()
  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3)

   Vec_FL_trans=copy.copy(Vec_FL)
   Vec_FL_trans.transpose()
   dot_vec=Vec_FL_trans*r
   dot_val=dot_vec.trace()
   r=r+(-1.00*dot_val*Vec_FL)



   for i in xrange(j+1):
    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 
  D_eig=eig_np(h)
  Lambda, index, Lambda_comp=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'r1', Lambda, Lambda_comp
  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m+=5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m+=5
     E1=Lambda
  p+=1

 Length=abs(E1/E1L)
 Length_val=-2.0*(1.00/math.log(Length))
 print "Length", Length,Length_val 
 return Length_val


def ED_up(c1,Tb1, Ta1,c2, a, b, c, d, Tb2, Ta2, Ta4, Tb4,Vec_uni):


 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D=",  D

 m=10
 W=2
 num=0
 E1=0
 p=0

 while p  <  (W+1):
  r=copy.copy(Vec_F)
  #r = r* (1.00/r.norm()) 
  
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()

  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_u(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4)
   #r.resize(D,1)
   for i in xrange(j+1):

    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 

  D_eig=eig_np(h)
  Lambda, index, Lambda_comp=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'u0', Lambda,Lambda_comp

  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  Vec_FL=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m=m+5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m=m+5
     E1=Lambda
  p+=1

 E1L=copy.copy(E1)
 Vec_FL=Vec_FL*(1.00/Vec_FL.norm())
 m=10
 W=2
 num=0
 E1=0
 p=0
 Vec_F.randomize()
 Vec_F=Vec_F*(1.00/Vec_F.norm())
 while p  <  (W+1):
  #print "norm", p, Vec_F.norm(),Vec_F[0], Vec_F[1], Vec_F[2], Vec_F[3] 
  r=copy.copy(Vec_F)
 
  Vec_FL_trans=copy.copy(Vec_FL)
  Vec_FL_trans.transpose()
  dot_vec=Vec_FL_trans*r
  dot_val=dot_vec.trace()
  r=r+(-1.00*dot_val*Vec_FL)
  
  
  #r = r* (1.00/r.norm()) 
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()
  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_u(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4)

   Vec_FL_trans=copy.copy(Vec_FL)
   Vec_FL_trans.transpose()
   dot_vec=Vec_FL_trans*r
   dot_val=dot_vec.trace()
   r=r+(-1.00*dot_val*Vec_FL)



   for i in xrange(j+1):
    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 
  D_eig=eig_np(h)
  Lambda, index,Lambda_comp=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'u1', Lambda,Lambda_comp
  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m+=5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m+=5
     E1=Lambda
  p+=1



 Length=abs(E1/E1L)
 Length_val=-2.0*(1.00/math.log(Length))
 print "Length", Length,Length_val 
 return Length_val

def Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3):
 CTM_1 = uni10.Network("Network/Right.net")
 CTM_1.putTensor('Vec_uni',Vec_uni)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 vec.permute([21,14,-14 , 7,-7,0],6)
 #print vec.printDiagram() 
 Vec_M=vec.getBlock()
 return Vec_M



def Multi_u(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4):
 CTM_1 = uni10.Network("Network/Up.net")
 CTM_1.putTensor('Vec_uni',Vec_uni)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 vec.permute([17, 18, -18, 19, -19,  20],6)
 Vec_M=vec.getBlock()
 return Vec_M









############################################################################################




















def  Corr_val_function(Iden,HH,H,H1,vec_right, vec_left,dis_val):

 vec_left.setLabel([10,23, 16, -16 , 9,-9, 2,-10])
 vec_right.setLabel([20,23,16,-16 , 9,-9,2,-20])
 Cor_norm=(vec_left*Iden)*vec_right
 Corr_val=(vec_left*HH)*vec_right
 Corr_val1=(vec_left*H)*vec_right
 Corr_val2=(vec_left*H1)*vec_right

 print  dis_val, Corr_val[0]/Cor_norm[0],Corr_val1[0]/Cor_norm[0],Corr_val2[0]/Cor_norm[0], Cor_norm[0] 

 print dis_val, (Corr_val[0]/Cor_norm[0]),  ((Corr_val1[0]*Corr_val2[0])/(Cor_norm[0]*Cor_norm[0]))  

 val=(Corr_val[0]/Cor_norm[0]) - ((Corr_val1[0]*Corr_val2[0])/(Cor_norm[0]*Cor_norm[0]))  

 return val


def  make_ap_openindex(a_u):
 a_u.setLabel([10,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.transpose()
 a_uc.setLabel([-3,-4,-10,-1,-2])
 result=a_uc*a_u
 result.permute([10,1,-1,2,-2,3,-3,4,-4,-10], 5)
 return result


def make_vleft(Tb4,Ta4,c1,c4):

 Tb4.setLabel([-5,3,-3,4])
 Ta4.setLabel([1,2,-2,-5])
 c1.setLabel([4,6])
 c4.setLabel([1,5])
 vec_left=(Tb4*Ta4)*(c1*c4)
 vec_left.permute([5,2,-2,3,-3,6],0)
 return vec_left

def  make_vright(Ta2,Tb2,c2,c3):

 Ta2.setLabel([-5,3,-3,4])
 Tb2.setLabel([1,2,-2,-5])
 c2.setLabel([6,4])
 c3.setLabel([5,1])
 vec_right=(Ta2*Tb2)*(c2*c3)
 vec_right.permute([5,2,-2,3,-3,6],6)
 return vec_right


def  Make_first_vecleft_a(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d):

 CTM_1 = uni10.Network("Network/LeftCorra.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_a(vec_right, Ta1, Ta3,Tb1,Tb3,ap,b,c,d):

 CTM_1 = uni10.Network("Network/RightCorra.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec

def  Make_first_vecleft_b(vec_left, Ta1, Ta3,Tb1,Tb3,a,bp,c,d):

 CTM_1 = uni10.Network("Network/LeftCorrb.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_b(vec_right, Ta1, Ta3,Tb1,Tb3,a,bp,c,d):

 CTM_1 = uni10.Network("Network/RightCorrb.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec


def  Make_first_vecleft_c(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,cp,d):

 CTM_1 = uni10.Network("Network/LeftCorrc.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_c(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,cp,d):

 CTM_1 = uni10.Network("Network/RightCorrc.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec


def  Make_first_vecleft_d(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,c,dp):

 CTM_1 = uni10.Network("Network/LeftCorrd.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_d(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,c,dp):

 CTM_1 = uni10.Network("Network/RightCorrd.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec

def  Make_midle_vecleft(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d):

 CTM_1 = uni10.Network("Network/LeftCorr1.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec



############################################################################################

def  make_ap_openindex(a_u):
 a_u.setLabel([10,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.transpose()
 a_uc.setLabel([-3,-4,-10,-1,-2])
 result=a_uc*a_u
 result.permute([10,1,-1,2,-2,3,-3,4,-4,-10], 5)
 return result


def make_down(c4,Ta3, Tb3,c3):
 Tb3.setLabel([-5,3,-3,4])
 Ta3.setLabel([1,2,-2,-5])
 c3.setLabel([4,6])
 c4.setLabel([5,1])
 vec_down=(Tb3*Ta3)*(c3*c4)
 vec_down.permute([5, 2, -2 , 3, -3 , 6],0)
 return vec_down

def  make_up(c1,Tb1, Ta1,c2):
 Ta1.setLabel([-5,3,-3,4])
 Tb1.setLabel([1,2,-2,-5])
 c2.setLabel([4,6])
 c1.setLabel([5,1])
 vec_up=(Ta1*Tb1)*(c1*c2)
 vec_up.permute([5,2,-2,3,-3,6],6)
 return vec_up




def  Make_first_vecdown_a(vec_down, ap, b, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/DownCorra.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec




def Make_first_vecup_a(vec_up, ap, b, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/UpCorra.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 return vec


def  Make_midle_vecdown(vec_down, a,b,c,d,Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/DownCorr1.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec=max_ten(vec)
 return vec



def  Make_first_vecdown_b(vec_down, a, bp, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/DownCorrb.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec




def Make_first_vecup_b(vec_up, a, bp, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/UpCorrb.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 return vec



def  Make_first_vecdown_c(vec_down, a, b, cp, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/DownCorrc.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec

def Make_first_vecup_c(vec_up, a, b, cp, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/UpCorrc.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 return vec




def  Make_first_vecdown_d(vec_down, a, b, c, dp, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/DownCorrd.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec

def Make_first_vecup_d(vec_up, a, b, c, dp, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network/UpCorrd.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 return vec















#################################################################################################





















def Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8):
 Gamma_a.save("Store/Gamma_a")
 Gamma_b.save("Store/Gamma_b")
 Gamma_c.save("Store/Gamma_c")
 Gamma_d.save("Store/Gamma_d")
 Landa_1.save("Store/Landa_1")
 Landa_2.save("Store/Landa_2")
 Landa_3.save("Store/Landa_3")
 Landa_4.save("Store/Landa_4")
 Landa_5.save("Store/Landa_5")
 Landa_6.save("Store/Landa_6")
 Landa_7.save("Store/Landa_7")
 Landa_8.save("Store/Landa_8")

def Reload_itebd():
 Gamma_a=uni10.UniTensor("Store/Gamma_a")
 Gamma_b=uni10.UniTensor("Store/Gamma_b")
 Gamma_c=uni10.UniTensor("Store/Gamma_c")
 Gamma_d=uni10.UniTensor("Store/Gamma_d")
 Landa_1=uni10.UniTensor("Store/Landa_1")
 Landa_2=uni10.UniTensor("Store/Landa_2")
 Landa_3=uni10.UniTensor("Store/Landa_3")
 Landa_4=uni10.UniTensor("Store/Landa_4")
 Landa_5=uni10.UniTensor("Store/Landa_5")
 Landa_6=uni10.UniTensor("Store/Landa_6")
 Landa_7=uni10.UniTensor("Store/Landa_7")
 Landa_8=uni10.UniTensor("Store/Landa_8")
 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8

def Store_Full(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.save("Store/a_u")
 b_u.save("Store/b_u")
 c_u.save("Store/c_u")
 d_u.save("Store/d_u")
 a.save("Store/a")
 b.save("Store/b")
 c.save("Store/c")
 d.save("Store/d")


def increase_norm(a_u,b_u,c_u,d_u,a,b,c,d, N):
 a_u=a_u*N
 b_u=b_u*N
 c_u=c_u*N
 d_u=d_u*N
 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)
 return a_u,b_u,c_u,d_u,a,b,c,d




def Store_Fullp(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.save("Store/ap_u")
 b_u.save("Store/bp_u")
 c_u.save("Store/cp_u")
 d_u.save("Store/dp_u")
 a.save("Store/ap")
 b.save("Store/bp")
 c.save("Store/cp")
 d.save("Store/dp")

def Reload_Fullp():
 ap_u=uni10.UniTensor("Store/ap_u")
 bp_u=uni10.UniTensor("Store/bp_u")
 cp_u=uni10.UniTensor("Store/cp_u")
 dp_u=uni10.UniTensor("Store/dp_u")
 ap=uni10.UniTensor("Store/ap")
 bp=uni10.UniTensor("Store/bp")
 cp=uni10.UniTensor("Store/cp")
 dp=uni10.UniTensor("Store/dp")
 return ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp


def Reload_Full():
 a_u=uni10.UniTensor("Store/a_u")
 b_u=uni10.UniTensor("Store/b_u")
 c_u=uni10.UniTensor("Store/c_u")
 d_u=uni10.UniTensor("Store/d_u")
 a=uni10.UniTensor("Store/a")
 b=uni10.UniTensor("Store/b")
 c=uni10.UniTensor("Store/c")
 d=uni10.UniTensor("Store/d")
 return a_u,b_u,c_u,d_u,a,b,c,d

def reconstruct_ab(a_u, a_uf):

 blk_qnums = a_u.blockQnum()
 blk_qnumsf = a_uf.blockQnum()

 for qnum in blk_qnums:
  if qnum in blk_qnumsf:
   mat_t=a_u.getBlock(qnum)
   dx=int(mat_t.row())
   dy=int(mat_t.col())
   mat_f=a_uf.getBlock(qnum)
   a_u.putBlock(qnum,mat_f.resize(dx,dy) )
 return a_u


def Reload_Full_previous(a_u, b_u, c_u, d_u):
 a_uf=uni10.UniTensor("Store/a_u")
 b_uf=uni10.UniTensor("Store/b_u")
 c_uf=uni10.UniTensor("Store/c_u")
 d_uf=uni10.UniTensor("Store/d_u")

 #print "a_uf",a_uf
 a_u=reconstruct_ab(a_u, a_uf)
 b_u=reconstruct_ab(b_u, b_uf)
 c_u=reconstruct_ab(c_u, c_uf)
 d_u=reconstruct_ab(d_u, d_uf)

 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)
# print b_uf
# print b_u

# print c_uf
# print c_u

# print d_uf
# print d_u

 return a_u,b_u,c_u,d_u,a,b,c,d



def Rand_env_total(Env):
 Env1=copy.copy(Env)
 for i in xrange(len(Env1)):
  Env1[i]=copy.copy(Env[i])
  Env1[i].randomize()
 return  Env1


def qr_parity(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
    
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([bd1,bd2, theta.bond(4),theta.bond(5)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]

    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).qr()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])

#    print LA
    return GA, LA

def lq_parity(theta):

    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())
    bd2=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([theta.bond(0),theta.bond(1),bd1,bd2])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA


def Decomposition(U):

 H=copy.copy(U)
 H.setLabel([0,1,2,3,4,5])
 H.permute([0,3,4,5,2,1],2)
 
 l,q=lq_parity(H)
 l.setLabel([0,3,-1,-2])

 q.setLabel([-1,-2,4,5,2,1])
 q.permute([-1,-2,4,1,5,2],4)
 qq, r=qr_parity(q)
 qq.setLabel([-1,-2,4,1,-3,-4])
 r.setLabel([-3,-4,5,2]) 
 
 l.permute([0,-1,-2,3],1)
 qq.permute([-1,-2,1,-4,-3,4],3)
 r.permute([-3,-4,2,5],3)
 
 MPO_list=[]
 MPO_list.append(l)
 MPO_list.append(qq)
 MPO_list.append(r)

 #H1=MPO_list[0]*MPO_list[1]*MPO_list[2]

 #H1.permute([0,1,2,3,4,5],3)
 #M=H1.getBlock()
 #D=M.eigh()
 #print D[0]
 #print "Test", H1.elemCmp(U), MPO_list[0],MPO_list[1],MPO_list[2]
 #print "hiihi",H1, D[0][0]
 return MPO_list

def slighty_random(a_u,b_u,c_u,d_u,a,b,c,d):
 rand=copy.copy(a_u)
 rand.randomize()
 a_u=a_u+(0.001)*rand

 rand=copy.copy(b_u)
 rand.randomize()
 b_u=b_u+(0.001)*rand
 
 rand=copy.copy(c_u)
 rand.randomize()
 c_u=c_u+(0.001)*rand
 
 rand=copy.copy(d_u)
 rand.randomize()
 d_u=d_u+(0.001)*rand


 a_u*=(1.00/MaxAbs(a_u)) 
 b_u*=(1.00/MaxAbs(b_u)) 
 c_u*=(1.00/MaxAbs(c_u)) 
 d_u*=(1.00/MaxAbs(d_u)) 
 
 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)
 
 return a_u,b_u,c_u,d_u,a,b,c,d

def choose_model(Model, h, d_phys):
 if Model is "Heisenberg":
   H0=Heisenberg0(h,d_phys)
   H00=Heisenberg00(h,d_phys)
   H1=Heisenberg1(h,d_phys)
   H2=threebody(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h,d_phys)
   H00=Heisenberg00_Z2(h,d_phys)
   H1=Heisenberg1_Z2(h,d_phys)
   H2=threebody_Z2(h,d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h,d_phys)
   H00=Heisenberg0_U1(h,d_phys)
   H1=H00
   H2=threebody_U1(h,d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h,d_phys)
   H00=Heisenberg0_U1Z2(h,d_phys)

 return H0, H00, H1, H2


def total_random(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.randomize()
 b_u.randomize()
 c_u.randomize()
 d_u.randomize()


 a_u=a_u*(1.00/MaxAbs(a_u)) 
 b_u=b_u*(1.00/MaxAbs(b_u)) 
 c_u=c_u*(1.00/MaxAbs(c_u)) 
 d_u=d_u*(1.00/MaxAbs(d_u)) 
 
 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)

# a=a*(1.00/MaxAbs(a)) 
# b=b*(1.00/MaxAbs(b)) 
# c=c*(1.00/MaxAbs(c)) 
# d=d*(1.00/MaxAbs(d)) 
 
 return a_u,b_u,c_u,d_u,a,b,c,d


def total_random1(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.orthoRand()
 b_u.orthoRand()
 c_u.orthoRand()
 d_u.orthoRand()


 a_u=a_u*(1.00/MaxAbs(a_u)) 
 b_u=b_u*(1.00/MaxAbs(b_u)) 
 c_u=c_u*(1.00/MaxAbs(c_u)) 
 d_u=d_u*(1.00/MaxAbs(d_u)) 
 
 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)

# a=a*(1.00/MaxAbs(a)) 
# b=b*(1.00/MaxAbs(b)) 
# c=c*(1.00/MaxAbs(c)) 
# d=d*(1.00/MaxAbs(d)) 
 
 return a_u,b_u,c_u,d_u,a,b,c,d



def max_ten(a):

 if ( MaxAbs(a) < 0.50e-1) or (MaxAbs(a) > 0.50e+1)   :
  a=a*(1.00/MaxAbs(a));
 else: a=a;
 
 return a

def Init_env(Env):
 c1=copy.copy(Env[0])
 c2=copy.copy(Env[1])
 c3=copy.copy(Env[2]) 
 c4=copy.copy(Env[3]) 
 Ta1=copy.copy(Env[4])
 Ta2=copy.copy(Env[5])
 Ta3=copy.copy(Env[6])
 Ta4=copy.copy(Env[7])
 Tb1=copy.copy(Env[8])
 Tb2=copy.copy(Env[9])
 Tb3=copy.copy(Env[10])
 Tb4=copy.copy(Env[11])
 return  c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4


def positve_Env(Env,a,b,c,d):
 c1=copy.copy(Env[0])
 c2=copy.copy(Env[1])
 c3=copy.copy(Env[2]) 
 c4=copy.copy(Env[3]) 
 Ta1=copy.copy(Env[4])
 Ta2=copy.copy(Env[5])
 Ta3=copy.copy(Env[6])
 Ta4=copy.copy(Env[7])
 Tb1=copy.copy(Env[8])
 Tb2=copy.copy(Env[9])
 Tb3=copy.copy(Env[10])
 Tb4=copy.copy(Env[11])

 a1=copy.copy(a)
 a1.setLabel([1,-1,2,-2,3,-3,4,-4])
 a1.partialTrace(1,-1)
 a1.partialTrace(4,-4)
 a1.combineBond([2,-2])
 a1.combineBond([3,-3])
 a1.permute([2,3],1)
 Env[0]=a1


 b1=copy.copy(b)
 b1.setLabel([1,-1,2,-2,3,-3,4,-4])
 b1.partialTrace(3,-3)
 b1.partialTrace(4,-4)
 b1.combineBond([1,-1])
 b1.combineBond([2,-2])
 b1.permute([1,2],2)
 Env[1]=b1


 c1=copy.copy(c)
 c1.setLabel([1,-1,2,-2,3,-3,4,-4])
 c1.partialTrace(1,-1)
 c1.partialTrace(2,-2)
 c1.combineBond([3,-3])
 c1.combineBond([4,-4])
 c1.permute([3,4],0)
 Env[2]=c1


 d1=copy.copy(d)
 d1.setLabel([1,-1,2,-2,3,-3,4,-4])
 d1.partialTrace(2,-2)
 d1.partialTrace(3,-3)
 d1.combineBond([1,-1])
 d1.combineBond([4,-4])
 d1.permute([1,4],1)
 Env[3]=d1


 Ta1=copy.copy(Env[4])
 Ta2=copy.copy(Env[5])
 Ta3=copy.copy(Env[6])
 Ta4=copy.copy(Env[7])
 Tb1=copy.copy(Env[8])
 Tb2=copy.copy(Env[9])
 Tb3=copy.copy(Env[10])
 Tb4=copy.copy(Env[11])

 a1=copy.copy(a)
 a1.setLabel([1,-1,2,-2,3,-3,4,-4])
 a1.partialTrace(4,-4)
 a1.combineBond([1,-1])
 a1.combineBond([3,-3])
 a1.permute([1,2,-2,3],3)
 Env[8]=a1
 Env[9]=copy.copy(a1)
 Env[4]=copy.copy(a1)
 Env[5]=copy.copy(a1)

 a1=copy.copy(a)
 a1.setLabel([1,-1,2,-2,3,-3,4,-4])
 a1.partialTrace(1,-1)
 a1.combineBond([2,-2])
 a1.combineBond([4,-4])
 a1.permute([2,3,-3,4],1)
 Env[6]=a1
 Env[7]=a1
 Env[10]=a1
 Env[11]=a1
 return Env


def Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4):
 c1.save("Store/c1")
 c2.save("Store/c2")
 c3.save("Store/c3")
 c4.save("Store/c4")
 Ta1.save("Store/Ta1")
 Ta2.save("Store/Ta2")
 Ta3.save("Store/Ta3")
 Ta4.save("Store/Ta4")
 Tb1.save("Store/Tb1")
 Tb2.save("Store/Tb2")
 Tb3.save("Store/Tb3")
 Tb4.save("Store/Tb4")

def Store_EnvEnv(Env,Env1,Env2,Env3):

 Env[0].save("Store/c11")
 Env[1].save("Store/c21")
 Env[2].save("Store/c31")
 Env[3].save("Store/c41")
 Env[4].save("Store/Ta11")
 Env[5].save("Store/Ta21")
 Env[6].save("Store/Ta31")
 Env[7].save("Store/Ta41")
 Env[8].save("Store/Tb11")
 Env[9].save("Store/Tb21")
 Env[10].save("Store/Tb31")
 Env[11].save("Store/Tb41")

 Env1[0].save("Store/c12")
 Env1[1].save("Store/c22")
 Env1[2].save("Store/c32")
 Env1[3].save("Store/c42")
 Env1[4].save("Store/Ta12")
 Env1[5].save("Store/Ta22")
 Env1[6].save("Store/Ta32")
 Env1[7].save("Store/Ta42")
 Env1[8].save("Store/Tb12")
 Env1[9].save("Store/Tb22")
 Env1[10].save("Store/Tb32")
 Env1[11].save("Store/Tb42")

 Env2[0].save("Store/c13")
 Env2[1].save("Store/c23")
 Env2[2].save("Store/c33")
 Env2[3].save("Store/c43")
 Env2[4].save("Store/Ta13")
 Env2[5].save("Store/Ta23")
 Env2[6].save("Store/Ta33")
 Env2[7].save("Store/Ta43")
 Env2[8].save("Store/Tb13")
 Env2[9].save("Store/Tb23")
 Env2[10].save("Store/Tb33")
 Env2[11].save("Store/Tb43")

 Env3[0].save("Store/c14")
 Env3[1].save("Store/c24")
 Env3[2].save("Store/c34")
 Env3[3].save("Store/c44")
 Env3[4].save("Store/Ta14")
 Env3[5].save("Store/Ta24")
 Env3[6].save("Store/Ta34")
 Env3[7].save("Store/Ta44")
 Env3[8].save("Store/Tb14")
 Env3[9].save("Store/Tb24")
 Env3[10].save("Store/Tb34")
 Env3[11].save("Store/Tb44")

def Reload_EnvEnv(Env,Env1,Env2,Env3):

 Env[0]=uni10.UniTensor("Store/c11")
 Env[1]=uni10.UniTensor("Store/c21")
 Env[2]=uni10.UniTensor("Store/c31")
 Env[3]=uni10.UniTensor("Store/c41")
 Env[4]=uni10.UniTensor("Store/Ta11")
 Env[5]=uni10.UniTensor("Store/Ta21")
 Env[6]=uni10.UniTensor("Store/Ta31")
 Env[7]=uni10.UniTensor("Store/Ta41")
 Env[8]=uni10.UniTensor("Store/Tb11")
 Env[9]=uni10.UniTensor("Store/Tb21")
 Env[10]=uni10.UniTensor("Store/Tb31")
 Env[11]=uni10.UniTensor("Store/Tb41")

 Env1[0]=uni10.UniTensor("Store/c12")
 Env1[1]=uni10.UniTensor("Store/c22")
 Env1[2]=uni10.UniTensor("Store/c32")
 Env1[3]=uni10.UniTensor("Store/c42")
 Env1[4]=uni10.UniTensor("Store/Ta12")
 Env1[5]=uni10.UniTensor("Store/Ta22")
 Env1[6]=uni10.UniTensor("Store/Ta32")
 Env1[7]=uni10.UniTensor("Store/Ta42")
 Env1[8]=uni10.UniTensor("Store/Tb12")
 Env1[9]=uni10.UniTensor("Store/Tb22")
 Env1[10]=uni10.UniTensor("Store/Tb32")
 Env1[11]=uni10.UniTensor("Store/Tb42")

 Env2[0]=uni10.UniTensor("Store/c13")
 Env2[1]=uni10.UniTensor("Store/c23")
 Env2[2]=uni10.UniTensor("Store/c33")
 Env2[3]=uni10.UniTensor("Store/c43")
 Env2[4]=uni10.UniTensor("Store/Ta13")
 Env2[5]=uni10.UniTensor("Store/Ta23")
 Env2[6]=uni10.UniTensor("Store/Ta33")
 Env2[7]=uni10.UniTensor("Store/Ta43")
 Env2[8]=uni10.UniTensor("Store/Tb13")
 Env2[9]=uni10.UniTensor("Store/Tb23")
 Env2[10]=uni10.UniTensor("Store/Tb33")
 Env2[11]=uni10.UniTensor("Store/Tb43")

 Env3[0]=uni10.UniTensor("Store/c14")
 Env3[1]=uni10.UniTensor("Store/c24")
 Env3[2]=uni10.UniTensor("Store/c34")
 Env3[3]=uni10.UniTensor("Store/c44")
 Env3[4]=uni10.UniTensor("Store/Ta14")
 Env3[5]=uni10.UniTensor("Store/Ta24")
 Env3[6]=uni10.UniTensor("Store/Ta34")
 Env3[7]=uni10.UniTensor("Store/Ta44")
 Env3[8]=uni10.UniTensor("Store/Tb14")
 Env3[9]=uni10.UniTensor("Store/Tb24")
 Env3[10]=uni10.UniTensor("Store/Tb34")
 Env3[11]=uni10.UniTensor("Store/Tb44")


def Reload_Env():
 c1=uni10.UniTensor("Store/c1")
 c2=uni10.UniTensor("Store/c2")
 c3=uni10.UniTensor("Store/c3")
 c4=uni10.UniTensor("Store/c4")

 Ta1=uni10.UniTensor("Store/Ta1")
 Ta2=uni10.UniTensor("Store/Ta2")
 Ta3=uni10.UniTensor("Store/Ta3")
 Ta4=uni10.UniTensor("Store/Ta4")

 Tb1=uni10.UniTensor("Store/Tb1")
 Tb2=uni10.UniTensor("Store/Tb2")
 Tb3=uni10.UniTensor("Store/Tb3")
 Tb4=uni10.UniTensor("Store/Tb4")
 return c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4





def reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env):
  
  Env[0]=copy.copy(c1)
  Env[1]=copy.copy(c2)
  Env[2]=copy.copy(c3) 
  Env[3]=copy.copy(c4) 
  Env[4]=copy.copy(Ta1)
  Env[5]=copy.copy(Ta2)
  Env[6]=copy.copy(Ta3)
  Env[7]=copy.copy(Ta4)
  Env[8]=copy.copy(Tb1)
  Env[9]=copy.copy(Tb2)
  Env[10]=copy.copy(Tb3)
  Env[11]=copy.copy(Tb4)
  
def Rand_env_slight(Env):
 
 for i in xrange(len(Env)):
  Env_tem=copy.copy(Env[i]) 
  Env_tem.randomize()
  Env[i]=Env[i]+0.01*Env_tem

def Rand_env_total(Env):
 Env1=copy.copy(Env)
 for i in xrange(len(Env1)):
  Env1[i]=copy.copy(Env[i])
  Env1[i].randomize()
 
 return  Env1

def make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1):

 q_list=[]
 q_list1=[]
 q_list2=[]

 degs = c2.bond(0).degeneracy()
 for qnum, dim in degs.iteritems():
  if qnum.U1() < 0:
   for i in xrange(dim):
    q_list.append(qnum)
    q_list1.append(-qnum)
  if qnum.U1() == 0:
   for i in xrange(dim):
    q_list2.append(qnum)

 q_list1.reverse()
 q_list=q_list+q_list2+q_list1
 bd=uni10.Bond(uni10.BD_IN,q_list)


 blk_qnums = c2.blockQnum()
 bd_list=[bd,c2.bond(1)]
 c2b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c2b.getBlock(qnum).row())
  dy=int(c2b.getBlock(qnum).col())
  c2b.putBlock(qnum,c2.getBlock(qnum).resize(dx,dy) )





 bd=uni10.Bond(uni10.BD_OUT,c2b.bond(0).Qlist())
 blk_qnums = Ta1.blockQnum()
 bd_list=[Ta1.bond(0),Ta1.bond(1),Ta1.bond(2), bd]
 Ta1b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Ta1b.getBlock(qnum).row())
  dy=int(Ta1b.getBlock(qnum).col())
  Ta1b.putBlock(qnum,Ta1.getBlock(qnum).resize(dx,dy) )




 q_list=[]
 q_list1=[]
 q_list2=[]

 
 degs = c1.bond(1).degeneracy()
 for qnum, dim in degs.iteritems():
  if qnum.U1() < 0:
   for i in xrange(dim):
    q_list.append(qnum)
    q_list1.append(-qnum)
  if qnum.U1() == 0:
   for i in xrange(dim):
    q_list2.append(qnum)

 q_list1.reverse()
 q_list=q_list+q_list2+q_list1
 bd=uni10.Bond(uni10.BD_OUT,q_list)
 
 
 blk_qnums = c1.blockQnum()
 bd_list=[c1.bond(0),bd]
 c1b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c1b.getBlock(qnum).row())
  dy=int(c1b.getBlock(qnum).col())
  c1b.putBlock(qnum,c1.getBlock(qnum).resize(dx,dy) )

 
 
 bd=uni10.Bond(uni10.BD_IN,c1b.bond(1).Qlist())
 blk_qnums = Tb1.blockQnum()
 bd_list=[bd,Tb1.bond(1),Tb1.bond(2), Tb1.bond(3)]
 Tb1b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Tb1b.getBlock(qnum).row())
  dy=int(Tb1b.getBlock(qnum).col())
  Tb1b.putBlock(qnum,Tb1.getBlock(qnum).resize(dx,dy) )
  
 
 blk_qnums = c3.blockQnum()
 bd_list=[c2b.bond(0),c3.bond(1)]
 c3b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c3b.getBlock(qnum).row())
  dy=int(c3b.getBlock(qnum).col())
  c3b.putBlock(qnum,c3.getBlock(qnum).resize(dx,dy) )
 
 blk_qnums = c4.blockQnum()
 bd_list=[c4.bond(0),c1b.bond(1)]
 c4b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c4b.getBlock(qnum).row())
  dy=int(c4b.getBlock(qnum).col())
  c4b.putBlock(qnum,c4.getBlock(qnum).resize(dx,dy) )



 bd=uni10.Bond(uni10.BD_OUT,c3b.bond(0).Qlist())
 blk_qnums = Tb3.blockQnum()
 bd_list=[Tb3.bond(0),Tb3.bond(1),Tb3.bond(2), bd]
 Tb3b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Tb3b.getBlock(qnum).row())
  dy=int(Tb3b.getBlock(qnum).col())
  Tb3b.putBlock(qnum,Tb3.getBlock(qnum).resize(dx,dy) )

 bd=uni10.Bond(uni10.BD_IN,c4b.bond(1).Qlist())
 blk_qnums = Ta3.blockQnum()
 bd_list=[bd,Ta3.bond(1),Ta3.bond(2), Ta3.bond(3)]
 Ta3b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Ta3b.getBlock(qnum).row())
  dy=int(Ta3b.getBlock(qnum).col())
  Ta3b.putBlock(qnum,Ta3.getBlock(qnum).resize(dx,dy) )



 return c1b, c2b, c3b, c4b, Tb3b, Ta3b, Ta1b, Tb1b

def rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):

 bd=uni10.Bond(uni10.BD_OUT,a.bond(0).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,a.bond(1).Qlist())

 tempo=copy.copy(Tb4)
 bd_list=[Tb4.bond(0),bd,bd1,Tb4.bond(3)]
 Tb4.assign(bd_list)
 blk_qnums = Tb4.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb4.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb4 dimension"      
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb4.putBlock(qnum,sv_mat.resize(dx,dy) )
 
 #print "Tb4",tempo.printDiagram(),Tb4.printDiagram(),a.printDiagram()
 #print "Tb4", tempo[10],Tb4[10],tempo.elemCmp(Tb4)#,Tb4.printDiagram(),
################################################################

 bd=uni10.Bond(uni10.BD_OUT,c.bond(0).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,c.bond(1).Qlist())
 
 tempo=copy.copy(Ta4)
 bd_list=[Ta4.bond(0),bd,bd1,Ta4.bond(3)]
 Ta4.assign(bd_list)
 blk_qnums = Ta4.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta4.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta4 dimension"      
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta4.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Ta4", Ta4.printDiagram(),Ta4[4]
 #print "Ta4", tempo.elemCmp(Ta4),Tb4.printDiagram(),

##################################################################

 bd=uni10.Bond(uni10.BD_IN,b.bond(4).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,b.bond(5).Qlist())

 tempo=copy.copy(Ta2)
 bd_list=[Ta2.bond(0),bd,bd1,Ta2.bond(3)]
 Ta2.assign(bd_list)
 blk_qnums = Ta2.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta2.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta2 dimension"      
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta2.putBlock(qnum,sv_mat.resize(dx,dy) )
# print "Ta2", Ta2.printDiagram(),Ta2[3]
# print "Ta2", tempo.elemCmp(Ta2),Tb4.printDiagram(),



 bd=uni10.Bond(uni10.BD_IN,d.bond(4).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,d.bond(5).Qlist())


 tempo=copy.copy(Tb2)
 bd_list=[Tb2.bond(0),bd,bd1,Tb2.bond(3)]
 Tb2.assign(bd_list)
 blk_qnums = Tb2.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb2.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb2 dimension"    
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb2.putBlock(qnum,sv_mat.resize(dx,dy) )
# print "Tb2", Tb2.printDiagram(),Tb2[4]
# print "Tb2", tempo.elemCmp(Tb2),Tb4.printDiagram(),

################################################################3
 bd=uni10.Bond(uni10.BD_IN,a.bond(6).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,a.bond(7).Qlist())

 tempo=copy.copy(Tb1)
 bd_list=[Tb1.bond(0),bd,bd1,Tb1.bond(3)]
 Tb1.assign(bd_list)
 blk_qnums = Tb1.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb1.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb1 dimension"  
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb1.putBlock(qnum,sv_mat.resize(dx,dy) )
# print "Tb1", Tb1.printDiagram(),Tb1[4]
# print "Tb1", tempo.elemCmp(Tb1),Tb4.printDiagram(),


 bd=uni10.Bond(uni10.BD_IN,b.bond(6).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,b.bond(7).Qlist())

 tempo=copy.copy(Ta1)
 bd_list=[Ta1.bond(0),bd,bd1,Ta1.bond(3)]
 Ta1.assign(bd_list)
 blk_qnums = Ta1.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta1.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta1 dimension"  
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta1.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Ta1",tempo.printDiagram(),Ta1.printDiagram(),b.printDiagram()
 #print "Ta1", tempo[10],Ta1[10],tempo.elemCmp(Ta1)#,Tb4.printDiagram(),

######################################################
 bd=uni10.Bond(uni10.BD_OUT,c.bond(2).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,c.bond(3).Qlist())

 tempo=copy.copy(Ta3)
 bd_list=[Ta3.bond(0),bd,bd1,Ta3.bond(3)]
 Ta3.assign(bd_list)
 blk_qnums = Ta3.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta3.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta3 dimension"  
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta3.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Ta3", Ta3.printDiagram(),Ta3[4]
 #print "Ta3", tempo.elemCmp(Ta3),Tb4.printDiagram(),
 
 
 
 bd=uni10.Bond(uni10.BD_OUT,d.bond(2).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,d.bond(3).Qlist())

 tempo=copy.copy(Tb3)
 bd_list=[Tb3.bond(0),bd,bd1,Tb3.bond(3)]
 Tb3.assign(bd_list)
 blk_qnums = Tb3.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb3.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb3 dimension"
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb3.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Tb3", Tb3.printDiagram(),Tb3[4]
 #print "Tb3", tempo.elemCmp(Tb3),Tb4.printDiagram(),

################################################ 

 return Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4

 
def Short_TrotterSteps1(N_iterF):
 List_delN=[]

 for i in xrange(5, 0, -2):
  Delta_N=(i*(1.0/100),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -2):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -2):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)

 return List_delN 
 
def Long_TrotterSteps(N_iterF):
 List_delN=[]


 for i in xrange(10, 1, -1):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -1):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)

 return List_delN 
 
 
def Long_TrotterSteps1(N_iterF):
 List_delN=[]


 for i in xrange(10, 1, -2):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -2):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)

 return List_delN 
 
#@profile
def Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U, Gauge):
 Gauge=copy.copy(Gauge)
 Gauge='Fixed'
 
 U.setLabel([-54,-55,-56,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,54,55,56])
 U_2=U*H

 U.setLabel([-54,-55,-56,54,55,56])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-54,-55,-56,54,55,56])

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---a---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 E_store=(((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp)))*(((E2*E3)*(b_up*b_dp)))
 A=E_store*((E1*E8)*(a_dp))
 A.permute([55,16,17,18,2],3)
 D_a=A

 if Gauge is "Fixed":
  A=E_store*((E1*E8)*(a_up))
  A.permute([-18,-2,55,-16,-17],2)
  A.transpose()
  D_a=D_a+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


 A=((((((E2*E3)*(b_up*b_d)))*((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_up*c_d)*U))*((E1*E8)*(a_d))
 A.permute([55,16,17,18,2],3)
 D_a=D_a+(-1.0)*A


 if Gauge is "Fixed":
  A=((((((E2*E3)*(b_u*b_dp)))*((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u*c_dp)*U))*((E1*E8)*(a_u))
  A.permute([-18,-2,-55,-16,-17],2)
  A.transpose()
  D_a=D_a+(-1.0)*A
#  D_a.transpose()
#  D_a.permute([55,16,17,18,2],3)


 D_a.transpose()
 D_a.permute([55,16,17,18,2],3)
 
#####################################----b----#############################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 E_store=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))
 A=E_store*(((E2*E3)*(b_dp)))
 A.permute([56,18,20,6,4],3)
 D_b=A

 if Gauge is "Fixed":
  A=E_store*(((E2*E3)*(b_up)))
  A.permute([-6,-4,56,-18,-20],2)
  A.transpose()
  D_b=D_b+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*((E4*E5)*(d_up*d_d)))*(((E2*E3)*(b_d)))
 A.permute([56,18,20,6,4],3)
 D_b=D_b+(-1.0)*A


 if Gauge is "Fixed":
  A=(((((E1*E8)*(a_u*a_dp)*U)*((E7*E6)*(c_u*c_dp))))*((E4*E5)*(d_u*d_dp)))*(((E2*E3)*(b_u)))
  A.permute([-6,-4,-56,-18,-20],2)
  A.transpose()
  D_b=D_b+(-1.0)*A
#  D_b.transpose()
#  D_b.permute([56,18,20,6,4],3)

 D_b.transpose()
 D_b.permute([56,18,20,6,4],3)

##################################---c---#################################################3


 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 E_store=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))

 A=E_store * ((E7*E6)*(c_dp))
 A.permute([54,14,12,19,17],3)
 D_c=A

 if Gauge is "Fixed":
  A=E_store * ((E7*E6)*(c_up))
  A.permute([-19,-17,54,-14,-12],2)
  A.transpose()
  D_c=D_c+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E2*E3)*(b_up*b_d)))*U)*(((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_d))
 A.permute([54,14,12,19,17],3)
 D_c=D_c+(-1.0)*A


 if Gauge is "Fixed":
  A=(((((E1*E8)*(a_u*a_dp))*((E2*E3)*(b_u*b_dp)))*U)*(((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u))
  A.permute([-19,-17,-54,-14,-12],2)
  A.transpose()
  D_c=D_c+(-1.0)*A
#  D_c.transpose()
#  D_c.permute([54,14,12,19,17],3)
 
 D_c.transpose()
 D_c.permute([54,14,12,19,17],3)
 
###################################----d----######################################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 E_store=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))

 A=E_store*((E4*E5)*(d_dp))
 A.permute([57,19,10,8,20],3)
 D_d=A

 if Gauge is "Fixed":
  A=E_store*((E4*E5)*(d_up))
  A.permute([-8,-20,57,-19,-10],2)
  A.transpose()
  D_d=D_d+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_d))
 A.permute([57,19,10,8,20],3)
 D_d=D_d+(-1.0)*A


 if Gauge is "Fixed":
  A=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u))
  A.permute([-8,-20,57,-19,-10],2)
  A.transpose()
  D_d=D_d+(-1.0)*A
#  D_d.transpose()
#  D_d.permute([57,19,10,8,20],3)
 # D_d=copy.copy(d_u)
 # D_d.set_zero()
 D_d.transpose()
 D_d.permute([57,19,10,8,20],3)
 
 
##################################################################################################
 return D_a, D_b, D_c, D_d 
 
#@profile 
def Obtain_grad_four1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Gauge):

 Gauge=copy.copy(Gauge)
 Gauge='Fixed'

 U.setLabel([-55,-56,-57,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,55,56,57])
 U_2=U*H

 U.setLabel([-55,-56,-57,55,56,57])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-55,-56,-57,55,56,57])

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---a---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


 E_store=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))


 A=E_store*((E1*E8)*(a_dp))
 A.permute([55,16,17,18,2],3)
 D_a=A


 if Gauge is "Fixed":
  A=E_store*((E1*E8)*(a_up))
  A.permute([-18,-2,55,-16,-17],2)
  A.transpose()
  D_a=D_a+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=((((((E2*E3)*(b_up*b_d))*U)*((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_up*c_d)))*((E1*E8)*(a_d))
 A.permute([55,16,17,18,2],3)
 D_a=D_a+(-1.0)*A

 if Gauge is "Fixed":
  A=((((((E2*E3)*(b_u*b_dp))*U)*((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u*c_dp)))*((E1*E8)*(a_u))
  A.permute([-18,-2,-55,-16,-17],2)
  A.transpose()
  D_a=D_a+(-1.0)*A
#  D_a.transpose()
#  D_a.permute([55,16,17,18,2],3)
 
 D_a.transpose()
 D_a.permute([55,16,17,18,2],3)
 
#####################################----b----#############################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 E_store=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))

 A=E_store*(((E2*E3)*(b_dp)))
 A.permute([56,18,20,6,4],3)
 D_b=A

 if Gauge is "Fixed":
  A=E_store*(((E2*E3)*(b_up)))
  A.permute([-6,-4,56,-18,-20],2)
  A.transpose()
  D_b=D_b+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d))))*((E4*E5)*(d_up*d_d))*U)*(((E2*E3)*(b_d)))
 A.permute([56,18,20,6,4],3)
 D_b=D_b+(-1.0)*A

 if Gauge is "Fixed":
  A=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp))))*((E4*E5)*(d_u*d_dp))*U)*(((E2*E3)*(b_u)))
  A.permute([-6,-4,-56,-18,-20],2)
  A.transpose()
  D_b=D_b+(-1.0)*A
#  D_b.transpose()
#  D_b.permute([56,18,20,6,4],3)

 D_b.transpose()
 D_b.permute([56,18,20,6,4],3)

##################################---c---#################################################3


 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 E_store=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))

 A=E_store * ((E7*E6)*(c_dp))
 A.permute([54,14,12,19,17],3)
 D_c=A

 if Gauge is "Fixed":
  A=E_store * ((E7*E6)*(c_up))
  A.permute([-19,-17,54,-14,-12],2)
  A.transpose()
  D_c=D_c+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E2*E3)*(b_up*b_d)))*U)*(((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_d))
 A.permute([54,14,12,19,17],3)
 D_c=D_c+(-1.0)*A

 if Gauge is "Fixed":
  A=(((((E1*E8)*(a_u*a_dp))*((E2*E3)*(b_u*b_dp)))*U)*(((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u))
  A.permute([-19,-17,54,-14,-12],2)
  A.transpose()
  D_c=D_c+(-1.0)*A
#  D_c.transpose()
#  D_c.permute([54,14,12,19,17],3)
 
 D_c.transpose()
 D_c.permute([54,14,12,19,17],3)
 
###################################----d----######################################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])



 E_store=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))

 A=E_store*((E4*E5)*(d_dp))
 A.permute([57,19,10,8,20],3)
 D_d=A

 if Gauge is "Fixed":
  A=E_store *((E4*E5)*(d_up))
  A.permute([-8,-20,57,-19,-10],2)
  A.transpose()
  D_d=D_d+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d))))*(((E2*E3)*(b_up*b_d))*U))*((E4*E5)*(d_d))
 A.permute([57,19,10,8,20],3)
 D_d=D_d+(-1.0)*A

 if Gauge is "Fixed":
  A=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp))))*(((E2*E3)*(b_u*b_dp))*U))*((E4*E5)*(d_u))
  A.permute([-8,-20,-57,-19,-10],2)
  A.transpose()
  D_d=D_d+(-1.0)*A
#  D_d.transpose()
#  D_d.permute([57,19,10,8,20],3)

 D_d.transpose()
 D_d.permute([57,19,10,8,20],3)

 
##################################################################################################
 return D_a, D_b, D_c, D_d 



def Obtain_grad_four_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):

 D_r=[0]*4
 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)

 a_d1=copy.copy(a_u)
 b_d1=copy.copy(b_u)
 c_d1=copy.copy(c_u)
 
###########################################################

 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])

 c_d1.setLabel([51,-14,-12,-19,-17])
 a_d1.setLabel([52,-16,-17,-18,-2])
 b_d1.setLabel([53,-18,-20,-6,-4])

 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])
##########################################################
 a_u.setLabel([55,16,64,66,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([52,-16,-64,-66,-2])

 b_u.setLabel([56,68,20,6,4])
 b_d=copy.copy(b_u)
 b_d.setLabel([53,-68,-20,-6,-4])

 c_u.setLabel([54,14,12,19,62])
 c_d=copy.copy(c_u)
 c_d.setLabel([51,-14,-12,-19,-62])
##################################  1  #####################################################
 c_ut=((c_u*(MPO_list[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,62,58,57],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-62,-58,-57])
 a_dt.setLabel([-55,-16,-17,-18,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])
 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*(((E1*E8)*(a_ut*a_dt))))*((E7*E6)*(c_ut*c_dt))

 A.permute([-62,-58,-57,-17,62,58,57,17],4)

 xt=copy.copy(plist[0])
 xt.setLabel([-62,-58,-57,-17])
 xt.permute([-62,-58,-57,-17],4)
 D_r[0]=xt*A
 D_r[0].permute([62,58,57,17],0)

 x=copy.copy(plist[0])
 x.permute([62,58,57,17],0)
 A.transpose()
 A=x*A
 A.permute([-62,-58,-57,-17],0)
 D_r[0]=D_r[0]+A
##########################################################################################
 
 
 Ap=((((E1*E8)*(a_u1*a_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-58,-57,-17],4)
 Ap.transpose()
 D_r[0]=D_r[0]+(-1.00)*Ap

 Ap=((((E1*E8)*(a_ut*a_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1)))*(((E7*E6)*(c_ut*c_d1)))
 Ap.permute([62,58,57,17],0)
 D_r[0]=D_r[0]+(-1.00)*Ap
 D_r[0].permute([62,58,57,17],3)
############################################################################################

##################################   2   #############################################
 c_ut=((c_u*(plist[0]*MPO_list[0])))
 a_ut=(a_u*(plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,64,58,57,18,2],5)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-64,-58,-57,-18,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])


 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))
 A.permute([-17,-64,-58,-57,17,64,58,57],4)

 xt=copy.copy(plist[1])
 xt.setLabel([-17,-64,-58,-57])
 xt.permute([-17,-64,-58,-57],4)
 D_r[1]=xt*A
 D_r[1].permute([17,64,58,57],0)

 x=copy.copy(plist[1])
 x.permute([17,64,58,57],0)
 A.transpose()
 A=x*A
 A.permute([-17,-64,-58,-57],0)
 D_r[1]=D_r[1]+A
##########################################################################################
 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-17,-64,-58,-57],4)
 Ap.transpose()
 D_r[1]=D_r[1]+(-1.00)*Ap

 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([17,64,58,57],0)
 D_r[1]=D_r[1]+(-1.00)*Ap
 D_r[1].permute([17,64,58,57],1)
############################################################################################


##################################  3  #####################################################
 a_ut=(a_u*(plist[1]*MPO_list[1]))
 c_ut=((c_u*(plist[0]*MPO_list[0])))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,66,59,60,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-17,-66,-59,-60,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))
 A.permute([-66,-59,-60,-18,66,59,60,18],4)

 xt=copy.copy(plist[2])
 xt.setLabel([-66,-59,-60,-18]) 
 xt.permute([-66,-59,-60,-18],4)
 D_r[2]=xt*A
 D_r[2].permute([66,59,60,18],0)

 x=copy.copy(plist[2])
 x.permute([66,59,60,18],0)
 A.transpose()
 A=x*A
 A.permute([-66,-59,-60,-18],0)
 D_r[2]=D_r[2]+A
##########################################################################################

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-66,-59,-60,-18],4)
 Ap.transpose()
 D_r[2]=D_r[2]+(-1.00)*Ap

 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([66,59,60,18],0)
 D_r[2]=D_r[2]+(-1.00)*Ap
 D_r[2].permute([66,59,60,18],3)
############################################################################################

##################################  4  #####################################################
 c_ut=((c_u*(plist[0]*MPO_list[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,68,59,60,20,6,4],5)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-17,-18,-2])
 b_dt.setLabel([-56,-68,-59,-60,-20,-6,-4])

 A=((((E1*E8)*(a_ut*a_dt)) *((E7*E6)*(c_ut*c_dt))) * ((((E4*E5)*d)))) * ((E2*E3)*(b_ut*b_dt))

 A.permute([-18,-68,-59,-60,18,68,59,60],4)

 xt=copy.copy(plist[3])
 xt.setLabel([-18,-68,-59,-60]) 
 xt.permute([-18,-68,-59,-60],4)
 D_r[3]=xt*A
 D_r[3].permute([18,68,59,60],0)

 x=copy.copy(plist[3])
 x.permute([18,68,59,60],0)
 A.transpose()
 A=x*A
 A.permute([-18,-68,-59,-60],0)
 D_r[3]=D_r[3]+A
##########################################################################################
 


 Ap=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-68,-59,-60],4)
 Ap.transpose()
 D_r[3]=D_r[3]+(-1.00)*Ap


 
 Ap=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)))*((E2*E3)*(b_ut*b_d1))
 Ap.permute([18,68,59,60],0)
 D_r[3]=D_r[3]+(-1.00)*Ap
 D_r[3].permute([18,68,59,60],1)
############################################################################################
 return D_r 





def Obtain_grad_four_MPO1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist):

 D_r=[0]*4
 c.setLabel([14,-14,12,-12,19,-19,17,-17])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 d_u1=copy.copy(d_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 b_u1.setLabel([55,18,20,6,4])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 b_u1=((b_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 b_u1.permute([-55,18,20,58,57,6,4,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 b_d1=copy.copy(b_u1)
 d_d1=copy.copy(d_u1)

 a_d1.transpose()
 b_d1.transpose()
 d_d1.transpose()
 

 a_d1.setLabel([-58,-57,-18,-2,-54,-16,-17])
 b_d1.setLabel([-6,-4,-59,-60,-55,-18,-20,-58,-57])
 d_d1.setLabel([-8,-20,-56,-19,-10,-59,-60])
 
 
 
 
##########################################################
 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(MPO_list[0])))
 b_ut=(b_u*(plist[1]*plist[3]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))

 a_ut.permute([-54,16,17,62,58,57,2],3)
 b_ut.permute([-55,18,20,6,4],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-62,-58,-57,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-18,-20])
 d_dt.setLabel([-8,-20,-56,-19,-10])
##################################  1  #####################################################

 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b_ut*b_dt)))*(((E7*E6)*(c))))*((E1*E8)*(a_ut*a_dt))

 A.permute([-62,-58,-57,-18,62,58,57,18],4)

 xt=copy.copy(plist[0])
 xt.setLabel([-62,-58,-57,-18])
 xt.permute([-62,-58,-57,-18],4)
 D_r[0]=xt*A
 D_r[0].permute([62,58,57,18],0)

 x=copy.copy(plist[0])
 x.permute([62,58,57,18],0)
 A.transpose()
 A=x*A
 A.permute([-62,-58,-57,-18],0)
 D_r[0]=D_r[0]+A
##########################################################################################
 
 Ap=((((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b_u1*b_dt)))*(((E7*E6)*(c))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-62,-58,-57,-18],4)
 Ap.transpose()
 D_r[0]=D_r[0]+(-1.00)*Ap

 Ap=((((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b_ut*b_d1)))*(((E7*E6)*(c))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([62,58,57,18],0)
 D_r[0]=D_r[0]+(-1.00)*Ap
 D_r[0].permute([62,58,57,18],3)
############################################################################################

##################################   2   #############################################

 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[3]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,64,58,57,20,6,4],5)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-64,-58,-57,-20])
 d_dt.setLabel([-8,-20,-56,-19,-10])


 A=(((((E4*E5)*(d_ut*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_dt))))*((E2*E3)*(b_ut*b_dt))
 A.permute([-18,-64,-58,-57,18,64,58,57],4)

 xt=copy.copy(plist[1])
 xt.setLabel([-18,-64,-58,-57])
 xt.permute([-18,-64,-58,-57],4)
 D_r[1]=xt*A
 D_r[1].permute([18,64,58,57],0)

 x=copy.copy(plist[1])
 x.permute([18,64,58,57],0)
 A.transpose()
 A=x*A
 A.permute([-18,-64,-58,-57],0)
 D_r[1]=D_r[1]+A
##########################################################################################
 Ap=(((((E4*E5)*(d_u1*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_u1*a_dt))))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-64,-58,-57],4)
 Ap.transpose()
 D_r[1]=D_r[1]+(-1.00)*Ap

 Ap=(((((E4*E5)*(d_ut*d_d1))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_d1))))*((E2*E3)*(b_ut*b_d1))
 Ap.permute([18,64,58,57],0)
 D_r[1]=D_r[1]+(-1.00)*Ap
 D_r[1].permute([18,64,58,57],1)
############################################################################################


##################################  3  #####################################################
 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[1]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,18,68,59,60,6,4],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-59,-60,-6,-4,-55,-18,-68])
 d_dt.setLabel([-8,-20,-56,-19,-10])

 A=(((((E4*E5)*(d_ut*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_dt))))*((E2*E3)*(b_ut*b_dt))
 A.permute([-60,-59,-20,-68,60,59,20,68],4)

 xt=copy.copy(plist[3])
 xt.setLabel([-60,-59,-20,-68]) 
 xt.permute([-60,-59,-20,-68],4)
 D_r[3]=xt*A
 D_r[3].permute([60,59,20,68],0)

 x=copy.copy(plist[3])
 x.permute([60,59,20,68],0)
 A.transpose()
 A=x*A
 A.permute([-60,-59,-20,-68],0)
 D_r[3]=D_r[3]+A
##########################################################################################

 Ap=(((((E4*E5)*(d_u1*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_u1*a_dt))))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-60,-59,-20,-68],4)
 Ap.transpose()
 D_r[3]=D_r[3]+(-1.00)*Ap

 Ap=(((((E4*E5)*(d_ut*d_d1))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_d1))))*((E2*E3)*(b_ut*b_d1))
 Ap.permute([60,59,20,68],0)
 D_r[3]=D_r[3]+(-1.00)*Ap
 D_r[3].permute([60,59,20,68],3)
############################################################################################

##################################  4  #####################################################
 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[1]*plist[3]*MPO_list[1]))
 d_ut=((d_u*(MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,18,20,6,4],3)
 d_ut.permute([-56,19,10,59,60,8,66],5)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-18,-20])
 d_dt.setLabel([-8,-66,-56,-19,-10,-59,-60])

 A=((((E2*E3)*(b_ut*b_dt)) )*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_dt))))* (((E4*E5)*(d_ut*d_dt)))

 A.permute([-66,-60,-59,-20,66,60,59,20],4)

 xt=copy.copy(plist[2])
 xt.setLabel([-66,-60,-59,-20]) 
 xt.permute([-66,-60,-59,-20],4)
 D_r[2]=xt*A
 D_r[2].permute([66,60,59,20],0)

 x=copy.copy(plist[2])
 x.permute([66,60,59,20],0)
 A.transpose()
 A=x*A
 A.permute([-66,-60,-59,-20],0)
 D_r[2]=D_r[2]+A
##########################################################################################
 


 Ap=((((E2*E3)*(b_u1*b_dt)) )*(((E7*E6)*(c))*((E1*E8)*(a_u1*a_dt))))* (((E4*E5)*(d_u1*d_dt)))
 Ap.permute([-66,-60,-59,-20],4)
 Ap.transpose()
 D_r[2]=D_r[2]+(-1.00)*Ap


 
 Ap=((((E2*E3)*(b_ut*b_d1)) )*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_d1))))* (((E4*E5)*(d_ut*d_d1)))
 Ap.permute([66,60,59,20],0)
 D_r[2]=D_r[2]+(-1.00)*Ap
 D_r[2].permute([66,60,59,20],1)
############################################################################################
 return D_r 



def make_equall_dis(a_u,b_u,c_u,d_u,a,b,c,d):

 for q in xrange(2):
  c_u, d_u=basicC.equall_dis_H(c_u, d_u)
  c_u, a_u=basicC.equall_dis_V(c_u, a_u)
  a_u, b_u=basicC.equall_dis_H(a_u, b_u)
  c_u, d_u=basicC.equall_dis_H(c_u, d_u)
  d_u, b_u=basicC.equall_dis_V(d_u, b_u)
  a_u, b_u=basicC.equall_dis_H(a_u, b_u)
  d_u, b_u=basicC.equall_dis_V(d_u, b_u)
  c_u, a_u=basicC.equall_dis_V(c_u, a_u)
  Maxa=MaxAbs(a_u)
  Maxb=MaxAbs(b_u)
  Maxc=MaxAbs(c_u)
  Maxd=MaxAbs(d_u)
  #print Maxa, Maxb, Maxc, Maxd

 a_u=max_ten(a_u)
 b_u=max_ten(b_u) 
 c_u=max_ten(c_u)
 d_u=max_ten(d_u)

 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)


 return a_u,b_u,c_u,d_u,a,b,c,d





