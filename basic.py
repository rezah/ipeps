import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time
import Move
import basic_FU


def Var_H(a_u,b_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi):
 Truncation=[0]
 c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Truncation=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D, Truncation)
 #print 'Truncation', Truncation[0]

 E1, E2, E3, E4, E5,E6=basic_FU.produce_Env_Hab(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 #basic_FU.test_env(E1, E2, E3, E4, E5,E6, a, b, c1,c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 E1, E2, E3, E4, E5,E6=basic_FU.proper_bond(E1, E2, E3, E4, E5,E6,D,d_phys)

 #basic_FU.test_energy(E1, E2, E3, E4, E5,E6, a, b, c1, c2, c3, c4, Ta1, Tb1, Ta2, Tb2, Ta3, Tb3, Ta4, Tb4, a_u, b_u, U)

 N_uni, l, l_d, r, r_d, q_u, qq_u = basic_FU.Qr_lQ_decom(a_u,b_u, E1, E2, E3, E4, E5,E6,D,d_phys)

 #basic_FU.test_energy_lr(N_uni, l, l_d, r, r_d, q_u,qq_u,U,E1, E2, E3, E4, E5,E6,a_u,b_u)

 lp, rp, lp_d, rp_d=basic_FU.initialize_lrprime(l, r, l_d, r_d, N_uni)
 lp, rp, lp_d, rp_d=basic_FU.initialize_SVD_lrprime(l, r, l_d, r_d, N_uni,U,D,d_phys)
 #lp, rp, lp_d, rp_d, N_uni, l, r, l_d, r_d,q_u, qq_u, a_u, b_u=basic_FU.initialize_Positiv_lrprime(l, r, l_d, r_d, N_uni, U, D, d_phys, q_u, qq_u,a_u,b_u)

 #basic_FU.test_energy_lr(N_uni, l, l_d, r, r_d, q_u,qq_u,U,E1, E2, E3, E4, E5,E6,a_u,b_u)

 #rp, rp_d, lp, lp_d=basic_FU.Do_optimization(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
 rp, rp_d, lp, lp_d=basic_FU.Do_optimization_Full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
 #rp, rp_d, lp, lp_d=basic_FU.Do_optimization_Grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
 
 
 lp, rp, lp_d, rp_d=basic_FU.Equall_Dist(lp, rp,D,d_phys)
 ap_u, bp_u=basic_FU.recover(lp, rp, q_u, qq_u)
 #Distance_val=basic_FU.Distance(l, r, lp, rp ,N_uni,U)
 #A=basic_FU.final_test_distance(ap_u, bp_u, a_u, b_u,E1, E2, E3, E4, E5,E6,U,N_uni)
 #Distance_val=basic_FU.Distance(l, r, lp, rp ,N_uni,U)
 #print 'Final', A[0], Distance_val[0]
 
 if ( (abs(bp_u.getBlock().absMax()) < 0.50e-1) or (abs(bp_u.getBlock().absMax()) > 0.50e+1)   ):
  bp_u=bp_u*(1.00/bp_u.getBlock().absMax());
  #print 'max', bp_u.getBlock().absMax()
 else: bp_u=bp_u;

 if ( (abs(ap_u.getBlock().absMax()) < 0.50e-1) or (abs(ap_u.getBlock().absMax()) > 0.50e+1)   ):
  ap_u=ap_u*(1.00/ap_u.getBlock().absMax()); 
  #print 'max', ap_u.getBlock().absMax()
 else: ap_u=ap_u;

 a=make_ab(ap_u)
 b=make_ab(bp_u)
 
 return ap_u, bp_u, a, b




def Var_V(c_u,a_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi):
 
 Truncation=[0]
 c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Truncation=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D, Truncation)
 
 #print 'Truncation', Truncation[0]
 #norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
 #print 'norm=', norm[0] 
 
 E1, E2, E3, E4, E5,E6=basic_FU.produce_Env_Hac(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E1, E2, E3, E4, E5,E6=reorder_env(E1, E2, E3, E4, E5,E6)

 c_u.permute([0,2,3,4,1],3)
 a_u.permute([0,2,3,4,1],3)
 c_u.setLabel([0,1,2,3,4])
 a_u.setLabel([0,1,2,3,4])
 c=make_ab(c_u)
 a=make_ab(a_u)
 
 #basic_FU.test_env(E1, E2, E3, E4, E5,E6, c, a, c1,c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 
 E1, E2, E3, E4, E5,E6=basic_FU.proper_bond(E1, E2, E3, E4, E5,E6,D,d_phys)
 
 #basic_FU.test_energy(E1, E2, E3, E4, E5,E6, c, a, c1, c2, c3, c4, Ta1, Tb1, Ta2, Tb2, Ta3, Tb3, Ta4, Tb4, c_u, a_u, U)
 
 N_uni, l, l_d, r, r_d, q_u, qq_u = basic_FU.Qr_lQ_decom(c_u,a_u, E1, E2, E3, E4, E5,E6,D,d_phys)

 #basic_FU.test_energy_lr(N_uni, l, l_d, r, r_d, q_u,qq_u,h,E1, E2, E3, E4, E5,E6,c_u,a_u)

 lp, rp, lp_d, rp_d=basic_FU.initialize_lrprime(l, r, l_d, r_d, N_uni)
 lp, rp, lp_d, rp_d=basic_FU.initialize_SVD_lrprime(l, r, l_d, r_d, N_uni,U,D,d_phys)
 #lp, rp, lp_d, rp_d, N_uni, l, r, l_d, r_d,q_u, qq_u, c_u, a_u=basic_FU.initialize_Positiv_lrprime(l, r, l_d, r_d, N_uni, U, D, d_phys, q_u, qq_u,c_u,a_u)

 #basic_FU.test_energy_lr(N_uni, l, l_d, r, r_d, q_u,qq_u,U,E1, E2, E3, E4, E5,E6,c_u,a_u)


 #rp, rp_d, lp, lp_d=basic_FU.Do_optimization(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
 rp, rp_d, lp, lp_d=basic_FU.Do_optimization_Full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
 #rp, rp_d, lp, lp_d=basic_FU.Do_optimization_Grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
 
 
 lp, rp, lp_d, rp_d=basic_FU.Equall_Dist(lp, rp,D,d_phys)
 cp_u, ap_u=basic_FU.recover(lp, rp, q_u, qq_u)
 #Distance_val=basic_FU.Distance(l, r, lp, rp ,N_uni,U)
 #A=basic_FU.final_test_distance(cp_u, ap_u, c_u, a_u,E1, E2, E3, E4, E5,E6,U,N_uni)
 #Distance_val=basic_FU.Distance(l, r, lp, rp ,N_uni,U)
 #print 'Final', A[0], Distance_val[0]
 
 cp_u.setLabel([0,1,2,3,4])
 ap_u.setLabel([0,1,2,3,4])

 cp_u.setLabel([0,2,3,4,1])
 ap_u.setLabel([0,2,3,4,1])

 cp_u.permute([0,1,2,3,4],3)
 ap_u.permute([0,1,2,3,4],3)

 
 if ( (abs(cp_u.getBlock().absMax()) < 0.50e-1) or (abs(cp_u.getBlock().absMax()) > 0.50e+1)   ):
  cp_u=cp_u*(1.00/cp_u.getBlock().absMax()); 
  #print 'max', cp_u.getBlock().absMax()
 else: cp_u=cp_u;

 if ( (abs(ap_u.getBlock().absMax()) < 0.50e-1) or (abs(ap_u.getBlock().absMax()) > 0.50e+1)   ):
  ap_u=ap_u*(1.00/ap_u.getBlock().absMax()); 
  #print 'max', ap_u.getBlock().absMax()
 else: ap_u=ap_u;
 
 
 cp=make_ab(cp_u)
 ap=make_ab(ap_u)
 
 

 return cp_u, ap_u, cp, ap


def reorder_env(E1, E2, E3, E4, E5,E6):

 E5p=copy.copy(E1)
 E6p=copy.copy(E2)
 E1p=copy.copy(E3)
 E2p=copy.copy(E4)
 E3p=copy.copy(E5)
 E4p=copy.copy(E6)
 return E1p,E2p,E3p,E4p,E5p,E6p



def Initialize_function(Gamma,Landa):
 for i in xrange(len(Gamma)):
  Gamma[i].randomize()
  Gamma[i]=Gamma[i]*(1.00/Gamma[i].norm())
  
 D=Gamma[0].bond(1).dim()
 matrix_Iden=uni10.Matrix(D, D)
 matrix_Iden.identity()
 for i in xrange(len(Landa)):
  Landa[i].putBlock(matrix_Iden)


 
def matSp():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0, 1, 0, 0]);

def matSm():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0, 0, 1, 0]);

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);

def transverseIsing(h):
    spin = 0.5
    sx = (matSp()+matSm())
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    ham =uni10.otimes(sz,sz)*(-1)+(-0.2500)*float(h)*(uni10.otimes(iden,sx)+uni10.otimes(sx,iden))
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "TFIM");
    H.putBlock(ham)
    return H


def Ham(h):
    spin = 0.5
    sx = (matSp()+matSm())
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    ham =uni10.otimes(sz,sz)*(-2.00)+(-0.500)*float(h)*(uni10.otimes(iden,sx)+uni10.otimes(sx,iden))
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "TFIM");
    H.putBlock(ham)
    return H

def makeTab(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor([bdi, bdi1, bdo])
 Tem0.randomize()
 Tem0*=(1.00/Tem0.norm())
 Tem1=uni10.UniTensor([bdi, bdi1, bdo])
 Tem1.randomize()
 Tem1*=(1.00/Tem1.norm())
 return Tem0, Tem1


def makec1(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 
 Tem0=uni10.UniTensor([bdi, bdo])
 Tem0.randomize()
 Tem0*=(1.00/Tem0.norm())
 Tem1=uni10.UniTensor([bdi, bdo])
 Tem1.randomize()
 Tem1*=(1.00/Tem1.norm())
 Tem2=uni10.UniTensor([bdi, bdo])
 Tem2.randomize()
 Tem2*=(1.00/Tem2.norm())
 Tem3=uni10.UniTensor([bdi, bdo])
 Tem3.randomize()
 Tem3*=(1.00/Tem3.norm())
 return Tem0,Tem1,Tem2,Tem3
 
def makeab(Landa,Gamma):

 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa)) ]
 Landa_sq=sqrt(Landa_cp)


 a_u=copy.copy(Gamma)
 a_d=copy.copy(Gamma)

 Landa_sq[0].setLabel([1,-1])
 Landa_sq[1].setLabel([2,-2])
 Landa_sq[2].setLabel([3,-3])
 Landa_sq[3].setLabel([4,-4])
 
 a_u.setLabel([0,1,2,3,4])
 a_d.setLabel([0,1,2,3,4])
 
 a_u=a_u*(Landa_sq[0]*Landa_sq[1]*Landa_sq[2]*Landa_sq[3])
 a_d=a_d*(Landa_sq[0]*Landa_sq[1]*Landa_sq[2]*Landa_sq[3])
 a_u.permute([0,-1,-2,-3,-4],3)
 a_d.permute([0,-1,-2,-3,-4],3)

 a_u.setLabel([0,1,2,3,4])
 a_d.setLabel([0,-1,-2,-3,-4])
 a=a_u*a_d
 a.combineBond([1,-1])
 a.combineBond([2,-2])
 a.combineBond([3,-3])
 a.combineBond([4,-4])
 a.permute([1,2,3,4],2)

 return a_u, a


 
def   sqrt(Landa):
 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa))   ]
 for q in xrange(len(Landa_cp)): 
  D=Landa_cp[q].bond()[0].dim()
  Landa_cpm=Landa_cp[q].getBlock()
  Landam=Landa[q].getBlock()
  for i in xrange(D):
   for j in xrange(D):
    Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
  Landa_cp[q].putBlock(Landa_cpm)
 return Landa_cp

 
def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 invL2 = uni10.Matrix(Landa2.bond()[0].dim(), Landa2.bond()[1].dim())
 D=Landa2.bond()[0].dim()
 for i in xrange(Landa2.bond()[0].dim()):
   for j in xrange(Landa2.bond()[0].dim()):
    invL2[i*D+j] = 0 if ((Landa2[i*D+j].real) < 1.0e-12) else (1.00 / (Landa2[i*D+j].real))
 invLanda2.putBlock(invL2)
 return invLanda2
def magnetization(a_u):
 a_u.setLabel([0,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.setLabel([5,-1,-2,-3,-4])
 sz = matSz()
 bdi = uni10.Bond(uni10.BD_IN, 2);
 bdo = uni10.Bond(uni10.BD_OUT, 2);
 Z =  uni10.UniTensor([bdi, bdo], "Z");
 Z.putBlock(sz)   
 Z.setLabel([5,0])
 result=a_uc*Z*a_u
 result.combineBond([1,-1])
 result.combineBond([2,-2])
 result.combineBond([3,-3])
 result.combineBond([4,-4])
 result.permute([1,2,3,4], 2)

 return result


def make_ab(a_u):
 a_u.setLabel([0,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.setLabel([0,-1,-2,-3,-4])
 result=a_uc*a_u
 result.combineBond([1,-1])
 result.combineBond([2,-2])
 result.combineBond([3,-3])
 result.combineBond([4,-4])
 result.permute([1,2,3,4], 2)

 return result




def corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation):

 z1=copy.copy(a)
 z1.randomize()
 z1.identity()
 Accuracy=1.00e-9
 Truncation=[0]
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 #print  '\n', '\n', 'CTM'
 while Loop_iter is 0: 
  c1, Ta4, Tb4, c4=Move.add_left(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D,Truncation)
  c1, Ta4, Tb4, c4=Move.add_left(c1,Tb4,Ta4,c4,Ta1,Tb3,b,d,chi,D,Truncation)

  c2, Ta2, Tb2, c3=Move.add_right(c2,Ta2,Tb2,c3,Ta1,Tb3,b,d,chi,D)
  c2, Ta2, Tb2, c3=Move.add_right(c2,Ta2,Tb2,c3,Tb1,Ta3,a,c,chi,D)

  c1, Ta1, Tb1, c2=Move.add_up(c1,Tb1,Ta1,c2,Tb4,Ta2,a,b,chi,D)
  c1, Ta1, Tb1, c2=Move.add_up(c1,Tb1,Ta1,c2,Ta4,Tb2,c,d,chi,D)

  c4, Ta3, Tb3, c3=Move.add_down(c4,Ta3,Tb3,c3,Ta4,Tb2,c,d,chi,D)
  c4, Ta3, Tb3, c3=Move.add_down(c4,Ta3,Tb3,c3,Tb4,Ta2,a,b,chi,D)

  norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  norm1=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
  E0=E1
  E1=abs(norm1[0])/abs(norm[0])
  if (abs((E0-E1)/E1) < Accuracy) and ( count > 2):Loop_iter=1;
  count+=1
  if (count > 20 ): print 'break! CTM'; break;
  #print E1,abs((E0-E1)/E1), count
  #print abs(norm1[0])/abs(norm[0]), Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()
 
 print 'CTM', norm[0], Truncation[0], '\n', '\n' 
 return c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, Truncation




def z_value(a,b,c,d,az,bz,cz,dz,chi,D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):

 z1=copy.copy(az)
 z2=copy.copy(bz)
 z3=copy.copy(cz)
 z4=copy.copy(dz)
 
 Truncation=[0]



 norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
 norm1=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
 norm2=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,z2,c,d)
 norm3=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,z3,d)
 norm4=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,z4)
 #print abs(norm1[0])/abs(1.00*norm[0]), abs(norm2[0])/abs(1.00*norm[0])
 
 return abs(norm1[0]+norm2[0]+norm3[0]+norm4[0])/abs(4.00*norm[0])
def E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi):
 E_ab=Energy_h(a_u,b_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 E_cd=Energy_h(c_u,d_u,c,d,a,b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 E_ba=Energy_h(b_u,a_u,b,a,d,c,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 E_dc=Energy_h(d_u,c_u,d,c,b,a,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)


 #print E_ab,E_ba,E_cd, E_dc, (E_ab+E_ba+E_cd+E_dc) / 4.00

 #print '\n','\n','\n'  

 E_ca=Energy_v(c_u,a_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 E_ac=Energy_v(a_u,c_u,c,d,a,b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 E_db=Energy_v(d_u,b_u,b,a,d,c,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 E_bd=Energy_v(b_u,d_u,d,c,b,a,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)

 #print E_ca,E_ac,E_db, E_bd, (E_ca+E_ac+E_db+E_bd) / 4.00

 return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc) / 4.00)

def Energy_v(c_u,a_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi):
 Truncation=[0]
 c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Truncation=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)

 E1, E2, E3, E4, E5,E6=produce_Env_Hac(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E1, E2, E3, E4, E5,E6=reorder_env(E1, E2, E3, E4, E5,E6)
 cp_u=copy.copy(c_u)
 ap_u=copy.copy(a_u)
 cp=copy.copy(c)
 ap=copy.copy(a)
 
 cp_u.setLabel([0,1,2,3,4])
 ap_u.setLabel([0,1,2,3,4])

 cp_u.permute([0,2,3,4,1],3)
 ap_u.permute([0,2,3,4,1],3)
 cp_u.setLabel([0,1,2,3,4])
 ap_u.setLabel([0,1,2,3,4])
 cp=make_ab(cp_u)
 ap=make_ab(ap_u)
 
 #test_env(E1, E2, E3, E4, E5,E6, cp, ap, c1,c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 
 E1, E2, E3, E4, E5,E6=proper_bond(E1, E2, E3, E4, E5,E6,D,d_phys)
 U=transverseIsing(h)
 #U=Heisenberg(h)
 
 
 E_ca=test_energy(E1, E2, E3, E4, E5,E6, cp, ap, c1, c2, c3, c4, Ta1, Tb1, Ta2, Tb2, Ta3, Tb3, Ta4, Tb4, cp_u, ap_u, U)

 return E_ca

def Energy_h(a_u,b_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi):
 Truncation=[0]
 c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Truncation=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)
 E1, E2, E3, E4, E5,E6=produce_Env_Hab(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 #test_env(E1, E2, E3, E4, E5,E6, a, b, c1,c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 E1, E2, E3, E4, E5,E6=proper_bond(E1, E2, E3, E4, E5,E6,D,d_phys)
 U=transverseIsing(h)
 #U=Heisenberg(h)
 
 E_ab=test_energy(E1, E2, E3, E4, E5,E6, a, b, c1, c2, c3, c4, Ta1, Tb1, Ta2, Tb2, Ta3, Tb3, Ta4, Tb4, a_u, b_u, U)
 
########################################################################3 
 return E_ab

def produce_Env_Hab(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys):

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 E1=c1*Tb1
 E1.permute([0,2,3],2)
 E1.setLabel([7,3,8])
 
 c2.setLabel([1,0])
 Ta1.setLabel([3,2,1])
 E6=c2*Ta1
 E6.permute([0,2,3],2)
 E6.setLabel([9,6,8])
 
 E5=Ta2
 E5.setLabel([10,5,9])
 
 E2=Tb4
 E2.setLabel([12,13,7])
 
 c3.setLabel([7,8])
 Tb2.setLabel([7,15,22])
 Tb3.setLabel([13,14,8])
 d.setLabel([16,14,15,23])
 E4=(((c3*Tb3)*Tb2)*d)
 E4.combineBond([13,16])
 E4.permute([13,23,22],2)
 E4.setLabel([11,4,10])
 
 c4.setLabel([11,10])
 Ta4.setLabel([11,17,18])
 Ta3.setLabel([10,12,13])
 c.setLabel([17,12,16,19])
 E3=(((c4*Ta4)*Ta3)*c)
 E3.combineBond([13,16])
 E3.permute([13,19,18],2)
 E3.setLabel([11,1,12])
 return E1, E2, E3, E4, E5,E6
def  proper_bond(E1, E2, E3, E4, E5,E6,D,d_phys):
 bd=uni10.Bond(E1.bond()[1].type(), D)
 A=uni10.UniTensor([E1.bond()[0],bd,bd,E1.bond()[2]])
 A.putBlock(E1.getBlock())
 E1=copy.copy(A)
 E1.setLabel([7,3,-3,8])
 
 bd=uni10.Bond(E2.bond()[1].type(), D)
 A=uni10.UniTensor([E2.bond()[0],bd,bd,E2.bond()[2]])
 A.putBlock(E2.getBlock())
 E2=copy.copy(A)
 E2.setLabel([12,13,-13,7])


 bd=uni10.Bond(E3.bond()[1].type(), D)
 A=uni10.UniTensor([E3.bond()[0],bd,bd,E3.bond()[2]])
 A.putBlock(E3.getBlock())
 E3=copy.copy(A)
 E3.setLabel([11,1,-1,12])


 bd=uni10.Bond(E4.bond()[1].type(), D)
 A=uni10.UniTensor([E4.bond()[0],bd,bd,E4.bond()[2]])
 A.putBlock(E4.getBlock())
 E4=copy.copy(A)
 E4.setLabel([11,4,-4,10])


 bd=uni10.Bond(E5.bond()[1].type(), D)
 A=uni10.UniTensor([E5.bond()[0],bd,bd,E5.bond()[2]])
 A.putBlock(E5.getBlock())
 E5=copy.copy(A)
 E5.setLabel([10,5,-5,9])

 bd=uni10.Bond(E6.bond()[1].type(), D)
 A=uni10.UniTensor([E6.bond()[0],bd,bd,E6.bond()[2]])
 A.putBlock(E6.getBlock())
 E6=copy.copy(A)
 E6.setLabel([9,6,-6,8])

 return E1, E2, E3, E4, E5,E6
 

def  test_energy(E1, E2, E3, E4, E5,E6, a, b, c1,c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, a_u,b_u,U):


 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.setLabel([20,-13,-1,-2,-3])
 #print a_u.printDiagram()
 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.setLabel([40,-2,-4,-5,-6])
 A=(((E2*(a_u*a_d))*E1)*E3)*(((((b_u*b_d)*E5)*E6)*E4))
 Norm_f=A
 #print 'Norm=', A[0]

 a_u.setLabel([20,13,1,2,3])
 a_d.setLabel([-20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d.setLabel([-40,-2,-4,-5,-6])


 U.setLabel([-20,-40,20,40])


 #print a_u.printDiagram(), b_u.printDiagram() 
 B=((((E2*(a_u*a_d))*E1)*E3)*U)*(((((b_u*b_d)*E5)*E6)*E4))
 return B[0]/A[0]

def test_env(E1, E2, E3, E4, E5,E6, a, b, c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 a.setLabel([13,1,2,3])
 b.setLabel([2,4,5,6])
 A=(((E2*a)*E1)*E3)*((((b*E5)*E6)*E4))
 print 'norm=', A[0]
def produce_Env_Hac(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys):


 c1.setLabel([0,1])
 Tb4.setLabel([3,2,0])
 E2=c1*Tb4
 E2.permute([3,2,1],2)
 E2.setLabel([10,2,12])
 E2.setLabel([8,6,9])
 E2.permute([9,6,8],1)
 
 c4.setLabel([0,1])
 Ta3.setLabel([1,2,3])
 E4=c4*Ta3
 E4.permute([0,2,3],2)
 E4.setLabel([9,6,8])
 E4.setLabel([7,13,12])
 E4.permute([12,13,7],1)
 
 E1=Tb1
 E1.setLabel([12,1,13])
 E1.setLabel([9,5,10])
 E1.permute([10,5,9],1)

 E3=Ta4
 E3.setLabel([9,5,10])
 E3.setLabel([7,3,8])
 
 
 c2.setLabel([5,6])
 Ta1.setLabel([21,20,5])
 Ta2.setLabel([22,19,6])
 b.setLabel([4,23,19,20])
 E6=(((c2*Ta2)*Ta1)*b)
 E6.combineBond([23,22])
 E6.permute([21,4,23],2)
 E6.setLabel([13,3,11])
 E6.setLabel([10,4,11])
 E6.permute([11,4,10],2)
 
 
 
 c3.setLabel([7,8])
 Tb2.setLabel([7,15,22])
 Tb3.setLabel([13,14,8])
 d.setLabel([16,14,15,23])
 E5=(((c3*Tb3)*Tb2)*d)
 E5.combineBond([23,22])
 E5.permute([13,16,23],2)
 E5.setLabel([8,7,11])
 E5.setLabel([12,1,11])
 E5.permute([11,1,12],2)



 return E1, E2, E3, E4, E5,E6
def reorder_env(E1, E2, E3, E4, E5,E6):

 E5p=copy.copy(E1)
 E6p=copy.copy(E2)
 E1p=copy.copy(E3)
 E2p=copy.copy(E4)
 E3p=copy.copy(E5)
 E4p=copy.copy(E6)
 return E1p,E2p,E3p,E4p,E5p,E6p
def make_ab(a_u):
 a_u.setLabel([0,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.setLabel([0,-1,-2,-3,-4])
 result=a_uc*a_u
 result.combineBond([1,-1])
 result.combineBond([2,-2])
 result.combineBond([3,-3])
 result.combineBond([4,-4])
 result.permute([1,2,3,4], 2)

 return result


