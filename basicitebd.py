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
import Move1

def corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D):

 z1=copy.copy(a)
 z1.randomize()
 Accuracy=1.00e-9
 Truncation=[0]
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 while Loop_iter is 0: 
  c1, Ta4, Tb4, c4=Move.add_left(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D)
  c1, Ta4, Tb4, c4=Move.add_left(c1,Tb4,Ta4,c4,Ta1,Tb3,b,d,chi,D)

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
  if (abs((E0-E1)/E1) < Accuracy) and ( count > 15):Loop_iter=1;
  count+=1
  #print E1,abs((E0-E1)/E1), count


  #print abs(norm1[0]+norm2[0]+norm3[0]+norm4[0])/abs(4.00*norm[0]), Truncation[0]
 #print 'hi', Tb1
 return c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4




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
 print abs(norm1[0])/abs(1.00*norm[0]), abs(norm2[0])/abs(1.00*norm[0])
 
 return abs(norm1[0]+norm2[0]+norm3[0]+norm4[0])/abs(4.00*norm[0])



def corner_transfer_matrix_onesite(a, az,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):

 z=copy.copy(az)
 for q in xrange(40):
  c1, c4, Ta4= Move1.left(c1,c4,Ta4,Ta1,Ta3,chi,a)
  c3, c4 , Ta3= Move1.down(c4,c3,Ta3,Ta4,Ta2,chi,a)
  c2, c3 , Ta2= Move1.right(c3,c2,Ta2,Ta1,Ta3,chi,a)
  c1, c2 , Ta1= Move1.up(c1,c2,Ta1,Ta4,Ta2,chi,a)

 norm0=Move1.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,a)
 norm1=Move1.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,z)
 #print 'norm',norm0[0]
 return norm1[0]/norm0[0]








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





def makeTab(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor([bdi, bdi1, bdo])
 Tem0.randomize()
 Tem1=uni10.UniTensor([bdi, bdi1, bdo])
 Tem1.randomize()
 
 return Tem0, Tem1


def makec1(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 
 Tem0=uni10.UniTensor([bdi, bdo])
 Tem0.randomize()
 Tem1=uni10.UniTensor([bdi, bdo])
 Tem1.randomize()
 Tem2=uni10.UniTensor([bdi, bdo])
 Tem2.randomize()
 Tem3=uni10.UniTensor([bdi, bdo])
 Tem3.randomize()
 
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
 

def Initialize_function(Gamma,Landa):
 for i in xrange(len(Gamma)):
  Gamma[i].randomize()
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



 
def update_rlink(Gamma,Landa,U,D,d_phys):
 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_a.setLabel([0,1,2,3,4])
 Gamma_b.setLabel([5,6,7,8,9])
 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([10,11,0,5])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa3p=copy.copy(Landa[2])


 Landa1.setLabel([3,6])
 Landa2.setLabel([2,-2])
 Landa3.setLabel([1,-1])
 Landa4.setLabel([4,-4])
 Landa8.setLabel([9,-9])
 Landa3p.setLabel([8,-8])
 Landa7.setLabel([7,-7])


 Tem_tensor=Gamma_a*(Landa1*Landa2*Landa3*Landa4)
 Tem_tensor2=Gamma_b*(Landa8*Landa3p*Landa7)
 
 Theta=(Hamiltonian*Tem_tensor)*Tem_tensor2
 Theta.permute([10,-1,-2,-4,11,-7,-8,-9],4)
 ##printTheta.printDiagram() 
 
 svd1=Theta.getBlock()
 svd = svd1.svd()
 #svd = Theta.getBlock().svd()
 # Truncation
 sv = svd[1]
 
 
 Sum=0
 for i in xrange(sv.row()):
   if (i>=D):
    Sum=sv[i]+Sum
 #print'truncation=', Sum
 norm = sv.resize(D, D)
 norm=norm.norm()
 sv *= 1.00 / norm
 #print'norm=', norm
 
 Landa[0].putBlock(sv) 
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 GsA=uni10.UniTensor([bdi_pys,bdi,bdi,bdi,bdo], "GsA")
 GsB=uni10.UniTensor([bdi,bdi_pyso,bdo,bdo,bdo], "GsB")
 
 GsA.setLabel([10,-1,-2,-4,3])
 GsB.setLabel([6,11,-7,-8,-9])
 
 GsA.putBlock(svd[0].resize(svd[0].row(), D))
 GsB.putBlock(svd[2].resize(D, svd[2].col()))
 ##printGsA.printDiagram(),GsB.printDiagram()
 
# GsA.permute([10,-1,-2,3,-4],3)
# GsB.permute([11,6,-7,-8,-9],3)

 
 invLanda2=inverse(Landa2)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda2.setLabel([-2,2])
 invLanda3.setLabel([-1,1])
 invLanda4.setLabel([-4,4])
 GsA=GsA*(invLanda2*invLanda3*invLanda4)

 invLanda8=inverse(Landa8)
 invLanda3=inverse(Landa3)
 invLanda7=inverse(Landa7)


 invLanda8.setLabel([-9,9])
 invLanda3.setLabel([-8,8])
 invLanda7.setLabel([-7,7])
 GsB=GsB*(invLanda8*invLanda3*invLanda7)
 ##print'identity', invLanda4.getBlock()*Landa4.getBlock()
 ##print'identity', invLanda3.getBlock()*Landa3.getBlock()
 ##print'identity', invLanda2.getBlock()*Landa2.getBlock()
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)
 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])


 
 Gamma[0].putBlock(GsA.getBlock())
 Gamma[1].putBlock(GsB.getBlock())
 Gamma[0]*=(1.00/Gamma[0].norm())
 Gamma[1]*=(1.00/Gamma[1].norm())

 
 ##printGamma[0].printDiagram(), Gamma[1].printDiagram() 

 #Gs[A].permute([-1, 3, 1], 1)
 #bondrm(Gs[A], Ls[B], 0)
 #bondrm(Gs[B], Ls[B], 1)  

def update_ulink(Gamma,Landa,U,D,d_phys):
 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_a.setLabel([0,1,2,3,4])
 Gamma_b.setLabel([5,6,7,8,9])
 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([10,11,0,5])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa4p=copy.copy(Landa[3])


 Landa1.setLabel([3,-3])
 Landa2.setLabel([2,9])
 Landa3.setLabel([1,-1])
 Landa4.setLabel([4,-4])
 Landa5.setLabel([6,-6])
 Landa4p.setLabel([7,-7])
 Landa6.setLabel([8,-8])


 Tem_tensor=Gamma_a*(Landa1*Landa2*Landa3*Landa4)
 Tem_tensor2=Gamma_b*(Landa5*Landa4p*Landa6)
 
 Theta=(Hamiltonian*Tem_tensor)*Tem_tensor2
 Theta.permute([10,-1,-3,-4,11,-6,-7,-8],4)
 ##printTheta.printDiagram() 
 
 svd1=Theta.getBlock()
 svd = svd1.svd()
 #svd = Theta.getBlock().svd()
 # Truncation
 sv = svd[1]
 
 
 Sum=0
 for i in xrange(sv.row()):
   if (i>=D):
    Sum=sv[i]+Sum
 #print'truncation=', Sum
 norm = sv.resize(D, D)
 norm=norm.norm()
 sv *= 1.00 / norm
 #print'norm=', norm
 
 Landa[1].putBlock(sv) 
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 GsA=uni10.UniTensor([bdi_pys,bdi,bdi,bdi,bdo], "GsA")
 GsB=uni10.UniTensor([bdi,bdi_pyso,bdo,bdo,bdo], "GsB")
 
 GsA.setLabel([10,-1,-3,-4,2])
 GsB.setLabel([9,11,-6,-7,-8])
 
 GsA.putBlock(svd[0].resize(svd[0].row(), D))
 GsB.putBlock(svd[2].resize(D, svd[2].col()))
 ##printGsA.printDiagram(),GsB.printDiagram()
 
# GsA.permute([10,-1,-2,3,-4],3)
# GsB.permute([11,6,-7,-8,-9],3)

 
 invLanda1=inverse(Landa1)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda1.setLabel([-3,3])
 invLanda3.setLabel([-1,1])
 invLanda4.setLabel([-4,4])
 GsA=GsA*(invLanda1*invLanda3*invLanda4)

 invLanda4=inverse(Landa4)
 invLanda6=inverse(Landa6)
 invLanda5=inverse(Landa5)


 invLanda4.setLabel([-7,7])
 invLanda5.setLabel([-6,6])
 invLanda6.setLabel([-8,8])
 GsB=GsB*(invLanda5*invLanda4*invLanda6)
 ##print'identity', invLanda4.getBlock()*Landa4.getBlock()
 ##print'identity', invLanda3.getBlock()*Landa3.getBlock()
 ##print'identity', invLanda2.getBlock()*Landa2.getBlock()
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)
 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])


 
 Gamma[0].putBlock(GsA.getBlock())
 Gamma[1].putBlock(GsB.getBlock())
 Gamma[0]*=(1.00/Gamma[0].norm())
 Gamma[1]*=(1.00/Gamma[1].norm())

 
 ##printGamma[0].printDiagram(), Gamma[1].printDiagram() 

 #Gs[A].permute([-1, 3, 1], 1)
 #bondrm(Gs[A], Ls[B], 0)
 #bondrm(Gs[B], Ls[B], 1)  
 
def update_rlink_eff(Gamma,Landa,U,D,d_phys):
 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_a.setLabel([0,1,2,3,4])
 Gamma_b.setLabel([5,6,7,8,9])
 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([10,11,0,5])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa3p=copy.copy(Landa[2])


 Landa1.setLabel([3,6])
 Landa2.setLabel([2,-2])
 Landa3.setLabel([1,-1])
 Landa4.setLabel([4,-4])
 Landa8.setLabel([9,-9])
 Landa3p.setLabel([8,-8])
 Landa7.setLabel([7,-7])

 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 q_uni=uni10.UniTensor([bdi,bdi,bdi,bdoB])
 r_uni=uni10.UniTensor([bdiB,bdi_pyso,bdo])

 q_uni.setLabel([-4,-1,-2,20])
 r_uni.setLabel([20,0,3])

 
 Left=Gamma_a*(Landa2*Landa3*Landa4)
 Left.permute([-4,-1,-2,0,3],3)
 Left_m=Left.getBlock()
 qr=Left_m.qr()

 q_uni.putBlock(qr[0])
 r_uni.putBlock(qr[1])

# q_unit=copy.copy(q_uni)
# q_unit.transpose()
# print '1', q_uni, q_uni, q_unit.getBlock()*q_uni.getBlock()

 l_uni=uni10.UniTensor([bdi_pys,bdi,bdoB])
 qq_uni=uni10.UniTensor([bdiB,bdo,bdo,bdo])
 l_uni.setLabel([5,6,40])
 qq_uni.setLabel([40,-7,-8,-9])


 
 Right=Gamma_b*(Landa8*Landa3p*Landa7)
 Right.permute([5,6,-7,-8,-9],2)
 Right_m=Right.getBlock()
 lq=Right_m.lq()

 l_uni.putBlock(lq[0])
 qq_uni.putBlock(lq[1])

 
 #qqt_uni=copy.copy( qq_uni)
 #qqt_uni.transpose()
 #print '2',  qq_uni,  l_uni, qq_uni.getBlock()*qqt_uni.getBlock()

 Theta=(r_uni*Landa1*Hamiltonian)*l_uni
 Theta.permute([10,20,11,40],2)
 ##printTheta.printDiagram() 
 
 svd=Theta.getBlock()
 svd = svd.svd()
 sv = svd[1]
 
 
 Sum=0
 for i in xrange(sv.row()):
   if (i>=D):
    Sum=sv[i]+Sum
 #print'truncation=', Sum
 norm = sv.resize(D, D)
 norm=norm.norm()
 sv *= 1.00 / norm
 #print'norm=', norm
 Landa[0].putBlock(sv) 
 
 U=uni10.UniTensor([bdi_pys,bdiB,bdo], "U")
 V=uni10.UniTensor([bdi,bdi_pyso,bdoB], "V")
 U.setLabel([10,20,3])
 V.setLabel([6,11,40])
 
 #print svd[0].row()
 U.putBlock(svd[0].resize(svd[0].row(), D))
 V.putBlock(svd[2].resize(D, svd[2].col()))
 GsA=U*q_uni
 GsB=V*qq_uni
 
 GsA.permute([10,-1,-2,3,-4],3)
 GsB.permute([11,6,-7,-8,-9],3)


 
 invLanda2=inverse(Landa2)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda2.setLabel([-2,2])
 invLanda3.setLabel([-1,1])
 invLanda4.setLabel([-4,4])
 GsA=GsA*(invLanda2*invLanda3*invLanda4)

 invLanda8=inverse(Landa8)
 invLanda3=inverse(Landa3)
 invLanda7=inverse(Landa7)


 invLanda8.setLabel([-9,9])
 invLanda3.setLabel([-8,8])
 invLanda7.setLabel([-7,7])
 GsB=GsB*(invLanda8*invLanda3*invLanda7)
 ##print'identity', invLanda4.getBlock()*Landa4.getBlock()
 ##print'identity', invLanda3.getBlock()*Landa3.getBlock()
 ##print'identity', invLanda2.getBlock()*Landa2.getBlock()
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)
 

 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])


 
 Gamma[0].putBlock(GsA.getBlock())
 Gamma[1].putBlock(GsB.getBlock())
 Gamma[0]*=(1.00/Gamma[0].norm())
 Gamma[1]*=(1.00/Gamma[1].norm())

 
 ##printGamma[0].printDiagram(), Gamma[1].printDiagram() 

 #Gs[A].permute([-1, 3, 1], 1)
 #bondrm(Gs[A], Ls[B], 0)
 #bondrm(Gs[B], Ls[B], 1) 

def update_ulink_eff(Gamma,Landa,U,D,d_phys):
 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_a.setLabel([0,1,2,3,4])
 Gamma_b.setLabel([5,6,7,8,9])
 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([10,11,0,5])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa4p=copy.copy(Landa[3])


 Landa1.setLabel([3,-3])
 Landa2.setLabel([2,9])
 Landa3.setLabel([1,-1])
 Landa4.setLabel([4,-4])
 Landa5.setLabel([6,-6])
 Landa4p.setLabel([7,-7])
 Landa6.setLabel([8,-8])


 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)


 q_uni=uni10.UniTensor([bdi,bdi,bdi,bdoB])
 r_uni=uni10.UniTensor([bdiB,bdi_pyso,bdo])

 q_uni.setLabel([-4,-1,-3,20])
 r_uni.setLabel([20,0,2])

 
 Left=Gamma_a*(Landa1*Landa3*Landa4)
 Left.permute([-4,-1,-3,0,2],3)
 Left_m=Left.getBlock()
 qr=Left_m.qr()

 q_uni.putBlock(qr[0])
 r_uni.putBlock(qr[1])






 l_uni=uni10.UniTensor([bdi_pys,bdi,bdoB])
 qq_uni=uni10.UniTensor([bdiB,bdo,bdo,bdo])
 l_uni.setLabel([5,9,40])
 qq_uni.setLabel([40,-7,-8,-6])

 
 Right=Gamma_b*(Landa5*Landa4p*Landa6)
 Right.permute([5,9,-7,-8,-6],2)
 Right_m=Right.getBlock()
 lq=Right_m.lq()

 l_uni.putBlock(lq[0])
 qq_uni.putBlock(lq[1])

 
 #qqt_uni=copy.copy( qq_uni)
 #qqt_uni.transpose()
 #print '2',  qq_uni,  l_uni, qq_uni.getBlock()*qqt_uni.getBlock()

 Theta=(r_uni*Landa2*Hamiltonian)*l_uni
 Theta.permute([10,20,11,40],2)

 
 svd1=Theta.getBlock()
 svd = svd1.svd()
 #svd = Theta.getBlock().svd()
 # Truncation
 sv = svd[1]
 
 
 Sum=0
 for i in xrange(sv.row()):
   if (i>=D):
    Sum=sv[i]+Sum
 #print'truncation=', Sum
 norm = sv.resize(D, D)
 norm=norm.norm()
 sv *= 1.00 / norm
 #print'norm=', norm
 
 Landa[1].putBlock(sv) 
 U=uni10.UniTensor([bdi_pys,bdiB,bdo], "U")
 V=uni10.UniTensor([bdi,bdi_pyso,bdoB], "V")
 U.setLabel([10,20,2])
 V.setLabel([9,11,40])
 
 #print svd[0].row()
 U.putBlock(svd[0].resize(svd[0].row(), D))
 V.putBlock(svd[2].resize(D, svd[2].col()))
 GsA=U*q_uni
 GsB=V*qq_uni


 
 invLanda1=inverse(Landa1)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda1.setLabel([-3,3])
 invLanda3.setLabel([-1,1])
 invLanda4.setLabel([-4,4])
 GsA=GsA*(invLanda1*invLanda3*invLanda4)

 invLanda4=inverse(Landa4)
 invLanda6=inverse(Landa6)
 invLanda5=inverse(Landa5)


 invLanda4.setLabel([-7,7])
 invLanda5.setLabel([-6,6])
 invLanda6.setLabel([-8,8])
 GsB=GsB*(invLanda5*invLanda4*invLanda6)
 ##print'identity', invLanda4.getBlock()*Landa4.getBlock()
 ##print'identity', invLanda3.getBlock()*Landa3.getBlock()
 ##print'identity', invLanda2.getBlock()*Landa2.getBlock()
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)
 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])


 
 Gamma[0].putBlock(GsA.getBlock())
 Gamma[1].putBlock(GsB.getBlock())
 Gamma[0]*=(1.00/Gamma[0].norm())
 Gamma[1]*=(1.00/Gamma[1].norm())
 
 ##printGamma[0].printDiagram(), Gamma[1].printDiagram() 

 #Gs[A].permute([-1, 3, 1], 1)
 #bondrm(Gs[A], Ls[B], 0)
 #bondrm(Gs[B], Ls[B], 1)  


def Energy_H(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys):
 
 
 
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

 a.setLabel([13,1,2,3])
 b.setLabel([2,4,5,6])

 A=((((E3*E4)*(E2*a))*(b*E5))*(E1*E6))
 print A[0]
 norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
 print norm[0]


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

 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.setLabel([20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.setLabel([40,-2,-4,-5,-6])



 A=((((E3*E4)*(E2*a_d*a_u))*(b_u*b_d*E5))*(E1*E6))
 print A[0]

 a_u.setLabel([20,13,1,2,3])
 a_d.setLabel([-20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d.setLabel([-40,-2,-4,-5,-6])

 H0=Ham(h)
 H0.setLabel([-20,-40,20,40])


 #print a_u.printDiagram(), b_u.printDiagram() 
 B=((((E3*E4)*(E2*a_d*a_u))*(H0)*(b_u*b_d*E5))*(E1*E6))
 print B[0]/A[0]




 a_uc=copy.copy(a_u)
 a_uc.permute([3,13,1,2,20],3)
 mat=a_uc.getBlock()
 qr=mat.qr()
 
 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)


 q_uni=uni10.UniTensor([bdi,bdi,bdi,bdoB])
 r_uni=uni10.UniTensor([bdiB,bdi_pyso,bdo])

 q_uni.setLabel([3,13,1,100])
 r_uni.setLabel([100,2,20])

 q_uni.putBlock(qr[0])
 r_uni.putBlock(qr[1])


 b_uc=copy.copy(b_u)
 b_uc.permute([40,2,4,5,6],2)
 mat=b_uc.getBlock()
 lq=mat.lq()

 l_uni=uni10.UniTensor([bdi_pys,bdi,bdoB])
 qq_uni=uni10.UniTensor([bdiB,bdo,bdo,bdo])
 l_uni.setLabel([40,2,200])
 qq_uni.setLabel([200,4,5,6])

 

 l_uni.putBlock(lq[0])
 qq_uni.putBlock(lq[1])

 r_uni_d=copy.copy(r_uni)
 r_uni_d.setLabel([-100,-2,-20])

 l_uni_d=copy.copy(l_uni)
 l_uni_d.setLabel([-40,-2,-200])


 q_uni_d=copy.copy(q_uni)
 q_uni_d.setLabel([-3,-13,-1,-100])
 qq_uni_d=copy.copy(qq_uni)
 qq_uni_d.setLabel([-200,-4,-5,-6])
 N=((((E3*E4)*(E2*q_uni*q_uni_d))*(qq_uni*qq_uni_d*E5))*(E1*E6))
 N.permute([100,-100,200,-200],2)

 A=(r_uni*r_uni_d)*H0*(l_uni*l_uni_d)



 #print N.printDiagram(), A.printDiagram()
 B=N*A
 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.setLabel([20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.setLabel([40,-2,-4,-5,-6])



 A=((((E3*E4)*(E2*a_d*a_u))*(b_u*b_d*E5))*(E1*E6))
 print B[0]/A[0]
 
def Energy_V(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys):
 
 
 
 c1.setLabel([0,1])
 Tb4.setLabel([3,2,0])
 E2=c1*Tb4
 E2.permute([3,2,1],2)
 E2.setLabel([10,2,12])
 
 
 c4.setLabel([0,1])
 Ta3.setLabel([1,2,3])
 E4=c4*Ta3
 E4.permute([0,2,3],2)
 E4.setLabel([9,6,8])



 E1=Tb1
 E1.setLabel([12,1,13])
 
 
 E3=Ta4
 E3.setLabel([9,5,10])
 
 

 c2.setLabel([5,6])
 Ta1.setLabel([21,20,5])
 Ta2.setLabel([22,19,6])
 b.setLabel([4,23,19,20])
 E6=(((c2*Ta2)*Ta1)*b)
 E6.combineBond([23,22])
 E6.permute([21,4,23],2)
 E6.setLabel([13,3,11])
 
 c3.setLabel([7,8])
 Tb2.setLabel([7,15,22])
 Tb3.setLabel([13,14,8])
 d.setLabel([16,14,15,23])
 E5=(((c3*Tb3)*Tb2)*d)
 E5.combineBond([23,22])
 E5.permute([13,16,23],2)
 E5.setLabel([8,7,11])

 a.setLabel([2,4,3,1])
 c.setLabel([5,6,7,4])

 A=((((E3*E4)*(E2*a))*(c*E5))*(E1*E6))
 print A[0]
 norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
 print norm[0]


 bd=uni10.Bond(E1.bond()[1].type(), D)
 A=uni10.UniTensor([E1.bond()[0],bd,bd,E1.bond()[2]])
 A.putBlock(E1.getBlock())
 E1=copy.copy(A)
 E1.setLabel([12,1,-1,13])
 
 bd=uni10.Bond(E2.bond()[1].type(), D)
 A=uni10.UniTensor([E2.bond()[0],bd,bd,E2.bond()[2]])
 A.putBlock(E2.getBlock())
 E2=copy.copy(A)
 E2.setLabel([10,2,-2,12])


 bd=uni10.Bond(E3.bond()[1].type(), D)
 A=uni10.UniTensor([E3.bond()[0],bd,bd,E3.bond()[2]])
 A.putBlock(E3.getBlock())
 E3=copy.copy(A)
 E3.setLabel([9,5,-5,10])


 bd=uni10.Bond(E4.bond()[1].type(), D)
 A=uni10.UniTensor([E4.bond()[0],bd,bd,E4.bond()[2]])
 A.putBlock(E4.getBlock())
 E4=copy.copy(A)
 E4.setLabel([9,6,-6,8])


 bd=uni10.Bond(E5.bond()[1].type(), D)
 A=uni10.UniTensor([E5.bond()[0],bd,bd,E5.bond()[2]])
 A.putBlock(E5.getBlock())
 E5=copy.copy(A)
 E5.setLabel([8,7,-7,11])

 bd=uni10.Bond(E6.bond()[1].type(), D)
 A=uni10.UniTensor([E6.bond()[0],bd,bd,E6.bond()[2]])
 A.putBlock(E6.getBlock())
 E6=copy.copy(A)
 E6.setLabel([13,3,-3,11])

 a_u.setLabel([20,2,4,3,1])
 a_d=copy.copy(a_u)
 a_d.setLabel([20,-2,-4,-3,-1])

 c_u.setLabel([40,5,6,7,4])
 c_d=copy.copy(c_u)
 c_d.setLabel([40,-5,-6,-7,-4])



 A=((((E3*E4)*(E2*a_d*a_u))*(c_u*c_d*E5))*(E1*E6))
 print A[0]

 a_u.setLabel([20,2,4,3,1])
 a_d.setLabel([-20,-2,-4,-3,-1])

 c_u.setLabel([40,5,6,7,4])
 c_d.setLabel([-40,-5,-6,-7,-4])

 H0=Ham(h)
 H0.setLabel([-20,-40,20,40])


 #print a_u.printDiagram(), b_u.printDiagram() 
 B=((((E3*E4)*(E2*a_d*a_u))*(H0)*(c_u*c_d*E5))*(E1*E6))
 print B[0]/A[0]




 a_uc=copy.copy(a_u)
 a_uc.permute([2,1,3,4,20],3)
 mat=a_uc.getBlock()
 qr=mat.qr()
 
 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)


 q_uni=uni10.UniTensor([bdi,bdi,bdi,bdoB])
 r_uni=uni10.UniTensor([bdiB,bdi_pyso,bdo])

 q_uni.setLabel([2,1,3,100])
 r_uni.setLabel([100,4,20])

 q_uni.putBlock(qr[0])
 r_uni.putBlock(qr[1])


 c_uc=copy.copy(c_u)
 c_uc.permute([40,4,5,6,7],2)
 mat=c_uc.getBlock()
 lq=mat.lq()

 l_uni=uni10.UniTensor([bdi_pys,bdi,bdoB])
 qq_uni=uni10.UniTensor([bdiB,bdo,bdo,bdo])
 l_uni.setLabel([40,4,200])
 qq_uni.setLabel([200,5,6,7])

 

 l_uni.putBlock(lq[0])
 qq_uni.putBlock(lq[1])

 r_uni_d=copy.copy(r_uni)
 r_uni_d.setLabel([-100,-4,-20])

 l_uni_d=copy.copy(l_uni)
 l_uni_d.setLabel([-40,-4,-200])


 q_uni_d=copy.copy(q_uni)
 q_uni_d.setLabel([-2,-1,-3,-100])
 qq_uni_d=copy.copy(qq_uni)
 qq_uni_d.setLabel([-200,-5,-6,-7])
 N=((((E3*E4)*(E2*q_uni*q_uni_d))*(qq_uni*qq_uni_d*E5))*(E1*E6))
 N.permute([100,-100,200,-200],2)

 A=(r_uni*r_uni_d)*H0*(l_uni*l_uni_d)



 #print N.printDiagram(), A.printDiagram()
 B=N*A
 a_u.setLabel([20,2,4,3,1])
 a_d=copy.copy(a_u)
 a_d.setLabel([20,-2,-4,-3,-1])

 c_u.setLabel([40,5,6,7,4])
 c_d=copy.copy(c_u)
 c_d.setLabel([40,-5,-6,-7,-4])





 A=((((E3*E4)*(E2*a_d*a_u))*(c_u*c_d*E5))*(E1*E6))
 print B[0]/A[0]


