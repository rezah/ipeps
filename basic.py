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

def corner_transfer_matrix_onesite(a, az,chi,D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):

 z=copy.copy(az)
 for q in xrange(20):
  c1, c4, Ta4= Move.left(c1,c4,Ta4,Ta1,Ta3,chi,a)
  
  c3, c4 , Ta3= Move.down(c4,c3,Ta3,Ta4,Ta2,chi,a)
  c2, c3 , Ta2= Move.right(c3,c2,Ta2,Ta1,Ta3,chi,a)

  c1, c2 , Ta1= Move.up(c1,c2,Ta1,Ta4,Ta2,chi,a)

  norm0=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,a)
  norm1=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,z)
  print abs(norm1[0])/abs(norm0[0]), norm1[0]/norm0[0]
 return abs(norm1[0])/abs(norm0[0])




def test_converg_CTM(a, b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,chi,Z):


 EnvCTM_net = uni10.Network("Full.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',Z)
 EnvCTM_net.putTensor('b',Z)
 EnvCTM_net.putTensor('a1',a)
 EnvCTM_net.putTensor('b1',b)
 Result_uni10=EnvCTM_net.launch()
 norm=Result_uni10.trace().real
 #print Result_uni10
 
 EnvCTM_net = uni10.Network("Full.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',a)
 EnvCTM_net.putTensor('b',b)
 EnvCTM_net.putTensor('a1',Z)
 EnvCTM_net.putTensor('b1',Z)
 Result_uni10=EnvCTM_net.launch()
 norm1=Result_uni10.trace().real

 print norm1/norm
 #print Z.norm(), a.norm(), b.norm()






def magnetization_value(a, b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,chi,az,bz):
 EnvCTM_net = uni10.Network("EnvCTM.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',a)
 EnvCTM_net.putTensor('b',b)
 EnvCTM_net.putTensor('a1',a)
 EnvCTM_net.putTensor('b1',b)
 Result_uni10=EnvCTM_net.launch()
 norm=Result_uni10.trace().real
 Res=copy.copy(Result_uni10)
 Res.partialTrace(-1,-3)
 Res.permute([-2,-4],1)
 #print'trace0', Res.trace().real
 #print'norm', norm

 

 EnvCTM_net = uni10.Network("EnvCTM.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',az)
 EnvCTM_net.putTensor('b',b)
 EnvCTM_net.putTensor('a1',a)
 EnvCTM_net.putTensor('b1',b)
 Result_uni10=EnvCTM_net.launch()
 z1=Result_uni10.trace().real
 ##print'z1', z1

 EnvCTM_net = uni10.Network("EnvCTM.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',a)
 EnvCTM_net.putTensor('b',bz)
 EnvCTM_net.putTensor('a1',a)
 EnvCTM_net.putTensor('b1',b)
 Result_uni10=EnvCTM_net.launch()
 z2=Result_uni10.trace().real
 ##print'z2', z2

 EnvCTM_net = uni10.Network("EnvCTM.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',a)
 EnvCTM_net.putTensor('b',b)
 EnvCTM_net.putTensor('a1',az)
 EnvCTM_net.putTensor('b1',b)
 Result_uni10=EnvCTM_net.launch()
 z3=Result_uni10.trace().real
 ##print'z3', z3
 
 
 EnvCTM_net = uni10.Network("EnvCTM.net")
 EnvCTM_net.putTensor('c1',c1)
 EnvCTM_net.putTensor('c2',c2)
 EnvCTM_net.putTensor('c3',c3)
 EnvCTM_net.putTensor('c4',c4)
 EnvCTM_net.putTensor('Ta1',Ta1)
 EnvCTM_net.putTensor('Ta2',Ta2)
 EnvCTM_net.putTensor('Ta3',Ta3)
 EnvCTM_net.putTensor('Ta4',Ta4)
 EnvCTM_net.putTensor('Tb1',Tb1)
 EnvCTM_net.putTensor('Tb2',Tb2)
 EnvCTM_net.putTensor('Tb3',Tb3)
 EnvCTM_net.putTensor('Tb4',Tb4)
 EnvCTM_net.putTensor('a',a)
 EnvCTM_net.putTensor('b',b)
 EnvCTM_net.putTensor('a1',a)
 EnvCTM_net.putTensor('b1',bz)
 Result_uni10=EnvCTM_net.launch()
 z4=Result_uni10.trace().real
 ##print'z4', z4
 z5=z1+z2+z3+z4
 z5*=(1.00/(4.00*norm))
 ##printz5
 return z5



def magnetization(a_u, b_u):
 a_u.setLabel([0,1,2,3,4])
 a_uc=copy.copy(a_u)
 M=a_uc.getBlock()
 M.conj()
 a_uc.putBlock(M)
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

 b_u.setLabel([0,1,2,3,4])
 b_uc=copy.copy(b_u)
 b_uc.setLabel([5,-1,-2,-3,-4])
 M=b_uc.getBlock()
 M.conj()
 b_uc.putBlock(M)

 sz = matSz()
 bdi = uni10.Bond(uni10.BD_IN, 2);
 bdo = uni10.Bond(uni10.BD_OUT, 2);
 Z =  uni10.UniTensor([bdi, bdo], "Z");
 Z.putBlock(sz)   
 Z.setLabel([5,0])
 result1=b_uc*Z*b_u
 result1.combineBond([1,-1])
 result1.combineBond([2,-2])
 result1.combineBond([3,-3])
 result1.combineBond([4,-4])
 result1.permute([1,2,3,4], 2)
 return result, result1

def permute(a, b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,1,2,3])
 a.permute([3,2,1,0],2)
 a.setLabel([0,1,2,3])
 ##print'a', a
 b.setLabel([0,1,2,3])
 b.permute([3,2,1,0],2)
 b.setLabel([0,1,2,3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([1,0],1)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([1,0],1)
 c4.setLabel([0,1])



 Ta1.setLabel([0,1,2])
 Ta1.permute([2,1,0],2)
 Ta1.setLabel([0,1,2])
 
 Ta2.setLabel([0,1,2])
 Ta2.permute([2,1,0],2)
 Ta2.setLabel([0,1,2])
 
 Ta3.setLabel([0,1,2])
 Ta3.permute([2,1,0],2)
 Ta3.setLabel([0,1,2])

 Ta4.setLabel([0,1,2])
 Ta4.permute([2,1,0],2)
 Ta4.setLabel([0,1,2])

 Tb1.setLabel([0,1,2])
 Tb1.permute([2,1,0],2)
 Tb1.setLabel([0,1,2])
 
 Tb2.setLabel([0,1,2])
 Tb2.permute([2,1,0],2)
 Tb2.setLabel([0,1,2])
 
 Tb3.setLabel([0,1,2])
 Tb3.permute([2,1,0],2)
 Tb3.setLabel([0,1,2])

 Tb4.setLabel([0,1,2])
 Tb4.permute([2,1,0],2)
 Tb4.setLabel([0,1,2])







def makeTab(a, b, chi, D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_IN, D*D)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor(uni10.CTYPE,[bdi, bdi1, bdo])
 Tem0.randomize()
 Tem1=uni10.UniTensor(uni10.CTYPE,[bdi, bdi1, bdo])
 Tem1.randomize()
 
 return Tem0, Tem1


def makec1(a, b, chi):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 
 Tem0=uni10.UniTensor(uni10.CTYPE,[bdi, bdo])
 Tem0.randomize()
 Tem1=uni10.UniTensor(uni10.CTYPE,[bdi, bdo])
 Tem1.randomize()
 Tem2=uni10.UniTensor(uni10.CTYPE,[bdi, bdo])
 Tem2.randomize()
 Tem3=uni10.UniTensor(uni10.CTYPE,[bdi, bdo])
 Tem3.randomize()
 
 return Tem0,Tem1,Tem2,Tem3
 
def makeab(Landa,Gamma_a, Gamma_b):

 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa)) ]
 Landa_sq=sqrt(Landa_cp)

 Landa_sq[0].setLabel([1,-1])
 Landa_sq[1].setLabel([2,-2])
 Landa_sq[2].setLabel([3,-3])
 Landa_sq[3].setLabel([4,-4])

 a_tempu=copy.copy(Gamma_a)
 a_tempd=copy.copy(Gamma_a)
 M=a_tempd.getBlock()
 M.conj()
 a_tempd.putBlock(M)
 a_tempu=a_tempu*(Landa_sq[0]*Landa_sq[1]*Landa_sq[2]*Landa_sq[3])
 a_tempd=a_tempd*(Landa_sq[0]*Landa_sq[1]*Landa_sq[2]*Landa_sq[3])
 a_tempu.permute([0,-1,-2,-3,-4],3)
 a_tempd.permute([0,-1,-2,-3,-4],3)

 a_tempu.setLabel([0,1,2,3,4])
 a_tempd.setLabel([0,-1,-2,-3,-4])
 a=a_tempu*a_tempd
 a.combineBond([1,-1])
 a.combineBond([2,-2])
 a.combineBond([3,-3])
 a.combineBond([4,-4])
 a.permute([1,2,3,4],2)

 b_tempu=copy.copy(Gamma_b)
 b_tempd=copy.copy(Gamma_b)
 M=b_tempd.getBlock()
 M.conj()
 b_tempd.putBlock(M)
 b_tempu=b_tempu*(Landa_sq[0]*Landa_sq[1]*Landa_sq[2]*Landa_sq[3])
 b_tempd=b_tempd*(Landa_sq[0]*Landa_sq[1]*Landa_sq[2]*Landa_sq[3])
 b_tempu.permute([0,-1,-2,-3,-4],3)
 b_tempd.permute([0,-1,-2,-3,-4],3)

 b_tempu.setLabel([0,1,2,3,4])
 b_tempd.setLabel([0,-1,-2,-3,-4])
 b=b_tempu*b_tempd
 b.combineBond([1,-1])
 b.combineBond([2,-2])
 b.combineBond([3,-3])
 b.combineBond([4,-4])
 b.permute([1,2,3,4],2)
 return a_tempu,b_tempu, a, b
 
 
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
 invLanda2=uni10.UniTensor(uni10.CTYPE,Landa2.bond())
 invL2 = uni10.CMatrix(Landa2.bond()[0].dim(), Landa2.bond()[1].dim())
 D=Landa2.bond()[0].dim()
 for i in xrange(Landa2.bond()[0].dim()):
   for j in xrange(Landa2.bond()[0].dim()):
    invL2[i*D+j] = 0 if ((Landa2[i*D+j].real) < 1.0e-12) else (1.00 / (Landa2[i*D+j].real))
 invLanda2.putBlock(invL2)
 return invLanda2
 

def Initialize_function(Gamma,Landa):
 for i in xrange(len(Gamma)):
  Gamma[i].randomize(uni10.CTYPE)
 D=Gamma[0].bond(1).dim()
 matrix_Iden=uni10.CMatrix(D, D)
 matrix_Iden.identity(uni10.CTYPE)
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

 
def update_rlink(Gamma,Landa,U,D,d):
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
 Landa1.setLabel([3,6])
 Landa2.setLabel([2,-2])
 Landa3.setLabel([1,-1])
 Landa4.setLabel([4,-4])
 Landa1p=copy.copy(Landa[0])
 Landa2p=copy.copy(Landa[1])
 Landa3p=copy.copy(Landa[2])
 Landa4p=copy.copy(Landa[3])
 Landa2p.setLabel([9,-9])
 Landa3p.setLabel([8,-8])
 Landa4p.setLabel([7,-7])
 Tem_tensor=Gamma_a*(Landa1*Landa2*Landa3*Landa4)
 Tem_tensor2=Gamma_b*(Landa2p*Landa3p*Landa4p)
 
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
 norm = sv.resize(D, D).norm()
 sv *= 1.00 / norm
 #print'norm=', norm
 
 Landa[0].putBlock(sv) 
 bdi_pys = uni10.Bond(uni10.BD_IN, d)
 bdi_pyso = uni10.Bond(uni10.BD_OUT, d)
 bdi = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, D)
 
 GsA=uni10.UniTensor(uni10.CTYPE,[bdi_pys,bdi,bdi,bdi,bdo], "GsA")
 GsB=uni10.UniTensor(uni10.CTYPE,[bdi,bdi_pyso,bdo,bdo,bdo], "GsB")
 
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
 invLanda2.setLabel([-9,9])
 invLanda3.setLabel([-8,8])
 invLanda4.setLabel([-7,7])
 GsB=GsB*(invLanda2*invLanda3*invLanda4)
 ##print'identity', invLanda4.getBlock()*Landa4.getBlock()
 ##print'identity', invLanda3.getBlock()*Landa3.getBlock()
 ##print'identity', invLanda2.getBlock()*Landa2.getBlock()
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)


 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])


 
 Gamma[0].putBlock(GsA.getBlock())
 Gamma[1].putBlock(GsB.getBlock())
 
 
 ##printGamma[0].printDiagram(), Gamma[1].printDiagram() 
 


 
 
 #Gs[A].permute([-1, 3, 1], 1)
 #bondrm(Gs[A], Ls[B], 0)
 #bondrm(Gs[B], Ls[B], 1)  


