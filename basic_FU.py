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
import basic
def  Do_optimization(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 Res=1
 Res1=2
 count=0
 Distance_val=[0]
 for q in xrange(20):
  #print 'Dis', Distance_val[0], abs(Res1-Res) / abs(Res)
  rp, rp_d=r_optimum(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
  lp, lp_d=l_optimum(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
  Distance_val=Distance(l, r, lp, rp ,N_uni,U)
  Res=Res1
  Res1=Distance_val[0]
  count+=1
  if count > 50: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val[0]) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    #print 'break, Dis', Distance_val[0], (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break

 return rp, rp_d, lp, lp_d

def  Do_optimization_Full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 
 Res=1
 Res1=2
 count=0
 Distance_val=Distance_val=Distance(l, r, lp, rp ,N_uni,U)
 for q in xrange(10):
  #print 'Dis', Distance_val[0], abs(Res1-Res) / abs(Res)
  rp, rp_d=r_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
  lp, lp_d=l_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
  Distance_val=Distance(l, r, lp, rp ,N_uni,U)
  Res=Res1
  Res1=Distance_val[0]
  count+=1
  if count > 50: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val[0]) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    #print 'break, Dis', Distance_val[0], (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break

 return rp, rp_d, lp, lp_d



def  Do_optimization_Grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 
  Es=Distance(l, r, lp, rp ,N_uni,U)
  Ef=[0]
  E2=[0.0]
  lp_first=copy.copy(lp)
  rp_first=copy.copy(rp)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=[0]
  count=0
  for i in xrange(60):
   count+=1
   E1=Distance(l, r, lp, rp ,N_uni,U)
   Ef=E1
   #print 'E=', E1[0], count
   #print lp, lp_d
   D_r, D_l=Obtain_grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
   D_r=(-1.00)*D_r
   D_l=(-1.00)*D_l
#   A=(D_r.norm()+D_l.norm())*(2.00)
#   D_r=D_r*(1.00/A)
#   D_l=D_l*(1.00/A)
   
   A=D_r*D_r
   B=D_l*D_l
   
   Norm_Z=A[0]+B[0]
   
   #print 'Norm', Norm_Z
   if (E1[0]<E_previous[0]) or (i is 0):
    if (abs(E1[0]) > 1.0e-10):
     if abs((E_previous[0]-E1[0])/E1[0]) < 1.0e-11:
      print 'Differnance Satisfied!', E_previous[0], E1[0], abs((E_previous[0]-E1[0])/E1[0]), i
      break
     else: 
      if abs((E_previous[0]-E1[0])) < 1.0e-11:
       print 'Differnance Satisfied!', E_previous[0], E1[0], abs((E_previous[0]-E1[0])), i
       break
      
   E_previous[0]=E1[0]
   
   if Norm_Z < 1.0e-8:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    count+=1
    rp_tem=rp+(2.00)*Gamma*D_r
    lp_tem=lp+(2.00)*Gamma*D_l
    E2=Distance(l, r, lp_tem, rp_tem ,N_uni,U)
    if E1[0]-E2[0] >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
     
   
   
   Break_loop=1
   while Break_loop is 1:
    count+=1
    rp_tem=rp+(1.00)*Gamma*D_r
    lp_tem=lp+(1.00)*Gamma*D_l
    E2=Distance(l, r, lp_tem, rp_tem ,N_uni,U)
    if abs((0.5)*Norm_Z*Gamma) <1.0e-15 or  abs(E1[0]-E2[0])<1.0e-15 :
     rp=rp+(1.00)*Gamma*D_r
     lp=lp+(1.00)*Gamma*D_l
     rp_d=copy.copy(rp)
     lp_d=copy.copy(lp)
     rp_d.setLabel([-3,-20,-1])
     lp_d.setLabel([-2,-40,-3])
     break
     
    if E1[0]-E2[0] < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   rp=rp+(1.00)*Gamma*D_r
   lp=lp+(1.00)*Gamma*D_l
   rp_d=copy.copy(rp)
   lp_d=copy.copy(lp)
   rp_d.setLabel([-3,-20,-1])
   lp_d.setLabel([-2,-40,-3])


  if(Ef[0] > Es[0]):
   print 'SD method, Fail, f<s', Ef[0], Es[0] 
   lp_first=copy.copy(lp_first)
   rp_first=copy.copy(rp_first)
  

  lp_d=copy.copy(lp)
  rp_d=copy.copy(rp)

  rp_d.setLabel([-3,-20,-1])
  lp_d.setLabel([-2,-40,-3])

 
  return rp, rp_d, lp, lp_d

 
def Obtain_grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 Iden=uni10.UniTensor(U.bond())
 d=U.bond()[0].dim()*U.bond()[1].dim()
 matrix_Iden=uni10.Matrix(d, d)
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])



 A2=(((lp*lp_d*rp)*N_uni)*Iden)
 A2.permute([-3,-20,-1],0)
 D_r=copy.copy(A2)
 A2=(((lp*lp_d*rp_d)*N_uni)*Iden)
 A2.permute([3,20,1],0)
 D_r=D_r+A2
 A3=((r)*U*(l*lp_d))*N_uni
 A3.permute([-3,-20,-1],0)
 D_r=D_r+(-1.00)*A3
 A3p=((r_d)*U*(l_d*lp))*N_uni
 A3p.permute([3,20,1],0)
 D_r=D_r+(-1.00)*A3p
 D_r.setLabel([3,20,1])
 D_r.permute([3,20,1],2)


 A2=(((rp*rp_d*lp)*N_uni)*Iden)
 A2.permute([-2,-40,-3],0)
 D_l=copy.copy(A2)
 A2=(((rp*rp_d*lp_d)*N_uni)*Iden)
 A2.permute([2,40,3],0)
 D_l=A2+D_l

 A3=((l)*U*(r*rp_d))*N_uni
 A3.permute([-2,-40,-3],0)
 D_l=D_l+(-1.00)*A3
 A3p=((l_d)*U*(rp*r_d))*N_uni
 A3p.permute([2,40,3],0)
 D_l=D_l+(-1.00)*A3p

 D_l.setLabel([2,40,3])
 D_l.permute([2,40,3],2)

 
 return D_r, D_l

def  r_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,30,50])
 #print H1, U
 Iden=uni10.UniTensor(U.bond())
 d=U.bond()[0].dim()*U.bond()[1].dim()
 matrix_Iden=uni10.Matrix(d, d)
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])
 A2=(((lp*lp_d)*N_uni)*Iden)
 A2.permute([3,20,1,-3,-20,-1],3)
 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 
 A3=((r)*U*(l*lp_d))*N_uni
 A3.permute([-3,-20,-1],0)
 A3p=((r_d)*U*(l_d*lp))*N_uni
 A3p.permute([3,20,1],0)
 A3=A3+A3p



 svd=A2.getBlock().svd()
 Landa=inv(svd[1])
 #print Landa.getBlock()*svd[1]
 v=copy.copy(svd[2])
 v.transpose()
 u=svd[0]
 u.transpose()
 s=Landa.getBlock()
 A2_inv=v*s*u
 A2_mat=A2.getBlock()
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv

 A=A3.getBlock()
 A=A*A2_inv
 A3.putBlock(A)
 A3.setLabel([3,20,1])
 A3.permute([3,20,1],2)

 rf=copy.copy(A3)
 rf_d=copy.copy(rf)
 rf_d.setLabel([-3,-20,-1])

 return rf, rf_d

def  l_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 Bond_val=U.bond()[0].dim()*U.bond()[1].dim()
 matrix_Iden=uni10.Matrix(Bond_val, Bond_val)
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])
 A2=(((rp*rp_d)*N_uni)*Iden)
 A2.permute([2,40,3,-2,-40,-3],3)
 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans




 A3=((l)*U*(r*rp_d))*N_uni
 A3.permute([-2,-40,-3],0)
 A3p=((l_d)*U*(rp*r_d))*N_uni
 A3p.permute([2,40,3],0)
 A3=A3+A3p

 
 
 
 svd=A2.getBlock().svd()
 Landa=inv(svd[1])
 #print Landa.getBlock()*svd[1]
 v=copy.copy(svd[2])
 v.transpose()
 u=svd[0]
 u.transpose()
 s=Landa.getBlock()
 A2_inv=v*s*u
 A2_mat=A2.getBlock()
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance=', distance_iden_val  
 A=A3.getBlock()
 A=A*A2_inv
 A3.putBlock(A)
 A3.setLabel([2,40,3])
 A3.permute([2,40,3],2)


 lf=copy.copy(A3)
 lf_d=copy.copy(lf)
 lf_d.setLabel([-2,-40,-3])

 return lf, lf_d



def  l_optimum(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 Bond_val=U.bond()[0].dim()*U.bond()[1].dim()
 matrix_Iden=uni10.Matrix(Bond_val, Bond_val)
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])
 A2=(((rp*rp_d)*N_uni)*Iden)
 A2.permute([2,40,3,-2,-40,-3],3)
 A3=((l)*U*(r*rp_d))*N_uni
 A3.permute([-2,-40,-3],0)
 svd=A2.getBlock().svd()
 Landa=inv(svd[1])
 #print Landa.getBlock()*svd[1]
 v=copy.copy(svd[2])
 v.transpose()
 u=svd[0]
 u.transpose()
 s=Landa.getBlock()
 A2_inv=v*s*u
 A2_mat=A2.getBlock()
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance=', distance_iden_val  
 A=A3.getBlock()
 A=A*A2_inv
 A3.putBlock(A)
 A3.setLabel([2,40,3])
 A3.permute([2,40,3],2)


 lf=copy.copy(A3)
 lf_d=copy.copy(lf)
 lf_d.setLabel([-2,-40,-3])

 return lf, lf_d


def  r_optimum(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,30,50])

 Iden=uni10.UniTensor(U.bond())
 d=U.bond()[0].dim()*U.bond()[1].dim()
 matrix_Iden=uni10.Matrix(d, d)
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])
 A2=(((lp*lp_d)*N_uni)*Iden)
 A2.permute([3,20,1,-3,-20,-1],3)


 A3=((r)*U*(l*lp_d))*N_uni
 A3.permute([-3,-20,-1],0)
 svd=A2.getBlock().svd()
 Landa=inv(svd[1])
 #print Landa.getBlock()*svd[1]
 v=copy.copy(svd[2])
 v.transpose()
 u=svd[0]
 u.transpose()
 s=Landa.getBlock()
 A2_inv=v*s*u
 A2_mat=A2.getBlock()
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv

 A=A3.getBlock()
 A=A*A2_inv
 A3.putBlock(A)
 A3.setLabel([3,20,1])
 A3.permute([3,20,1],2)

 rf=copy.copy(A3)
 rf_d=copy.copy(rf)
 rf_d.setLabel([-3,-20,-1])

 return rf, rf_d

def inv(svd):
 bdi=uni10.Bond(uni10.BD_IN, svd.col())
 bdo=uni10.Bond(uni10.BD_OUT, svd.row())
 Landa=uni10.UniTensor([bdi,bdo])
 Landa.putBlock(svd)
 Landa=basic.inverse(Landa)
 return Landa
 
def distance_iden(A2,A2_inv):
 A=A2*A2_inv
 d=int(A.row())
 matrixI=uni10.Matrix(d, d)
 for i in xrange(d):
  for j in xrange(d):
   if(i==j):
    matrixI[i*d+j]=1.0
   else:
    matrixI[i*d+j]=0.0
 sum=0 
 for i in xrange(d):
  for j in xrange(d):
   sum=abs(A[i*d+j]+(-1.00)*matrixI[i*d+j])
 return sum
 
 
def  recover(l, r, q,qq):
 a_u=q*r
 a_u.permute([20,4,5,3,2],3)
 a_u.setLabel([0,1,2,3,4])
 #print a_ut.similar(a_u)#,a_ut.printDiagram() 
 
 b_u=qq*l
 b_u.permute([40,3,4,5,6],3)
 b_u.setLabel([0,1,2,3,4])
 #print b_ut.similar(b_u) 
 return a_u, b_u





def  test_energy_lr(N_uni, l, l_d, r, r_d, q,qq,U,E1, E2, E3, E4, E5,E6,a_u,b_u):
 U.setLabel([-20,-40,20,40])
 Iden=uni10.UniTensor(U.bond())
 Bond_val=U.bond()[0].dim()*U.bond()[1].dim()
 matrix_Iden=uni10.Matrix(Bond_val, Bond_val)
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])

 
 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.setLabel([-20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.setLabel([-40,-2,-4,-5,-6])

 B=(((r*r_d)*U)*(l*l_d))*(N_uni)
 A=(((r*r_d)*Iden)*(l*l_d))*(N_uni)
 print 'E=', B[0]/A[0]
 
 
 a_ut=q*r
 a_ut.permute([20,4,5,3,2],3)
 print a_ut.similar(a_u),distance_two(a_u.getBlock(),a_ut.getBlock()) 
 
 b_ut=qq*l
 b_ut.permute([40,3,4,5,6],3)
 
 print b_ut.similar(b_u), distance_two(b_u.getBlock(),b_ut.getBlock())






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
 print 'Norm=', A[0]

 a_u.setLabel([20,13,1,2,3])
 a_d.setLabel([-20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d.setLabel([-40,-2,-4,-5,-6])


 U.setLabel([-20,-40,20,40])


 #print a_u.printDiagram(), b_u.printDiagram() 
 B=((((E2*(a_u*a_d))*E1)*E3)*U)*(((((b_u*b_d)*E5)*E6)*E4))
 print 'E=', B[0]/A[0]


def test_env(E1, E2, E3, E4, E5,E6, a, b, c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 a.setLabel([13,1,2,3])
 b.setLabel([2,4,5,6])
 A=(((E2*a)*E1)*E3)*((((b*E5)*E6)*E4))
 print 'norm=', A[0]


def Distance(l, r, lp, rp ,N_uni, U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,30,50])
  
 U.setLabel([30,50,20,40])
 H=U*H1
 H.permute([-20,-40,20,40],2)
 #H.setLabel([-20,-40,20,40]) 
 
 H1.setLabel([-20,-40,20,40]) 
 U.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 matrix_Iden=uni10.Matrix(U.bond()[0].dim()*U.bond()[1].dim(), U.bond()[2].dim()*U.bond()[3].dim())
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])
 
 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)

 rp_d.setLabel([-3,-20,-1])
 lp_d.setLabel([-2,-40,-3])

 
 l_d=copy.copy(l)
 r_d=copy.copy(r)

 r_d.setLabel([-3,-20,-1])
 l_d.setLabel([-2,-40,-3])
 
 
 
 A1=((((r*r_d)*H)*(N_uni))*l*l_d)
 A2=((((rp*rp_d)*Iden)*(N_uni))*lp*lp_d)
 A3=((((r*rp_d)*U)*(N_uni))*l*lp_d)
 A4=((((rp*r_d)*H1)*(N_uni))*lp*l_d)
 A=A1+A2+(-1.00)*A3+(-1.00)*A4
 return A

def initialize_lrprime(l, r, l_d, r_d, N_uni):
 lp=copy.copy(l)
 rp=copy.copy(r)
 lp_d=copy.copy(l_d)
 rp_d=copy.copy(r_d)

 rp_d.setLabel([-3,-20,-1])
 lp_d.setLabel([-2,-40,-3])
 
 return  lp, rp, lp_d, rp_d

def initialize_Positiv_lrprime(l, r, l_d, r_d, N_uni, U, D, d_phys, q_u, qq_u,a_u,b_u,Positive):
 U.setLabel([-20,-40,20,40])
 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdiBB=uni10.Bond(uni10.BD_IN, d_phys*D*d_phys*D)
 bdoBB=uni10.Bond(uni10.BD_OUT, d_phys*D*d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 #A=(r*r_d)*U*(l*l_d)*N_uni
 #print 'test', A[0]
 #print N_uni
 N=copy.copy(N_uni)
 N.setLabel([1,-1,2,-2])
 N.permute([-1,-2,1,2], 2)
 N1=copy.copy(N)
 N1.transpose()
 N=(N+N1)*(1.00/2.00)
 if Positive is 'Restrict':
  M=N.getBlock()
  N_2=M*M 
  eig=N_2.eigh()
  #print eig[0]
  e=Sqrt(eig[0])
  #print e
  U_trans=copy.copy(eig[1])
  U_trans.transpose()
  M=U_trans*e*eig[1]
  N.putBlock(M)

 eig=N.getBlock()
 eig=eig.eigh() 
 e=Positiv(eig[0])
 
# print '\n','\n'
# for i1 in xrange(int(e.row())):
#  print e[i1]
# print '\n','\n'
 for q in xrange(int(eig[0].row())):
    e[q]=(e[q]**(1.00/2.00)) 
 X=e*eig[1]
 eig[1].transpose()
 X_tran=eig[1]*e
 X_u=uni10.UniTensor([bdiBB,bdoB,bdoB])
 X_u.putBlock(X)
 
 X_tran_u=uni10.UniTensor([bdiB,bdiB,bdoBB])
 X_tran_u.putBlock(X_tran)

#################################################
 X_u.permute([1,0,2],1)
 #print X_u
 lq1=X_u.getBlock().lq()
 l1=lq1[0]
 #print l1, l1[0] 
 l1_inv=l1.inverse()
 #print l1*l1_inv
 
 X_u.permute([1,0,2],2)
 qr1=X_u.getBlock().qr()
 r1=qr1[1]
 r1_inv=r1.inverse()

 r1_uni_inv=uni10.UniTensor([bdiB,bdoB])
 r1_uni_inv.putBlock(r1_inv)
 r1_uni_inv.setLabel([2,-2])
 qq_u=qq_u*r1_uni_inv
 qq_u.permute([-2,4,5,6],3)
 qq_u.setLabel([2,4,5,6])

 
 
 l1_uni_inv=uni10.UniTensor([bdiB,bdoB])
 l1_uni_inv.putBlock(l1_inv)
 l1_uni_inv.setLabel([-1,1])
 #print q_u.printDiagram()
 q_u=q_u*l1_uni_inv
 q_u.permute([2,4,5,-1],3)
 q_u.setLabel([2,4,5,1])
 
 
 X_u=X_u*l1_uni_inv*r1_uni_inv
 X_u.permute([0,-1,-2],1)
 X_u.setLabel([0,1,2])
 
##############################################
# X_tran_u.permute([0,2,1],1)
# lq2=X_tran_u.getBlock().lq()
# l2=lq2[0]
# l2_inv=l2.inverse()
# #print '1', distance_two(l2_inv, l1_inv)
# X_tran_u.permute([0,2,1],2)
# qr2=X_tran_u.getBlock().qr()
# r2=qr2[1]
# r2_inv=r2.inverse()
# #print '2', distance_two(r2_inv, r1_inv)

# r2_uni_inv=uni10.UniTensor([bdiB,bdoB])
# r2_uni_inv.putBlock(r2_inv)
# r2_uni_inv.setLabel([1,-2])
# 
# l2_uni_inv=uni10.UniTensor([bdiB,bdoB])
# l2_uni_inv.putBlock(l2_inv)
# l2_uni_inv.setLabel([-1,0])
# 
# X_tran_u=X_tran_u*l2_uni_inv*r2_uni_inv
# X_tran_u.permute([-1,-2,2],2)
# X_tran_u.setLabel([-1,-2,0])

 X_tran_u=copy.copy(X_u)
 X_tran_u.transpose()
 X_tran_u.setLabel([-1,-2,0])
 
#################################################
 
 #print l2.similar(l1)  
 #print r2.similar(r1)
# A=copy.copy(X_tran_u)
# A.transpose()
# dis=distance_two(X_u.getBlock(), A.getBlock())
# print 'hi', X_u.similar(A), dis
# 
 N=X_tran_u*X_u
 N.permute([1,-1,2,-2],2)
 r1_uni=uni10.UniTensor([bdiB,bdoB])
 r1_uni.putBlock(r1)

 l1_uni=uni10.UniTensor([bdiB,bdoB])
 l1_uni.putBlock(l1)

 #print l, r
 
 
 l1_uni.setLabel([1,0])
 r=r*l1_uni
 r.permute([3,20,0],2)
 r.setLabel([3,20,1])
 r_d=copy.copy(r)
 r_d.setLabel([-3,-20,-1])
 
 
 r1_uni.setLabel([0,2])
 l=r1_uni*l
 l.permute([0,40,3],2)
 l.setLabel([2,40,3])
 l_d=copy.copy(l)
 l_d.setLabel([-2,-40,-3])
 #A=(r*r_d)*U*(l*l_d)*N
 #print 'test', A[0]
 
# A=X_u*r*l*X_tran_u*r_d*l_d
# #print A, Tes1
# dis=distance_two(A.getBlock(),Tes1.getBlock() )

 lp=copy.copy(l)
 rp=copy.copy(r)
 lp_d=copy.copy(l_d)
 rp_d=copy.copy(r_d)

 rp_d.setLabel([-3,-20,-1])
 lp_d.setLabel([-2,-40,-3])

#######################################################
 
 U.setLabel([-20,-40,20,40])
 lp=copy.copy(l)
 rp=copy.copy(r)
 #lp.setLabel([1,40,-3])
 Teta=U*lp*rp
 Teta.permute([1,-20,2,-40],2)
 svd=Teta.getBlock().svd()
 s=Sqrt(svd[1])
 U=svd[0]*s
 V=s*svd[2]
 U.resize(U.row(),rp.bond()[0].dim() )
 V.resize(rp.bond()[0].dim(), V.col())
 
 rp=uni10.UniTensor([bdiB, bdi_pys, bdo ])
 rp.putBlock(U)
 rp.setLabel([1,20,3])
 rp.permute([3,20,1],2)
 
 lp=uni10.UniTensor([ bdi, bdoB, bdi_pyso])
 lp.putBlock(V)
 lp.setLabel([3,2,40])
 lp.permute([2,40,3],2)
 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)
 rp_d.setLabel([-3,-20,-1])
 lp_d.setLabel([-2,-40,-3])
##########################################################

 a_u=q_u*r
 a_u.permute([20,4,5,3,2],3)
 
 b_u=qq_u*l
 b_u.permute([40,3,4,5,6],3)
 
 return  lp, rp, lp_d, rp_d, N, l, r, l_d, r_d, q_u, qq_u, a_u, b_u 



def Positiv(e):
 d=int(e.row())
 if (e[0] > 0) and (e[d-1] > 0):
  return e

 if (e[0] < 0) and (e[d-1] < 0):
  return (-1.00)*e

 if abs(e[0]) >  abs(e[d-1]) :
   print 'e_small > e_large',  e, e[0],  e[d-1]
   e=(-1.00)*e


 for q in xrange(d):
  if e[q] < 0: e[q]=0;
 return e






def Sqrt(e):
 d=int(e.row())
 
 for q in xrange(d):
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e



def  distance_two(A,B ):
 d=int(A.row()*A.col())
 sum=0.0
 for q in xrange(d):
  sum+=abs(A[q]-B[q])
 return sum

def Equall_Dist(l, r,D,d_phys):
 
 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 lp=copy.copy(l)
 rp=copy.copy(r)
 rp.permute([3,20,1],1)
 lq_r=rp.getBlock().lq()
 
 
 
 lp.permute([2,40,3],2)
 qr_l=lp.getBlock().qr()
 
 
 
 
 #lp.setLabel([1,40,-3])
 Teta=qr_l[1]*lq_r[0]
 svd=Teta.svd()
 s=Sqrt(svd[1])
 U=svd[0]*s
 V=s*svd[2]
 #U.resize(U.row(),rp.bond()[2].dim() )
 #V.resize(rp.bond()[2].dim(), V.col())
 #rp=uni10.UniTensor([bdi, bdi_pys,bdoB])
 
 #print lq_r[1].col(),  U.row()
 A=V*lq_r[1]
 B=qr_l[0]*U
 rp.putBlock(A)
 lp.putBlock(B)
 
 rp.permute([3,20,1],2)
 lp.permute([2,40,3],2)
 
  
 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)
 rp_d.setLabel([-3,-20,-1])
 lp_d.setLabel([-2,-40,-3])
 
 return  lp, rp, lp_d, rp_d




def initialize_SVD_lrprime(l, r, l_d, r_d, N_uni,U,D,d_phys):
 
 bdiB=uni10.Bond(uni10.BD_IN, d_phys*D)
 bdoB=uni10.Bond(uni10.BD_OUT, d_phys*D)
 bdi_pys=uni10.Bond(uni10.BD_IN, d_phys)
 bdi_pyso=uni10.Bond(uni10.BD_OUT, d_phys)
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 U.setLabel([-20,-40,20,40])
 lp=copy.copy(l)
 rp=copy.copy(r)
 #lp.setLabel([1,40,-3])
 Teta=U*lp*rp
 Teta.permute([1,-20,2,-40],2)
 svd=Teta.getBlock().svd()
 s=Sqrt(svd[1])
 U=svd[0]*s
 V=s*svd[2]
 U.resize(U.row(),rp.bond()[0].dim() )
 V.resize(rp.bond()[0].dim(), V.col())
 
 rp=uni10.UniTensor([bdiB, bdi_pys, bdo ])
 rp.putBlock(U)
 rp.setLabel([1,20,3])
 rp.permute([3,20,1],2)
 
 lp=uni10.UniTensor([ bdi, bdoB, bdi_pyso])
 lp.putBlock(V)
 lp.setLabel([3,2,40])
 lp.permute([2,40,3],2)
 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)
 rp_d.setLabel([-3,-20,-1])
 lp_d.setLabel([-2,-40,-3])
 
 return  lp, rp, lp_d, rp_d


 
def Qr_lQ_decom(a_u,b_u, E1, E2, E3, E4, E5,E6,D,d_phys):
 a_uc=copy.copy(a_u)
 a_uc.setLabel([20,13,1,2,3])
 a_uc.permute([3,13,1,20,2],3)
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
 r_uni.setLabel([100,20,2])

 q_uni.putBlock(qr[0])
 r_uni.putBlock(qr[1])


 b_uc=copy.copy(b_u)
 b_uc.setLabel([40,2,4,5,6])
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
 r_uni_d.setLabel([-100,-20,-2])

 l_uni_d=copy.copy(l_uni)
 l_uni_d.setLabel([-40,-2,-200])


 q_uni_d=copy.copy(q_uni)
 q_uni_d.setLabel([-3,-13,-1,-100])
 qq_uni_d=copy.copy(qq_uni)
 qq_uni_d.setLabel([-200,-4,-5,-6])
 N=((((E1*q_uni)*E2)*q_uni_d)*E3)*((((E6*qq_uni)*E5)*qq_uni_d)*E4)
 N.permute([100,-100,200,-200],2)

 N.setLabel([1,-1,2,-2])
 #r_uni.setLabel([100,20,2])
 r_uni.setLabel([1,20,3])
 r_uni.permute([3,20,1],2)
 
 r_uni_d=copy.copy(r_uni)
 r_uni_d.setLabel([-3,-20,-1])

 #l_uni.setLabel([40,2,200])
 l_uni.setLabel([40,3,2])
 l_uni.permute([2,40,3],2)
 
 l_uni_d=copy.copy(l_uni)
 l_uni_d.setLabel([-2,-40,-3])

 #q_uni.setLabel([3,13,1,100])
 q_uni.setLabel([2,4,5,1])
 #qq_uni.setLabel([200,4,5,6])
 qq_uni.setLabel([2,4,5,6])
 
 
 return   N, l_uni, l_uni_d, r_uni, r_uni_d, q_uni,qq_uni
 
 
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
 
 
 
 
def final_test_distance(ap_u, bp_u, a_u, b_u,E1, E2, E3, E4, E5,E6,U,N_uni):


 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,30,50])
  
 U.setLabel([30,50,20,40])
 H=U*H1
 H.permute([-20,-40,20,40],2)
 #H.setLabel([-20,-40,20,40]) 
 
 H1.setLabel([-20,-40,20,40]) 
 U.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 matrix_Iden=uni10.Matrix(U.bond()[0].dim()*U.bond()[1].dim(), U.bond()[2].dim()*U.bond()[3].dim())
 matrix_Iden.identity()
 Iden.putBlock(matrix_Iden)
 Iden.setLabel([-20,-40,20,40])
 
 
 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.setLabel([-20,-13,-1,-2,-3])

 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.setLabel([-40,-2,-4,-5,-6])


 ap_u.setLabel([20,13,1,2,3])
 ap_d=copy.copy(ap_u)
 ap_d.setLabel([-20,-13,-1,-2,-3])

 bp_u.setLabel([40,2,4,5,6])
 bp_d=copy.copy(bp_u)
 bp_d.setLabel([-40,-2,-4,-5,-6])
 #N=(((((E2*)*E1)*E3)*Iden)*((((E5)*E6)*E4)))

# a_u=q*r
# a_u.permute([20,4,5,3,2],3)
# a_u.setLabel([0,1,2,3,4])
# #print a_ut.similar(a_u)#,a_ut.printDiagram() 
# 
# b_u=qq*l
# b_u.permute([40,3,4,5,6],3)
# b_u.setLabel([0,1,2,3,4])
# #print b_ut.similar(b_u) 
 
 A1=(((((E2*(a_u*a_d))*E1)*E3)*H)*(((((b_u*b_d)*E5)*E6)*E4))) 
 A2=((((E2*(ap_u*ap_d))*E1)*E3)*Iden)*(((((bp_u*bp_d)*E5)*E6)*E4))
 A3=((((E2*(a_u*ap_d))*E1)*E3)*U)*(((((b_u*bp_d)*E5)*E6)*E4))
 A4=((((E2*(ap_u*a_d))*E1)*E3)*H1)*(((((bp_u*b_d)*E5)*E6)*E4))

 A=A1+A2+(-1.00)*A3+(-1.00)*A4
 return A
 
 
 
 
 

