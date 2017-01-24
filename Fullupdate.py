import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time
import basic
import Move


def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,N_iter,delta,h):

 c1, c2,c3,c4=basic.makec1(chi,D*D)
 Ta1, Tb1=basic.makeTab(chi,D*D)
 Ta2, Tb2=basic.makeTab(chi,D*D)
 Ta3, Tb3=basic.makeTab(chi,D*D)
 Ta4, Tb4=basic.makeTab( chi,D*D)
 az=basic.magnetization(a_u)
 bz=basic.magnetization(b_u)
 cz=basic.magnetization(c_u)
 dz=basic.magnetization(d_u)
 Truncation=[0]
 
 
 c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, Truncation=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)
 E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
 z_value=basic.z_value(a,b,c,d,az,bz,cz,dz,chi,D*D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 print 'E_toal=', E_value
 print 'z_value=', z_value


 print 'Ham field=', h
 H0=basic.transverseIsing(h);
 U = uni10.UniTensor(H0.bond(), "U");
 U.putBlock(uni10.takeExp(-delta, H0.getBlock()))
 #print U

 for q in xrange(0):

  #print '1'
  a_u, b_u, a, b=basic.Var_H(a_u,b_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
  #print a, b, c, d
  #print '2'
  b_u, a_u, b, a=basic.Var_H(b_u,a_u,b,a,d,c,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
  #print '3'
  c_u, d_u, c, d=basic.Var_H(c_u,d_u,c,d,a,b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
  #print '4'
  d_u, c_u, d, c=basic.Var_H(d_u,c_u,d,c,b,a,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)

  #print '5'
  c_u, a_u, c, a=basic.Var_V(c_u,a_u,a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)

  #print '6'
  a_u, c_u, a, c=basic.Var_V(a_u,c_u,c,d,a,b,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
  #print '7'
  d_u, b_u, d, b=basic.Var_V(d_u,b_u,b,a,d,c,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
  #print '8'
  b_u, d_u, b, d=basic.Var_V(b_u,d_u,d,c,b,a,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)

  az=basic.magnetization(a_u)
  bz=basic.magnetization(b_u)
  cz=basic.magnetization(c_u)
  dz=basic.magnetization(d_u)
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, Truncation=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)

  z_value=basic.z_value(a,b,c,d,az,bz,cz,dz,chi,D*D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
  print 'z_value=', z_value
  E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
  print 'E_toal=', E_value

 E_0=1.0
 E_1=2.0
 for i in xrange(1,60):
  delta=1.00/pow(10.00,i)
  print 'delta =', delta
  if delta>1.0e-1:
   N_iter=40
  if delta<=1.0e-1 and delta>1.0e-3:
   N_iter=40
  if delta<1.0e-3  and delta>1.0e-4:
   N_iter=40
  if delta<1.0e-4:
   break

  print "N_iter=", N_iter

  U.putBlock(uni10.takeExp(-delta, H0.getBlock()))

  for q in xrange(N_iter):
   #print '1'
   a_u, b_u, a, b=basic.Var_H(a_u,b_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '2'
   b_u, a_u, b, a=basic.Var_H(b_u,a_u,b,a,d,c,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '3'
   c_u, d_u, c, d=basic.Var_H(c_u,d_u,c,d,a,b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '4'
   d_u, c_u, d, c=basic.Var_H(d_u,c_u,d,c,b,a,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '5'
   c_u, a_u, c, a=basic.Var_V(c_u,a_u,a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '6'
   a_u, c_u, a, c=basic.Var_V(a_u,c_u,c,d,a,b,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '7'
   d_u, b_u, d, b=basic.Var_V(d_u,b_u,b,a,d,c,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)
   #print '8'

   b_u, d_u, b, d=basic.Var_V(b_u,d_u,d,c,b,a,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi)

   az=basic.magnetization(a_u)
   bz=basic.magnetization(b_u)
   cz=basic.magnetization(c_u)
   dz=basic.magnetization(d_u)
   c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, Truncation=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)

   print  a.norm(),b.norm(), c.norm(), d.norm()
   print  c1.norm(), c2.norm(), c3.norm(), c4.norm()
   print  Ta1.norm(), Ta2.norm(),Ta3.norm(), Ta4.norm()
   print  Tb1.norm(), Tb2.norm(),Tb3.norm(), Tb4.norm()
   print  Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)

   z_value=basic.z_value(a,b,c,d,az,bz,cz,dz,chi,D*D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
   print 'z_value=', z_value
   E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
   print 'E_toal=', E_value
   
   E_0=E_1
   E_1=E_value
   print 'E_diff=', abs((E_0-E_1) / E_0) 
   if ( abs((E_0-E_1) / E_0) ) < 1.00e-6: 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0) 
    break;

 return a_u,b_u,c_u,d_u,a,b,c,d
 
 
