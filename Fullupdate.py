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

def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,N_iter,delta,h,Env,Gauge,Positive):

 Truncation=[0]
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)


 print 'Ham field=', h
 H0=basic.transverseIsing(h);
 U = uni10.UniTensor(H0.bond(), "U");
 U.putBlock(uni10.takeExp(-delta, H0.getBlock()))
 #print U

 E_0=1.0
 E_1=2.0
 for i in xrange(1,800):

  delta, N_iter=basic.Def_deltaNiter(i)
  delta, N_iter=basic.Def_deltaNiter_less(i)
  if delta is 0: break;
  print 'delta =', delta
  print "N_iter=", N_iter

  U.putBlock(uni10.takeExp(-delta, H0.getBlock()))

  for q in xrange(N_iter):

   a_u, b_u, a, b=basic.Var_H(a_u,b_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)

   c_u, a_u, c, a=basic.Var_V(c_u,a_u,a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)


   c_u, d_u, c, d=basic.Var_H(c_u,d_u,c,d,a,b,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)

   a_u, c_u, a, c=basic.Var_V(a_u,c_u,c,d,a,b,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)


   b_u, a_u, b, a=basic.Var_H(b_u,a_u,b,a,d,c,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)
   

   d_u, b_u, d, b=basic.Var_V(d_u,b_u,b,a,d,c,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)
   

   d_u, c_u, d, c=basic.Var_H(d_u,c_u,d,c,b,a,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)

   


   b_u, d_u, b, d=basic.Var_V(b_u,d_u,d,c,b,a,c1,c2,c3,c4,Ta1,Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,U,d_phys,chi,Gauge,Positive)

   az=basic.magnetization(a_u)
   bz=basic.magnetization(b_u)
   cz=basic.magnetization(c_u)
   dz=basic.magnetization(d_u)
   c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, Truncation=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)

   print  a.norm(), b.norm(), c.norm(), d.norm()
   print  c1.norm(), c2.norm(), c3.norm(), c4.norm()
   print  Ta1.norm(), Ta2.norm(), Ta3.norm(), Ta4.norm()
   print  Tb1.norm(), Tb2.norm(), Tb3.norm(), Tb4.norm()
   norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
   print  norm[0]
   z_value=basic.z_value(a,b,c,d,az,bz,cz,dz,chi,D*D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
   print 'z_value=', z_value
   E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
   print 'E_toal=', E_value
   
   E_0=E_1
   E_1=E_value
   print 'E_diff=', abs((E_0-E_1) / E_0) , 'Num_iter=', q
   if ( abs((E_0-E_1) / E_0) ) < 1.00e-6: 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0) 
    break;
 
 Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]
 return a_u,b_u,c_u,d_u,a,b,c,d, Env

def Init_env(Env):
 rand_env=copy.copy(Env[0])
 rand_env.randomize()
 c1=Env[0]+(1.00e-5)*rand_env
 c2=Env[1]+(1.00e-5)*rand_env
 c3=Env[2]+(1.00e-5)*rand_env 
 c4=Env[3]+(1.00e-5)*rand_env 
 rand_env=copy.copy(Env[4])
 rand_env.randomize() 
 Ta1=Env[4]+(1.00e-5)*rand_env
 Ta2=Env[5]+(1.00e-5)*rand_env
 Ta3=Env[6]+(1.00e-5)*rand_env
 Ta4=Env[7]+(1.00e-5)*rand_env
 rand_env=copy.copy(Env[8])
 rand_env.randomize() 
 Tb1=Env[8]+(1.00e-5)*rand_env
 Tb2=Env[9]+(1.00e-5)*rand_env
 Tb3=Env[10]+(1.00e-5)*rand_env
 Tb4=Env[11]+(1.00e-5)*rand_env
 return  c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4


