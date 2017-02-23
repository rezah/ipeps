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

def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,delta,h,Env,Gauge,Positive,Corner_method,N_iterF,Acc_E,Steps,Model):
 Steps_copy=copy.copy(Steps)
 Truncation=[0]


 if Model is "Ising":
  H0=basic.transverseIsing_Z2(h,d_phys)
 if Model is "Heisenberg":
  H0=basic.Heisenberg_Z2(h,d_phys)
 U = uni10.UniTensor(H0.bond(), "U");

 E_0=1.0
 E_1=2.0
 for i in xrange(1,800):

  delta=1.00/pow(10.00,i)
  if delta>=1.0e-1:
   N_iter=N_iterF
  if delta<=1.0e-1 and delta>1.0e-3:
   N_iter=N_iterF
  if delta<1.0e-3  and delta>1.0e-5:
   N_iter=N_iterF
  if delta<1.0e-5:
   break
  
  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
      U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))


  for q in xrange(N_iter):

   #t0=time.time()
   a_u, b_u, a, b=basic.Var_H(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"


   #t0=time.time()
   c_u, a_u, c, a=basic.Var_V(c_u,a_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"


   #t0=time.time()
   c_u, d_u, c, d=basic.Var_H(c_u,d_u,c,d,a,b,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   #t0=time.time()
   a_u, c_u, a, c=basic.Var_V(a_u,c_u,c,d,a,b,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   #t0=time.time()
   b_u, a_u, b, a=basic.Var_H(b_u,a_u,b,a,d,c,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   #t0=time.time()
   d_u, b_u, d, b=basic.Var_V(d_u,b_u,b,a,d,c,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"
   
   #t0=time.time()
   d_u, c_u, d, c=basic.Var_H(d_u,c_u,d,c,b,a,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"


   #t0=time.time()
   b_u, d_u, b, d=basic.Var_V(b_u,d_u,d,c,b,a,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"



   E_value=basic.E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
   #print 'E_toal=', E_value
   
   E_0=E_1
   E_1=E_value
   print 'E_diff=', abs((E_0-E_1) / E_0) , 'Num_iter=', q 
   if (( abs((E_0-E_1) / E_0) ) < Acc_E) or ( q is int(N_iter-1)): 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0)
    E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
    print 'E_toal=', E_value   
    break;
 
 Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]
 return a_u,b_u,c_u,d_u,a,b,c,d, Env

def Init_env(Env):
 c1=Env[0]
 c2=Env[1]
 c3=Env[2] 
 c4=Env[3] 
 Ta1=Env[4]
 Ta2=Env[5]
 Ta3=Env[6]
 Ta4=Env[7]
 Tb1=Env[8]
 Tb2=Env[9]
 Tb3=Env[10]
 Tb4=Env[11]
 return  c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4


 


