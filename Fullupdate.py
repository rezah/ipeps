import pyUni10 as uni10
import sys
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random
import copy
import time
import basic
import Move
import basicA
import basicB
import basicC

def Full_Update(a_u,b_u,a,b,chi,d_phys,D,delta,h,Env,Env1,Gauge,Corner_method,N_iterF,Acc_E,Steps,Model):
 Steps_copy=copy.copy(Steps)

 basic.Reload_EnvEnv(Env,Env1)

 if Model is "Heisenberg":
   H0=basic.Heisenberg0(h[0],h[1])
   H00=basic.Heisenberg00(h[0],h[1])
   H1=basic.Heisenberg1(h[2])
 if Model is "Heisenberg_Z2":
   H0=basic.Heisenberg0_Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg00_Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)

  
 U = uni10.UniTensor(H0.bond(), "U");
 blk_qnums = H0.blockQnum()
 for qnum in blk_qnums:
  U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))

 U0 = uni10.UniTensor(H00.bond(), "U");
 blk_qnums = H00.blockQnum()
 for qnum in blk_qnums:
  U0.putBlock(qnum, uni10.takeExp(-delta, H00.getBlock(qnum)))


 U1 = uni10.UniTensor(H1.bond(), "U");
 blk_qnums = H1.blockQnum()
 for qnum in blk_qnums:
  U1.putBlock(qnum, uni10.takeExp(-delta, H1.getBlock(qnum)))
 
 
 
 E_0=1.0
 E_1=2.0
 E_min=1.e14
 
 List_delN=basic.Short_TrotterSteps(N_iterF)
 #List_delN=basic.Short_TrotterSteps1(N_iterF)
 #List_delN=basic.Long_TrotterSteps(N_iterF)
 #List_delN=basic.Long_TrotterSteps1(N_iterF)
 print List_delN
 for delta, N_iter in List_delN:

  print delta, N_iter

  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
      U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))
  blk_qnums = H00.blockQnum()
  for qnum in blk_qnums:
      U0.putBlock(qnum, uni10.takeExp(-delta, H00.getBlock(qnum)))


  for q in xrange(N_iter):
   t0=time.time()

#################################################################################################
#   #print U
#   #t0=time.time()
#   if delta>=0.05: break;
#   if delta>=0.04: break;
   c=copy.copy(b)
   d=copy.copy(a)

   a_u, b_u, a, b=basicB.Var_ab(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0)
   #print time.time() - t0, "Seconds, Left"

   c=copy.copy(b)
   d=copy.copy(a)

   b_u, a_u, b, a=basicB.Var_ab(b_u,a_u,b,a,d,c,Env1,D,U,d_phys,chi,Gauge,Corner_method,H0)

#########
   c=copy.copy(b)
   d=copy.copy(a)

   a_u, b_u, a, b=basicB.Var_ca(a_u,b_u,b,a,d,c,Env1,D,U0,d_phys,chi,Gauge,Corner_method,H00)

   c=copy.copy(b)
   d=copy.copy(a)

   b_u, a_u, b, a=basicB.Var_ca(b_u,a_u,a,b,c,d,Env,D,U0,d_phys,chi,Gauge,Corner_method,H00)


####################################################################################################


   
   E_value=basic.E_total_conv(a_u,b_u,a,b,Env,Env1,D,h,d_phys,chi,Corner_method,Model)
   print 'E_s=', E_value, time.time() - t0
   #M_value=basic.M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   #print 'M_s=', M_value

   basic.Store_Full(a_u,b_u,a,b)
   basic.Store_EnvEnv(Env,Env1)


   E_0=E_1
   E_1=E_value
   if E_1 < E_min:
    basic.Store_Full(a_u,b_u,a,b)
    E_min=E_1
    print "E_m=", E_min
   print 'E_d=', abs((E_0-E_1) / E_0) , 'Num_iter=', q 


   if (( abs((E_0-E_1) / E_0) ) < Acc_E) or ( q is int(N_iter-1)): 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0)
    a_u,b_u,a,b=basic.Reload_Full()
    E_value=basic.E_total(a_u,b_u,a,b,Env,Env1,D,h,d_phys,chi,Corner_method,Model)
    basic.Store_Full(a_u,b_u,a,b)
    E_0=10.0
    E_1=20.0
    print 'E_t=', E_value   
    break;

 
 return a_u,b_u,a,b, Env

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


 


