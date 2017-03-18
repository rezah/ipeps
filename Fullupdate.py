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

def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,delta,h,Env,Env1,Env2,Env3,Gauge,Positive,Corner_method,N_iterF,Acc_E,Steps,Model):
 Steps_copy=copy.copy(Steps)
 
 EnvP=[ Env[i] for i in xrange(len(Env)) ]
 EnvP1=[ Env[i] for i in xrange(len(Env)) ]
 EnvP2=[ Env[i] for i in xrange(len(Env)) ]
 EnvP3=[ Env[i] for i in xrange(len(Env)) ]

 if Model is "Heisenberg":
   H0=basic.Heisenberg0(h[0],h[1])
   H1=basic.Heisenberg1(h[2])
 if Model is "Heisenberg_Z2":
   H0=basic.Heisenberg0_Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)

  
 U = uni10.UniTensor(H0.bond(), "U");
 blk_qnums = H0.blockQnum()
 for qnum in blk_qnums:
  U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))
 U1 = uni10.UniTensor(H1.bond(), "U");
 blk_qnums = H1.blockQnum()
 for qnum in blk_qnums:
  U1.putBlock(qnum, uni10.takeExp(-delta, H1.getBlock(qnum)))

 
 
 E_0=1.0
 E_1=2.0
 E_min=1.e14
 for i in xrange(1,600):

  if i is 1:
   delta=1.0e-1
   N_iter=N_iterF
  if i is 2:
   delta=0.50e-1
   N_iter=N_iterF
  if i is 3:
   delta=1.0e-2
   N_iter=N_iterF
  if i is 4:
   delta=0.50e-2
  if i is 5:
   delta=1.0e-3
  if i is 6:
   delta=0.50e-3
   N_iter=N_iterF
  if i is 7:
   delta=1.0e-4
   N_iter=N_iterF
  if i is 8:
   break
  print delta, N_iter
  
  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
      U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))
  blk_qnums = H1.blockQnum()
  for qnum in blk_qnums:
      U1.putBlock(qnum, uni10.takeExp(-delta, H1.getBlock(qnum)))
  
  #MPO_list=basic.Decomposition(U)

  for q in xrange(N_iter):


   E_value=basic.E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   print 'E_s=', E_value
   E_0=E_1
   E_1=E_value
   if E_1 < E_min:
    basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
    E_min=E_1
    print "E_m=", E_min
   print 'E_d=', abs((E_0-E_1) / E_0) , 'Num_iter=', q 


   if (( abs((E_0-E_1) / E_0) ) < Acc_E) or ( q is int(N_iter-1)): 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0)
    a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
    E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
    basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
    E_0=10.0
    E_1=20.0
    print 'E_t=', E_value   
    break;


#################################################################################################
   print "1"
#   #print U
#   #t0=time.time()
   a_u, b_u, a, b=basicB.Var_ab(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   c_u, a_u, c, a=basicB.Var_ca(c_u,a_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   
   c_u, a_u, b_u , c, a, b=basicA.Var_cb(c_u,a_u,b_u,a,b,c,d,Env,D,U1,d_phys,chi,Gauge,Positive,Corner_method)

   a_u, b_u, d_u, a, b, d=basicA.Var_ad(a_u, b_u, d_u,a,b,c,d,Env,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
#################################################################################################


###############################################################################################
   print "2"
   c_u, d_u, c, d=basicB.Var_ab(c_u,d_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   a_u, c_u, a, c=basicB.Var_ca(a_u,c_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Positive,Corner_method)


   a_u, c_u, d_u, a, c, d=basicA.Var_cb(a_u,c_u,d_u,c,d,a,b,Env1,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
   c_u, d_u, b_u, c, d, b=basicA.Var_ad(c_u, d_u, b_u,c,d,a,b,Env1,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
#################################################################################################
   print "3"

   b_u, a_u, b, a=basicB.Var_ab(b_u,a_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   d_u, b_u, d, b=basicB.Var_ca(d_u,b_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Positive,Corner_method)


   d_u, b_u, a_u, d, b, a=basicA.Var_cb(d_u,b_u,a_u,b,a,d,c,Env2,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
   b_u, a_u, c_u, b, a, c=basicA.Var_ad(b_u, a_u, c_u,b,a,d,c,Env2,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
###################################################################################################
   print "4"

   d_u, c_u, d, c=basicB.Var_ab(d_u,c_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   b_u, d_u, b, d=basicB.Var_ca(b_u,d_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Positive,Corner_method)

   b_u, d_u, c_u, b, d, c=basicA.Var_cb(b_u,d_u,c_u,d,c,b,a,Env3,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
   d_u, c_u, a_u, d, c, a=basicA.Var_ad(d_u, c_u, a_u,d,c,b,a,Env3,D,U1,d_phys,chi,Gauge,Positive,Corner_method)

##################################################################################################



#####################################################################################

 
 #Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]
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


 


