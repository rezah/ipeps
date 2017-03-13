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
 Truncation=[0]

 if Model is "Ising":
  H0=basic.transverseIsing_Z2(h,d_phys)
 if Model is "Heisenberg":
  H0=basic.Heisenberg_Z2(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=basic.Heisenberg_Z2(h,d_phys)
 if Model is "Heisenberg_U1":
   H0=basic.Heisenberg_U1(h,d_phys)
 if Model is "threebody_Z2":
   H0=basic.threebody_Z2(h,d_phys)
 if Model is "threebody_U1":
   H0=basic.threebody_U1(h,d_phys)
 if Model is "threebody":
   H0=basic.threebody(h,d_phys)
  
 U = uni10.UniTensor(H0.bond(), "U");
 blk_qnums = H0.blockQnum()
 for qnum in blk_qnums:
  U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))
# MPO_list=basic.Decomposition(U)
# plist, plistd=basicA.initialize_plist(a_u, b_u, c_u, MPO_list)
# plist1, plistd1=basicB.initialize_plist(a_u, c_u, d_u, MPO_list)
 
 
 E_0=1.0
 E_1=2.0
 
 for i in xrange(1,2):

  if i is 1:
   delta=1.0e-1
   N_iter=N_iterF
  if i is 2:
   delta=0.50e-1
   N_iter=N_iterF
  if i is 1:
   delta=1.0e-2
   N_iter=N_iterF
  if i is 2:
   delta=0.50e-2
  if i is 3:
   delta=1.0e-3
  if i is 4:
   delta=0.50e-3
   N_iter=N_iterF
  if i is 5:
   delta=1.0e-4
   N_iter=N_iterF
  if i is 6:
   break
  #delta=0.0
  N_iter=4
  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
      U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))
  
  #MPO_list=basic.Decomposition(U)

  for q in xrange(N_iter):


#################################################################################################
#   a_u, b_u, c_u, a, b, c=basicA.Var_cab(a_u, b_u, c_u,a,b,c,d,Env,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist, plistd)

   


#   a_u, c_u, d_u, a, c, d=basicB.Var_acd(a_u,c_u,d_u,a,b,c,d,Env,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist1, plistd1)

####################################################################################################
#   b_u, a_u, d_u, b, a, d=basicA.Var_cab(b_u,a_u,d_u,b,a,d,c,Env1,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist, plistd)


#   b_u, d_u, c_u, b, d, c=basicB.Var_acd(b_u,d_u,c_u,b,a,d,c,Env1,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist1, plistd1)

#####################################################################################################
#   c_u, d_u, a_u, c, d, a=basicA.Var_cab(c_u,d_u,a_u,c,d,a,b,Env2,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist, plistd)

#   c_u, a_u, b_u, c, a, b=basicB.Var_acd(c_u,a_u,b_u,c,d,a,b,Env2,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist1, plistd1)

######################################################################################################
#   d_u, c_u, b_u, d, c, b=basicA.Var_cab(d_u,c_u,b_u,d,c,b,a,Env3,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist, plistd)


#   d_u, b_u, a_u, d, b, a=basicB.Var_acd(d_u,b_u,a_u,d,c,b,a,Env3,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist1, plistd1)

######################################################################################################









   #t0=time.time()
   a_u, b_u, a, b=basic.Var_H(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"


   #t0=time.time()
   c_u, a_u, c, a=basic.Var_V(c_u,a_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"


   #t0=time.time()
   c_u, d_u, c, d=basic.Var_H(c_u,d_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   #t0=time.time()
   a_u, c_u, a, c=basic.Var_V(a_u,c_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   #t0=time.time()
   b_u, a_u, b, a=basic.Var_H(b_u,a_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   #t0=time.time()
   d_u, b_u, d, b=basic.Var_V(d_u,b_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"
   
   #t0=time.time()
   d_u, c_u, d, c=basic.Var_H(d_u,c_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"


   #t0=time.time()
   b_u, d_u, b, d=basic.Var_V(b_u,d_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"



   E_value=basic.E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   #print 'E_toal=', E_value
   E_value=1.00*q
   E_0=E_1+20
   E_1=E_value+10
   print 'E_diff=', abs((E_0-E_1) / E_0) , 'Num_iter=', q 
   if (( abs((E_0-E_1) / E_0) ) < Acc_E) or ( q is int(N_iter-1)): 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0)
    #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
    print 'E_toal=', E_value   
    break;
 
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


 


