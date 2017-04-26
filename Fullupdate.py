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

def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,delta,h,Env,Env1,Env2,Env3,Gauge,Corner_method,N_iterF,Acc_E,Steps,Model,N_grad, Opt_method):

 Steps_copy=copy.copy(Steps)

 #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)

 if Model is "Heisenberg":
   H0=basic.Heisenberg0(h[0],h[1])
   H00=basic.Heisenberg00(h[0],h[1])
   H1=basic.Heisenberg1(h[2])
   H2=basic.threebody(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=basic.Heisenberg0_Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg00_Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)
   H2=basic.threebody_U1(h,d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)

  
 U = uni10.UniTensor(H0.bond(), "U");
 blk_qnums = H0.blockQnum()
 for qnum in blk_qnums:
  U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))




 U0 = uni10.UniTensor(H00.bond(), "U0");
 U1 = uni10.UniTensor(H1.bond(), "U1");
 U2 = uni10.UniTensor(H2.bond(), "U2");
 
 
 E_0=1.0
 E_1=2.0
 E_min=1.e+14
 
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

  blk_qnums = H2.blockQnum()
  for qnum in blk_qnums:
      U2.putBlock(qnum, uni10.takeExp(-delta, H2.getBlock(qnum)))

  blk_qnums = H1.blockQnum()
  for qnum in blk_qnums:
      U1.putBlock(qnum, uni10.takeExp(-delta, H1.getBlock(qnum)))


  for q in xrange(N_iter):
   t0=time.time()

#################################################################################################
   
#   a_u, b_u, a, b=basicB.Var_ab(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0)
#   #print time.time() - t0, "Seconds, Left"

#   c_u, d_u, c, d=basicB.Var_ab(c_u,d_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Corner_method,H0)
########

#   b_u, a_u, b, a=basicB.Var_ab(b_u,a_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Corner_method,H0)

#   d_u, c_u, d, c=basicB.Var_ab(d_u,c_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Corner_method,H0)

##########

#   c_u, a_u, c, a=basicB.Var_ca(c_u,a_u,a,b,c,d,Env,D,U0,d_phys,chi,Gauge,Corner_method,H00)

#   d_u, b_u, d, b=basicB.Var_ca(d_u,b_u,b,a,d,c,Env2,D,U0,d_phys,chi,Gauge,Corner_method,H00)

#########
#   a_u, c_u, a, c=basicB.Var_ca(a_u,c_u,c,d,a,b,Env1,D,U0,d_phys,chi,Gauge,Corner_method,H00)

#   b_u, d_u, b, d=basicB.Var_ca(b_u,d_u,d,c,b,a,Env3,D,U0,d_phys,chi,Gauge,Corner_method,H00)
###################################################################################################


#######################################################################################
   print "\n"
   a_u, b_u,c_u,d_u, a,b,c,d=basicC.Var_cab(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

   print "\n"
   b_u, a_u,d_u,c_u, b,a,d,c=basicC.Var_cab(b_u, a_u,d_u,c_u,b,a,d,c,Env1,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

   print "\n"
   c_u, d_u,a_u,b_u, c,d,a,b=basicC.Var_cab(c_u, d_u,a_u,b_u,c,d,a,b,Env2,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

   print "\n"
   d_u, c_u,b_u,a_u, d,c,b,a=basicC.Var_cab(d_u, c_u,b_u,a_u,d,c,b,a,Env3,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

###########################

   print "\n"
   a_u, b_u,c_u,d_u, a,b,c,d=basicC.Var_abd(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

   print "\n"
   b_u, a_u,d_u,c_u, b,a,d,c=basicC.Var_abd(b_u, a_u,d_u,c_u,b,a,d,c,Env1,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

   print "\n"
   c_u, d_u,a_u,b_u, c,d,a,b=basicC.Var_abd(c_u, d_u,a_u,b_u,c,d,a,b,Env2,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)

   print "\n"
   d_u, c_u,b_u,a_u, d,c,b,a=basicC.Var_abd(d_u, c_u,b_u,a_u,d,c,b,a,Env3,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method)



#############################################################################################
   
   E_value=basic.E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   print 'E_s=', E_value, time.time() - t0
   #M_value=basic.M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   #print 'M_s=', M_value
   #break


   E_0=E_1
   E_1=E_value
   if E_1 < E_min:
    E_min=E_1
    print "E_m=", E_min
    basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
    basic.Store_EnvEnv(Env,Env1,Env2,Env3)
   else:
    a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
    break; 

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

 
 return a_u,b_u,c_u,d_u,a,b,c,d, Env


 


