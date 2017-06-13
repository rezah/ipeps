import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
#import random
import copy
import time
import basic
#import basicA
import basicB
import basicC

def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,h,Env,Env1,Env2,Env3,Gauge,Corner_method,N_iterF,Acc_E,Model,N_grad, Opt_method,Inv_method,N_svd,N_env,method,check_step):

 H0, H00, H1, H2=basic.choose_model(Model, h, d_phys)

 U = uni10.UniTensor(H0.bond(), "U");
 U0 = uni10.UniTensor(H00.bond(), "U0");
 U1 = uni10.UniTensor(H1.bond(), "U1");
 U2 = uni10.UniTensor(H2.bond(), "U2");
 

 E_0=1.0
 E_1=2.0
 E_min=1.e+14


 MPO_list=basic.Decomposition(H2)
 plist=basicC.initialize_plist(a_u, b_u, c_u, MPO_list)
 #basicC.Reload_plist(plist)
 plist1=basicC.initialize_plist(b_u, a_u, d_u, MPO_list)
 plist2=basicC.initialize_plist(c_u, d_u, a_u, MPO_list)
 plist3=basicC.initialize_plist(d_u, c_u, b_u, MPO_list)

 plist00=basicC.initialize_plist1(a_u, b_u, d_u, MPO_list)
 plist11=basicC.initialize_plist1(b_u, a_u, c_u, MPO_list)
 plist22=basicC.initialize_plist1(c_u, d_u, b_u, MPO_list)
 plist33=basicC.initialize_plist1(d_u, c_u, a_u, MPO_list)
 
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
  MPO_list=basic.Decomposition(U2)

  blk_qnums = H1.blockQnum()
  for qnum in blk_qnums:
      U1.putBlock(qnum, uni10.takeExp(-delta, H1.getBlock(qnum)))


  for q in xrange(N_iter):
   t0=time.time()

#################################################################################################
##   
#   a_u, b_u, a, b=basicB.Var_ab(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0,N_env,N_svd)
#   #print time.time() - t0, "Seconds, Left"

#   c_u, d_u, c, d=basicB.Var_ab(c_u,d_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Corner_method,H0,N_env,N_svd)
########

#   b_u, a_u, b, a=basicB.Var_ab(b_u,a_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Corner_method,H0,N_env,N_svd)

#   d_u, c_u, d, c=basicB.Var_ab(d_u,c_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Corner_method,H0,N_env,N_svd)

##########

#   c_u, a_u, c, a=basicB.Var_ca(c_u,a_u,a,b,c,d,Env,D,U0,d_phys,chi,Gauge,Corner_method,H00,N_env,N_svd)

#   d_u, b_u, d, b=basicB.Var_ca(d_u,b_u,b,a,d,c,Env2,D,U0,d_phys,chi,Gauge,Corner_method,H00,N_env,N_svd)

#########
#   a_u, c_u, a, c=basicB.Var_ca(a_u,c_u,c,d,a,b,Env1,D,U0,d_phys,chi,Gauge,Corner_method,H00,N_env,N_svd)

#   b_u, d_u, b, d=basicB.Var_ca(b_u,d_u,d,c,b,a,Env3,D,U0,d_phys,chi,Gauge,Corner_method,H00,N_env,N_svd)
#####################################################################################################

######################################################################################

   print "\n"
   a_u, b_u,c_u,d_u, a,b,c,d=basicC.Var_cab(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist,MPO_list,Inv_method,N_svd,N_env,method,check_step)
   print "\n"
   b_u, a_u,d_u,c_u, b,a,d,c=basicC.Var_cab(b_u, a_u,d_u,c_u,b,a,d,c,Env1,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist1,MPO_list,Inv_method,N_svd,N_env,method,check_step)
   print "\n"
   c_u, d_u,a_u,b_u, c,d,a,b=basicC.Var_cab(c_u, d_u,a_u,b_u,c,d,a,b,Env2,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist2,MPO_list,Inv_method,N_svd,N_env,method,check_step)
   print "\n"
   d_u, c_u,b_u,a_u, d,c,b,a=basicC.Var_cab(d_u, c_u,b_u,a_u,d,c,b,a,Env3,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist3,MPO_list,Inv_method,N_svd,N_env,method,check_step)

###########################

   print "\n"
   a_u, b_u,c_u,d_u, a,b,c,d=basicC.Var_abd(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist00,MPO_list,Inv_method,N_svd,N_env,method,check_step)

   print "\n"
   b_u, a_u,d_u,c_u, b,a,d,c=basicC.Var_abd(b_u, a_u,d_u,c_u,b,a,d,c,Env1,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist11,MPO_list,Inv_method,N_svd,N_env,method,check_step)

   print "\n"
   c_u, d_u,a_u,b_u, c,d,a,b=basicC.Var_abd(c_u, d_u,a_u,b_u,c,d,a,b,Env2,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist22,MPO_list,Inv_method,N_svd,N_env,method,check_step)

   print "\n"
   d_u, c_u,b_u,a_u, d,c,b,a=basicC.Var_abd(d_u, c_u,b_u,a_u,d,c,b,a,Env3,D,U2,d_phys,chi,Gauge,Corner_method,H2,N_grad, Opt_method,plist33,MPO_list,Inv_method,N_svd,N_env,method,check_step)

#############################################################################################
   
   E_value=basic.E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
   print 'E_s=', E_value, time.time() - t0
   M_value=basic.M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   print 'M_s=', M_value
   basic.Translational_sym(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)

   #break


   E_0=E_1
   E_1=E_value
   if E_1 < E_min:
    E_min=E_1
    print "E_m=", E_min
    basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
    basic.Store_EnvEnv(Env,Env1,Env2,Env3)
    basicC.Store_plist(plist)
    basicC.Store_plist1(plist00)
   elif q>2:
    a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
    basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
    basicC.Reload_plist(plist)
    basicC.Reload_plist1(plist00)
    E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
    print "E_min=", E_value
    break; 

   print 'E_d=', abs((E_0-E_1) / E_0) , 'Num_iter=', q 


   if (( abs((E_0-E_1) / E_0) ) < Acc_E) or ( q is int(N_iter-1)): 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0)
    a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
    basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
    basicC.Reload_plist(plist)
    basicC.Reload_plist1(plist00)
    E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model,N_env)
    #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
    E_0=10.0
    E_1=20.0
    print 'E_t=', E_value
    break;

 
 return a_u, b_u, c_u, d_u, a, b, c, d, Env, Env1, Env2, Env3


 


