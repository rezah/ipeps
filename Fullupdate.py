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

def Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,delta,h,Env,Env1,Env2,Env3,Gauge,Positive,Corner_method,N_iterF,Acc_E,Steps,Model):
 Steps_copy=copy.copy(Steps)

 #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 EnvP=[ Env[i] for i in xrange(len(Env)) ]
 EnvP1=[ Env[i] for i in xrange(len(Env)) ]
 EnvP2=[ Env[i] for i in xrange(len(Env)) ]
 EnvP3=[ Env[i] for i in xrange(len(Env)) ]

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
 
# H2=basic.threebody(h,d_phys)
# U2 = uni10.UniTensor(H2.bond(), "U");
# blk_qnums = H2.blockQnum()
# for qnum in blk_qnums:
#  U2.putBlock(qnum, uni10.takeExp(-delta, H2.getBlock(qnum)))
# MPO_list=basic.Decomposition(U2)

# plist=basicC.initialize_plist(a_u, b_u, c_u, MPO_list)
# plist1=[copy.copy(plist[i])  for i in xrange(len(plist)) ]
# plist2=[copy.copy(plist[i])  for i in xrange(len(plist)) ]
# plist3=[copy.copy(plist[i])  for i in xrange(len(plist)) ]
# 
# plist00=basicC.initialize_plist1(a_u, c_u, d_u, MPO_list)
# plist11=[copy.copy(plist00[i])  for i in xrange(len(plist00)) ]
# plist22=[copy.copy(plist00[i])  for i in xrange(len(plist00)) ]
# plist33=[copy.copy(plist00[i])  for i in xrange(len(plist00)) ]
 
 
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

#  blk_qnums = H1.blockQnum()
#  for qnum in blk_qnums:
#      U1.putBlock(qnum, uni10.takeExp(-delta, H1.getBlock(qnum)))
#  for qnum in blk_qnums:
#      U2.putBlock(qnum, uni10.takeExp(-delta, H2.getBlock(qnum)))
#  
#  MPO_list=basic.Decomposition(U2)

  for q in xrange(N_iter):


#################################################################################################
#   #print U
#   #t0=time.time()
#   if delta>=0.05: break;
#   if delta>=0.04: break;

   a_u, b_u, a, b=basicB.Var_ab(a_u,b_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method)
   #print time.time() - t0, "Seconds, Left"

   c_u, d_u, c, d=basicB.Var_ab(c_u,d_u,c,d,a,b,Env1,D,U,d_phys,chi,Gauge,Positive,Corner_method)


#######


   b_u, a_u, b, a=basicB.Var_ab(b_u,a_u,b,a,d,c,Env2,D,U,d_phys,chi,Gauge,Positive,Corner_method)

   d_u, c_u, d, c=basicB.Var_ab(d_u,c_u,d,c,b,a,Env3,D,U,d_phys,chi,Gauge,Positive,Corner_method)


#########

   c_u, a_u, c, a=basicB.Var_ca(c_u,a_u,a,b,c,d,Env,D,U0,d_phys,chi,Gauge,Positive,Corner_method)



   d_u, b_u, d, b=basicB.Var_ca(d_u,b_u,b,a,d,c,Env2,D,U0,d_phys,chi,Gauge,Positive,Corner_method)


########



   a_u, c_u, a, c=basicB.Var_ca(a_u,c_u,c,d,a,b,Env1,D,U0,d_phys,chi,Gauge,Positive,Corner_method)



   b_u, d_u, b, d=basicB.Var_ca(b_u,d_u,d,c,b,a,Env3,D,U0,d_phys,chi,Gauge,Positive,Corner_method)



################################################################################################
################################################################################################
####################################################################################################
#   c_u, a_u, b_u , c, a, b=basicA.Var_cb(c_u,a_u,b_u,a,b,c,d,Env,D,U1,d_phys,chi,Gauge,Positive,Corner_method)

#   a_u, c_u, d_u, a, c, d=basicA.Var_cb(a_u,c_u,d_u,c,d,a,b,Env1,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
######

#   d_u, b_u, a_u, d, b, a=basicA.Var_cb(d_u,b_u,a_u,b,a,d,c,Env2,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
#   
#   b_u, d_u, c_u, b, d, c=basicA.Var_cb(b_u,d_u,c_u,d,c,b,a,Env3,D,U1,d_phys,chi,Gauge,Positive,Corner_method)

######

#   a_u, b_u, d_u, a, b, d=basicA.Var_ad(a_u, b_u, d_u,a,b,c,d,Env,D,U1,d_phys,chi,Gauge,Positive,Corner_method)
#   
#   c_u, d_u, b_u, c, d, b=basicA.Var_ad(c_u, d_u, b_u,c,d,a,b,Env1,D,U1,d_phys,chi,Gauge,Positive,Corner_method)

######
#   b_u, a_u, c_u, b, a, c=basicA.Var_ad(b_u, a_u, c_u,b,a,d,c,Env2,D,U1,d_phys,chi,Gauge,Positive,Corner_method)


#   d_u, c_u, a_u, d, c, a=basicA.Var_ad(d_u, c_u, a_u,d,c,b,a,Env3,D,U1,d_phys,chi,Gauge,Positive,Corner_method)

######################################################################################


#######################################################################################
#   print "\n"
#   c_u, a_u, b_u, c, a, b=basicC.Var_cab(c_u, a_u, b_u,a,b,c,d,Env,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist)

#   print "\n"

#   b_u, d_u, c_u, b, d, c=basicC.Var_cab(b_u, d_u, c_u,d,c,b,a,Env3,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist3)

#   print "\n"

#   d_u, b_u, a_u, d, b, a=basicC.Var_cab(d_u, b_u, a_u,b,a,d,c,Env1,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist1)

#   print "\n"

#   a_u, c_u, d_u, a, c, d=basicC.Var_cab(a_u, c_u, d_u,c,d,a,b,Env2,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist2)

###########################
#   print "\n"
#   a_u, c_u, d_u, a, c, d=basicC.Var_acd(a_u,c_u,d_u,a,b,c,d,Env,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist00)

#   print "\n"

#   d_u, b_u, a_u, d, b, a=basicC.Var_acd(d_u,b_u,a_u,d,c,b,a,Env3,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist33)
#   print "\n"

#   b_u, d_u, c_u, b, d, c=basicC.Var_acd(b_u,d_u,c_u,b,a,d,c,Env1,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist11)
#   print "\n"

#   c_u, a_u, b_u, c, a, b=basicC.Var_acd(c_u,a_u,b_u,c,d,a,b,Env2,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist22)

#############################################################################################


   
   E_value=basic.E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   print 'E_s=', E_value
   #M_value=basic.M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model)
   #print 'M_s=', M_value

   basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
   basic.Store_EnvEnv(Env,Env1,Env2,Env3)


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


 


