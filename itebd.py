import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time
import basicitebd
import basic


def itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,
Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iterF, h,Model,q_D):
  
  
  if Model is "Ising":
   H0=basic.transverseIsing_Z2(h,d_phys)
  if Model is "Heisenberg":
   H0=basic.Heisenberg_Z2(h,d_phys)
  U = uni10.UniTensor(H0.bond(), "U");
  
  for i in xrange(1,200):

   delta=1.00/pow(5.0,i) 

   if delta>1.0e-1:
    N_iter=N_iterF
    #delta=1.0e-1
   if delta<1.0e-1 and delta>1.0e-3:
    N_iter=N_iterF
   if delta<1.0e-3  and delta>1.0e-5:
    N_iter=N_iterF
   if delta<1.0e-5:
    break


   print 'delta =', delta
   print "N_iterF=", N_iterF

   for q in xrange(N_iterF):
    blk_qnums = H0.blockQnum()
    for qnum in blk_qnums:
        U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))
        
    #print "U", U
#############################################################################
    Gamma=[Gamma_a,Gamma_b]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
  #rlink
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)


  #rlink
    Gamma=[Gamma_b, Gamma_a]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_5,Landa_6,Landa_2,Landa_4]
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)

    Gamma=[Gamma_c,Gamma_d]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
  #rlink
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)
  
  #rlink
    Gamma=[Gamma_d, Gamma_c]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_3,Landa_1,Landa_4,Landa_2]
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)

############################################################################

 #ulink
    Gamma=[Gamma_a, Gamma_c]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys,q_D)
    
 #ulink
    Gamma=[ Gamma_c,Gamma_a]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_7,Landa_8]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys,q_D)

 #ulink
    Gamma=[ Gamma_b,Gamma_d]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys,q_D)

 #ulink
    Gamma=[ Gamma_d,Gamma_b]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_2,Landa_4]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys,q_D)


################################################################################



  return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8  



