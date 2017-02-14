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


def itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,
Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iterF,
delta, h,Steps,Model):
  
  
  if Model is "Ising":
   H0=basicitebd.transverseIsing(h)
  if Model is "Heisenberg":
   H0=basicitebd.Heisenberg(h)

  U = uni10.UniTensor(H0.bond(), "U");
  Steps_copy=copy.copy(Steps)
  for i in xrange(1,1000):
   delta, N_iterF=basicitebd.Def_deltaNiter(i,N_iterF,Steps_copy)
   if delta is 0: break;
   print 'delta =', delta
   print "N_iterF=", N_iterF


   for q in xrange(N_iterF):
   
    U.putBlock(uni10.takeExp(-delta, H0.getBlock()))

    Gamma=[Gamma_a,Gamma_b]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
  #rlink
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys)
  
  #rlink
    Gamma=[Gamma_b, Gamma_a]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_5,Landa_6,Landa_2,Landa_4]
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys)

    Gamma=[Gamma_c,Gamma_d]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
  #rlink
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys)
  
  #rlink
    Gamma=[Gamma_d, Gamma_c]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_3,Landa_1,Landa_4,Landa_2]
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys)
################################################################################

 #ulink
    Gamma=[Gamma_a, Gamma_c]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys)

 #ulink
    Gamma=[ Gamma_c,Gamma_a]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_7,Landa_8]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys)

 #ulink
    Gamma=[ Gamma_b,Gamma_d]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys)

 #ulink
    Gamma=[ Gamma_d,Gamma_b]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_2,Landa_4]
    basicitebd.update_ulink_eff(Gamma,Landa,U,D,d_phys)


  return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8  


def itebd_standard(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,
Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iterF,
delta, h):

  H0=basicitebd.transverseIsing(h)
  U = uni10.UniTensor(H0.bond(), "U");

  for i in xrange(1,100):
   delta=1.00/pow(2.00,i)
   print 'delta =', delta
   if delta>1.0e-1:
    N_iter=N_iterF
   if delta<1.0e-1 and delta>1.0e-4:
    N_iter=N_iterF
   if delta<1.0e-4  and delta>1.0e-5:
    N_iter=N_iterF
   if delta<1.0e-5:
    break

   #print "N_iter=", N_iter


   for q in xrange(N_iter):
   
    U.putBlock(uni10.takeExp(-delta, H0.getBlock()))

    Gamma=[Gamma_a,Gamma_b]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
  #rlink
    basicitebd.update_rlink(Gamma,Landa,U,D,d_phys)
  
  #rlink
    Gamma=[Gamma_b, Gamma_a]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_5,Landa_6,Landa_2,Landa_4]
    basicitebd.update_rlink(Gamma,Landa,U,D,d_phys)

    Gamma=[Gamma_c,Gamma_d]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
  #rlink
    basicitebd.update_rlink(Gamma,Landa,U,D,d_phys)
  
  #rlink
    Gamma=[Gamma_d, Gamma_c]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_3,Landa_1,Landa_4,Landa_2]
    basicitebd.update_rlink(Gamma,Landa,U,D,d_phys)


 #ulink
    Gamma=[Gamma_a, Gamma_c]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
    basicitebd.update_ulink(Gamma,Landa,U,D,d_phys)

 #ulink
    Gamma=[ Gamma_c,Gamma_a]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_7,Landa_8]
    basicitebd.update_ulink(Gamma,Landa,U,D,d_phys)

 #ulink
    Gamma=[ Gamma_b,Gamma_d]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
    basicitebd.update_ulink(Gamma,Landa,U,D,d_phys)

 #ulink
    Gamma=[ Gamma_d,Gamma_b]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_2,Landa_4]
    basicitebd.update_ulink(Gamma,Landa,U,D,d_phys)


  return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8  
