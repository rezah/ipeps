import pyUni10 as uni10
import sys
import numpy as np
#import matplotlib
#matplotlib.use('pdf')
#import matplotlib.pyplot as plt
#import pylab
import random
import copy
import time
import basic
import itebd
import Fullupdate
import Move
###################### Initialize parameters ###########################
#Model="Heisenberg_Z2"         #Heisenberg, Ising
Model="Heisenberg_U1"         #Heisenberg, Ising
#Model="Heisenberg_U1Z2"         #Heisenberg, Ising
#Model="Heisenberg"         #Heisenberg, Ising
#Model="threebody"         #Heisenberg, Ising
#Model="threebody_Z2"         #Heisenberg, Ising
#Model="threebody_U1"         #Heisenberg, Ising
#Model="threebody"         #Heisenberg, Ising

#D=[4]
#chi=[20]
#d_phys=[8]
##d_phys=[2]

#D=[2,1]
#chi=[10,10]
#d_phys=[4,4]

D=[2,2,1]
chi=[10,20,20,20,10]
d_phys=[1,1]

#D=[2,2,2,2,2,2]
#chi=[10,10,10,10,10,10]
#d_phys=[1,1]


N_iteritebd=200
N_iterF=30
Gauge='Fixed'
Corner_method='CTMFull'   #CTM, CTMRG, CTMFull
Acc_E=1.00e-7
Steps=[1.0e-1,4.0e-2,4.0e-3,4.0e-4,5.0e-5,1.0e-6] #,[Start,steps,End] 
delta=0.001
###################################################################
######################### No-symmetry #############################################

#q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
#q_list=[q0_even]
#qchi_list=[q0_even]
#q_phys=[q0_even]*d_phys[0]

###################### Z(2) ######################################

#q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
#q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
#q_list=[q0_even,q0_odd]
#qchi_list=[q0_even,q0_odd]
#q_phys=[q0_even,q0_even,q0_even,q0_even,q0_odd,q0_odd,q0_odd,q0_odd]
##q_phys=[q0_even,q0_odd]

##########################  U(1)  ################################

q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
q1_even = uni10.Qnum(1,uni10.PRT_EVEN);
q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
q4_even = uni10.Qnum(4,uni10.PRT_EVEN);
q5_even = uni10.Qnum(5,uni10.PRT_EVEN);

q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);
q_5_even = uni10.Qnum(-5,uni10.PRT_EVEN);

#qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
#qchi_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even]
#qchi_list=[q_1_even,q0_even,q1_even]
#qchi_list=[q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even]
#qchi_list=[q_5_even,q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even,q5_even]
#qchi_list=[q_1_even,q1_even]

q_list=[q_1_even,q0_even,q1_even]
#q_list=[q_1_even,q0_even]
#q_list=[q_2_even,q_1_even,q0_even,q1_even, q2_even]
#q_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even, q2_even,q3_even]
#q_list=[q_1_even,q1_even]

q_phys=[q_1_even,q1_even]

##########################  Z2*U(1)  ################################
#q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
#q1_even = uni10.Qnum(1,uni10.PRT_EVEN);
#q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
#q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
#q4_even = uni10.Qnum(4,uni10.PRT_EVEN);

#q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
#q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
#q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
#q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);


#q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
#q1_odd = uni10.Qnum(1,uni10.PRT_ODD);
#q2_odd = uni10.Qnum(2,uni10.PRT_ODD);
#q3_odd = uni10.Qnum(3,uni10.PRT_ODD);
#q4_odd = uni10.Qnum(4,uni10.PRT_ODD);

#q_1_odd = uni10.Qnum(-1,uni10.PRT_ODD);
#q_2_odd = uni10.Qnum(-2,uni10.PRT_ODD);
#q_3_odd = uni10.Qnum(-3,uni10.PRT_ODD);
#q_4_odd = uni10.Qnum(-4,uni10.PRT_ODD);

#qchi_list=[q_1_even,q_1_odd,q0_even,q0_odd,q1_even,q1_odd]
#q_list=[q_1_even,q_1_odd,q0_even,q0_odd,q1_even,q1_odd]
#q_phys=[q_1_odd,q1_even]
#############################################################################



q_D=[]
q_D2=[]
q_chi=[]

for i in xrange(len(D)):
 for q in xrange(D[i]):
  q_D.append(q_list[i])


for i in xrange(len(chi)):
 for q in xrange(chi[i]):
  q_chi.append(qchi_list[i])

Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.produce_GammaLanda(q_D, q_phys)


Gamma=[Gamma_a,Gamma_b,Gamma_c,Gamma_d]
Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
basic.Initialize_function(Gamma,Landa)

Landa1=[Landa[2],Landa[1],Landa[0],Landa[3]]
a_u,a=basic.makeab(Landa1,Gamma_a)
Landa1=[Landa[0],Landa[6],Landa[2],Landa[7]]
b_u,b=basic.makeab(Landa1,Gamma_b)
Landa1=[Landa[4],Landa[3],Landa[5],Landa[1]]
c_u,c=basic.makeab(Landa1,Gamma_c)
Landa1=[Landa[5],Landa[7],Landa[4],Landa[6]]
d_u,d=basic.makeab(Landa1,Gamma_d)




Env=basic.produce_env_init(q_chi,q_D)
Env1=basic.Rand_env_total(Env)
Env2=basic.Rand_env_total(Env)
Env3=basic.Rand_env_total(Env)

zlist=[]
J1_list=[h*0.0100 for h in range(270,400)]
J1_list=[1.0]
J2_list=[0.50]
h_list=[1.0]

for h , J1, J2 in zip( h_list, J1_list, J2_list):
 print h, J1, J2 
 h=[h, J1, J2]
#########################################################################################
 
 Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.Reload_itebd()
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,q_phys,D,N_iteritebd,h,Model,q_D)
 
 basic.Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8)


 print Landa_1#,Landa_1[0],Landa_1[3],Landa_1[4],Landa_1[7]
 print Landa_2#,Landa_2[0],Landa_2[3],Landa_2[4],Landa_2[7]
 print Landa_3#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]
 print Landa_4#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]

 print Landa_5#,Landa_1[0],Landa_1[3],Landa_1[4],Landa_1[7]
 print Landa_6#,Landa_2[0],Landa_2[3],Landa_2[4],Landa_2[7]
 print Landa_7#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]
 print Landa_8#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]



 Landa=[Landa_3,Landa_2,Landa_1,Landa_4]
 a_u,a=basic.makeab(Landa,Gamma_a)
 Landa=[Landa_1,Landa_7,Landa_3,Landa_8]
 b_u,b=basic.makeab(Landa,Gamma_b)
 Landa=[Landa_5,Landa_4,Landa_6,Landa_2]
 c_u,c=basic.makeab(Landa,Gamma_c)
 Landa=[Landa_6,Landa_8,Landa_5,Landa_7]
 d_u,d=basic.makeab(Landa,Gamma_d)
 #print "a",a.printDiagram()
 #Env=basic.positve_Env(Env,a,b,c,d)
 #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model)

 #print 'E_toal=', E_value
 #break
#########################################################################################

############################################################################

 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.total_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full_previous(a_u, b_u, c_u, d_u)

 a_u,b_u,c_u,d_u,a,b,c,d,Env=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,q_phys,D,delta,h,Env,Env1,Env2,Env3,Gauge,Corner_method,N_iterF,Acc_E,Steps,Model)

 #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model)
# print "E_final", E_value

 basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)



