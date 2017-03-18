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
Model="Heisenberg"         #Heisenberg, Ising
#Model="threebody"         #Heisenberg, Ising
#Model="threebody_Z2"         #Heisenberg, Ising
#Model="threebody_U1"         #Heisenberg, Ising
#Model="threebody"         #Heisenberg, Ising

D=[4]
chi=[40]
d_phys=[1]



#D=[2,2]
#chi=[10,10]
#d_phys=[1,1]

#D=[2,2,2]
#chi=[4,4,4,4,4]
#d_phys=[1,1]


N_iteritebd=50
N_iterF=30
Gauge='Fixed'
Positive='Restrict'
Corner_method='CTM'   #CTM, CTMRG, CTMFull
Acc_E=1.00e-7
Steps=[1.0e-1,4.0e-2,4.0e-3,4.0e-4,5.0e-5,1.0e-6] #,[Start,steps,End] 
delta=0.001
###################################################################
zlist=[]
Elist=[]
zlist1=[]
Elist1=[]
zlist2=[]
Elist2=[]

######################### No-symmetry #############################################3

q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
q_list=[q0_even,q0_even]
qchi_list=[q0_even,q0_even]
q_phys=[q0_even,q0_even]


###################### Z(2) ######################################

#q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
#q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
#q_list=[q0_even,q0_odd]
#qchi_list=[q0_even,q0_odd]
#q_phys=[q0_even,q0_odd]
#q_phys=[q0_even,q0_odd]


##########################  U(1)  ################################

#q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
#q1_even = uni10.Qnum(1,uni10.PRT_EVEN);
#q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
#q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
#q4_even = uni10.Qnum(4,uni10.PRT_EVEN);

#q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
#q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
#q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
#q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);

##q_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
##qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
#qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]

##q_list=[q_1_even,q0_even,q1_even]
##qchi_list=[q_1_even,q0_even,q1_even]
#q_list=[q_1_even,q0_even,q1_even]
##q_list=[q_1_even,q0_even]
#q_phys=[q_1_even,q1_even]


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
##############################################################################



q_D=[]
q_D2=[]
q_chi=[]

for i in xrange(len(D)):
 for q in xrange(D[i]):
  q_D.append(q_list[i])




for i in xrange(len(chi)):
 for q in xrange(chi[i]):
  q_chi.append(qchi_list[i])

file = open("Data/varianceAll.txt", "w")

bdi = uni10.Bond(uni10.BD_IN, q_D)
bdo = uni10.Bond(uni10.BD_OUT, q_D)
bdi_pys = uni10.Bond(uni10.BD_IN, q_phys)

#print '\n','\n',bdi_pys.combine(bdi_pys.combine(bdi_pys))

Gamma_a=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_a")
Gamma_b=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_b")
Gamma_c=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_c")
Gamma_d=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_d")

Landa_1=uni10.UniTensor([bdi,bdo],"Landa_1")
Landa_2=uni10.UniTensor([bdi,bdo],"Landa_2")
Landa_3=uni10.UniTensor([bdi,bdo],"Landa_3")
Landa_4=uni10.UniTensor([bdi,bdo],"Landa_4")
Landa_5=uni10.UniTensor([bdi,bdo],"Landa_5")
Landa_6=uni10.UniTensor([bdi,bdo],"Landa_6")
Landa_7=uni10.UniTensor([bdi,bdo],"Landa_7")
Landa_8=uni10.UniTensor([bdi,bdo],"Landa_8")

Gamma=[Gamma_a,Gamma_b,Gamma_c,Gamma_d]
Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
basic.Initialize_function(Gamma,Landa)

#print  Landa_1
Landa1=[Landa[2],Landa[1],Landa[0],Landa[3]]
a_u,a=basic.makeab(Landa1,Gamma_a)
Landa1=[Landa[0],Landa[6],Landa[2],Landa[7]]
b_u,b=basic.makeab(Landa1,Gamma_b)
Landa1=[Landa[4],Landa[3],Landa[5],Landa[1]]
c_u,c=basic.makeab(Landa1,Gamma_c)
Landa1=[Landa[5],Landa[7],Landa[4],Landa[6]]
d_u,d=basic.makeab(Landa1,Gamma_d)





c1, c2,c3,c4=basic.makec1(q_chi,q_D)
Ta1, Tb1=basic.makeTab(q_chi,q_D)
Ta2, Tb2=basic.makeTab(q_chi,q_D)
Ta3, Tb3=basic.makeTab1(q_chi,q_D)
Ta4, Tb4=basic.makeTab1(q_chi,q_D)
Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]
Env1=basic.Rand_env_total(Env)
Env2=basic.Rand_env_total(Env)
Env3=basic.Rand_env_total(Env)

zlist=[]
J1_list=[h*0.0100 for h in range(270,400)]
J1_list=[0.0]
J2_list=[1.0]
h_list=[0.0]

for h, J2, J1 in zip( h_list, J2_list, J1_list):
 print h, J2, J1 
 h=[h, J1, J2]
#########################################################################################
 
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.Reload_itebd()
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,q_phys,D,N_iteritebd,h,Model,q_D)
 
 #basic.Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8)


# print Landa_1#,Landa_1[0],Landa_1[3],Landa_1[4],Landa_1[7]
# print Landa_2#,Landa_2[0],Landa_2[3],Landa_2[4],Landa_2[7]
# print Landa_3#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]
# print Landa_4#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]

# print Landa_5#,Landa_1[0],Landa_1[3],Landa_1[4],Landa_1[7]
# print Landa_6#,Landa_2[0],Landa_2[3],Landa_2[4],Landa_2[7]
# print Landa_7#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]
# print Landa_8#,Landa_8[0],Landa_8[3],Landa_8[4],Landa_8[7]



 Landa=[Landa_3,Landa_2,Landa_1,Landa_4]
 a_u,a=basic.makeab(Landa,Gamma_a)
 Landa=[Landa_1,Landa_7,Landa_3,Landa_8]
 b_u,b=basic.makeab(Landa,Gamma_b)
 Landa=[Landa_5,Landa_4,Landa_6,Landa_2]
 c_u,c=basic.makeab(Landa,Gamma_c)
 Landa=[Landa_6,Landa_8,Landa_5,Landa_7]
 d_u,d=basic.makeab(Landa,Gamma_d)
 #print "a",a.printDiagram()
 
 #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model)

 #print 'E_toal=', E_value
 #break
#########################################################################################

############################################################################
 Gauge='nFixed'
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.total_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()

 #print a_u.printDiagram(), b_u.printDiagram() 
 a_u,b_u,c_u,d_u,a,b,c,d,Env=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,q_phys,D,delta,h,Env,Env1,Env2,Env3,Gauge,Positive,Corner_method,N_iterF,Acc_E,Steps,Model)

 E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model)

# zlist1.append(z_value)
# Elist1.append(E_value)
 basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 print "E_final", E_value

###########################################################################################



##########################################################################################

##plt.plot( hlist, zlist,'b*',label='Ising, D=2, simple',markersize=np.sqrt(200.))
#plt.plot( hlist, zlist1,'g>',label='Ising, D=2, FU',markersize=np.sqrt(200.))
#plt.plot( hlist, zlist2,'r<',label='Ising, D=2, FU, Guage',markersize=np.sqrt(200.))
##plt.plot( hlist, zlist4,'m^',label='Ising, D=5',markersize=np.sqrt(200.))

#plt.xlabel('h', fontsize=25)
#plt.ylabel('Z', fontsize=25)
#plt.legend(loc='upper right')
#plt.savefig('Z.pdf')
#plt.show()
#plt.clf()

##plt.plot( hlist, Elist,'b*',label='Ising, D=2, simple',markersize=np.sqrt(200.))
#plt.plot( hlist, Elist1,'g>',label='Ising, D=2, FU',markersize=np.sqrt(200.))
#plt.plot( hlist, Elist2,'r<',label='Ising, D=2, FU, Guage',markersize=np.sqrt(200.))
##plt.plot( hlist, Elist4,'m^',label='Ising, D=5',markersize=np.sqrt(200.))
#plt.xlabel('h', fontsize=25)
#plt.ylabel('E', fontsize=25)
#plt.legend(loc='upper right')
#plt.savefig('E.pdf')
#plt.show()
#plt.clf()



