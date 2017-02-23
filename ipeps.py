import pyUni10 as uni10
import sys
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pylab
import random
import copy
import time
import basic
import itebd
import Fullupdate
import Move
###################### Initialize parameters ###########################
Model="Heisenberg"         #Heisenberg, Ising
D=[2,2]
chi=[10,10]
d_phys=[1,1]
N_iteritebd=200
N_iterF=40
Gauge='Non-Fixed'
Positive='Restrict'
Corner_method='CTM'   #CTM, CTMRG, CTMFull
Acc_E=1.00e-6
Steps=[1.0e-1,4.0e-2,4.0e-3,4.0e-4,5.0e-5,1.0e-6] #,[Start,steps,End] 
delta=0.001
###################################################################
zlist=[]
Elist=[]
zlist1=[]
Elist1=[]
zlist2=[]
Elist2=[]
q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
q_list=[q0_even,q0_odd]
qchi_list=[q0_even,q0_odd]

q_D=[]
q_D2=[]
q_chi=[]
q_phys=[q0_even,q0_odd]

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

#print  Gamma_a
Landa1=[Landa[2],Landa[1],Landa[0],Landa[3]]
a_u,a=basic.makeab(Landa1,Gamma_a)
Landa1=[Landa[0],Landa[6],Landa[2],Landa[7]]
b_u,b=basic.makeab(Landa1,Gamma_b)
Landa1=[Landa[4],Landa[3],Landa[5],Landa[1]]
c_u,c=basic.makeab(Landa1,Gamma_c)
Landa1=[Landa[5],Landa[7],Landa[4],Landa[6]]
d_u,d=basic.makeab(Landa1,Gamma_d)


Landa1=[Landa[2],Landa[1],Landa[0],Landa[3]]
ap_u,ap=basic.makeab(Landa1,Gamma_a)
Landa1=[Landa[0],Landa[6],Landa[2],Landa[7]]
bp_u,bp=basic.makeab(Landa1,Gamma_b)
Landa1=[Landa[4],Landa[3],Landa[5],Landa[1]]
cp_u,cp=basic.makeab(Landa1,Gamma_c)
Landa1=[Landa[5],Landa[7],Landa[4],Landa[6]]
dp_u,dp=basic.makeab(Landa1,Gamma_d)



c1, c2,c3,c4=basic.makec1(q_chi,q_D)
Ta1, Tb1=basic.makeTab(q_chi,q_D)
Ta2, Tb2=basic.makeTab(q_chi,q_D)
Ta3, Tb3=basic.makeTab1(q_chi,q_D)
Ta4, Tb4=basic.makeTab1(q_chi,q_D)
Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]

zlist=[]
hlist=[h*0.0100 for h in range(270,400)]
hlist=[0.00]

for h in hlist:
 print h

#########################################################################################
 
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.Reload_itebd()
 print  Landa_1
 Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,
Landa_8,chi,q_phys,D,N_iteritebd,h,Model,q_D)
 
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

 E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,q_phys,chi,Corner_method,Model)

 #print 'E_toal=', E_value

#########################################################################################

############################################################################
 Gauge='Non-Fixed'
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()

 a_u,b_u,c_u,d_u,a,b,c,d,Env=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,q_phys,D,delta,h,Env,Gauge,Positive,Corner_method,N_iterF,Acc_E,Steps,Model)

 E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,q_phys,chi,Corner_method,Model)
 z_value=basic.z_value(a,b,c,d,a_u,b_u,c_u,d_u,chi,D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Corner_method)
 zlist1.append(z_value)
 Elist1.append(E_value)
 basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
###########################################################################################



###########################################################################################
# ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp=basic.Reload_Fullp()
# basic.Store_Fullp(ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp)

# print 'E_toal=', E_value
# print 'z_value=', z_value
# Gauge='Non-Fixed'

# ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp,Env=Fullupdate.Full_Update(ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp,chi,d_phys,D,delta,h,Env,Gauge,Positive,Corner_method,N_iterF,Acc_E,Steps,Model)

# E_value=basic.E_total(ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi,Corner_method,Model)
# z_value=basic.z_value(ap,bp,cp,dp,ap_u,bp_u,cp_u,dp_u,chi,D*D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Corner_method)

# print 'E_toal=', E_value
# print 'z_value=', z_value

# zlist2.append(z_value)
# Elist2.append(E_value)
# basic.Store_Fullp(ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp)

# basic.Store(hlist,zlist1, zlist1,zlist2,Elist1, Elist1 , Elist2 , file)

##########################################################################################

#plt.plot( hlist, zlist,'b*',label='Ising, D=2, simple',markersize=np.sqrt(200.))
plt.plot( hlist, zlist1,'g>',label='Ising, D=2, FU',markersize=np.sqrt(200.))
plt.plot( hlist, zlist2,'r<',label='Ising, D=2, FU, Guage',markersize=np.sqrt(200.))
#plt.plot( hlist, zlist4,'m^',label='Ising, D=5',markersize=np.sqrt(200.))

plt.xlabel('h', fontsize=25)
plt.ylabel('Z', fontsize=25)
plt.legend(loc='upper right')
plt.savefig('Z.pdf')
plt.show()
plt.clf()

#plt.plot( hlist, Elist,'b*',label='Ising, D=2, simple',markersize=np.sqrt(200.))
plt.plot( hlist, Elist1,'g>',label='Ising, D=2, FU',markersize=np.sqrt(200.))
plt.plot( hlist, Elist2,'r<',label='Ising, D=2, FU, Guage',markersize=np.sqrt(200.))
#plt.plot( hlist, Elist4,'m^',label='Ising, D=5',markersize=np.sqrt(200.))
plt.xlabel('h', fontsize=25)
plt.ylabel('E', fontsize=25)
plt.legend(loc='upper right')
plt.savefig('E.pdf')
plt.show()
plt.clf()



