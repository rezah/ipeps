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
chi=20
d_phys=2
D=2
N_iter=200
delta=0.001
###################################################################

bdi = uni10.Bond(uni10.BD_IN, D)
bdo = uni10.Bond(uni10.BD_OUT, D)
bdi_pys = uni10.Bond(uni10.BD_IN, d_phys)

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

Landa=[Landa_3,Landa_2,Landa_1,Landa_4]

a_u,a=basic.makeab(Landa,Gamma_a)
Landa=[Landa_1,Landa_7,Landa_3,Landa_8]
b_u,b=basic.makeab(Landa,Gamma_b)
Landa=[Landa_5,Landa_4,Landa_6,Landa_2]
c_u,c=basic.makeab(Landa,Gamma_c)
Landa=[Landa_5,Landa_4,Landa_6,Landa_2]
d_u,d=basic.makeab(Landa,Gamma_d)


zlist=[]
hlist=[h*0.05 for h in range(80)]
hlist=[3.0]

for h in hlist:
 print h

 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa6, Landa7,Landa8=itebd.itebd_standard(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iter,delta,h)
 
 Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa6, Landa7,Landa8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iter,delta,h)
 
 

 Landa=[Landa_3,Landa_2,Landa_1,Landa_4]
 a_u,a=basic.makeab(Landa,Gamma_a)
 Landa=[Landa_1,Landa_7,Landa_3,Landa_8]
 b_u,b=basic.makeab(Landa,Gamma_b)
 Landa=[Landa_5,Landa_4,Landa_6,Landa_2]
 c_u,c=basic.makeab(Landa,Gamma_c)
 Landa=[Landa_5,Landa_4,Landa_6,Landa_2]
 d_u,d=basic.makeab(Landa,Gamma_d)

a_u,b_u,c_u,d_u,a,b,c,d=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,d_phys,D,N_iter,delta,h)

c1, c2,c3,c4=basic.makec1(chi,D*D)
Ta1, Tb1=basic.makeTab(chi,D*D)
Ta2, Tb2=basic.makeTab(chi,D*D)
Ta3, Tb3=basic.makeTab(chi,D*D)
Ta4, Tb4=basic.makeTab( chi,D*D)
az=basic.magnetization(a_u)
bz=basic.magnetization(b_u)
cz=basic.magnetization(c_u)
dz=basic.magnetization(d_u)
Truncation=[0]
c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, Truncation=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,Truncation)
E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,h,d_phys,chi)
z_value=basic.z_value(a,b,c,d,az,bz,cz,dz,chi,D*D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

print 'E_toal=', E_value
print 'z_value=', z_value



