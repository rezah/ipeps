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
#import env
#import optimize 
#import plotdata


###################### Initialize parameters ###########################
chi=40
d=2
D=2
N_iter=100
delta=0.1
###################################################################


bdi = uni10.Bond(uni10.BD_IN, D)
bdo = uni10.Bond(uni10.BD_OUT, D)
bdi_pys = uni10.Bond(uni10.BD_IN, d)

Gamma_a=uni10.UniTensor(uni10.CTYPE,[bdi_pys,bdi,bdi,bdo,bdo], "Gamma_a")
Gamma_b=uni10.UniTensor(uni10.CTYPE,[bdi_pys,bdi,bdi,bdo,bdo], "Gamma_b")
Gamma=[Gamma_a,Gamma_b]


Landa_1=uni10.UniTensor(uni10.CTYPE,[bdi,bdo],"Landa_1")
Landa_2=uni10.UniTensor(uni10.CTYPE,[bdi,bdo],"Landa_2")
Landa_3=uni10.UniTensor(uni10.CTYPE,[bdi,bdo],"Landa_3")
Landa_4=uni10.UniTensor(uni10.CTYPE,[bdi,bdo],"Landa_4")
Landa=[Landa_1,Landa_2,Landa_3,Landa_4]


basic.Initialize_function(Gamma,Landa)
a_u,b_u, a, b=basic.makeab(Landa,Gamma_a, Gamma_b)
c1, c2,c3,c4=basic.makec1(a, b, chi)
Ta1, Tb1=basic.makeTab(a, b, chi,D)
Ta2, Tb2=basic.makeTab(a, b, chi,D)
Ta3, Tb3=basic.makeTab(a, b, chi,D)
Ta4, Tb4=basic.makeTab(a, b, chi,D)
#print Landa[0], Landa[1],Landa[2],Landa[3]
#print Gamma[0]
zlist=[]
hlist=[h*0.05 for h in range(80)]

for h in hlist:
 print h
 H0=basic.transverseIsing(h)
 U = uni10.UniTensor(H0.bond(), "U");

 for i in xrange(1,100):
  delta=1.00/pow(2.00,i)
  print 'delta =', delta
  if delta>1.0e-1:
   N_iter=300
  if delta<1.0e-1 and delta>1.0e-5:
   N_iter=300
  if delta<1.0e-5  and delta>1.0e-6:
   N_iter=300
  if delta<1.0e-6:
   break

  print "N_iter=", N_iter


  for q in xrange(N_iter):

   Gamma=[Gamma_a,Gamma_b]
   Landa=[Landa_1,Landa_2,Landa_3,Landa_4]
   U.putBlock(uni10.takeExp(-delta, H0.getBlock()))
 #rlink
   basic.update_rlink(Gamma,Landa,U,D,d)

 #llink
   Gammap=[Gamma_b, Gamma_a]
   Landap=[Landa_3,Landa_4,Landa_1,Landa_2]
   basic.update_rlink(Gammap,Landap,U,D,d)
   #basic.update_llink(Gamma,Landa,U,D,d)

 #ulink
   Gamma_ap=copy.copy(Gamma_a)
   Gamma_bp=copy.copy(Gamma_b)
   Gamma_ap.permute([0,4,1,2,3],3)
   Gamma_bp.permute([0,4,1,2,3],3)
   Gamma_ap.setLabel([0,1,2,3,4])
   Gamma_bp.setLabel([0,1,2,3,4])
   Gammap=[Gamma_ap, Gamma_bp]

   Landap=[Landa_2,Landa_3,Landa_4,Landa_1]
   basic.update_rlink(Gammap,Landap,U,D,d)

 #dlink
   Gammap=[Gamma_bp, Gamma_ap]
   Landap=[Landa_4,Landa_1,Landa_2,Landa_3]
   basic.update_rlink(Gammap,Landap,U,D,d)
   Gamma_ap.permute([0,2,3,4,1],3)
   Gamma_bp.permute([0,2,3,4,1],3)
   Gamma_ap.setLabel([0,1,2,3,4])
   Gamma_bp.setLabel([0,1,2,3,4])
   Gamma_a=Gamma_ap
   Gamma_b=Gamma_bp
 #print Landa[0],Landa[1],Landa[2] ,Landa[3], Landa[0][0], Landa[1][0], Landa[2][0], Landa[3][0]
 #print Landa[0][5*4+4], Landa[1][5*4+4], Landa[2][5*4+4], Landa[3][5*4+4]
 #print Landa[0][5*2+2], Landa[1][5*2+2], Landa[2][5*2+2], Landa[3][5*2+2]
 #print Landa[0][5*3+3], Landa[1][5*3+3], Landa[2][5*3+3], Landa[3][5*3+3]
 #print U
 #print U

 a_u,b_u, a, b=basic.makeab(Landa,Gamma_a, Gamma_b)
 ############################## CTM ######################################
 ########################### initialize ############################
 print a.similar(b)#,a,b,a[10], b[10]
 az, bz=basic.magnetization(a_u, b_u)
 ########################################################################################
 z_value=basic.corner_transfer_matrix_onesite(a, az,chi,D,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 zlist.append(z_value.real)
 #basic.corner_transfer_matrix_onesite(a, bz,chi,D)
print 'hi'
plt.plot( hlist, zlist,'g.')
print 'hi'
plt.xlabel('h', fontsize=20)
plt.ylabel('Z', fontsize=20)
#plt.legend(loc='upper right')
plt.savefig('Z.pdf')
plt.show()
plt.clf()


 
