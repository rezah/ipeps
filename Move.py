import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time
import TruncateU

def norm_CTM(c):

 if ( (abs(c.getBlock().absMax()) < 0.50e-1) or (abs(c.getBlock().absMax()) > 0.50e+1)   ):
  c*=(1.00/c.getBlock().absMax()); 
 return c;


def  add_down1(c4,Ta3,Tb3,c3,Ta4,Tb2,c,d,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]
#################################################
 c4.setLabel([0,1])
 Ta4.setLabel([0,2,3])
 c4bar=c4*Ta4
 c4bar.permute([1,2,3],2)
 
 c3.setLabel([1,0])
 Tb2.setLabel([0,2,3])
 c3bar=c3*Tb2
 c3bar.permute([3,1,2],1)
 
 c_f=copy.copy(c4bar)
 c_f_trans=copy.copy(c4bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(c3bar)
 c_s_trans=copy.copy(c3bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)
 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()
 
 
#############################
 
 c4.setLabel([0,1])
 Ta4.setLabel([0,4,5])
 Ta3.setLabel([1,2,3])
 c.setLabel([4,2,6,7])
 Q=((c4*Ta4)*(Ta3))*c
 Q.permute([3,6,5,7],2)
 Q.combineBond([5,7])
 
 c3.setLabel([1,0])
 Tb2.setLabel([0,2,3])
 Tb3.setLabel([7,4,1])
 d.setLabel([5,4,2,6])
 Q1=((c3*Tb2)*(Tb3))*d
 Q1.permute([3,6,7,5],2)
 Q1.combineBond([3,6])
 
 
 c_f=copy.copy(Q)
 c_f_trans=copy.copy(Q)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(Q1)
 c_s_trans=copy.copy(Q1)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 W_trans=copy.copy(W)
 W_trans.transpose()

 ##############################
 Z_trans.setLabel([4,1,2])
 c4bar=c4bar*Z_trans
 c4bar.permute([3,4],2)
 ############################
 Z.setLabel([1,2,4])
 c3bar=c3bar*Z
 c3bar.permute([4,3],1)
 #############################
 Ta3.setLabel([0,1,2])
 c.setLabel([3,1,4,5])
 Z.setLabel([0,3,7])
 W_trans.setLabel([6,2,4])
 Ta3bar=((Ta3*c)*Z)*W_trans
 Ta3bar.permute([7,5,6],1)
 ###########################
 Tb3.setLabel([0,1,2])
 d.setLabel([3,1,4,5])
 W.setLabel([0,3,6])
 Z_trans.setLabel([7,2,4])
 Tb3bar=((Tb3*d)*Z_trans)*W
 Tb3bar.permute([6,5,7],1)
 ###########################
 
 
 c3=norm_CTM(c3bar)
 c4=norm_CTM(c4bar)
 Ta3=norm_CTM(Ta3bar)
 Tb3=norm_CTM(Tb3bar)

 return c4, Ta3, Tb3, c3




def add_up1(c1,Tb1,Ta1,c2,Tb4,Ta2,a,b,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]


 c1.setLabel([0,1])
 Tb4.setLabel([3,2,0])
 c1bar=c1*Tb4
 c1bar.permute([1,2,3],2)

 c2.setLabel([0,1])
 Ta2.setLabel([3,2,1])
 c2bar=c2*Ta2
 c2bar.permute([3,0,2],1)

 c_f=copy.copy(c1bar)
 c_f_trans=copy.copy(c1bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(c2bar)
 c_s_trans=copy.copy(c2bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)
 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()
######################################

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 Tb4.setLabel([4,5,0])
 a.setLabel([5,6,7,2])
 Q=((c1*Tb1)*(Tb4))*a
 Q.permute([3,7,4,6],2)
 Q.combineBond([4,6])



 c2.setLabel([0,1])
 Ta2.setLabel([3,2,1])
 Ta1.setLabel([4,5,0])
 b.setLabel([6,7,2,5])
 Q1=((c2*Ta2)*(Ta1))*b
 Q1.permute([3,7,4,6],2)
 Q1.combineBond([3,7])



 c_f=copy.copy(Q)
 c_f_trans=copy.copy(Q)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(Q1)
 c_s_trans=copy.copy(Q1)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 W_trans=copy.copy(W)
 W_trans.transpose()


 ##############################
 Z_trans.setLabel([4,1,2])
 c1bar=c1bar*Z_trans
 c1bar.permute([3,4],1)
 ############################
 Z.setLabel([0,2,4])
 c2bar=c2bar*Z
 c2bar.permute([4,3],0)
 #############################
 Tb1.setLabel([0,1,2])
 a.setLabel([3,4,5,1])
 Z.setLabel([0,3,7])
 W_trans.setLabel([6,2,5])
 Tb1bar=((Tb1*a)*Z)*W_trans
 Tb1bar.permute([7,4,6],2)
#####################################
 Ta1.setLabel([0,1,2])
 b.setLabel([3,4,5,1])
 Z_trans.setLabel([7,2,5])
 W.setLabel([0,3,6])
 Ta1bar=((Ta1*b)*Z_trans)*W
 Ta1bar.permute([6,4,7],2)
##############################################3

 c1=norm_CTM(c1bar)
 c2=norm_CTM(c2bar)
 Ta1=norm_CTM(Ta1bar)
 Tb1=norm_CTM(Tb1bar)

 return  c1, Ta1, Tb1, c2





def  add_right1(c2,Ta2,Tb2,c3,Ta1,Tb3,b,d,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]
 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,3],2)
 
 c3.setLabel([1,0])
 Tb3.setLabel([3,2,1])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2],1)

 c_f=copy.copy(c2bar)
 c_f_trans=copy.copy(c2bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(c3bar)
 c_s_trans=copy.copy(c3bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)
 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

################################################################
 c2.setLabel([0,1])
 Ta2.setLabel([3,2,1])
 Ta1.setLabel([4,5,0])
 b.setLabel([6,7,2,5])
 Q=((c2*Ta2)*(Ta1))*b
 Q.permute([3,7,4,6],2)
 Q.combineBond([4,6])
 
 c3.setLabel([1,0])
 Tb2.setLabel([0,2,3])
 Tb3.setLabel([5,4,1])
 d.setLabel([7,4,2,6])
 Q1=((c3*Tb2)*(Tb3))*d
 Q1.permute([5,7,3,6],2)
 Q1.combineBond([5,7])
 

 c_f=copy.copy(Q)
 c_f_trans=copy.copy(Q)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(Q1)
 c_s_trans=copy.copy(Q1)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 W_trans=copy.copy(W)
 W_trans.transpose()


 ##############################
 Z_trans.setLabel([4,1,2])
 c2bar=c2bar*Z_trans
 c2bar.permute([3,4],2)
 ############################
 Z.setLabel([0,2,4])
 c3bar=c3bar*Z
 c3bar.permute([3,4],1)
 #############################
 Ta2.setLabel([0,1,2])
 b.setLabel([4,5,1,3])
 Z.setLabel([2,3,7])
 W_trans.setLabel([6,0,5])
 Ta2bar=((Ta2*b)*Z)*W_trans
 Ta2bar.permute([6,4,7],2)
 ###########################
 Tb2.setLabel([0,1,2])
 d.setLabel([4,5,1,3])
 W.setLabel([2,3,7])
 Z_trans.setLabel([6,0,5])
 Tb2bar=((Tb2*d)*Z_trans)*W
 Tb2bar.permute([6,4,7],2)
 ###########################
 ###########################

 c3=norm_CTM(c3bar)
 c2=norm_CTM(c2bar)
 Ta2=norm_CTM(Ta2bar)
 Tb2=norm_CTM(Tb2bar)


 return c2, Ta2, Tb2, c3









def  add_right(c2,Ta2,Tb2,c3,Ta1,Tb3,b,d,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]
 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,3],2)
 
 c3.setLabel([0,1])
 Tb3.setLabel([1,2,3])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2],1)

 c_f=copy.copy(c2bar)
 c_f_trans=copy.copy(c2bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(c3bar)
 c_s_trans=copy.copy(c3bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)
 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

 Z.setLabel([0,1,2])
 Z.permute([1,2,0],2)
 
 Z_trans.setLabel([0,1,2])
 Z_trans.permute([1,2,0],1)
 
################################################################
 c2.setLabel([0,1])
 Ta2.setLabel([1,2,3])
 Ta1.setLabel([4,5,0])
 b.setLabel([6,7,2,5])
 Q=((c2*Ta2)*(Ta1))*b
 Q.permute([3,7,4,6],2)
# Q.combineBond([4,6])

 
 c3.setLabel([-1,-4])
 Tb2.setLabel([3,-2,-1])
 Tb3.setLabel([-4,-3,6])
 d.setLabel([4,-3,-2,7])
 Q1=((c3*Tb2)*(Tb3))*d
 Q1.permute([4,6,3,7],2)
 #Q1.combineBond([5,7])
 

 Q_trans=copy.copy(Q)
 Q_trans.transpose()
 Q.setLabel([3,7,4,6])
 Q_trans.setLabel([4,6,-3,-7])
 Q.permute([4,6,7,3],3)
 Q_trans.permute([-3,4,6,-7],1)
 Q=Q*Q_trans
 Q.permute([3,7,-3,-7],2)

 Q1_trans=copy.copy(Q1)
 Q1_trans.transpose()
 Q1.setLabel([4,6,3,7])
 Q1_trans.setLabel([-3,-7,4,6])
 Q1.permute([4,3,6,7],2)
 Q1_trans.permute([6,-7,4,-3],2)
 Q1=Q1_trans*Q1
 Q1.permute([-3,-7,3,7],2)
 Q1.setLabel([3,7,-3,-7])

 theta=Q+Q1
 theta.setLabel([0,1,2,3])
 theta.permute([0,1,2,3],2)

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)


 W_trans=copy.copy(W)
 W_trans.transpose()

 W.setLabel([0,1,2])
 W.permute([1,2,0],2)

 W_trans.setLabel([0,1,2])
 W_trans.permute([1,2,0],1)


 ##############################
 Z_trans.setLabel([1,2,4])
 c2bar=c2bar*Z_trans
 c2bar.permute([3,4],1)
 ############################
 Z.setLabel([2,4,0])
 c3bar=c3bar*Z
 c3bar.permute([4,3],1)
 #############################
 Ta2.setLabel([2,1,0])
 b.setLabel([4,5,1,3])
 Z.setLabel([3,7,2])
 W_trans.setLabel([0,5,6])
 Ta2bar=((Ta2*b)*Z)*W_trans
 Ta2bar.permute([7,4,6],2)
 ###########################
 Tb2.setLabel([2,1,0])
 d.setLabel([4,5,1,3])
 W.setLabel([3,7,2])
 Z_trans.setLabel([0,5,6])
 Tb2bar=((Tb2*d)*Z_trans)*W
 Tb2bar.permute([7,4,6],2)
 ###########################
 ###########################


 c3=norm_CTM(c3bar)
 c2=norm_CTM(c2bar)
 Ta2=norm_CTM(Ta2bar)
 Tb2=norm_CTM(Tb2bar)

 return c2, Ta2, Tb2, c3


def  add_left(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]


 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,3],2)

 c4.setLabel([0,1])
 Ta3.setLabel([3,2,0])
 c4bar=c4*Ta3
 c4bar.permute([3,1,2],1)

 c_f=copy.copy(c1bar)
 c_f_trans=copy.copy(c1bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(c4bar)
 c_s_trans=copy.copy(c4bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])

 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)

 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

####################################

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 Tb4.setLabel([4,5,0])
 a.setLabel([5,6,7,2])
 Q=((c1*Tb1)*(Tb4))*a
 Q.permute([4,6,3,7],2)
 Q.combineBond([3,7])

 c4.setLabel([-1,-2])
 Ta3.setLabel([1,-3,-1])
 Ta4.setLabel([-2,-4,4])
 c.setLabel([-4,-3,0,6])
 Q1=((c4*Ta3)*(Ta4))*c
 Q1.permute([0,1,4,6],2)

 #Q1.combineBond([5,6])



 Q_trans=copy.copy(Q)
 Q_trans.transpose()
 Q.setLabel([4,6,2])
 Q_trans.setLabel([2,-4,-6])
 Q=Q*Q_trans
 Q.permute([4,6,-4,-6],2)

 Q1_trans=copy.copy(Q1)
 Q1_trans.transpose()
 Q1_trans.setLabel([-4,-6,0,1])
 Q1.permute([1,0,4,6],1)
 Q1_trans.permute([0,-4,-6,1],3)
 Q1=Q1*Q1_trans
 Q1.permute([-4,-6,4,6],2)
 Q1.setLabel([4,6,-4,-6])
 theta=Q+Q1
 theta.setLabel([0,1,2,3])

##################################################

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 W_trans=copy.copy(W)
 W_trans.transpose()

 ##############################
 Z_trans.setLabel([4,0,2])
 c1bar=c1bar*Z_trans
 c1bar.permute([4,3],1)
# print 'c1',c1bar.printDiagram()
 ############################
 Z.setLabel([1,2,4])
 c4bar=c4bar*Z
 c4bar.permute([3,4],1)
 #############################
 Tb4.setLabel([0,1,2])
 a.setLabel([1,3,4,5])
 Z.setLabel([2,5,6])
 W_trans.setLabel([7,0,3])
 Tb4bar=((Tb4*a)*Z)*W_trans
 Tb4bar.permute([7,4,6],2)
 ###########################
 Ta4.setLabel([0,1,2])
 c.setLabel([1,3,4,5])
 W.setLabel([2,5,6])
 Z_trans.setLabel([7,0,3])
 Ta4bar=((Ta4*c)*Z_trans)*W
 Ta4bar.permute([7,4,6],2)
#############################


 c1=norm_CTM(c1bar)
 c4=norm_CTM(c4bar)
 Ta4=norm_CTM(Ta4bar)
 Tb4=norm_CTM(Tb4bar)


 return c1, Ta4, Tb4, c4
 










def  add_left1(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]


 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,3],2)

 c4.setLabel([0,1])
 Ta3.setLabel([1,2,3])
 c4bar=c4*Ta3
 c4bar.permute([3,0,2],1)

 c_f=copy.copy(c1bar)
 c_f_trans=copy.copy(c1bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(c4bar)
 c_s_trans=copy.copy(c4bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])

 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)

 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

####################################

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 Tb4.setLabel([4,5,0])
 a.setLabel([5,6,7,2])
 Q=((c1*Tb1)*(Tb4))*a
 Q.permute([4,6,3,7],2)
 Q.combineBond([3,7])

 c4.setLabel([0,1])
 Ta3.setLabel([1,4,5])
 Ta4.setLabel([0,2,3])
 c.setLabel([2,4,6,7])
 Q1=((c4*Ta3)*(Ta4))*c
 Q1.permute([5,6,3,7],2)
 Q1.combineBond([5,6])


 c_f=copy.copy(Q)
 c_f_trans=copy.copy(Q)
 c_f_trans.transpose()
 c_f.setLabel([0,1,2])
 c_f_trans.setLabel([2,3,4])
 c_s=copy.copy(Q1)
 c_s_trans=copy.copy(Q1)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4])
 c_s_trans.setLabel([0,1,2])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,3,4],2)

##################################################

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 W_trans=copy.copy(W)
 W_trans.transpose()

 ##############################
 Z_trans.setLabel([4,0,2])
 c1bar=c1bar*Z_trans
 c1bar.permute([4,3],1)
 ############################
 Z.setLabel([0,2,4])
 c4bar=c4bar*Z
 c4bar.permute([4,3],0)
 #############################
 Tb4.setLabel([0,1,2])
 a.setLabel([1,3,4,5])
 Z.setLabel([2,5,6])
 W_trans.setLabel([7,0,3])
 Tb4bar=((Tb4*a)*Z)*W_trans
 Tb4bar.permute([7,4,6],1)
 ###########################
 Ta4.setLabel([0,1,2])
 c.setLabel([1,3,4,5])
 W.setLabel([2,5,6])
 Z_trans.setLabel([7,0,3])
 Ta4bar=((Ta4*c)*Z_trans)*W
 Ta4bar.permute([7,4,6],1)
#############################
 c1=norm_CTM(c1bar)
 c4=norm_CTM(c4bar)
 Ta4=norm_CTM(Ta4bar)
 Tb4=norm_CTM(Tb4bar)
###############################



 return c1, Ta4, Tb4, c4


 
def  magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):


 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([4,24])
 c4.setLabel([22,18])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([7,10,14]) 
 Ta3.setLabel([23,19,22]) 
 Ta4.setLabel([18,15,11]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([14,17,4]) 
 Tb3.setLabel([24,20,23]) 
 Tb4.setLabel([11,8,4])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([15,19,16,12])
 d.setLabel([16,20,17,13])
 norm=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*Ta1)*b))*(((c1*Tb1)*Tb4)*a)
 
 #print 'hi1', '\n',norm,
 return norm
 
def  magnetization_value1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):


 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([24,4])
 c4.setLabel([18,22])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta3.setLabel([22,19,23]) 
 Ta4.setLabel([18,15,11]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([4,17,14]) 
 Tb3.setLabel([23,20,24]) 
 Tb4.setLabel([11,8,4])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([15,19,16,12])
 d.setLabel([16,20,17,13])
 norm=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*Ta1)*b))*(((c1*Tb1)*Tb4)*a) 
 #print 'hi1', '\n',norm,
 return norm
 
 
 
def permute(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,1,2,3])
 a.permute([3,2,1,0],2)
 a.setLabel([0,1,2,3])
 ##print'a', a
 b.setLabel([0,1,2,3])
 b.permute([3,2,1,0],2)
 b.setLabel([0,1,2,3])

 c.setLabel([0,1,2,3])
 c.permute([3,2,1,0],2)
 c.setLabel([0,1,2,3])

 d.setLabel([0,1,2,3])
 d.permute([3,2,1,0],2)
 d.setLabel([0,1,2,3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([1,0],1)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([1,0],1)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,2])
 Ta1.permute([2,1,0],1)
 Ta1.setLabel([0,1,2])
 
 Ta2.setLabel([0,1,2])
 Ta2.permute([2,1,0],1)
 Ta2.setLabel([0,1,2])
 
 Ta3.setLabel([0,1,2])
 Ta3.permute([2,1,0],2)
 Ta3.setLabel([0,1,2])

 Ta4.setLabel([0,1,2])
 Ta4.permute([2,1,0],2)
 Ta4.setLabel([0,1,2])

 Tb1.setLabel([0,1,2])
 Tb1.permute([2,1,0],1)
 Tb1.setLabel([0,1,2])
 
 Tb2.setLabel([0,1,2])
 Tb2.permute([2,1,0],1)
 Tb2.setLabel([0,1,2])
 
 Tb3.setLabel([0,1,2])
 Tb3.permute([2,1,0],2)
 Tb3.setLabel([0,1,2])

 Tb4.setLabel([0,1,2])
 Tb4.permute([2,1,0],2)
 Tb4.setLabel([0,1,2])

 
def permute1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,1,2,3])
 a.permute([3,2,1,0],2)
 a.setLabel([0,1,2,3])
 ##print'a', a
 b.setLabel([0,1,2,3])
 b.permute([3,2,1,0],2)
 b.setLabel([0,1,2,3])

 c.setLabel([0,1,2,3])
 c.permute([3,2,1,0],2)
 c.setLabel([0,1,2,3])

 d.setLabel([0,1,2,3])
 d.permute([3,2,1,0],2)
 d.setLabel([0,1,2,3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([1,0],1)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([1,0],1)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,2])
 Ta1.permute([2,1,0],2)
 Ta1.setLabel([0,1,2])
 
 Ta2.setLabel([0,1,2])
 Ta2.permute([2,1,0],2)
 Ta2.setLabel([0,1,2])
 
 Ta3.setLabel([0,1,2])
 Ta3.permute([2,1,0],1)
 Ta3.setLabel([0,1,2])

 Ta4.setLabel([0,1,2])
 Ta4.permute([2,1,0],1)
 Ta4.setLabel([0,1,2])

 Tb1.setLabel([0,1,2])
 Tb1.permute([2,1,0],2)
 Tb1.setLabel([0,1,2])
 
 Tb2.setLabel([0,1,2])
 Tb2.permute([2,1,0],2)
 Tb2.setLabel([0,1,2])
 
 Tb3.setLabel([0,1,2])
 Tb3.permute([2,1,0],1)
 Tb3.setLabel([0,1,2])

 Tb4.setLabel([0,1,2])
 Tb4.permute([2,1,0],1)
 Tb4.setLabel([0,1,2])
 
 
def permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,1,2,3])
 a.permute([3,2,1,0],2)
 a.setLabel([0,1,2,3])
 ##print'a', a
 b.setLabel([0,1,2,3])
 b.permute([3,2,1,0],2)
 b.setLabel([0,1,2,3])

 c.setLabel([0,1,2,3])
 c.permute([3,2,1,0],2)
 c.setLabel([0,1,2,3])

 d.setLabel([0,1,2,3])
 d.permute([3,2,1,0],2)
 d.setLabel([0,1,2,3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([0,1],2)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([0,1],0)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,2])
 Ta1.permute([2,1,0],1)
 Ta1.setLabel([0,1,2])
 
 Ta2.setLabel([0,1,2])
 Ta2.permute([2,1,0],1)
 Ta2.setLabel([0,1,2])
 
 Ta3.setLabel([0,1,2])
 Ta3.permute([2,1,0],2)
 Ta3.setLabel([0,1,2])

 Ta4.setLabel([0,1,2])
 Ta4.permute([2,1,0],2)
 Ta4.setLabel([0,1,2])

 Tb1.setLabel([0,1,2])
 Tb1.permute([2,1,0],1)
 Tb1.setLabel([0,1,2])
 
 Tb2.setLabel([0,1,2])
 Tb2.permute([2,1,0],1)
 Tb2.setLabel([0,1,2])
 
 Tb3.setLabel([0,1,2])
 Tb3.permute([2,1,0],2)
 Tb3.setLabel([0,1,2])

 Tb4.setLabel([0,1,2])
 Tb4.permute([2,1,0],2)
 Tb4.setLabel([0,1,2])

 
def permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,1,2,3])
 a.permute([3,2,1,0],2)
 a.setLabel([0,1,2,3])
 ##print'a', a
 b.setLabel([0,1,2,3])
 b.permute([3,2,1,0],2)
 b.setLabel([0,1,2,3])

 c.setLabel([0,1,2,3])
 c.permute([3,2,1,0],2)
 c.setLabel([0,1,2,3])

 d.setLabel([0,1,2,3])
 d.permute([3,2,1,0],2)
 d.setLabel([0,1,2,3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([0,1],2)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([0,1],0)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,2])
 Ta1.permute([2,1,0],2)
 Ta1.setLabel([0,1,2])
 
 Ta2.setLabel([0,1,2])
 Ta2.permute([2,1,0],2)
 Ta2.setLabel([0,1,2])
 
 Ta3.setLabel([0,1,2])
 Ta3.permute([2,1,0],1)
 Ta3.setLabel([0,1,2])

 Ta4.setLabel([0,1,2])
 Ta4.permute([2,1,0],1)
 Ta4.setLabel([0,1,2])

 Tb1.setLabel([0,1,2])
 Tb1.permute([2,1,0],2)
 Tb1.setLabel([0,1,2])
 
 Tb2.setLabel([0,1,2])
 Tb2.permute([2,1,0],2)
 Tb2.setLabel([0,1,2])
 
 Tb3.setLabel([0,1,2])
 Tb3.permute([2,1,0],1)
 Tb3.setLabel([0,1,2])

 Tb4.setLabel([0,1,2])
 Tb4.permute([2,1,0],1)
 Tb4.setLabel([0,1,2])
 


def  distance(c1,c2,c3,c4,c1_f,c2_f,c3_f,c4_f):

 s1=c1.getBlock().svd()
 s1_f=c1_f.getBlock().svd()
 dis_val1=dis(s1[1],s1_f[1] )
 #print dis_val1
 s1=c2.getBlock().svd()
 s1_f=c2_f.getBlock().svd()
 dis_val2=dis(s1[1],s1_f[1] )
 #print dis_val2
 s1=c3.getBlock().svd()
 s1_f=c3_f.getBlock().svd()
 dis_val3=dis(s1[1],s1_f[1] )
 #print dis_val3
 s1=c4.getBlock().svd()
 s1_f=c4_f.getBlock().svd()
 dis_val4=dis(s1[1],s1_f[1] )
 #print dis_val4
 return (dis_val1+dis_val2+dis_val3+dis_val4) / 4.00



def dis(s1,s2):
 sum=0
 for q in xrange(int(s1.row())):
  sum=sum+abs(s1[q]-s2[q])
 return sum


#print '\n','\n','\n'
#Move.corner_transfer_matrix_onesite(a, z,chi,D)
 
 
 
 
 
 
 
 
