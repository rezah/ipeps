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

def MaxAbs(c):
 blk_qnums = c.blockQnum()
 max_list=[]
 for qnum in blk_qnums:
    c_mat=c.getBlock(qnum)
    max_list.append(c_mat.absMax())
 #sv_mat = uni10.Matrix( len(max_list), len(max_list), max_list, True)
 #return sv_mat.absMax()
 max_list_f=[abs(x) for x in max_list]
 #print max_list_f, max(max_list_f)
 return max(max_list_f)

def norm_CTM(c):

 if ( (abs(MaxAbs(c)) < 0.50e-1) or (abs(MaxAbs(c)) > 0.50e+1)   ):
  c*=(1.00/MaxAbs(c)); 
 return c;




def  add_right1(c2,Ta2,Tb2,c3,Ta1,Tb3,b,d,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]
 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,-2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,-2,3],3)
 
 c3.setLabel([1,0])
 Tb3.setLabel([3,2,-2,1])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2,-2],1)

 c_f=copy.copy(c2bar)
 c_f_trans=copy.copy(c2bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,-1,2])
 c_f_trans.setLabel([2,3,4,-4])
 c_s=copy.copy(c3bar)
 c_s_trans=copy.copy(c3bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4,-4])
 c_s_trans.setLabel([0,1,-1,2])

 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,-1,3,4,-4],3)

 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

################################################################
 c2.setLabel([0,1])
 Ta2.setLabel([3,2,-2,1])
 Ta1.setLabel([4,5,-5,0])
 b.setLabel([6,-6,7,-7,2,-2,5,-5])
 Q=((c2*Ta2)*(Ta1))*b
 Q.permute([3,7,-7,4,6,-6],3)
 #Q.combineBond([4,6])
 
 c3.setLabel([1,0])
 Tb2.setLabel([0,2,-2,3])
 Tb3.setLabel([5,4,-4,1])
 d.setLabel([7,-7,4,-4,2,-2,6,-6])
 Q1=((c3*Tb2)*(Tb3))*d
 Q1.permute([5,7,-7,3,6,-6],3)
 #Q1.combineBond([5,7])
 

 c_f=copy.copy(Q)
 c_f_trans=copy.copy(Q)
 c_f_trans.transpose()
 c_f.setLabel([0,1,-1,2,3,-3])
 c_f_trans.setLabel([2,3,-3,4,5,-5])
 c_s=copy.copy(Q1)
 c_s_trans=copy.copy(Q1)
 c_s_trans.transpose()
 c_s.setLabel([2,3,-3,4,5,-5])
 c_s_trans.setLabel([0,1,-1,2,3,-3])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,-1,4,5,-5],3)

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 W_trans=copy.copy(W)
 W_trans.transpose()



 ##############################
 Z_trans.setLabel([4,1,2,-2])
 c2bar=c2bar*Z_trans
 c2bar.permute([3,4],2)
 ############################
 Z.setLabel([0,2,-2,4])
 c3bar=c3bar*Z
 c3bar.permute([3,4],1)
 #############################
 Ta2.setLabel([0,1,-1,2])
 b.setLabel([4,-4,5,-5,1,-1,3,-3])
 Z.setLabel([2,3,-3,7])
 W_trans.setLabel([6,0,5,-5])
 Ta2bar=((Ta2*b)*Z)*W_trans
 Ta2bar.permute([6,4,-4,7],3)
 ###########################
 Tb2.setLabel([0,1,-1,2])
 d.setLabel([4,-4,5,-5,1,-1,3,-3])
 W.setLabel([2,3,-3,7])
 Z_trans.setLabel([6,0,5,-5])
 Tb2bar=((Tb2*d)*Z_trans)*W
 Tb2bar.permute([6,4,-4,7],3)
 ###########################
 ###########################

 c3=norm_CTM(c3bar)
 c2=norm_CTM(c2bar)
 Ta2=norm_CTM(Ta2bar)
 Tb2=norm_CTM(Tb2bar)


 return c2, Ta2, Tb2, c3




def  add_left1(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]


 c1.setLabel([0,1])
 Tb1.setLabel([1,2,-2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,-2,3],3)

 c4.setLabel([0,1])
 Ta3.setLabel([1,2,-2,3])
 c4bar=c4*Ta3
 c4bar.permute([3,0,2,-2],1)

 c_f=copy.copy(c1bar)
 c_f_trans=copy.copy(c1bar)
 c_f_trans.transpose()
 c_f.setLabel([0,1,-1,2])
 c_f_trans.setLabel([2,3,4,-4])
 c_s=copy.copy(c4bar)
 c_s_trans=copy.copy(c4bar)
 c_s_trans.transpose()
 c_s.setLabel([2,3,4,-4])
 c_s_trans.setLabel([0,1,-1,2])

 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,-1,3,4,-4],3)

 Z,  V,  s = TruncateU.setTruncation(theta, chi_dim)

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

####################################

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,-2,3])
 Tb4.setLabel([4,5,-5,0])
 a.setLabel([5,-5,6,-6,7,-7,2,-2])
 Q=((c1*Tb1)*(Tb4))*a
 Q.permute([4,6,-6,3,7,-7],3)
 #Q.combineBond([3,7])

 c4.setLabel([0,1])
 Ta3.setLabel([1,4,-4,5])
 Ta4.setLabel([0,2,-2,3])
 c.setLabel([2,-2,4,-4,6,-6,7,-7])
 Q1=((c4*Ta3)*(Ta4))*c
 Q1.permute([5,6,-6,3,7,-7],3)
 #Q1.combineBond([5,6])


 c_f=copy.copy(Q)
 c_f_trans=copy.copy(Q)
 c_f_trans.transpose()
 c_f.setLabel([0,1,-1,2,3,-3])
 c_f_trans.setLabel([2,3,-3,4,5,-5])
 c_s=copy.copy(Q1)
 c_s_trans=copy.copy(Q1)
 c_s_trans.transpose()
 c_s.setLabel([2,3,-3,4,5,-5])
 c_s_trans.setLabel([0,1,-1,2,3,-3])
 theta=c_f*c_f_trans+c_s_trans*c_s
 theta.permute([0,1,-1,4,5,-5],3)

##################################################

 W,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 W_trans=copy.copy(W)
 W_trans.transpose()

 ##############################
 Z_trans.setLabel([4,0,2,-2])
 c1bar=c1bar*Z_trans
 c1bar.permute([4,3],1)
 ############################
 Z.setLabel([0,2,-2,4])
 c4bar=c4bar*Z
 c4bar.permute([4,3],0)
 #############################
 Tb4.setLabel([0,1,-1,2])
 a.setLabel([1,-1,3,-3,4,-4,5,-5])
 Z.setLabel([2,5,-5,6])
 W_trans.setLabel([7,0,3,-3])
 Tb4bar=((Tb4*a)*Z)*W_trans
 Tb4bar.permute([7,4,-4,6],1)
 ###########################
 Ta4.setLabel([0,1,-1,2])
 c.setLabel([1,-1,3,-3,4,-4,5,-5])
 W.setLabel([2,5,-5,6])
 Z_trans.setLabel([7,0,3,-3])
 Ta4bar=((Ta4*c)*Z_trans)*W
 Ta4bar.permute([7,4,-4,6],1)
#############################
 c1=norm_CTM(c1bar)
 c4=norm_CTM(c4bar)
 Ta4=norm_CTM(Ta4bar)
 Tb4=norm_CTM(Tb4bar)
###############################



 return c1, Ta4, Tb4, c4


 

 
def  magnetization_value1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):


 CTM_net = uni10.Network("Network/CTM.net")
 CTM_net.putTensor('c1',c1)
 CTM_net.putTensor('c2',c2)
 CTM_net.putTensor('c3',c3)
 CTM_net.putTensor('c4',c4)
 CTM_net.putTensor('Ta1',Ta1)
 CTM_net.putTensor('Ta2',Ta2)
 CTM_net.putTensor('Ta3',Ta3)
 CTM_net.putTensor('Ta4',Ta4)
 CTM_net.putTensor('Tb1',Tb1)
 CTM_net.putTensor('Tb2',Tb2)
 CTM_net.putTensor('Tb3',Tb3)
 CTM_net.putTensor('Tb4',Tb4)
 CTM_net.putTensor('a',a)
 CTM_net.putTensor('b',b)
 CTM_net.putTensor('c',c)
 CTM_net.putTensor('d',d)
 norm=CTM_net.launch()
 
# c1.setLabel([4,1])
# c2.setLabel([3,7])
# c3.setLabel([24,4])
# c4.setLabel([18,22])
# Ta1.setLabel([2,6,-6,3]) 
# Ta2.setLabel([14,10,-10,7]) 
# Ta3.setLabel([22,19,-19,23]) 
# Ta4.setLabel([18,15,-15,11]) 
# Tb1.setLabel([1,5,-5,2]) 
# Tb2.setLabel([4,17,-17,14]) 
# Tb3.setLabel([23,20,-20,24]) 
# Tb4.setLabel([11,8,-8,4])
# a.setLabel([8,-8,12,-12,9,-9,5,-5])
# b.setLabel([9,-9,13,-13,10,-10,6,-6])
# c.setLabel([15,-15,19,-19,16,-16,12,-12])
# d.setLabel([16,-16,20,-20,17,-17,13,-13])
# norm=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*Ta1)*b))*(((c1*Tb1)*Tb4)*a) 
 #print 'hi1', '\n',norm,
 return norm
 

def permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 a.permute([3,-3,2,-2,1,-1,0,10],4)
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 ##print'a', a
 b.setLabel([0,10,1,-1,2,-2,3,-3])
 b.permute([3,-3,2,-2,1,-1,0,10],4)
 b.setLabel([0,10,1,-1,2,-2,3,-3])

 c.setLabel([0,10,1,-1,2,-2,3,-3])
 c.permute([3,-3,2,-2,1,-1,0,10],4)
 c.setLabel([0,10,1,-1,2,-2,3,-3])

 d.setLabel([0,10,1,-1,2,-2,3,-3])
 d.permute([3,-3,2,-2,1,-1,0,10],4)
 d.setLabel([0,10,1,-1,2,-2,3,-3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([0,1],0)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([0,1],2)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,-1,2])
 Ta1.permute([2,1,-1,0],1)
 Ta1.setLabel([0,1,-1,2])
 
 Ta2.setLabel([0,1,-1,2])
 Ta2.permute([2,1,-1,0],1)
 Ta2.setLabel([0,1,-1,2])
 
 Ta3.setLabel([0,1,-1,2])
 Ta3.permute([2,1,-1,0],3)
 Ta3.setLabel([0,1,-1,2])

 Ta4.setLabel([0,1,-1,2])
 Ta4.permute([2,1,-1,0],3)
 Ta4.setLabel([0,1,-1,2])

 Tb1.setLabel([0,1,-1,2])
 Tb1.permute([2,1,-1,0],1)
 Tb1.setLabel([0,1,-1,2])
 
 Tb2.setLabel([0,1,-1,2])
 Tb2.permute([2,1,-1,0],1)
 Tb2.setLabel([0,1,-1,2])
 
 Tb3.setLabel([0,1,-1,2])
 Tb3.permute([2,1,-1,0],3)
 Tb3.setLabel([0,1,-1,2])

 Tb4.setLabel([0,1,-1,2])
 Tb4.permute([2,1,-1,0],3)
 Tb4.setLabel([0,1,-1,2])

 
def permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 a.permute([3,-3,2,-2,1,-1,0,10],4)
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 ##print'a', a
 b.setLabel([0,10,1,-1,2,-2,3,-3])
 b.permute([3,-3,2,-2,1,-1,0,10],4)
 b.setLabel([0,10,1,-1,2,-2,3,-3])

 c.setLabel([0,10,1,-1,2,-2,3,-3])
 c.permute([3,-3,2,-2,1,-1,0,10],4)
 c.setLabel([0,10,1,-1,2,-2,3,-3])

 d.setLabel([0,10,1,-1,2,-2,3,-3])
 d.permute([3,-3,2,-2,1,-1,0,10],4)
 d.setLabel([0,10,1,-1,2,-2,3,-3])

 
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


 Ta1.setLabel([0,1,-1,2])
 Ta1.permute([2,1,-1,0],3)
 Ta1.setLabel([0,1,-1,2])
 
 Ta2.setLabel([0,1,-1,2])
 Ta2.permute([2,1,-1,0],3)
 Ta2.setLabel([0,1,-1,2])
 
 Ta3.setLabel([0,1,-1,2])
 Ta3.permute([2,1,-1,0],1)
 Ta3.setLabel([0,1,-1,2])

 Ta4.setLabel([0,1,-1,2])
 Ta4.permute([2,1,-1,0],1)
 Ta4.setLabel([0,1,-1,2])

 Tb1.setLabel([0,1,-1,2])
 Tb1.permute([2,1,-1,0],3)
 Tb1.setLabel([0,1,-1,2])
 
 Tb2.setLabel([0,1,-1,2])
 Tb2.permute([2,1,-1,0],3)
 Tb2.setLabel([0,1,-1,2])
 
 Tb3.setLabel([0,1,-1,2])
 Tb3.permute([2,1,-1,0],1)
 Tb3.setLabel([0,1,-1,2])

 Tb4.setLabel([0,1,-1,2])
 Tb4.permute([2,1,-1,0],1)
 Tb4.setLabel([0,1,-1,2])
 



 
 
 
 
 
 
 
