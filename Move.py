import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time

def  add_right(c2,Ta2,Tb2,c3,Ta1,Tb3,b,d,chi,D):
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,3],2)
 
 c3.setLabel([0,1])
 Tb3.setLabel([3,2,1])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2],1)
 
 Mc2=c2bar.getBlock()
 Mc2_trans=copy.copy(Mc2)
 Mc2_trans.transpose()

 Mc3=c3bar.getBlock()
 Mc3_trans=copy.copy(Mc3)
 Mc3_trans.transpose()
 Mfinal=Mc2*Mc2_trans+Mc3_trans*Mc3
 svd=Mfinal.svd()
 Z=uni10.UniTensor([c2bar.bond()[0],c2bar.bond()[1] ,bdo] )
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.transpose()

 c2.setLabel([0,1])
 Ta2.setLabel([3,2,1])
 Ta1.setLabel([4,5,0])
 b.setLabel([6,7,2,5])
 Q=((c2*Ta2)*(Ta1))*b
 Q.permute([3,7,4,6],2)
 
 c3.setLabel([0,1])
 Tb2.setLabel([0,2,3])
 Tb3.setLabel([5,4,1])
 d.setLabel([7,4,2,6])
 Q1=((c3*Tb2)*(Tb3))*d
 Q1.permute([5,7,3,6],2)
 
 
 Mc1=Q.getBlock()
 Mc1_trans=copy.copy(Mc1)
 Mc1_trans.transpose()

 Mc4=Q1.getBlock()
 Mc4_trans=copy.copy(Mc4)
 Mc4_trans.transpose()
 Mfinal=Mc1*Mc1_trans+Mc4_trans*Mc4

 svd=Mfinal.svd()
 W=uni10.UniTensor([Q.bond()[0],Q.bond()[1] ,bdo] )
 W.putBlock(svd[0].resize(svd[0].row(), chi))

 W_trans=copy.copy(W)
 W_trans.transpose()
 Sum=0

 ##############################
 Z_trans.setLabel([4,1,2])
 c2bar=c2bar*Z_trans
 c2bar.permute([3,4],1)
 ############################
 Z.setLabel([0,2,4])
 c3bar=c3bar*Z
 c3bar.permute([4,3],1)
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
 if ( (abs(c3bar.getBlock().absMax()) < 0.50e-1) or (abs(c3bar.getBlock().absMax()) > 0.50e+1)   ):
  c3=c3bar*(1.00/c3bar.getBlock().absMax()); 
 else: c3=c3bar;
 #print 'norm000', c3.norm(), c3.getBlock().absMax()
 
 if ( (abs(c2bar.getBlock().absMax()) < 0.50e-1) or (abs(c2bar.getBlock().absMax()) > 0.50e+1) ):
  c2=c2bar*(1.00/c2bar.getBlock().absMax()); 
 else: c2=c2bar;
 #print 'norm111', c2.norm(), c2.getBlock().absMax()
 
 if ( (abs(Ta2bar.getBlock().absMax()) < 0.50e-1) or (abs(Ta2bar.getBlock().absMax()) > 0.50e+1) ):
  Ta2=Ta2bar*(1.00/Ta2bar.getBlock().absMax()); 
 else: Ta2=Ta2bar;
 #print 'norm222', Ta2.norm(), Ta2.getBlock().absMax()

 if ( (abs(Tb2bar.getBlock().absMax()) < 0.50e-1) or (abs(Tb2bar.getBlock().absMax()) > 0.50e+1) ):
  Tb2=Tb2bar*(1.00/Tb2bar.getBlock().absMax()); 
 else: Tb2=Tb2bar;
 #print 'norm333', Tb2.norm(), Tb2.getBlock().absMax()
 
 return c2, Ta2, Tb2, c3



def  add_left(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D,Truncation):


 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,3],2)

 c4.setLabel([0,1])
 Ta3.setLabel([1,2,3])
 c4bar=c4*Ta3
 c4bar.permute([3,0,2],1)

 Mc1=c1bar.getBlock()
 Mc1_trans=copy.copy(Mc1)
 Mc1_trans.transpose()

 Mc4=c4bar.getBlock()
 Mc4_trans=copy.copy(Mc4)
 Mc4_trans.transpose()
 
 Mfinal=Mc1*Mc1_trans+Mc4_trans*Mc4
 svd=Mfinal.svd()

 Z=uni10.UniTensor([c1bar.bond()[0],c1bar.bond()[1] ,bdo])
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.transpose()
 Sum=0
 norm_trace=svd[1].trace()
 
 #s=svd[1]*(1.00/norm_trace)
 
 for i in xrange(svd[1].row()):
   if (i>=chi):
    Sum=svd[1][i]+Sum
 #print'truncation0=', Sum
 
 Truncation[0]=Sum
####################################
 
 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 Tb4.setLabel([4,5,0])
 a.setLabel([5,6,7,2])
 Q=((c1*Tb1)*(Tb4))*a
 Q.permute([4,6,3,7],2)
 
 c4.setLabel([0,1])
 Ta3.setLabel([1,4,5])
 Ta4.setLabel([0,2,3])
 c.setLabel([2,4,6,7])
 Q1=((c4*Ta3)*(Ta4))*c
 Q1.permute([5,6,3,7],2)
 
 
 Mc1=Q.getBlock()
 Mc1_trans=copy.copy(Mc1)
 Mc1_trans.transpose()

 Mc4=Q1.getBlock()
 Mc4_trans=copy.copy(Mc4)
 Mc4_trans.transpose()
 Mfinal=Mc1*Mc1_trans+Mc4_trans*Mc4

 svd=Mfinal.svd()
 W=uni10.UniTensor([Q.bond()[0],Q.bond()[1] ,bdo] )
 W.putBlock(svd[0].resize(svd[0].row(), chi))

 W_trans=copy.copy(W)
 W_trans.transpose()
 Sum=0
 for i in xrange(svd[1].row()):
   if (i>=chi):
    Sum=svd[1][i]+Sum
# print'truncation1=', Sum

 ##############################
 Z_trans.setLabel([4,0,2])
 c1bar=c1bar*Z_trans
 c1bar.permute([4,3],1)
 ############################
 Z.setLabel([0,2,4])
 c4bar=c4bar*Z
 c4bar.permute([4,3],1)
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
 if ( (abs(c1bar.getBlock().absMax()) < 0.50e-1) or (abs(c1bar.getBlock().absMax()) > 0.50e+1)   ):
  c1=c1bar*(1.00/c1bar.getBlock().absMax()); 
 else: c1=c1bar;
 #print 'norm0', c1.norm(), c1.getBlock().absMax()
 
 if ( (abs(c4bar.getBlock().absMax()) < 0.50e-1) or (abs(c4bar.getBlock().absMax()) > 0.50e+1) ):
  c4=c4bar*(1.00/c4bar.getBlock().absMax()); 
 else: c4=c4bar;
 #print 'norm1', c4.norm(), c4.getBlock().absMax()
 
 if ( (abs(Ta4bar.getBlock().absMax()) < 0.50e-1) or (abs(Ta4bar.getBlock().absMax()) > 0.50e+1) ):
  Ta4=Ta4bar*(1.00/Ta4bar.getBlock().absMax()); 
 else: Ta4=Ta4bar;
 #print 'norm2', Ta4.norm(), Ta4.getBlock().absMax()

 if ( (abs(Tb4bar.getBlock().absMax()) < 0.50e-1) or (abs(Tb4bar.getBlock().absMax()) > 0.50e+1) ):
  #print 'norm3', Tb4bar.norm(), Tb4bar.getBlock().absMax()
  Tb4=Tb4bar*(1.00/Tb4bar.getBlock().absMax()); 
 else: Tb4=Tb4bar;
 #print 'norm3', Tb4.norm(), Tb4.getBlock().absMax()
 
 
###############################
 return c1, Ta4, Tb4, c4
 

 
 
 
def  magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):
 c1.setLabel([0,1])
 c2.setLabel([5,6])
 c3.setLabel([7,8])
 c4.setLabel([11,10])

 Ta1.setLabel([21,20,5])
 Ta2.setLabel([22,19,6])
 Ta3.setLabel([10,12,13])
 Ta4.setLabel([11,17,18])

 Tb1.setLabel([1,2,21])
 Tb2.setLabel([7,15,22])
 Tb3.setLabel([13,14,8])
 Tb4.setLabel([18,3,0])


 a.setLabel([3,9,4,2])
 b.setLabel([4,23,19,20])
 d.setLabel([16,14,15,23])
 c.setLabel([17,12,16,9])
 
 norm=((((((c4*Ta4)*Ta3)*c)*(((c3*Tb3)*Tb2)*d)))*((((c2*Ta2)*Ta1)*b)*(((c1*Tb1)*Tb4)*a)))
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

# c2.setLabel([0,1])
# c2.permute([1,0],1)
# c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

# c4.setLabel([0,1])
# c4.permute([1,0],1)
# c4.setLabel([0,1])


 Ta1.setLabel([0,1,2])
 Ta1.permute([2,1,0],2)
 Ta1.setLabel([0,1,2])
 
 Ta2.setLabel([0,1,2])
 Ta2.permute([2,1,0],2)
 Ta2.setLabel([0,1,2])
 
 Ta3.setLabel([0,1,2])
 Ta3.permute([2,1,0],2)
 Ta3.setLabel([0,1,2])

 Ta4.setLabel([0,1,2])
 Ta4.permute([2,1,0],2)
 Ta4.setLabel([0,1,2])

 Tb1.setLabel([0,1,2])
 Tb1.permute([2,1,0],2)
 Tb1.setLabel([0,1,2])
 
 Tb2.setLabel([0,1,2])
 Tb2.permute([2,1,0],2)
 Tb2.setLabel([0,1,2])
 
 Tb3.setLabel([0,1,2])
 Tb3.permute([2,1,0],2)
 Tb3.setLabel([0,1,2])

 Tb4.setLabel([0,1,2])
 Tb4.permute([2,1,0],2)
 Tb4.setLabel([0,1,2])

def test_env_Ten(c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4):
 Rand_ten=copy.copy(c1)
 Rand_ten.randomize()
 if ( c1.getBlock().norm()  < 1.0e-5 ) or ( c1.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, c1', c1.getBlock().norm()
  c1=c1*(1.00/c1.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( c2.getBlock().norm()  < 1.0e-5 ) or ( c2.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, c2', c2.getBlock().norm()
  c2=c2*(1.00/c2.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( c3.getBlock().norm()  < 1.0e-5 ) or ( c3.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, c3', c3.getBlock().norm()
  c3=c3*(1.00/c3.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( c4.getBlock().norm()  < 1.0e-5 ) or ( c4.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, c4', c4.getBlock().norm()
  c4=c4*(1.00/c4.getBlock().norm())+(0.001)*Rand_ten(); 
 Rand_ten=copy.copy(Ta1)
 Rand_ten.randomize()
 if ( Ta1.getBlock().norm()  < 1.0e-5 ) or ( Ta1.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Ta1', Ta1.getBlock().norm()
  Ta1=Ta1*(1.00/Ta1.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( Ta2.getBlock().norm()  < 1.0e-5 ) or ( Ta2.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Ta2', Ta2.getBlock().norm()
  Ta2=Ta2*(1.00/Ta2.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( Ta3.getBlock().norm()  < 1.0e-5 ) or ( Ta3.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Ta3', Ta3.getBlock().norm()
  Ta3=Ta3*(1.00/Ta3.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( Ta4.getBlock().norm()  < 1.0e-5 ) or ( Ta4.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Ta4', Ta4.getBlock().norm()
  Ta4=Ta4*(1.00/Ta4.getBlock().norm())+(0.001)*Rand_ten(); 


 if ( Tb1.getBlock().norm()  < 1.0e-5 ) or ( Tb1.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Tb1', Tb1.getBlock().norm()
  Tb1=Tb1*(1.00/Tb1.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( Tb2.getBlock().norm()  < 1.0e-5 ) or ( Tb2.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Tb2', Tb2.getBlock().norm()
  Tb2=Tb2*(1.00/Tb2.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( Tb3.getBlock().norm()  < 1.0e-5 ) or ( Tb3.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Tb3', Tb3.getBlock().norm()
  Tb3=Tb3*(1.00/Tb3.getBlock().norm())+(0.001)*Rand_ten(); 
 if ( Tb4.getBlock().norm()  < 1.0e-5 ) or ( Tb4.getBlock().norm()  > 1.0e+5 )   :
  print 'Warning, norm, Tb4', Tb4.getBlock().norm()
  Tb4=Tb4*(1.00/Tb4.getBlock().norm())+(0.001)*Rand_ten(); 

 return c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4




def  distance(c1,c2,c3,c4,c1_f,c2_f,c3_f,c4_f):

 c1=c1*(1.00/c1.norm())
 c2=c2*(1.00/c2.norm())
 c3=c3*(1.00/c3.norm())
 c4=c4*(1.00/c4.norm())

 c1_f=c1_f*(1.00/c1_f.norm())
 c2_f=c2_f*(1.00/c2_f.norm())
 c3_f=c3_f*(1.00/c3_f.norm())
 c4_f=c4_f*(1.00/c4_f.norm())

 sum=0.0
 s1=c1.getBlock().svd()
 s1_f=c1_f.getBlock().svd()
 dis_val1=dis(s1[1],s1_f[1] )
 sum+=s1[1].trace()+s1_f[1].trace()

 #print dis_val1, s1[1].trace(),s1_f[1].trace()
 s1=c2.getBlock().svd()
 s1_f=c2_f.getBlock().svd()
 dis_val2=dis(s1[1],s1_f[1] )
 sum+=s1[1].trace()+s1_f[1].trace()

 #print dis_val2, s1[1].trace(),s1_f[1].trace()
 s1=c3.getBlock().svd()
 s1_f=c3_f.getBlock().svd()
 dis_val3=dis(s1[1],s1_f[1] )
 sum+=s1[1].trace()+s1_f[1].trace()

 #print dis_val3, s1[1].trace(),s1_f[1].trace()
 s1=c4.getBlock().svd()
 s1_f=c4_f.getBlock().svd()
 dis_val4=dis(s1[1],s1_f[1] )
 sum+=s1[1].trace()+s1_f[1].trace()
 #print dis_val4, s1[1].trace(),s1_f[1].trace()
 return (dis_val1+dis_val2+dis_val3+dis_val4) / (4.00)



def dis(s1,s2):
 sum=0
 for q in xrange(int(s1.row())):
  sum=sum+abs(s1[q]-s2[q])
 return sum




 
