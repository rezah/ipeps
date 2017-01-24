import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time


def add_up(c1,Tb1,Ta1,c2,Tb4,Ta2,a,b,chi,D):
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c1.setLabel([0,1])
 Tb4.setLabel([3,2,0])
 c1bar=c1*Tb4

def magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,a):
 c1.setLabel([0,1])
 c2.setLabel([5,6])
 c3.setLabel([7,8])
 c4.setLabel([11,10])
 Ta1.setLabel([1,2,5])
 Ta2.setLabel([7,4,6])
 Ta3.setLabel([10,9,8])
 Ta4.setLabel([11,3,0])
 a.setLabel([3,9,4,2])
 norm=((c1*Ta1*Ta4)*(a)*((c3*Ta3*Ta2))*(c2*c4))
 return norm 

def up(c1,c2,Ta1,Ta4,Ta2,chi,a):
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c1.setLabel([0,1])
 Ta4.setLabel([3,2,0])
 c1bar=c1*Ta4
 c1bar.permute([1,2,3],2)

 c2.setLabel([0,1])
 Ta2.setLabel([3,2,1])
 c2bar=c2*Ta2
 c2bar.permute([3,0,2],1)

 Mc1=c1bar.getBlock()
 Mc1_trans=copy.copy(Mc1)

 Mc1_trans.transpose()

 Mc2=c2bar.getBlock()
 Mc2_trans=copy.copy(Mc2)
 Mc2_trans.transpose()
 
 Mfinal=Mc1*Mc1_trans+Mc2_trans*Mc2
 svd=Mfinal.svd()
 Z=uni10.UniTensor([c1bar.bond()[0],c1bar.bond()[1] ,bdo] )
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.transpose()


 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 Tb4.setLabel([4,5,0])
 a.setLabel([5,6,7,2])
 Q=((c1*Tb1)*(Tb4))*a
 Q.permute([3,7,4,6],2)



 c2.setLabel([0,1])
 Ta2.setLabel([3,2,1])
 Ta1.setLabel([4,5,0])
 b.setLabel([6,7,2,5])
 Q1=((c2*Ta2)*(Ta1))*b
 Q1.permute([3,7,4,6],2)

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

 Mc1_trans.cTranspose()


 Mc2=c2bar.getBlock()
 Mc2_trans=copy.copy(Mc2)
 Mc2_trans.cTranspose()

 Mfinal=Mc1*Mc1_trans+Mc2_trans*Mc2
 #Mfinal=Mc1_trans*Mc1+Mc4_trans*Mc4
 #print Mfinal
 #D=Mfinal.eigh()
 svd=Mfinal.svd()
 #D_trans=copy.copy(D[1])
 #D_trans.transpose()
 #print 'Start',Mfinal, D_trans*D[0]*D[1],D[0],D[1], 'End'
 Z=uni10.UniTensor(uni10.CTYPE,[c1bar.bond()[0],c1bar.bond()[1] ,bdo] )
 #Z.putBlock(D[1].resize(D[1].row(), chi))
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.cTranspose()


 ##############################
 Z_trans.setLabel([4,1,2])
 c1bar=c1bar*Z_trans
 c1bar.permute([3,4],1)
 ############################
 Z.setLabel([0,2,4])
 c2bar=c2bar*Z
 c2bar.permute([4,3],1)
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
 if ( (abs(c1bar.getBlock().absMax()) < 0.50e-1) or (abs(c1bar.getBlock().absMax()) > 0.50e+1)   ):
  c1=c1bar*(1.00/c1bar.getBlock().absMax()); 
 else: c1=c1bar;
 #print 'norm0000', c1.norm(), c1.getBlock().absMax()
 
 if ( (abs(c2bar.getBlock().absMax()) < 0.50e-1) or (abs(c2bar.getBlock().absMax()) > 0.50e+1) ):
  c2=c2bar*(1.00/c2bar.getBlock().absMax()); 
 else: c2=c2bar;
 #print 'norm1111', c2.norm(), c2.getBlock().absMax()
 
 if ( (abs(Ta1bar.getBlock().absMax()) < 0.50e-1) or (abs(Ta1bar.getBlock().absMax()) > 0.50e+1) ):
  Ta1=Ta1bar*(1.00/Ta1bar.getBlock().absMax()); 
 else: Ta1=Ta1bar;
 #print 'norm2222', Ta1.norm(), Ta1.getBlock().absMax()

 if ( (abs(Tb1bar.getBlock().absMax()) < 0.50e-1) or (abs(Tb1bar.getBlock().absMax()) > 0.50e+1) ):
  Tb1=Tb1bar*(1.00/Tb1bar.getBlock().absMax()); 
 else: Tb1=Tb1bar;
 #print 'norm3333', Tb1.norm(), Tb1.getBlock().absMax()
 
 
 
 return  c1, Ta1, Tb1, c2
 

def  add_right(c2,Ta2,Tb2,c3,Ta1,Tb3,b,d,chi,D):
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 #################3################################

 Ta1.setLabel([0,1,2])
 a.setLabel([3,4,5,1])
 Z_trans.setLabel([7,2,5])
 Z.setLabel([0,3,6])
 Ta1bar=((Ta1*a)*Z)*Z_trans
 Ta1bar.permute([6,4,7],2)
#####################################
 c1=c1bar*(1.00/c1bar.norm())
 c2=c2bar*(1.00/c2bar.norm())
 Ta1=Ta1bar*(1.00/Ta1bar.norm())
 return c1, c2 , Ta1
 



def  right(c3,c2,Ta2,Ta1,Ta3,chi,a):
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,3],2)
 
 c3.setLabel([0,1])

 Tb3.setLabel([3,2,1])
 c3bar=c3*Tb3

 Ta3.setLabel([3,2,1])
 c3bar=c3*Ta3
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

 Mc2_trans.cTranspose()
 
 
 Mc3=c3bar.getBlock()
 Mc3_trans=copy.copy(Mc3)
 Mc3_trans.cTranspose()

 Mfinal=Mc2*Mc2_trans+Mc3_trans*Mc3
 svd=Mfinal.svd()
 Z=uni10.UniTensor(uni10.CTYPE,[c2bar.bond()[0],c2bar.bond()[1] ,bdo] )
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.cTranspose()

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

def  add_down(c4,Ta3,Tb3,c3,Ta4,Tb2,c,d,chi,D):
#################################################

 a.setLabel([4,5,1,3])
 Z.setLabel([2,3,6])
 Z_trans.setLabel([7,0,5])
 Ta2bar=((Ta2*a)*Z)*Z_trans
 Ta2bar.permute([7,4,6],2)
 #print Ta4bar
 ###########################
 c3=c3bar*(1.00/c3bar.norm())
 c2=c2bar*(1.00/c2bar.norm())
 Ta2=Ta2bar*(1.00/Ta2bar.norm())
 
 return  c2, c3 , Ta2


def  down(c4,c3,Ta3,Ta4,Ta2,chi,a):
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c4.setLabel([0,1])
 Ta4.setLabel([0,2,3])
 c4bar=c4*Ta4
 c4bar.permute([1,2,3],2)
 
 c3.setLabel([0,1])

 Tb2.setLabel([0,2,3])
 c3bar=c3*Tb2
 c3bar.permute([3,1,2],1)
 
 Mc4=c4bar.getBlock()
 Mc4_trans=copy.copy(Mc4)
 Mc4_trans.transpose()

 Mc3=c3bar.getBlock()
 Mc3_trans=copy.copy(Mc3)
 Mc3_trans.transpose()
 
 
 Mfinal=Mc4*Mc4_trans+Mc3_trans*Mc3
 svd=Mfinal.svd()
 Z=uni10.UniTensor([c4bar.bond()[0],c4bar.bond()[1] ,bdo] )
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.transpose()
 Sum=0
 for i in xrange(svd[1].row()):
   if (i>=chi):
    Sum=svd[1][i]+Sum
 #print 'truncation1=', Sum
#############################
 
 c4.setLabel([0,1])
 Ta4.setLabel([0,4,5])
 Ta3.setLabel([1,2,3])
 c.setLabel([4,2,6,7])
 Q=((c4*Ta4)*(Ta3))*c
 Q.permute([3,6,5,7],2)
 
 c3.setLabel([0,1])
 Tb2.setLabel([0,2,3])
 Tb3.setLabel([7,4,1])
 d.setLabel([5,4,2,6])
 Q1=((c3*Tb2)*(Tb3))*d
 Q1.permute([3,6,7,5],2)
 
 
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


 Ta2.setLabel([0,2,3])
 c3bar=c3*Ta2
 c3bar.permute([3,1,2],1)
 
 Mc4=c4bar.getBlock()
 Mc4_ctrans=copy.copy(Mc4)
 Mc4_ctrans.cTranspose()
 #print Mc4, Mc4_ctrans
 Mc3=c3bar.getBlock()
 Mc3_ctrans=copy.copy(Mc3)
 Mc3_ctrans.cTranspose()
 
 Mfinal=Mc4*Mc4_ctrans+Mc3_ctrans*Mc3
 svd=Mfinal.svd()
 #print svd[0] 
# bdi1 = uni10.Bond(uni10.BD_OUT, D)
# bdi2 = uni10.Bond(uni10.BD_OUT, chi)
# bdi1=c4bar.bond()[1]
# bdi2=c4bar.bond()[2]
# bdi2.change(uni10.BD_IN)
# bdi1.change(uni10.BD_IN)
 Z=uni10.UniTensor(uni10.CTYPE,[c4bar.bond()[0],c4bar.bond()[1] ,bdo] )
 Z.putBlock(svd[0].resize(svd[0].row(), chi))
 #print Z
 Z_trans=copy.copy(Z)
 Z_trans.cTranspose()

 #print Z_trans
 ##############################
 Z_trans.setLabel([4,1,2])
 c4bar=c4bar*Z_trans
 c4bar.permute([3,4],1)
 ############################
 Z.setLabel([1,2,4])
 c3bar=c3bar*Z
 c3bar.permute([3,4],1)
 #############################
 Ta3.setLabel([0,1,2])

 c.setLabel([3,1,4,5])
 Z.setLabel([0,3,7])
 W_trans.setLabel([6,2,4])
 Ta3bar=((Ta3*c)*Z)*W_trans
 Ta3bar.permute([7,5,6],2)
 ###########################
 Tb3.setLabel([0,1,2])
 d.setLabel([3,1,4,5])
 W.setLabel([0,3,6])
 Z_trans.setLabel([7,2,4])
 Tb3bar=((Tb3*d)*Z_trans)*W
 Tb3bar.permute([6,5,7],2)
 ###########################
 if ( (abs(c3bar.getBlock().absMax()) < 0.50e-1) or (abs(c3bar.getBlock().absMax()) > 0.50e+1)   ):
  c3=c3bar*(1.00/c3bar.getBlock().absMax()); 
 else: c3=c3bar;
 #print 'norm00', c3.norm(), c3.getBlock().absMax()
 
 if ( (abs(c4bar.getBlock().absMax()) < 0.50e-1) or (abs(c4bar.getBlock().absMax()) > 0.50e+1) ):
  c4=c4bar*(1.00/c4bar.getBlock().absMax()); 
 else: c4=c4bar;
 #print 'norm11', c4.norm(), c4.getBlock().absMax()
 
 if ( (abs(Ta3bar.getBlock().absMax()) < 0.50e-1) or (abs(Ta3bar.getBlock().absMax()) > 0.50e+1) ):
  Ta3=Ta3bar*(1.00/Ta3bar.getBlock().absMax()); 
 else: Ta3=Ta3bar;
 #print 'norm22', Ta3.norm(), Ta3.getBlock().absMax()

 if ( (abs(Tb3bar.getBlock().absMax()) < 0.50e-1) or (abs(Tb3bar.getBlock().absMax()) > 0.50e+1) ):
  Tb3=Tb3bar*(1.00/Tb3bar.getBlock().absMax()); 
 else: Tb3=Tb3bar;
 #print 'norm33', Tb3.norm(), Tb3.getBlock().absMax()
 

 return c4, Ta3, Tb3, c3


def  add_left(c1,Tb4,Ta4,c4,Tb1,Ta3,a,c,chi,D,Truncation):


 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 c1bar=c1*Tb1

 a.setLabel([5,1,3,4])
 Z.setLabel([0,5,6])
 Z_trans.setLabel([7,2,3])
 Ta3bar=((Ta3*a)*Z)*Z_trans
 Ta3bar.permute([6,4,7],2)
 ###########################
 c3=c3bar*(1.00/c3bar.norm())
 c4=c4bar*(1.00/c4bar.norm())
 Ta3=Ta3bar*(1.00/Ta3bar.norm())
 #################3################################
 
 return c3, c4 , Ta3

 
def  left(c1,c4,Ta4,Ta1,Ta3,chi,a):
 bdo = uni10.Bond(uni10.BD_OUT, chi)

 c1.setLabel([0,1])
 Ta1.setLabel([1,2,3])
 c1bar=c1*Ta1
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

 
 
 Mc1_trans=copy.copy(Mc1)
 Mc1_trans.cTranspose()
 
 Mc4=c4bar.getBlock()
 Mc4_trans=copy.copy(Mc4)
 Mc4_trans.cTranspose()
 
 Mfinal=Mc1*Mc1_trans+Mc4_trans*Mc4

 svd=Mfinal.svd()

 Z=uni10.UniTensor(uni10.CTYPE,[c1bar.bond()[0],c1bar.bond()[1] ,bdo] )
 Z.putBlock(svd[0].resize(svd[0].row(), chi))

 Z_trans=copy.copy(Z)
 Z_trans.cTranspose()

 Sum=0
 for i in xrange(svd[1].row()):
   if (i>=chi):
    Sum=svd[1][i]+Sum

# print'truncation1=', Sum


 print'truncation1=', Sum
 
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
 
 

 Ta4.setLabel([0,1,2])
 a.setLabel([1,3,4,5])
 Z.setLabel([2,5,6])
 Z_trans.setLabel([7,0,3])
 Ta4bar=((Ta4*a)*Z)*Z_trans
 Ta4bar.permute([7,4,6],2)
 #print Ta4bar
 ###########################
 c1=c1bar*(1.00/c1bar.norm())
 c4=c4bar*(1.00/c4bar.norm())
 Ta4=Ta4bar*(1.00/Ta4bar.norm())

 return c1, c4, Ta4
 
 
