import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time
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
 Ta3.setLabel([3,2,1])
 c3bar=c3*Ta3
 c3bar.permute([3,0,2],1)
 
 Mc2=c2bar.getBlock()
 Mc2_trans=copy.copy(Mc2)
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
 
 
