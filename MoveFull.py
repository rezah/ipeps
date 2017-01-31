import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time

def  add_left(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D):

 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([4,24])
 c4.setLabel([18,22])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta3.setLabel([22,19,23]) 
 Ta4.setLabel([18,15,-3]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([4,17,14]) 
 Tb3.setLabel([23,20,24]) 
 Tb4.setLabel([-1,8,4])
 a.setLabel([8,-2,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([15,19,16,-4])
 d.setLabel([16,20,17,13])
 Contract=((((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*(Ta1))*b))*(((c1*Tb1)*(Tb4))*a))
 Contract.permute([-1,-2,-3,-4],2)
 #print (Contract.trace())
 #print Contract.printDiagram(),Contract1.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U1x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U1x")
 U1x.putBlock(svd[0].resize(svd[0].row(), chi))
 U1x_trans=copy.copy(U1x)
 U1x_trans.transpose()
##############################################################################################


 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([4,24])
 c4.setLabel([18,22])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([-1,10,7]) 
 Ta3.setLabel([22,19,23]) 
 Ta4.setLabel([18,15,11]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([4,17,-3]) 
 Tb3.setLabel([23,20,24]) 
 Tb4.setLabel([11,8,4])
 a.setLabel([8,12,9,5])
 b.setLabel([9,-2,10,6])
 c.setLabel([15,19,16,12])
 d.setLabel([16,20,17,-4])
 Contract=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c1*Tb1)*(Tb4))*a))*(((c2*Ta2)*(Ta1))*b)
 Contract.permute([-1,-2,-3,-4],2)

 #print (Contract.trace())
 #print Contract.printDiagram(),Contract1.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U1=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U1")
 U1.putBlock(svd[0].resize(svd[0].row(), chi))
 U1_trans=copy.copy(U1)
 U1_trans.transpose()

##########################################################################################

 c1.setLabel([-1,1])
 c2.setLabel([3,7])
 c3.setLabel([4,24])
 c4.setLabel([18,22])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta3.setLabel([22,19,23]) 
 Ta4.setLabel([18,15,11]) 
 Tb1.setLabel([1,-2,2]) 
 Tb2.setLabel([4,17,14]) 
 Tb3.setLabel([23,20,24]) 
 Tb4.setLabel([11,8,-3])
 a.setLabel([8,12,9,-4])
 b.setLabel([9,13,10,6])
 c.setLabel([15,19,16,12])
 d.setLabel([16,20,17,13])
 Contract=((((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*(Ta1))*b))* (Tb4*a))*(c1*Tb1)
 Contract.permute([-1,-2,-3,-4],2)
 Q=copy.copy(Contract)
 #print (Contract.trace())
 #print Contract.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U2x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U2x")
 U2x.putBlock(svd[0].resize(svd[0].row(), chi))
 U2x_trans=copy.copy(U2x)
 U2x_trans.transpose()
##############################################################################################

 c1.setLabel([4,1])
 c2.setLabel([3,-1])
 c3.setLabel([4,24])
 c4.setLabel([18,22])
 Ta1.setLabel([2,-2,3]) 
 Ta2.setLabel([14,10,-3]) 
 Ta3.setLabel([22,19,23]) 
 Ta4.setLabel([18,15,11]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([4,17,14]) 
 Tb3.setLabel([23,20,24]) 
 Tb4.setLabel([11,8,4])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,-4])
 c.setLabel([15,19,16,12])
 d.setLabel([16,20,17,13])
 Contract=((((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c1*Tb1)*(Tb4))*a))*(Ta2*b))*(c2*Ta1)
 Contract.permute([-1,-2,-3,-4],2)
 Q1=copy.copy(Contract)

 #print (Contract.trace())
 #print Contract.printDiagram(),Contract1.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U2=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U2")
 U2.putBlock(svd[0].resize(svd[0].row(), chi))
 U2_trans=copy.copy(U2)
 U2_trans.transpose()

##########################################################################################


 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([4,24])
 c4.setLabel([-3,22])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta3.setLabel([22,-4,23]) 
 Ta4.setLabel([-1,15,11]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([4,17,14]) 
 Tb3.setLabel([23,20,24]) 
 Tb4.setLabel([11,8,4])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([15,-2,16,12])
 d.setLabel([16,20,17,13])
 Contract=((((((c1*Tb1)*(Tb4))*a)  * (((c2*Ta2)*(Ta1))*b)) * (((c3*Tb2)*Tb3)*(d))) * (c*Ta4))*(c4*Ta3)
 Contract.permute([-1,-2,-3,-4],2)
 Q_b=copy.copy(Contract)

 #print (Contract.trace())
 #print Contract.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U3x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U3x")
 U3x.putBlock(svd[0].resize(svd[0].row(), chi))
 U3x_trans=copy.copy(U3x)
 U3x_trans.transpose()

###########################################################################################



 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([-3,24])
 c4.setLabel([18,22])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta3.setLabel([22,19,23]) 
 Ta4.setLabel([18,15,11]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([-1,17,14]) 
 Tb3.setLabel([23,-4,24]) 
 Tb4.setLabel([11,8,4])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([15,19,16,12])
 d.setLabel([16,-2,17,13])
 Contract=((((((c1*Tb1)*(Tb4))*a)  * (((c2*Ta2)*(Ta1))*b)) * (((c4*Ta3)*Ta4)*c)) * (d*Tb2))*(c3*Tb3)
 Contract.permute([-1,-2,-3,-4],2)
 Q1_b=copy.copy(Contract)

 #print (Contract.trace())
 #print Contract.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U3=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U3")
 U3.putBlock(svd[0].resize(svd[0].row(), chi))
 U3_trans=copy.copy(U3)
 U3_trans.transpose()
####################################################################

 M_1=Q.getBlock()
 M_2=Q_b.getBlock()

 M_final=M_1+M_2
 svd = M_final
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Uup=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "Uup")
 Uup.putBlock(svd[0].resize(svd[0].row(), chi))
 Uup_trans=copy.copy(Uup)
 Uup_trans.transpose()



 M_1=Q1.getBlock()
 M_2=Q1_b.getBlock()

 M_final=M_1+M_2
 svd = M_final
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Udo=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "Udo")
 Udo.putBlock(svd[0].resize(svd[0].row(), chi))
 Udo_trans=copy.copy(Udo)
 Udo_trans.transpose()


 U2=copy.copy(Udo)
 U2_trans=copy.copy(Udo_trans)
 U3=copy.copy(Udo)
 U3_trans=copy.copy(Udo_trans)


 U2x=copy.copy(Uup)
 U2x_trans=copy.copy(Uup_trans)
 U3x=copy.copy(Uup)
 U3x_trans=copy.copy(Uup_trans)


###############################################################################
 c1.setLabel([0,1])
 Tb1.setLabel([1,2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,3],2)

 c4.setLabel([0,1])
 Ta3.setLabel([1,2,3])
 c4bar=c4*Ta3
 c4bar.permute([0,2,3],1)

 ##############################
 U2x_trans.setLabel([4,0,2])
 c1bar=c1bar*U2x_trans
 c1bar.permute([4,3],1)
 ############################
 U3x.setLabel([0,2,4])
 c4bar=c4bar*U3x
 c4bar.permute([4,3],1)
 #############################
 Tb4.setLabel([0,1,2])
 a.setLabel([1,3,4,5])
 U2x.setLabel([2,5,6])
 U1x_trans.setLabel([7,0,3])
 Tb4bar=((Tb4*a)*U2x)*U1x_trans
 Tb4bar.permute([7,4,6],2)
 ###########################
 Ta4.setLabel([0,1,2])
 c.setLabel([1,3,4,5])
 U1x.setLabel([2,5,6])
 U3x_trans.setLabel([7,0,3])
 Ta4bar=((Ta4*c)*U3x_trans)*U1x
 Ta4bar.permute([7,4,6],2)
#############################


 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,3],2)
 
 c3.setLabel([0,1])
 Tb3.setLabel([3,2,1])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2],1)

 ##############################
 U2_trans.setLabel([4,1,2])
 c2bar=c2bar*U2_trans
 c2bar.permute([3,4],1)
 ############################
 U3.setLabel([0,2,4])
 c3bar=c3bar*U3
 c3bar.permute([4,3],1)
 #############################
 Ta2.setLabel([0,1,2])
 b.setLabel([4,5,1,3])
 U2.setLabel([2,3,7])
 U1_trans.setLabel([6,0,5])
 Ta2bar=((Ta2*b)*U2)*U1_trans
 Ta2bar.permute([6,4,7],2)
 ###########################
 Tb2.setLabel([0,1,2])
 d.setLabel([4,5,1,3])
 U1.setLabel([2,3,7])
 U3_trans.setLabel([6,0,5])
 Tb2bar=((Tb2*d)*U3_trans)*U1
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


 return c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3
 
 
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
 ##print 'hi1', '\n',norm,
 return norm
 
 
def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 invL2 = uni10.Matrix(Landa2.bond()[0].dim(), Landa2.bond()[1].dim())
 D=Landa2.bond()[0].dim()
 for i in xrange(Landa2.bond()[0].dim()):
   for j in xrange(Landa2.bond()[0].dim()):
    invL2[i*D+j] = 0 if ((Landa2[i*D+j].real) < 1.0e-10) else (1.00 / (Landa2[i*D+j].real))
 invLanda2.putBlock(invL2)
 return invLanda2
def sqt(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 invL2 = uni10.Matrix(Landa2.bond()[0].dim(), Landa2.bond()[1].dim())
 D=Landa2.bond()[0].dim()
 for i in xrange(Landa2.bond()[0].dim()):
   for j in xrange(Landa2.bond()[0].dim()):
    invL2[i*D+j] = ((Landa2[i*D+j])**(1.00/2.00))
 invLanda2.putBlock(invL2)
 return invLanda2


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

