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
 Contract=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*Ta1)*b))*(((c1*Tb1)*Tb4)*a)
 Contract.permute([-1,-2,-3,-4],2)
 #print Contract.trace()#, Contract.printDiagram()
 svd = Contract.getBlock()
 svd=svd.svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U1x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U1x")
 U1x.putBlock(svd[0].resize(svd[0].row(), chi))
 U1x_trans=copy.copy(U1x)
 U1x_trans.transpose()


 







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
 Contract=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c1*Tb1)*Tb4)*a))*(((c2*Ta2)*Ta1)*b)
 Contract.permute([-1,-2,-3,-4],2)
 #print Contract.trace()#, Contract.printDiagram()
 svd = Contract.getBlock().svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U2x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U2x")
 U2x.putBlock(svd[0].resize(svd[0].row(), chi))
 U2x_trans=copy.copy(U2x)
 U2x_trans.transpose()




 Ta4p=copy.copy(Ta4)
 Tb4p=copy.copy(Tb4)
 Ta2p=copy.copy(Ta2)
 Tb2p=copy.copy(Tb2)
 U1xb=copy.copy(U1x)
 U1xb_trans=copy.copy(U1x_trans)
 U2xb=copy.copy(U2x)
 U2xb_trans=copy.copy(U2x_trans)
 ap=copy.copy(a)
 bp=copy.copy(b)
 cp=copy.copy(c)
 dp=copy.copy(d)
 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([47,50])
 c4.setLabel([44,48])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta2p.setLabel([34,30,27]) 
 Ta3.setLabel([48,45,49]) 
 Ta4.setLabel([-1,21,17]) 
 Ta4p.setLabel([44,41,37]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([27,23,20])
 Tb2p.setLabel([47,43,40]) 
 Tb3.setLabel([49,46,50]) 
 Tb4.setLabel([11,8,4])
 Tb4p.setLabel([31,28,-3])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([21,-2,22,18])
 d.setLabel([22,26,23,19])
 ap.setLabel([28,32,29,-4])
 bp.setLabel([29,33,30,26])
 cp.setLabel([41,45,42,38])
 dp.setLabel([42,46,43,39])
 U1x.setLabel([17,18,15])
 U1x_trans.setLabel([15,11,12])
 U2x.setLabel([20,19,16])
 U2x_trans.setLabel([16,14,13])
 U1xb.setLabel([37,38,35])
 U1xb_trans.setLabel([35,31,32])
 U2xb.setLabel([40,39,36])
 U2xb_trans.setLabel([36,34,33])

 Contract=( ( (((((c1*Tb1 )*Tb4)*a) * U1x_trans) * ((((c2*Ta1)*Ta2)*b)*  U2x_trans))
 *  ( (((((c4*Ta3)*Ta4p)*cp)*U1xb) * ((((c3*Tb3)*Tb2p)*dp)*U2xb)) *
( ((Tb2*U2x)*d)   *   ((Ta2p*U2xb_trans)*bp) ) ) ) * ((Ta4 *U1x) *c) ) * ((Tb4p* ap)* U1xb_trans) 

 Contract.permute([-1,-2,-3,-4],2)
 #print Contract.trace()
 #print Contract.printDiagram()
 svd = Contract.getBlock().svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U3x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U3x")
 U3x.putBlock(svd[0].resize(svd[0].row(), chi))
 U3x_trans=copy.copy(U3x)
 U3x_trans.transpose()



 c1.setLabel([4,1])
 c2.setLabel([3,7])
 c3.setLabel([47,50])
 c4.setLabel([44,48])
 Ta1.setLabel([2,6,3]) 
 Ta2.setLabel([14,10,7]) 
 Ta2p.setLabel([34,30,-3]) 
 Ta3.setLabel([48,45,49]) 
 Ta4.setLabel([24,21,17]) 
 Ta4p.setLabel([44,41,37]) 
 Tb1.setLabel([1,5,2]) 
 Tb2.setLabel([-1,23,20])
 Tb2p.setLabel([47,43,40]) 
 Tb3.setLabel([49,46,50]) 
 Tb4.setLabel([11,8,4])
 Tb4p.setLabel([31,28,24])
 a.setLabel([8,12,9,5])
 b.setLabel([9,13,10,6])
 c.setLabel([21,25,22,18])
 d.setLabel([22,-2,23,19])
 ap.setLabel([28,32,29,25])
 bp.setLabel([29,33,30,-4])
 cp.setLabel([41,45,42,38])
 dp.setLabel([42,46,43,39])
 U1x.setLabel([17,18,15])
 U1x_trans.setLabel([15,11,12])
 U2x.setLabel([20,19,16])
 U2x_trans.setLabel([16,14,13])
 U1xb.setLabel([37,38,35])
 U1xb_trans.setLabel([35,31,32])
 U2xb.setLabel([40,39,36])
 U2xb_trans.setLabel([36,34,33])



 Contract=( ( (((((c1*Tb1 )*Tb4)*a) * U1x_trans) * ((((c2*Ta1)*Ta2)*b)*  U2x_trans))
 *  ( (((((c4*Ta3)*Ta4p)*cp)*U1xb) * ((((c3*Tb3)*Tb2p)*dp)*U2xb)) *
( ((Ta4 *U1x) *c)    *  ((Tb4p* ap)* U1xb_trans)  ) ) ) * ((Tb2*U2x)*d) ) *((Ta2p*U2xb_trans)*bp)

 Contract.permute([-1,-2,-3,-4],2)
 #print Contract.trace()
 #print Contract.printDiagram()
 svd = Contract.getBlock().svd()
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 U4x=uni10.UniTensor([Contract.bond()[0],Contract.bond()[1], bdo], "U4x")
 U4x.putBlock(svd[0].resize(svd[0].row(), chi))
 U4x_trans=copy.copy(U4x)
 U4x_trans.transpose()



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
 U3x_trans.setLabel([4,0,2])
 c1bar=c1bar*U3x_trans
 c1bar.permute([4,3],1)
 ############################
 U3x.setLabel([0,2,4])
 c4bar=c4bar*U3x
 c4bar.permute([4,3],1)
 #############################
 Tb4.setLabel([0,1,2])
 a.setLabel([1,3,4,5])
 U3x.setLabel([2,5,6])
 U1x_trans.setLabel([7,0,3])
 Tb4bar=((Tb4*a)*U3x)*U1x_trans
 Tb4bar.permute([7,4,6],2)
 ###########################
 Ta4.setLabel([0,1,2])
 c.setLabel([1,3,4,5])
 U1x.setLabel([2,5,6])
 U3x_trans.setLabel([7,0,3])
 Ta4bar=((Ta4*c)*U3x_trans)*U1x
 Ta4bar.permute([7,4,6],2)
#############################
 c1bar=c1bar*(1.00/c1bar.norm())
 c4bar=c4bar*(1.00/c4bar.norm())
 Ta4bar=Ta4bar*(1.00/Ta4bar.norm())
 Tb4bar=Tb4bar*(1.00/Tb4bar.norm())
###############################

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
 U4x_trans.setLabel([4,1,2])
 c2bar=c2bar*U4x_trans
 c2bar.permute([3,4],1)
 ############################
 U4x.setLabel([0,2,4])
 c3bar=c3bar*U4x
 c3bar.permute([4,3],1)
 #############################
 Ta2.setLabel([0,1,2])
 b.setLabel([4,5,1,3])
 U4x.setLabel([2,3,7])
 U2x_trans.setLabel([6,0,5])
 Ta2bar=((Ta2*b)*U4x)*U2x_trans
 Ta2bar.permute([6,4,7],2)
 ###########################
 Tb2.setLabel([0,1,2])
 d.setLabel([4,5,1,3])
 U2x.setLabel([2,3,7])
 U4x_trans.setLabel([6,0,5])
 Tb2bar=((Tb2*d)*U4x_trans)*U2x
 Tb2bar.permute([6,4,7],2)
 ###########################
 ###########################
 c3bar=c3bar*(1.00/c3bar.norm())
 c2bar=c2bar*(1.00/c2bar.norm())
 Ta2bar=Ta2bar*(1.00/Ta2bar.norm())
 Tb2bar=Tb2bar*(1.00/Tb2bar.norm())


 return c1bar, Ta4bar, Tb4bar, c4bar, c2bar, Ta2bar, Tb2bar, c3bar
 
 
 
 
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

 
