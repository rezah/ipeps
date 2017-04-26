import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
#import random
import copy
#import line_profiler
import TruncateU
import basicB
#import time


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
 Max_val=abs(MaxAbs(c))
 if ( Max_val < 0.50e-1) or (Max_val > 0.50e+1) and abs(Max_val)>1.0e-12:
  c*=(1.00/Max_val) 
 return c


def distance(theta,A):
   blk_qnums = theta.blockQnum()
   val=0
   for qnum in blk_qnums:
    T1=theta.getBlock(qnum)
    T2=A.getBlock(qnum)
    print "col", theta.getBlock(qnum).col(), qnum 
    for  i  in xrange( int(theta.getBlock(qnum).col()) ):
     if abs(T1[i]) > 1.00e-11:  
      val=val+abs((T1[i]-T2[i]) / T1[i])
      #if abs((T1[i]-T2[i]) / T1[i]) > 1.00e+1: print "hi", T1[i], T2[i], i 
     else: val=val+(T1[i]-T2[i]); #print T1[i]-T2[i]
   return val 

#@profile 
def produce_projectives(theta,theta1,chi_dim):
 theta=copy.copy(theta)
 theta1=copy.copy(theta1)
 
 theta.setLabel([1,2,20,3,4,40])
 theta.permute([1,2,20,3,4,40],3) 

 
 U, s, V=TruncateU.svd_parity1(theta)
 

 U.setLabel([1,2,20,-1,-2,-3])
 s.setLabel([-1,-2,-3,3,4,5])
 R=U*s
 R.permute([1,2,20,3,4,5],3)



 theta1.setLabel([1,2,20,3,4,40])
 theta1.permute([1,2,20,3,4,40],3) 

 
 U, s, V=TruncateU.svd_parity1(theta1)



 U.setLabel([1,2,20,-1,-2,-3])
 s.setLabel([-1,-2,-3,6,7,8])
 Rb=U*s
 Rb.permute([1,2,20,6,7,8],3)
 Rb.permute([6,7,8,1,2,20],3)
 
 
  
 A=R*Rb
 A.permute([6,7,8,3,4,5],3)

# print A
 V, U, s=TruncateU.setTruncation(A, chi_dim) 



 U.setLabel([-1,3,4,5])
 V.setLabel([6,7,8,-1])
 if MaxAbs(s) > 1.0e-8:
  s=s*(1.00/MaxAbs(s)) 
 s=TruncateU.inverse(s)
 s=TruncateU.Sqrt(s)
 
 s.setLabel([6,-1])
 U=s*U
 U.permute([6,3,4,5],1)

 s.setLabel([-1,9])
 V=V*s
 V.permute([6,7,8,9],3) 


 R.permute([1,2,20,3,4,5],3)
 U.transpose()
 U1x=R*U
 U1x.permute([1,2,20,6],3)



 Rb.permute([6,7,8,1,2,20],3)
 V.transpose()
 U1x_trans=Rb*V
 U1x_trans.permute([9,1,2,20],1)

 return U1x, U1x_trans 









#@profile 
def  add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D):

 
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]

# t0=time.time()

 CTM_1 = uni10.Network("Network/CTM1.net")
 CTM_1.putTensor('c1',c1)
 CTM_1.putTensor('c2',c2)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 theta=CTM_1.launch()
 theta.permute([100, 300, -300 , 400, 200 ,-200],3)

# print time.time() - t0, "up"



 
 CTM_2 = uni10.Network("Network/CTM2.net")
 CTM_2.putTensor('c3',c3)
 CTM_2.putTensor('c4',c4)
 CTM_2.putTensor('Ta3',Ta3)
 CTM_2.putTensor('Ta4',Ta4)
 CTM_2.putTensor('Tb2',Tb2)
 CTM_2.putTensor('Tb3',Tb3)
 CTM_2.putTensor('c',c)
 CTM_2.putTensor('d',d)
 theta1=CTM_2.launch()
 theta1.permute([100, 300, -300 , 400, 200 ,-200],3)


 
 U1x, U1x_trans=produce_projectives(theta,theta1, chi_dim)

 
 theta.permute([  400, 200 ,-200, 100, 300, -300], 3)
 theta1.permute([  400, 200 ,-200, 100, 300, -300], 3)

 U2x, U2x_trans=produce_projectives(theta,theta1, chi_dim)




 
# t0=time.time()

 Ta4p=copy.copy(Ta4)
 Ta4p.setName("Ta4p")
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

 CTM_3 = uni10.Network("Network/CTM3.net")
 CTM_3.putTensor('c1',c1)
 CTM_3.putTensor('c2',c2)
 CTM_3.putTensor('Ta1',Ta1)
 CTM_3.putTensor('Ta2',Ta2)
 CTM_3.putTensor('Ta4',Ta4)
 CTM_3.putTensor('Tb1',Tb1)
 CTM_3.putTensor('Tb2',Tb2)
 CTM_3.putTensor('Tb4',Tb4)
 CTM_3.putTensor('a',a)
 CTM_3.putTensor('b',b)
 CTM_3.putTensor('c',c)
 CTM_3.putTensor('d',d)
 CTM_3.putTensor('U1x',U1x)
 CTM_3.putTensor('U1x_trans',U1x_trans)
 CTM_3.putTensor('U2x',U2x)
 CTM_3.putTensor('U2x_trans',U2x_trans)
 theta=CTM_3.launch() 
 theta.permute([100, 300, -300 , 400, 200 ,-200],3)

# print time.time() - t0, "upup"

 

 
 CTM_4 = uni10.Network("Network/CTM4.net")
 CTM_4.putTensor('c3',c3)
 CTM_4.putTensor('c4',c4)
 CTM_4.putTensor('Ta2p',Ta2p)
 CTM_4.putTensor('Ta3',Ta3)
 CTM_4.putTensor('Ta4p',Ta4p)
 CTM_4.putTensor('Tb2p',Tb2p)
 CTM_4.putTensor('Tb3',Tb3)
 CTM_4.putTensor('Tb4p',Tb4p)
 CTM_4.putTensor('ap',ap)
 CTM_4.putTensor('bp',bp)
 CTM_4.putTensor('cp',cp)
 CTM_4.putTensor('dp',dp)
 CTM_4.putTensor('U1xb',U1xb)
 CTM_4.putTensor('U1xb_trans',U1xb_trans)
 CTM_4.putTensor('U2xb',U2xb)
 CTM_4.putTensor('U2xb_trans',U2xb_trans)
 theta1=CTM_4.launch() 
 theta1.permute([100, 300, -300 , 400, 200 ,-200],3)
 #print "hi", Ta4.printDiagram(), Tb4p.printDiagram()
 

 U3x, U3x_trans=produce_projectives(theta,theta1, chi_dim)

 theta.permute([  400, 200 ,-200, 100, 300, -300], 3)
 theta1.permute([  400, 200 ,-200, 100, 300, -300], 3)

 U4x, U4x_trans=produce_projectives(theta,theta1, chi_dim)


###############################################################################
 c1.setLabel([0,1])
 Tb1.setLabel([1,2,-2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,-2,3],3)

 c4.setLabel([0,1])
 Ta3.setLabel([1,2,-2,3])
 c4bar=c4*Ta3
 c4bar.permute([0,2,-2,3],1)

 ##############################
 U3x_trans.setLabel([4,0,2,-2])
 c1bar=c1bar*U3x_trans
 c1bar.permute([4,3],1)
 ############################
 U3x.setLabel([0,2,-2,4])
 c4bar=c4bar*U3x
 c4bar.permute([4,3],0)

 
 #############################
 Tb4.setLabel([0,1,-1,2])
 a.setLabel([1,-1,3,-3,4,-4,5,-5])
 U3x.setLabel([2,5,-5,6])
 U1x_trans.setLabel([7,0,3,-3])
 Tb4bar=((Tb4*U3x)*a)*U1x_trans
 Tb4bar.permute([7,4,-4,6],1)
# absorb_1 = uni10.Network("Network/absorb1.net")
# absorb_1.putTensor('Tb4',Tb4)
# absorb_1.putTensor('a',a)
# absorb_1.putTensor('U3x',U3x)
# absorb_1.putTensor('U1x_trans',U1x_trans)
# Tb4bar=absorb_1.launch() 

 ###########################
 Ta4.setLabel([0,1,-1,2])
 c.setLabel([1,-1,3,-3,4,-4,5,-5])
 U1x.setLabel([2,5,-5,6])
 U3x_trans.setLabel([7,0,3,-3])
 Ta4bar=((Ta4*U3x_trans)*c)*U1x
 Ta4bar.permute([7,4,-4,6],1)
#############################
###############################
 c1=norm_CTM(c1bar)
 c4=norm_CTM(c4bar)
 Ta4=norm_CTM(Ta4bar)
 Tb4=norm_CTM(Tb4bar)
 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,-2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,-2,3],3)
# c3.setLabel([0,1])
 c3.setLabel([1,0])
 Tb3.setLabel([3,2,-2,1])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2,-2],1)
 ##############################
 U4x_trans.setLabel([4,1,2,-2])
 c2bar=c2bar*U4x_trans
 c2bar.permute([3,4],2)
 ############################
 U4x.setLabel([0,2,-2,4])
 c3bar=c3bar*U4x
# c3bar.permute([4,3],1)
 c3bar.permute([3,4],1)
 #############################
 Ta2.setLabel([0,1,-1,2])
 b.setLabel([4,-4,5,-5,1,-1,3,-3])
 U4x.setLabel([2,3,-3,7])
 U2x_trans.setLabel([6,0,5,-5])
 Ta2bar=((Ta2*U4x)*b)*U2x_trans
 Ta2bar.permute([6,4,-4,7],3)
 ###########################
 Tb2.setLabel([0,1,-1,2])
 d.setLabel([4,-4,5,-5,1,-1,3,-3])
 U2x.setLabel([2,3,-3,7])
 U4x_trans.setLabel([6,0,5,-5])
 Tb2bar=((Tb2*U4x_trans)*d)*U2x
 Tb2bar.permute([6,4,-4,7],3)
 ###########################
 ###########################

 c3=norm_CTM(c3bar)
 c2=norm_CTM(c2bar)
 Ta2=norm_CTM(Ta2bar)
 Tb2=norm_CTM(Tb2bar)



 return c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3
 
 
 
 
 
#@profile 
 
def  magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):

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


def  Env_energy_h(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4, Tb4,D,d_phys)
 
 E_ab=basicB.Energy_ab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0, a_u, b_u)

 return E_ab

def  Env_energy_v(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,c_u,a_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E_ca=basicB.Energy_ca(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0,c_u,a_u)

 return E_ca





#@profile 
 
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
 



 
 
 
 
 
 
 
