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
 if ( Max_val < 0.50e-1) or (Max_val > 0.50e+1):
  c*=(1.00/Max_val) 
 return c
 
 
def check_eigenvalues(s):
 M_s=s.getBlock()
 p=0
 for i in xrange(int(M_s.row())):
    for j in xrange(int(M_s.col())):
     if i==j:
      if M_s[i*int(M_s.col())+j]<1.e-8:
       p=i
       #print "hiiiiiiiii", p, M_s[i*int(M_s.col())+j]
       break
    else:
        continue  # executed if the loop ended normally (no break)
    break  # executed if 'continue' was skipped (break)

 return p
 
 
 
def pick_vec(p,U):
 x=int(U.row())
 y=int(U.col())
 #print x, y, p
 Vec_up=uni10.Matrix(x,1)
 for i in xrange(x):
      Vec_up[i]=U[i*y+p]
 #print Vec_up
 return Vec_up


def add_vec_to_mat(U,vec):
 x=int(U.row())
 y=int(U.col())
 U_new=copy.copy(U)
 U_new.resize(x, y+1)
 p=y
 for i in xrange(x):
      U_new[i*(y+1)+p]=vec[i*(1)+0]
 return U_new



def Add_onevector(V_mat, U_mat, p, Vp, Up):
 Up_mat=Up.getBlock()
 Vp_mat=Vp.getBlock()
 Vp_mat.transpose()

 Vec_up=pick_vec(p,Up_mat)
 Vec_vp=pick_vec(p,Vp_mat)

# Vec_up.randomize()
# Vec_vp.randomize()
 #print "hi", U_mat.col(), Vec_vp.row() 
 Vec_vp=U_mat*Vec_vp
 V_mat.transpose()
 Vec_up=V_mat*Vec_up
 U_mat=add_vec_to_mat(U_mat,Vec_up)
 V_mat=add_vec_to_mat(V_mat,Vec_vp)
 V_mat.transpose()

 A=V_mat*U_mat 
 #print A.col(), A.row()  
 return A, V_mat,  U_mat

#@profile 

def produce_projectives(theta,chi_dim):
 theta=copy.copy(theta)
 
 theta.setLabel([1,2,20,3,4,40])
 theta.permute([1,2,20,3,4,40],3) 
 #theta.transpose()
 
 U, V, s=TruncateU.setTruncation(theta,chi_dim)

 if MaxAbs(s) > 1.0e+3 or MaxAbs(s) < 1.0e-1:
  s=s*(1.00/MaxAbs(s))
# s=TruncateU.inverse(s)
 s=TruncateU.Sqrt(s)

# M_s=s.getBlock()
# for i in xrange(int(M_s.row())):
#  for j in xrange(int(M_s.col())):
#   if i==j:
#    print "1", M_s[i*int(M_s.col())+j]
  
 U.setLabel([1,2,20,0])
 V.setLabel([0,1,2,20])
 s.setLabel([0,-1])
 U=U*s
 s.setLabel([-2,0])
 V=s*V
 U.setLabel([1,2,20,-1])
 V.setLabel([-2,1,2,20])
 

 V.transpose()
 A=V*U
 A.permute([-2,-1],1)

 Up, Vp, sp=TruncateU.setTruncation3(A, chi_dim) 

# p=check_eigenvalues(sp)
# #print "p", p
# p=14
# A_mat, V_mat,  U_mat=Add_onevector(V.getBlock(),U.getBlock(),p,Up, Vp)
# 
# print "\n"
# M_s=sp.getBlock()
# for i in xrange(int(M_s.row())):
#  for j in xrange(int(M_s.col())):
#   if i==j:
#    print "0", M_s[i*int(M_s.col())+j]

# x=int(A.bond(0).dim())+1
# y=int(A.bond(1).dim())+1

# bdi=uni10.Bond(uni10.BD_IN,x)
# bdo=uni10.Bond(uni10.BD_OUT,y)
# A_new=uni10.UniTensor([bdi,bdo])

# A_new.putBlock(A_mat)

# U=uni10.UniTensor([U.bond(0),U.bond(1), U.bond(2),bdo])
# V=uni10.UniTensor([bdi,V.bond(1),V.bond(2), V.bond(3)])
# #print Rb.printDiagram(), Rb_mat.col(), Rb_mat.row()
# U.putBlock(U_mat)
# V.putBlock(V_mat)
# U.setLabel([1,2,20,-1])
# V.setLabel([-2,1,2,20])

# Up, Vp, sp=TruncateU.setTruncation3(A_new, chi_dim+1) 
## print  A_new.printDiagram() 
## print "After1", "\n"
## M_s=sp.getBlock()
## for i in xrange(int(M_s.row())):
##  for j in xrange(int(M_s.col())):
##   if i==j:
##    print i, M_s[i*int(M_s.col())+j]


# p=15
# A_mat, V_mat,  U_mat=Add_onevector(V.getBlock(),U.getBlock(),p,Up, Vp)
# 
## print "After2"
## M_s=sp.getBlock()
## for i in xrange(int(M_s.row())):
##  for j in xrange(int(M_s.col())):
##   if i==j:
##    print "1", M_s[i*int(M_s.col())+j]

# x=int(A.bond(0).dim())+2
# y=int(A.bond(1).dim())+2

# bdi=uni10.Bond(uni10.BD_IN,x)
# bdo=uni10.Bond(uni10.BD_OUT,y)
# A_new=uni10.UniTensor([bdi,bdo])

# A_new.putBlock(A_mat)

# U=uni10.UniTensor([U.bond(0),U.bond(1), U.bond(2),bdo])
# V=uni10.UniTensor([bdi,V.bond(1),V.bond(2), V.bond(3)])
# #print Rb.printDiagram(), Rb_mat.col(), Rb_mat.row()
# U.putBlock(U_mat)
# V.putBlock(V_mat)
# U.setLabel([1,2,20,-1])
# V.setLabel([-2,1,2,20])

# Up, Vp, sp=TruncateU.setTruncation3(A_new, chi_dim+2) 
## print "AfterFinal", "\n"
## M_s=sp.getBlock()
## for i in xrange(int(M_s.row())):
##  for j in xrange(int(M_s.col())):
##   if i==j:
##    print i, M_s[i*int(M_s.col())+j]


 if MaxAbs(sp) > 1.0e+3 or MaxAbs(sp) < 1.0e-1:
  sp=sp*(1.00/MaxAbs(s))

 sp=TruncateU.inverse(sp)
 sp=TruncateU.Sqrt(sp)

 Vp.transpose()
 Vp.setLabel([-1,0])
 sp.setLabel([0,3])
 #print Vp.printDiagram(), sp.printDiagram() 
 U1x=U*Vp*sp
 U1x.permute([1,2,20,3],3)

 Up.transpose()
 Up.setLabel([0,-2])
 sp.setLabel([3,0])
 U1x_trans=sp*Up*V
 U1x_trans.permute([3,1,2,20],1)



 #U1x.setLabel([1,2,20,-3])
 #A=U1x_trans*U1x
 #print A[0]


 return U1x, U1x_trans  

# U1x_trans.transpose()
# U1x.transpose()
# return U1x_trans, U1x 
 
 
 
 
 

#@profile  
def  add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D):
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]

 
 #t0=time.time()

 CTM_1 = uni10.Network("Network/CTM_1.net")
 CTM_1.putTensor('c1',c1)
 CTM_1.putTensor('c2',c2)
 CTM_1.putTensor('c3',c3)
 CTM_1.putTensor('c4',c4)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 theta=CTM_1.launch()
 theta.permute([100,300,-300,200,400,-400],3)
 
 #print time.time() - t0, "up"


 #print theta.trace()#, Contract.printDiagram()
 #print theta.printDiagram(),MaxAbs(theta)
 #print theta.printDiagram() 

 #t0=time.time()

 U1x,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 U1x_trans=copy.copy(U1x)
 U1x_trans.transpose()



 #U1x, U1x_trans   = produce_projectives(theta, chi_dim)
 #print time.time() - t0, "svd"





# print "U1x_trans", U1x_trans.printDiagram()
# print "U1x", U1x.printDiagram()
 
# t0=time.time()

 CTM_2 = uni10.Network("Network/CTM_2.net")
 CTM_2.putTensor('c1',c1)
 CTM_2.putTensor('c2',c2)
 CTM_2.putTensor('c3',c3)
 CTM_2.putTensor('c4',c4)
 CTM_2.putTensor('Ta1',Ta1)
 CTM_2.putTensor('Ta2',Ta2)
 CTM_2.putTensor('Ta3',Ta3)
 CTM_2.putTensor('Ta4',Ta4)
 CTM_2.putTensor('Tb1',Tb1)
 CTM_2.putTensor('Tb2',Tb2)
 CTM_2.putTensor('Tb3',Tb3)
 CTM_2.putTensor('Tb4',Tb4)
 CTM_2.putTensor('a',a)
 CTM_2.putTensor('b',b)
 CTM_2.putTensor('c',c)
 CTM_2.putTensor('d',d)
 theta=CTM_2.launch()

 theta.permute([100,300,-300,200,400,-400],3)
 #print theta.trace()#, Contract.printDiagram()
# print time.time() - t0, "du"

 
 U2x,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 U2x_trans=copy.copy(U2x)
 U2x_trans.transpose()


# U2x, U2x_trans   = produce_projectives(theta, chi_dim)

# print "U2x_trans", U2x_trans.printDiagram()
# print "U2x", U2x.printDiagram()

 
 #t0=time.time()

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

 CTM_3 = uni10.Network("Network/CTM_3.net")
 CTM_3.putTensor('c1',c1)
 CTM_3.putTensor('c2',c2)
 CTM_3.putTensor('c3',c3)
 CTM_3.putTensor('c4',c4)
 CTM_3.putTensor('Ta1',Ta1)
 CTM_3.putTensor('Ta2',Ta2)
 CTM_3.putTensor('Ta2p',Ta2p)
 CTM_3.putTensor('Ta3',Ta3)
 CTM_3.putTensor('Ta4',Ta4)
 CTM_3.putTensor('Ta4p',Ta4p)
 CTM_3.putTensor('Tb1',Tb1)
 CTM_3.putTensor('Tb2',Tb2)
 CTM_3.putTensor('Tb2p',Tb2p)
 CTM_3.putTensor('Tb3',Tb3)
 CTM_3.putTensor('Tb4',Tb4)
 CTM_3.putTensor('Tb4p',Tb4p)
 CTM_3.putTensor('a',a)
 CTM_3.putTensor('b',b)
 CTM_3.putTensor('c',c)
 CTM_3.putTensor('d',d)
 CTM_3.putTensor('ap',ap)
 CTM_3.putTensor('bp',bp)
 CTM_3.putTensor('cp',cp)
 CTM_3.putTensor('dp',dp)
 CTM_3.putTensor('U1x',U1x)
 CTM_3.putTensor('U1x_trans',U1x_trans)
 CTM_3.putTensor('U1xb',U1xb)
 CTM_3.putTensor('U1xb_trans',U1xb_trans)
 CTM_3.putTensor('U2x',U2x)
 CTM_3.putTensor('U2x_trans',U2x_trans)
 CTM_3.putTensor('U2xb',U2xb)
 CTM_3.putTensor('U2xb_trans',U2xb_trans)
 theta=CTM_3.launch() 
 #print CTM_3, CTM_3.profile()
 #print theta.printDiagram()
 theta.permute([100,300,-300,200,400,-400],3)
 #print time.time() - t0, "upup"
 

 #print theta.trace()#, Contract.printDiagram()
 U3x,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 U3x_trans=copy.copy(U3x)
 U3x_trans.transpose()

# U3x, U3x_trans   = produce_projectives(theta, chi_dim)
 

# t0=time.time()

 CTM_4 = uni10.Network("Network/CTM_4.net")
 CTM_4.putTensor('c1',c1)
 CTM_4.putTensor('c2',c2)
 CTM_4.putTensor('c3',c3)
 CTM_4.putTensor('c4',c4)
 CTM_4.putTensor('Ta1',Ta1)
 CTM_4.putTensor('Ta2',Ta2)
 CTM_4.putTensor('Ta2p',Ta2p)
 CTM_4.putTensor('Ta3',Ta3)
 CTM_4.putTensor('Ta4',Ta4)
 CTM_4.putTensor('Ta4p',Ta4p)
 CTM_4.putTensor('Tb1',Tb1)
 CTM_4.putTensor('Tb2',Tb2)
 CTM_4.putTensor('Tb2p',Tb2p)
 CTM_4.putTensor('Tb3',Tb3)
 CTM_4.putTensor('Tb4',Tb4)
 CTM_4.putTensor('Tb4p',Tb4p)
 CTM_4.putTensor('a',a)
 CTM_4.putTensor('b',b)
 CTM_4.putTensor('c',c)
 CTM_4.putTensor('d',d)
 CTM_4.putTensor('ap',ap)
 CTM_4.putTensor('bp',bp)
 CTM_4.putTensor('cp',cp)
 CTM_4.putTensor('dp',dp)
 CTM_4.putTensor('U1x',U1x)
 CTM_4.putTensor('U1x_trans',U1x_trans)
 CTM_4.putTensor('U1xb',U1xb)
 CTM_4.putTensor('U1xb_trans',U1xb_trans)
 CTM_4.putTensor('U2x',U2x)
 CTM_4.putTensor('U2x_trans',U2x_trans)
 CTM_4.putTensor('U2xb',U2xb)
 CTM_4.putTensor('U2xb_trans',U2xb_trans)
 theta=CTM_4.launch() 
 #print CTM_4, CTM_4.profile()
 #print theta.printDiagram()
 theta.permute([100,300,-300,200,400,-400],3)
# print time.time() - t0, "dudu"

 
 #print theta.trace()#, Contract.printDiagram()
 U4x,  V,  s = TruncateU.setTruncation(theta, chi_dim)
 U4x_trans=copy.copy(U4x)
 U4x_trans.transpose()

# U4x, U4x_trans   = produce_projectives(theta, chi_dim)

 #print "U4x_trans", U4x_trans.printDiagram()



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
 return norm[0]
 





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

#@profile 
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
