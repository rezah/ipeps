import pyUni10 as uni10
import sys
import numpy as np
from numpy import linalg as LA
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


def Mat_np_to_Uni(Mat_np):
 d0=np.size(Mat_np,0)
 d1=np.size(Mat_np,1)
 Mat_uni=uni10.Matrix(d0,d1)
 for i in xrange(d0):
  for j in xrange(d1):
   Mat_uni[i*d1+j]=Mat_np[i,j]
 return  Mat_uni


def Mat_nptoUni(Mat_np):
 d0=np.size(Mat_np,0)
 Mat_uni=uni10.Matrix(d0,d0, True)
 for i in xrange(d0):
   Mat_uni[i]=Mat_np[i]
 return  Mat_uni

 
def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np

def eig_np(A):
 D_eig=[A]*2
 A_np=Mat_uni_to_np(A)
 w, v = LA.eig(A_np)
 D_eig[0]=Mat_nptoUni(w)
 D_eig[1]=Mat_np_to_Uni(v)
 return D_eig

def  make_Q(q_vec): 
 D=int(q_vec[0].row())
 m=len(q_vec)
 Q=uni10.Matrix(D, m)
 for i in xrange(m):
  for j in xrange(D):
    Q[j*m+i]=q_vec[i][j]
 return Q

def return_vec(A, index ):
 D=int(A.row())
 vec_tem=uni10.Matrix(D,1)
 for i in xrange(D): 
  vec_tem[i]=A[i*D+index].real
 return vec_tem


def find_maxindex(A):
 D=int(A.row())
 max_val=0
 index=0
 for i in xrange(D):
  if (i == 0) or ( max_val < abs(A[i]) ):
   max_val = abs(A[i])
   index=i
 return max_val, index



def Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3):
 CTM_1 = uni10.Network("Network/Right.net")
 CTM_1.putTensor('Vec_uni',Vec_uni)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 vec.permute([21,14,-14 , 7,-7,0],6)
 #print vec.printDiagram() 
 Vec_M=vec.getBlock()
 return Vec_M


def Multi_u(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4):
 CTM_1 = uni10.Network("Network/Up.net")
 CTM_1.putTensor('Vec_uni',Vec_uni)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 vec.permute([17, 18, -18, 19, -19,  20],6)
 Vec_M=vec.getBlock()
 return Vec_M


def Multi_d(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4):

 Vec_uni1=copy.copy(Vec_uni)
 Vec_uni1.transpose()

 CTM_1 = uni10.Network("Network/Down.net")
 CTM_1.putTensor('Vec_uni1',Vec_uni1)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 vec.permute([3, 4, -4 , 5, -5 , 6],0)

 vec.transpose()
 Vec_M=vec.getBlock()

 return Vec_M






def Multi_l(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3):

 Vec_uni1=copy.copy(Vec_uni)
 Vec_uni1.transpose()
 #print Vec_uni1.printDiagram()
 CTM_1 = uni10.Network("Network/Left.net")
 CTM_1.putTensor('Vec_uni1',Vec_uni1)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([23, 16, -16 , 9,-9, 2],0)
 vec.transpose()
 Vec_M=vec.getBlock()
 return Vec_M


def make_vec_right(c2,Ta2,Tb2,c3):

 Ta2.setLabel([-5,3,-3,4])
 Tb2.setLabel([1,2,-2,-5])
 c2.setLabel([6,4])
 c3.setLabel([5,1])
 vec_right=(Ta2*Tb2)*(c2*c3)
 vec_right.permute([5,2,-2,3,-3,6],6)
 return vec_right

def make_vec_left(c1,Tb4,Ta4,c4):

 Tb4.setLabel([-5,3,-3,4])
 Ta4.setLabel([1,2,-2,-5])
 c1.setLabel([4,6])
 c4.setLabel([1,5])
 vec_left=(Tb4*Ta4)*(c1*c4)
 vec_left.permute([5,2,-2,3,-3,6],0)
 vec_left.transpose()
 return vec_left

def make_vec_down(c4,Ta3, Tb3,c3):

 Tb3.setLabel([-5,3,-3,4])
 Ta3.setLabel([1,2,-2,-5])
 c3.setLabel([4,6])
 c4.setLabel([5,1])
 vec_down=(Tb3*Ta3)*(c3*c4)
 vec_down.permute([5, 2, -2 , 3, -3 , 6],0)
 vec_down.transpose()
 return vec_down



def make_vec_up(c1,Tb1, Ta1,c2):

 Ta1.setLabel([-5,3,-3,4])
 Tb1.setLabel([1,2,-2,-5])
 c2.setLabel([4,6])
 c1.setLabel([5,1])
 vec_up=(Ta1*Tb1)*(c1*c2)
 vec_up.permute([5,2,-2,3,-3,6],6)
 return vec_up

def distance(theta,A):
   blk_qnums = theta.blockQnum()
   val=0
   for qnum in blk_qnums:
    T1=theta.getBlock(qnum)
    T2=A.getBlock(qnum)
    print "row", theta.getBlock(qnum).row(), qnum 
    for  i  in xrange( int(theta.getBlock(qnum).row()) ):
     if abs(T1[i]) > 1.00e-11:  
      val=val+abs((T1[i]-T2[i]) / T1[i])
      #if abs((T1[i]-T2[i]) / T1[i]) > 0.00004: print "hi", T1[i], T2[i], i 
     else: val=val+(T1[i]-T2[i]); #print "hi",  T1[i]-T2[i] 
   return val 
 
def distance1(theta,A):
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
 
 
def decompos_r(Vec_uni):
 Vec_uni.setLabel([1,2,-2,3,-3,4])
 Vec_uni.permute([1,2,-2,3,-3,4],1)
 chi_dim=Vec_uni.bond(0).dim()


 c3, V, s=TruncateU.setTruncation1(Vec_uni,chi_dim )
 s.setLabel([1,0])
 V.setLabel([0,2,-2,3,-3,4])
 c3.setLabel([1,-10])
 c3.permute([1,-10],1)

 V=s*V
 V.permute([1,2,-2,3,-3,4],5)
 chi_dim=Vec_uni.bond(5).dim()
 V1, c2, s =TruncateU.setTruncation2(V,chi_dim )
 c2.setLabel([-20,4])
 c2.permute([4,-20],2)
 #print c2.printDiagram()
 s.setLabel([0,4])
 V1.setLabel([1,2,-2,3,-3,0])
 V1=V1*s
 V1.permute([1,2,-2,3,-3,4],3)

# chi_dim1=(Vec_uni.bond(0).dim())*(Vec_uni.bond(1).dim())*(Vec_uni.bond(2).dim())
# chi_dim2=(Vec_uni.bond(5).dim())*(Vec_uni.bond(4).dim())*(Vec_uni.bond(3).dim())
# if chi_dim1 >= chi_dim2:
#   chi_dim=chi_dim2
# else:  
#   chi_dim=chi_dim1

 chi_dim1=(Vec_uni.bond(0).dim())
 chi_dim2=(Vec_uni.bond(5).dim())
 if chi_dim1 >= chi_dim2:
   chi_dim=chi_dim2
 else:  
   chi_dim=chi_dim1



 Tb2, Ta2, s =TruncateU.setTruncation(V1,chi_dim )
 Tb2.setLabel([-10,2,-2,10])
 s.setLabel([10,11])
 Tb2=Tb2*s
 Tb2.permute([-10,2,-2,11],3)
 Ta2.setLabel([11,3,-3,-20])
 Ta2.permute([11,3,-3,-20],3)
 
# vec_tem=Ta2*Tb2*c2*c3
# vec_tem.permute([1,2,-2,3,-3,4],6)
# Vec_uni.permute([1,2,-2,3,-3,4],6)
# print "hi", vec_tem[10], Vec_uni[10], distance(vec_tem, Vec_uni), "end" 
 #print Ta2.printDiagram(), Tb2.printDiagram() 
 return c2,Ta2,Tb2,c3

def decompos_u(Vec_uni):
 Vec_uni.setLabel([1,2,-2,3,-3,4])
 Vec_uni.permute([1,2,-2,3,-3,4],1)
 chi_dim=Vec_uni.bond(0).dim()


 c1, V, s=TruncateU.setTruncation1(Vec_uni,chi_dim )
 s.setLabel([1,0])
 V.setLabel([0,2,-2,3,-3,4])
 c1.setLabel([1,-10])
 c1.permute([1,-10],1)
 #print "c1", c1.printDiagram()
 V=s*V
 V.permute([1,2,-2,3,-3,4],5)
 chi_dim=Vec_uni.bond(5).dim()
 V1,c2,s =TruncateU.setTruncation2(V,chi_dim )
 c2.setLabel([-20,4])
 c2.permute([-20,4],2)

 #print "c2", c2.printDiagram()

 s.setLabel([0,4])
 V1.setLabel([1,2,-2,3,-3,0])
 V1=V1*s
 V1.permute([1,2,-2,3,-3,4],3)

# chi_dim1=(Vec_uni.bond(0).dim())*(Vec_uni.bond(1).dim())*(Vec_uni.bond(2).dim())
# chi_dim2=(Vec_uni.bond(5).dim())*(Vec_uni.bond(4).dim())*(Vec_uni.bond(3).dim())
# if chi_dim1 >= chi_dim2:
#   chi_dim=chi_dim2
# else:  
#   chi_dim=chi_dim1

 chi_dim1=(Vec_uni.bond(0).dim())
 chi_dim2=(Vec_uni.bond(5).dim())
 if chi_dim1 >= chi_dim2:
   chi_dim=chi_dim2
 else:  
   chi_dim=chi_dim1

 Tb1, Ta1, s =TruncateU.setTruncation(V1,chi_dim )
 Tb1.setLabel([-10,2,-2,10])
 s.setLabel([10,11])
 Tb1=Tb1*s
 Tb1.permute([-10,2,-2,11],3)
 Ta1.setLabel([11,3,-3,-20])
 Ta1.permute([11,3,-3,-20],3)
 
# vec_tem=Ta1*Tb1*c1*c2
# vec_tem.permute([1,2,-2,3,-3,4],6)
# Vec_uni.permute([1,2,-2,3,-3,4],6)
# print "hi", vec_tem[0], Vec_uni[0],  distance(vec_tem, Vec_uni), "end" 
 return c1,Tb1,Ta1,c2


def decompos_l(Vec_uni):
 Vec_uni.transpose()
 Vec_uni.setLabel([1,2,-2,3,-3,4])
 Vec_uni.permute([1,2,-2,3,-3,4],1)
 chi_dim=Vec_uni.bond(0).dim()

 c4, V, s=TruncateU.setTruncation1(Vec_uni,chi_dim )
 s.setLabel([1,0])
 V.setLabel([0,2,-2,3,-3,4])
 c4.setLabel([1,-10])
 c4.permute([-10,1],0)
 #print "c4", c4.printDiagram()
 V=s*V
 V.permute([1,2,-2,3,-3,4],5)
 chi_dim=Vec_uni.bond(5).dim()
 V1, c1, s =TruncateU.setTruncation2(V,chi_dim )
 c1.setLabel([-20,4])
 c1.permute([-20,4],1)
 
 #print "c1", c1.printDiagram()
 
 
 s.setLabel([0,4])
 V1.setLabel([1,2,-2,3,-3,0])
 V1=V1*s
 V1.permute([1,2,-2,3,-3,4],3)

# chi_dim1=(Vec_uni.bond(0).dim())*(Vec_uni.bond(1).dim())*(Vec_uni.bond(2).dim())
# chi_dim2=(Vec_uni.bond(5).dim())*(Vec_uni.bond(4).dim())*(Vec_uni.bond(3).dim())
# if chi_dim1 >= chi_dim2:
#   chi_dim=chi_dim2
# else:  
#   chi_dim=chi_dim1

 chi_dim1=(Vec_uni.bond(0).dim())
 chi_dim2=(Vec_uni.bond(5).dim())
 if chi_dim1 >= chi_dim2:
   chi_dim=chi_dim2
 else:  
   chi_dim=chi_dim1


 Ta4, Tb4, s =TruncateU.setTruncation(V1,chi_dim )
 Ta4.setLabel([-10,2,-2,10])
 s.setLabel([10,11])
 Ta4=Ta4*s
 #print "hi0", Ta4.printDiagram(), s.printDiagram()

 Ta4.permute([-10,2,-2,11],1)
 Tb4.setLabel([11,3,-3,-20])
 Tb4.permute([11,3,-3,-20],1)
 #print "hi", Ta4.printDiagram(), Tb4.printDiagram()
# vec_tem=Ta4*Tb4*c4*c1
# vec_tem.permute([1,2,-2,3,-3,4],0)
# Vec_uni.permute([1,2,-2,3,-3,4],0)
# print "hi", vec_tem[10], Vec_uni[10], distance1(vec_tem, Vec_uni), "end" 
 return c1,Tb4,Ta4,c4



def decompos_d(Vec_uni):
 Vec_uni.transpose()
 Vec_uni.setLabel([1,2,-2,3,-3,4])
 Vec_uni.permute([1,2,-2,3,-3,4],1)
 chi_dim=Vec_uni.bond(0).dim()


 c4, V, s=TruncateU.setTruncation1(Vec_uni,chi_dim )
 s.setLabel([1,0])
 V.setLabel([0,2,-2,3,-3,4])
 c4.setLabel([1,-10])
 c4.permute([1,-10],0)

 V=s*V
 V.permute([1,2,-2,3,-3,4],5)
 chi_dim=Vec_uni.bond(5).dim()
 V1,c3, s =TruncateU.setTruncation2(V,chi_dim )
 c3.setLabel([-20,4])
 c3.permute([-20,4],1)
 
 
 
 s.setLabel([0,4])
 V1.setLabel([1,2,-2,3,-3,0])
 V1=V1*s
 V1.permute([1,2,-2,3,-3,4],3)


# chi_dim1=(Vec_uni.bond(0).dim())*(Vec_uni.bond(1).dim())*(Vec_uni.bond(2).dim())
# chi_dim2=(Vec_uni.bond(5).dim())*(Vec_uni.bond(4).dim())*(Vec_uni.bond(3).dim())
# if chi_dim1 >= chi_dim2:
#   chi_dim=chi_dim2
# else:  
#   chi_dim=chi_dim1


 chi_dim1=(Vec_uni.bond(0).dim())
 chi_dim2=(Vec_uni.bond(5).dim())
 if chi_dim1 >= chi_dim2:
   chi_dim=chi_dim2
 else:  
   chi_dim=chi_dim1



 Ta3, Tb3, s =TruncateU.setTruncation(V1,chi_dim )
 Ta3.setLabel([-10,2,-2,10])
 s.setLabel([10,11])
 Ta3=Ta3*s
 Ta3.permute([-10,2,-2,11],1)
 Tb3.setLabel([11,3,-3,-20])
 Tb3.permute([11,3,-3,-20],1)
 
# vec_tem=Ta3*Tb3*c4*c3
# vec_tem.permute([1,2,-2,3,-3,4],0)
# Vec_uni.permute([1,2,-2,3,-3,4],0)
# print "hi", vec_tem[10], Vec_uni[10], distance1(vec_tem, Vec_uni), "end" 
 return c4,Ta3,Tb3,c3






def ED_right(c2,Ta2,Tb2,c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3):

 Vec_uni=make_vec_right(c2,Ta2,Tb2,c3)

 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D=",  D

 m=5
 W=1
 num=0
 E1=0
 p=0

 while p  <  (W+1):
  #print "norm", p, Vec_F.norm(),Vec_F[0], Vec_F[1], Vec_F[2], Vec_F[3] 
  r=copy.copy(Vec_F)
  #r = r* (1.00/r.norm()) 
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()
  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3)
   for i in xrange(j+1):
    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 
  D_eig=eig_np(h)
  Lambda, index=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  #print 'r', Lambda
  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m+=5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m+=5
     E1=Lambda
  p+=1
 #Vec_F=Vec_F*(1.00/Vec_F.norm())
 Vec_uni.putBlock(Vec_F)
 c2,Ta2,Tb2,c3=decompos_r(Vec_uni)
 return c2,Ta2,Tb2,c3

def ED_left(c1,Tb4,Ta4,c4, a, b, c, d, Tb1, Ta1, Ta3, Tb3):

 Vec_uni=make_vec_left(c1,Tb4,Ta4,c4)

 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D", D

 m=5
 W=1
 num=0
 E1=0
 p=0

 while p  <  (W+1):
  r=copy.copy(Vec_F)
  #r = r* (1.00/r.norm()) 
  
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()

  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_l(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3)
   #r.resize(D,1)
   for i in xrange(j+1):
    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 

  D_eig=eig_np(h)
  Lambda, index=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  #print "l", Lambda

  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m=m+5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m=m+5
     E1=Lambda
  p+=1
 #Vec_F=Vec_F*(1.00/Vec_F.norm())
 Vec_uni.putBlock(Vec_F)
 c1,Tb4,Ta4,c4=decompos_l(Vec_uni)
 return c1,Tb4,Ta4,c4

def ED_up(c1,Tb1, Ta1,c2, a, b, c, d, Tb2, Ta2, Ta4, Tb4):

 Vec_uni=make_vec_up(c1,Tb1, Ta1,c2)

 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D=",  D

 m=5
 W=1
 num=0
 E1=0
 p=0

 while p  <  (W+1):
  r=copy.copy(Vec_F)
  #r = r* (1.00/r.norm()) 
  
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()

  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_u(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4)
   #r.resize(D,1)
   for i in xrange(j+1):

    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 

  D_eig=eig_np(h)
  Lambda, index=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  #print "u", Lambda

  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m=m+5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m=m+5
     E1=Lambda
  p+=1
 #Vec_F=Vec_F*(1.00/Vec_F.norm())
 Vec_uni.putBlock(Vec_F)
 c1,Tb1,Ta1,c2=decompos_u(Vec_uni)
 return c1,Tb1,Ta1,c2


def ED_down(c4,Ta3, Tb3,c3, a, b, c, d, Tb2, Ta2, Ta4, Tb4):

 Vec_uni=make_vec_down(c4,Ta3,Tb3,c3)

 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D", D

 m=5
 W=1
 num=0
 E1=0
 p=0

 while p  <  (W+1):
  r=copy.copy(Vec_F)
  #r = r* (1.00/r.norm()) 
  
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()

  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_d(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4)

   for i in xrange(j+1):
    q_vec_trans=copy.copy(q_vec[i])
    q_vec_trans.transpose()
    dot_vec=q_vec_trans*r
    h[i*m+j]=dot_vec.trace()
    
    r=r+((-1.00)*(h[i*m+j]*q_vec[i]))
   if j<(m-1):
    h[((j+1)*m)+j]=r.norm()
    if r.norm() > 1.0e-8:
     q_vec.append(r*(1.00/r.norm()))
    else:  break; 

  D_eig=eig_np(h)
  Lambda, index=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  #print "d", Lambda

  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  if p==W and num==0:
   p=-1
   m+=5
   E1=copy.copy(Lambda)
   num+=1
  elif p==W:
   num+=1
   if abs(Lambda) > 1.e-9: 
    if  (((abs(Lambda-E1))/(abs(Lambda)))< 1.e-9): num+=1
    elif m<=20:
     p=-1
     m=m+5
     E1=Lambda
   else:
    if  (abs(Lambda-E1))< 1.e-9:
     num+=1
    elif m<=20: 
     p=-1
     m=m+5
     E1=Lambda
  p+=1
 #Vec_F=Vec_F*(1.00/Vec_F.norm())
 Vec_uni.putBlock(Vec_F)
 c4,Ta3,Tb3,c3=decompos_d(Vec_uni)
 
 
 
 return c4,Ta3,Tb3,c3










def  add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D):


 c2,Ta2,Tb2,c3=ED_right(c2,Ta2,Tb2,c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3)
 c1,Tb4,Ta4,c4=ED_left(c1,Tb4,Ta4,c4, a, b, c, d, Tb1, Ta1, Ta3, Tb3)
 c1,Tb1,Ta1,c2=ED_up(c1,Tb1, Ta1,c2, a, b, c, d, Tb2, Ta2, Ta4, Tb4)
 c4,Ta3,Tb3,c3=ED_down(c4,Ta3,Tb3,c3, a, b, c, d, Tb2, Ta2, Ta4, Tb4)


# Vec_uni=make_vec_down(c4,Ta3,Tb3,c3)
# Vec_uni_u=make_vec_up(c1,Tb1, Ta1,c2)
# A=Vec_uni*Vec_uni_u
# print "Norm", A[0]

# Vec_uni=make_vec_down(c4,Ta3,Tb3,c3)
# Vec_uni_u=make_vec_up(c1,Tb1, Ta1,c2)
# A=Vec_uni*Vec_uni_u
# print "Norm", A[0]

# Vec_uni=make_vec_right(c2,Ta2,Tb2,c3)
# Vec_uni_l=make_vec_left(c1,Tb4,Ta4,c4)
# A=Vec_uni*Vec_uni_l
# print "Norm", A[0]

# Vec_uni=make_vec_right(c2,Ta2,Tb2,c3)
# Vec_uni_l=make_vec_left(c1,Tb4,Ta4,c4)
# A=Vec_uni*Vec_uni_l
# print "Norm", A[0]

# print "c1", c1.printDiagram()
# print "c2", c2.printDiagram()
# print "c3", c3.printDiagram()
# print "c4", c4.printDiagram()

# print "Tb1", Tb1.printDiagram()
# print "Ta1", Ta1.printDiagram()

# print "Tb2", Tb2.printDiagram()
# print "Ta2", Ta2.printDiagram()

# print "Ta3", Ta3.printDiagram()
# print "Tb3", Tb3.printDiagram()

# print "Ta4", Ta4.printDiagram()
# print "Tb4", Tb4.printDiagram()

 return c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4
 
 
 
 
 
 
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
 

 
 
 
 
 
 
 
