import pyUni10 as uni10
import sys
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random
import copy
import time
import Move
 
def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1,True)
  invLt = uni10.Matrix(D, D1,True)
  invLt=Landa2.getBlock(qnum,True)
  #print invLt[0], invLt[1], invLt[2], invLt[3]
  for i in xrange(D):
      invL2[i] = 0 if ((invLt[i].real) < 1.0e-12) else (1.00 / (invLt[i].real))
  invLanda2.putBlock(qnum,invL2)
 return invLanda2

def Renew_dim(dims,dims_val,chi,dim_svd):
   free_par=chi-dims_val
   #print "free_par", free_par
   #print dim_svd
   cnt=0;
   cnt1=-10;
   while cnt < free_par:
    for i in xrange(len(dims)):
     if  ((int(len(dims)/2)+i)<=(len(dims)-1)):
         if (dims[int(len(dims)/2)+i]<dim_svd[int(len(dims)/2)+i]) and  (dims[int(len(dims)/2)+i] > 0):
          dims[int(len(dims)/2)+i]+=1
          cnt+=1
          if cnt >= free_par:
           break
         if (int(len(dims)/2)-i) >= 0 and ((int(len(dims)/2)-i)<=(len(dims)-1)):
          if (dims[int(len(dims)/2)-i]<dim_svd[int(len(dims)/2)-i]) and (dims[int(len(dims)/2)-i] > 0):
           dims[int(len(dims)/2)-i]+=1
           cnt+=1
           if cnt >= free_par:
            break
    
    if (cnt-cnt1) is 0:
      break;
    cnt1=cnt;

   return dims


def setTruncation_long(theta, chi):

    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    #print blk_qnums,theta.printDiagram()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
        #print svs
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#     
#    #print dims_val
#    if dims_val < chi:
#      dims=Renew_dim(dims,dims_val,chi,dim_svd) 

#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#    print svs#, dims_val
    
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),theta.bond(2),theta.bond(3), bdo_mid])
    GB.assign([bdi_mid, theta.bond(4), theta.bond(5),theta.bond(6),theta.bond(7),theta.bond(8),theta.bond(9),theta.bond(10)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    sv_mat = uni10.Matrix(len(svs), len(svs), svs, True)
    sv_mat.resize(int(bdi_mid.dim()), int(bdo_mid.dim()))
    norm = sv_mat.norm()
#    print bdi_mid
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim) * (1.00/norm)  )
#    print LA
    return GA, GB, LA


def setTruncation_long1(theta, chi):

    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    #print blk_qnums,theta.printDiagram()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
        #print svs
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#     
#    #print dims_val
#    if dims_val < chi:
#      dims=Renew_dim(dims,dims_val,chi,dim_svd) 

#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#    print dims#, dims_val
    
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GB.assign([bdi_mid,theta.bond(7), theta.bond(8),theta.bond(9),theta.bond(10)])
    GA.assign([theta.bond(0), theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5),theta.bond(6),bdo_mid])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    sv_mat = uni10.Matrix(len(svs), len(svs), svs, True)
    sv_mat.resize(int(bdi_mid.dim()), int(bdo_mid.dim()))
    norm = sv_mat.norm()
#    print bdi_mid
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim) * (1.00/norm)  )
#    print LA
    return GA, GB, LA




def setTruncation_short(theta, chi):

    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    #print blk_qnums,theta.printDiagram()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
        #print svs
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#     
#    #print dims_val
#    if dims_val < chi:
#      dims=Renew_dim(dims,dims_val,chi,dim_svd) 

#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#    print dims#, dims_val
    
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),theta.bond(2),theta.bond(3), bdo_mid])
    GB.assign([bdi_mid, theta.bond(4), theta.bond(5),theta.bond(6)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    sv_mat = uni10.Matrix(len(svs), len(svs), svs, True)
    sv_mat.resize(int(bdi_mid.dim()), int(bdo_mid.dim()))
    norm = sv_mat.norm()
#    print bdi_mid
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim) * (1.00/norm)  )
#    print LA
    return GA, GB, LA






def setTruncation_short1(theta, chi):

    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    #print blk_qnums,theta.printDiagram()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
        #print svs
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#     
#    #print dims_val
#    if dims_val < chi:
#      dims=Renew_dim(dims,dims_val,chi,dim_svd) 

#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#    print dims#, dims_val
    
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),theta.bond(2), bdo_mid])
    GB.assign([bdi_mid, theta.bond(3), theta.bond(4),theta.bond(5),theta.bond(6)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    sv_mat = uni10.Matrix(len(svs), len(svs), svs, True)
    sv_mat.resize(int(bdi_mid.dim()), int(bdo_mid.dim()))
    norm = sv_mat.norm()
#    print bdi_mid
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim) * (1.00/norm)  )
#    print LA
    return GA, GB, LA




def setTruncation(theta, chi):

    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    #print blk_qnums,theta.printDiagram()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
        #print svs
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#     
#    #print dims_val
#    if dims_val < chi:
#      dims=Renew_dim(dims,dims_val,chi,dim_svd) 

#    dims_val=0
#    for i in xrange(len(dims)):
#     dims_val+=dims[i]
#    print dims#, dims_val
    
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),theta.bond(2), bdo_mid])
    GB.assign([bdi_mid, theta.bond(3), theta.bond(4),theta.bond(5)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    sv_mat = uni10.Matrix(len(svs), len(svs), svs, True)
    sv_mat.resize(int(bdi_mid.dim()), int(bdo_mid.dim()))
    norm = sv_mat.norm()
#    print bdi_mid
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim) * (1.00/norm)  )
#    print LA
    return GA, GB, LA




def sv_merge(svs, bidxs, bidx, sv_mat, chi, len_qn):
#Return the length (the number of items) of an object. The argument may be a sequence (such as a string, bytes, tuple, list, or range) or a collection (such as a dictionary, set, or frozen set).
    #print 'helo'
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        #print 'length', length
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > 0.0e-12):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > 0.0e-12):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > 0.0e-12):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > 0.0e-12):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    #print svs, bidxs,'\n','\n'
    return svs, bidxs

 

def lq_parity(theta):
#    bd1=copy.copy(theta.bond(0))
#    bd2=copy.copy(theta.bond(1))
#    bd1.change(uni10.BD_OUT)
#    bd2.change(uni10.BD_OUT)
    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())
    bd2=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())


    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4)])
    LA=uni10.UniTensor([theta.bond(0),theta.bond(1),bd1,bd2])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    #print theta
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA


def MaxAbs(c):
 blk_qnums = c.blockQnum()
 max_list=[]
 for qnum in blk_qnums:
    c_mat=c.getBlock(qnum)
    max_list.append(c_mat.absMax())
 max_list_f=[abs(x) for x in max_list]
 return max(max_list_f)

def max_ten(a):
 Max_val=MaxAbs(a)
 if ( Max_val < 0.1e-1) or (Max_val > 0.1e+1)   :

  if Max_val >= 1:
   #print ">1",Max_val
   a=a*(1.00/Max_val)
  if Max_val < 1: 
   #print "<1",Max_val
   a=a*(1.00/Max_val)

#  if Max_val >= 1:
#   print ">1",Max_val, Max_val**(1./2.)
#   a=a*(1.00/(Max_val**(1./2.)))
#  if Max_val < 1: 
#   print "<1",Max_val
#   a=a*(1.00/Max_val)

 else: a=a;
 return a



 
def qr_parity(theta):

#    bd1=copy.copy(theta.bond(3))
#    bd2=copy.copy(theta.bond(4))
#    bd1.change(uni10.BD_IN)
#    bd2.change(uni10.BD_IN)
    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4)])
    LA=uni10.UniTensor([bd1,bd2, theta.bond(3),theta.bond(4)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]

    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).qr()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])

#    print LA
    return GA, LA

def norm_Symmetry(LA):
 norm=0
 blk_qnums = LA.blockQnum()
 for qnum in blk_qnums:
  M=LA.getBlock(qnum)  
  norm=norm+(M.norm()*M.norm())
 norm=norm**(1.00/2.00)
 return norm

 
def update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd):

 if fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[0].bond(3).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)



# print q_D
# print  D


 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]
 #print D_dim


 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_a.setLabel([0,1,2,3,4])
 Gamma_b.setLabel([5,6,7,8,9])
 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([10,11,0,5])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa3p=copy.copy(Landa[2])


 Landa1.setLabel([3,6])
 Landa2.setLabel([-2,2])
 Landa3.setLabel([-1,1])
 Landa4.setLabel([4,-4])
 Landa8.setLabel([9,-9])
 Landa3p.setLabel([8,-8])
 Landa7.setLabel([-7,7])


 
 Left=Gamma_a*(Landa2*Landa3*Landa4)
 Left.permute([-4,-1,-2,0,3],3)
 
 q_uni,r_uni=qr_parity(Left)

 q_uni.setLabel([-4,-1,-2,20,200])
 r_uni.setLabel([20,200,0,3])
 r_uni.permute([20,200,0,3],3)


 
# A=r_uni*q_uni
# A.permute([-4,-1,-2,0,3],3)
# print A.elemCmp(Left),A[0],Left[0]

 
 Right=Gamma_b*(Landa8*Landa3p*Landa7)
 Right.permute([5,6,-7,-8,-9],2)

 l_uni,qq_uni=lq_parity(Right)
 
 l_uni.setLabel([5,6,40,400])
 l_uni.permute([5,6,40,400],2)
 qq_uni.setLabel([40,400,-7,-8,-9])
 
 
# qqt_uni=copy.copy( qq_uni)
# qqt_uni.transpose()
# A=l_uni*qq_uni
# A.permute([5,6,-7,-8,-9],2)
 #print A.elemCmp(Right),A,Right

 Theta=(r_uni*Landa1*Hamiltonian)*l_uni
 Theta.permute([10,20,200,11,40,400],3)
 ##printTheta.printDiagram() 
 #blk_qnums=Theta.blockQnum()
 #print Theta.printDiagram(),Theta.getBlock(blk_qnums[0]),Theta.getBlock(blk_qnums[1])


 if fixbond_itebd is 'off':
  U,V,LA=setTruncation(Theta,D_dim) 
#########################################################################################
 elif fixbond_itebd is 'on':
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(3), Theta.bond(4),Theta.bond(5)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
##########################################################################################


 blk_qnums = LA.blockQnum()
 Landa[0].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[0].putBlock(qnum,LA.getBlock(qnum))


 U.setLabel([10,20,200,3])
 U.permute([20,200,10,3],2)

 V.setLabel([6,11,40,400])
 V.permute([11,6,40,400],2)
 
 GsA=U*q_uni
 GsB=V*qq_uni
 
 GsA.permute([10,-1,-2,3,-4],3)
 GsB.permute([11,6,-7,-8,-9],3)


 invLanda2=inverse(Landa2)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda2.setLabel([2,-2])
 invLanda3.setLabel([1,-1])
 invLanda4.setLabel([-4,4])
 
 GsA=GsA*(invLanda2*invLanda3*invLanda4)

 invLanda8=inverse(Landa8)
 invLanda3=inverse(Landa3)
 invLanda7=inverse(Landa7)


 invLanda8.setLabel([-9,9])
 invLanda3.setLabel([-8,8])
 invLanda7.setLabel([7,-7])
 GsB=GsB*(invLanda8*invLanda3*invLanda7)
 ##print'identity', invLanda4.getBlock()*Landa4.getBlock()
 ##print'identity', invLanda3.getBlock()*Landa3.getBlock()
 ##print'identity', invLanda2.getBlock()*Landa2.getBlock()
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)
 

 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])
 #Gamma[0]=uni10.UniTensor(GsA.bond())
 #Gamma[1]=uni10.UniTensor(GsB.bond())


 blk_qnums = GsA.blockQnum()
 Gamma[0].assign(GsA.bond())
 for qnum in blk_qnums:
  Gamma[0].putBlock(qnum,GsA.getBlock(qnum))
  
 blk_qnums = GsB.blockQnum()
 Gamma[1].assign(GsB.bond())
 for qnum in blk_qnums:
  Gamma[1].putBlock(qnum,GsB.getBlock(qnum))
 
#@profile 
###########################################################################################################################
def update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd):
 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]
 #print D_dim


 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_c=copy.copy(Gamma[2])

 Gamma_a.setLabel([55,1,2,3,4])
 Gamma_b.setLabel([56,6,7,8,9])
 Gamma_c.setLabel([54,10,11,12,-2])

 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([51,52,53,54,55,56])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa3p=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa4p=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])


 Landa1.setLabel([3,6])
 Landa2.setLabel([-2,2])
 Landa3.setLabel([-1,1])
 Landa3p.setLabel([8,-8])
 Landa4.setLabel([4,-4])
 Landa4p.setLabel([-11,11])
 Landa8.setLabel([9,-9])
 Landa7.setLabel([-7,7])
 Landa5.setLabel([-10,10])
 Landa6.setLabel([12,-12])


 
 Gamma_c=(((Gamma_c*Landa5)*Landa6)*Landa4p)
 #Gamma_c.permute([-10,-11,-12,54,-2],3)
 
# q_uni,r_uni=qr_parity(Gamma_c)

# q_uni.setLabel([-10,-11,-12,80,81])
# r_uni.setLabel([80,81,54,-2])
# r_uni.permute([80,81,54,-2],3)

 
 Gamma_b=(((Gamma_b*Landa8)*Landa3p)*Landa7)
 #Gamma_b.permute([56,6,-7,-8,-9],2)

# l_uni,qq_uni=lq_parity(Gamma_b)
 
# l_uni.setLabel([56,6,82,83])
# qq_uni.setLabel([82,83,-7,-8,-9])


 Theta=(Gamma_c*Hamiltonian)*(((((Gamma_a*Landa3)*Landa2)*Landa4)*Landa1)*Gamma_b)
 Theta.permute([51,-10,-11,-12,-1,-4,52,-7,-8,-9,53],4)

 if fixbond_itebd is 'off':
  U,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[2].bond(4).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
##########################################################################################



 blk_qnums = LA.blockQnum()
 Landa[1].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[1].putBlock(qnum,LA.getBlock(qnum))
 #print "1", Landa[1]


 U.setLabel([51,-10,-11,-12,84])

 GsC=copy.copy(U)
 GsC.permute([51,-10,-11,-12,84],3)


 invLanda5=inverse(Landa5)
 invLanda4=inverse(Landa4)
 invLanda6=inverse(Landa6)
 invLanda5.setLabel([10,-10])
 invLanda4.setLabel([11,-11])
 invLanda6.setLabel([-12,12])
 
 GsC=(((GsC*invLanda5)*invLanda4)*invLanda6)
 GsC.permute([51,10,11,12,84],3)


 Theta.permute([-7,-8,-9,53,51,-10,-11,-12,-1,-4,52],4)


 if fixbond_itebd is 'off':
  U_1,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[1].bond(1).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U_1=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U_1.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)



 
 blk_qnums = LA.blockQnum()
 Landa[0].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[0].putBlock(qnum,LA.getBlock(qnum))
 Landa[0].setLabel([0,1])
 Landa[0].permute([1,0],1)
 #print "0", Landa[0]



 U_1.setLabel([-7,-8,-9,53,85])
 GsB=copy.copy(U_1)
 GsB.permute([53,85,-7,-8,-9],3)
 invLanda7=inverse(Landa7)
 invLanda3=inverse(Landa3)
 invLanda8=inverse(Landa8)
 invLanda7.setLabel([7,-7])
 invLanda3.setLabel([-8,8])
 invLanda8.setLabel([-9,9])
 GsB=(((GsB*invLanda8)*invLanda3)*invLanda7)
 GsB.permute([53,85,7,8,9],3)
 
 
 U.transpose()
 U_1.transpose()
 U=(Theta*U)*U_1
 U.permute([52,-1,84,85,-4],3)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda3.setLabel([1,-1])
 invLanda4.setLabel([-4,4])

 invLanda1=inverse(Landa[0])
 invLanda2=inverse(Landa[1])
 invLanda1.setLabel([85,-85])
 invLanda2.setLabel([-84,84])


 GsA=(U*invLanda3)*invLanda4*(invLanda1*invLanda2)
 GsA.permute([52,1,-84,-85,4],3)
 #GsA*=(1.0/GsA.norm())
 GsA=max_ten(GsA) 
 #print GsA.printDiagram() 
 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])
 GsC.setLabel([0,1,2,3,4])


 blk_qnums = GsA.blockQnum()
 Gamma[0].assign(GsA.bond())
 for qnum in blk_qnums:
  Gamma[0].putBlock(qnum,GsA.getBlock(qnum))
  
 blk_qnums = GsB.blockQnum()
 Gamma[1].assign(GsB.bond())
 for qnum in blk_qnums:
  Gamma[1].putBlock(qnum,GsB.getBlock(qnum))
 
 blk_qnums = GsC.blockQnum()
 Gamma[2].assign(GsC.bond())
 for qnum in blk_qnums:
  Gamma[2].putBlock(qnum,GsC.getBlock(qnum))


def update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd):
 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]
 #print D_dim


 Gamma_c=copy.copy(Gamma[0])
 Gamma_d=copy.copy(Gamma[1])
 Gamma_b=copy.copy(Gamma[2])

 Gamma_c.setLabel([54,10,11,12,-2])
 Gamma_d.setLabel([55,-12,14,13,-7])
 Gamma_b.setLabel([56,6,7,-1,9])

 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([51,52,53,54,55,56])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa3p=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa5p=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa8p=copy.copy(Landa[7])


 Landa1.setLabel([3,6])
 Landa2.setLabel([-2,2])
 Landa3.setLabel([-1,1])
 Landa4.setLabel([-11,11])
 Landa8.setLabel([9,-9])
 Landa8p.setLabel([-14,14])
 Landa7.setLabel([-7,7])
 Landa5.setLabel([-10,10])
 Landa5p.setLabel([13,-13])
 Landa6.setLabel([12,-12])


 
 Gamma_c=(((Gamma_c*Landa5)*Landa2)*Landa4)
 Gamma_b=(((Gamma_b*Landa8)*Landa3)*Landa1)

 Theta=(Gamma_c*Hamiltonian)*(((((Gamma_d*Landa6)*Landa7)*Landa8p)*Landa5p)*Gamma_b)
 #print  Theta.printDiagram()
 Theta.permute([51,-11,-10,2,52,-14,-13,1,53,-9,3],4)

 if fixbond_itebd is 'off':
  U,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[0].bond(3).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
##########################################################################################



 blk_qnums = LA.blockQnum()
 Landa[5].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[5].putBlock(qnum,LA.getBlock(qnum))


 U.setLabel([51,-11,-10,2,84])
 

 GsC=copy.copy(U)
 GsC.permute([51,-10,-11,84,2],3)
 invLanda5=inverse(Landa5)
 invLanda4=inverse(Landa4)
 invLanda2=inverse(Landa2)
 invLanda5.setLabel([10,-10])
 invLanda4.setLabel([11,-11])
 invLanda2.setLabel([2,-2])
 
 GsC=(((GsC*invLanda5)*invLanda4)*invLanda2)
 GsC.permute([51,10,11,84,-2],3)



 Theta.permute([1,53,-9,3, 51,-11,-10,2,52,-14,-13],4)

 if fixbond_itebd is 'off':
  U_1,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[2].bond(2).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U_1=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U_1.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)


 
 blk_qnums = LA.blockQnum()
 Landa[6].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[6].putBlock(qnum,LA.getBlock(qnum))

 Landa[6].setLabel([0,1])
 Landa[6].permute([1,0],1)
 Landa[6].setLabel([0,1])

 U_1.setLabel([1,53,-9,3,85])
 
 
 GsB=copy.copy(U_1)
 GsB.permute([53,3,85,1,-9],3)
 invLanda1=inverse(Landa1)
 invLanda3=inverse(Landa3)
 invLanda8=inverse(Landa8)
 invLanda1.setLabel([-3,3])
 invLanda3.setLabel([1,-1])
 invLanda8.setLabel([-9,9])
 GsB=(((GsB*invLanda8)*invLanda3)*invLanda1)
 GsB.permute([53,-3,85,-1,9],3)
 
 



 U.transpose()
 U_1.transpose()
 U=(Theta*U)*U_1
 U.permute([52,84,-14,-13,85],3)
 invLanda8=inverse(Landa8)
 invLanda5=inverse(Landa5)
 invLanda8.setLabel([14,-14])
 invLanda5.setLabel([-13,13])

 invLanda6=inverse(Landa[5])
 invLanda7=inverse(Landa[6])
 invLanda6.setLabel([-84,84])
 invLanda7.setLabel([85,-85])


 GsD=(U*invLanda8)*invLanda5*(invLanda6*invLanda7)
 GsD.permute([52,-84,14,13,-85],3)
 GsD=max_ten(GsD) 
 

 
 GsD.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])
 GsC.setLabel([0,1,2,3,4])

 blk_qnums = GsC.blockQnum()
 Gamma[0].assign(GsC.bond())
 for qnum in blk_qnums:
  Gamma[0].putBlock(qnum,GsC.getBlock(qnum))


 blk_qnums = GsD.blockQnum()
 Gamma[1].assign(GsD.bond())
 for qnum in blk_qnums:
  Gamma[1].putBlock(qnum,GsD.getBlock(qnum))
  
 blk_qnums = GsB.blockQnum()
 Gamma[2].assign(GsB.bond())
 for qnum in blk_qnums:
  Gamma[2].putBlock(qnum,GsB.getBlock(qnum))
 
 
 #print "hi_inner", GsA[6], GsB[4]
 #print Landa_tem*(1.0/Landa_tem.norm())
 #return GsA, GsB, Landa_tem*(1.0/Landa_tem.norm()) 







def update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd):
 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]
 #print D_dim

 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_d=copy.copy(Gamma[2])

 Gamma_a.setLabel([54,1,2,3,4])
 Gamma_b.setLabel([55,6,7,8,9])
 Gamma_d.setLabel([56,12,10,11,-7])

 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([51,52,53,54,55,56])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa3p=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa8p=copy.copy(Landa[7])


 Landa1.setLabel([3,6])
 Landa2.setLabel([-2,2])
 Landa3.setLabel([-1,1])
 Landa4.setLabel([4,-4])
 Landa8.setLabel([9,-9])
 Landa3p.setLabel([8,-8])
 Landa7.setLabel([-7,7])
 Landa8p.setLabel([-10,10])
 Landa5.setLabel([11,-11])
 Landa6.setLabel([-12,12])


 
 Gamma_d=(((Gamma_d*Landa5)*Landa6)*Landa8p)
# Gamma_d.permute([-12,-10,-11,56,-7],3)
# 
# q_uni,r_uni=qr_parity(Gamma_d)

# q_uni.setLabel([-12,-10,-11,80,81])
# r_uni.setLabel([80,81,56,-7])
# r_uni.permute([80,81,56,-7],3)

 

 
 Gamma_a=(((Gamma_a*Landa3)*Landa2)*Landa4)
# Gamma_a.permute([-1,-4,-2,54,3],3)

# qq_uni, l_uni=qr_parity(Gamma_a)
# 
# l_uni.setLabel([82,83,54,3])
# l_uni.permute([82,83,54,3],3)
# qq_uni.setLabel([-1,-4,-2,82,83])
# 
# 
## qqt_uni=copy.copy( qq_uni)
## qqt_uni.transpose()
## A=l_uni*qq_uni
## A.permute([5,6,-7,-8,-9],2)
# #print A.elemCmp(Right),A,Right

 Theta=(Gamma_a*Hamiltonian)*(((((Gamma_b*Landa3p)*Landa8)*Landa7)*Landa1)*Gamma_d)

 Theta.permute([53,-12,-10,-11,52,-9,-8,51,-1,-2,-4],4)

 if fixbond_itebd is 'off':
  U,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[2].bond(4).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)



 blk_qnums = LA.blockQnum()
 Landa[6].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[6].putBlock(qnum,LA.getBlock(qnum))


 U.setLabel([53,-12,-10,-11,84])
 
 
 GsD=copy.copy(U)


 invLanda5=inverse(Landa5)
 invLanda8=inverse(Landa8)
 invLanda6=inverse(Landa6)
 invLanda5.setLabel([-11,11])
 invLanda8.setLabel([10,-10])
 invLanda6.setLabel([12,-12])
 
 GsD=(((GsD*invLanda5)*invLanda8)*invLanda6)
 GsD.permute([53,12,10,11,84],3)




 Theta=(Gamma_a*Hamiltonian)*(((((Gamma_b*Landa3p)*Landa8)*Landa7)*Landa1)*Gamma_d)
 Theta.permute([51,-1,-2,-4,52,-9,-8,53,-12,-10,-11],4)


 if fixbond_itebd is 'off':
  U_1,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  #print Gamma_a.printDiagram() 
  q_D=Gamma[0].bond(3).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U_1=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U_1.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)


 blk_qnums = LA.blockQnum()
 Landa[0].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[0].putBlock(qnum,LA.getBlock(qnum))

 U_1.setLabel([51,-1,-2,-4,85])
 GsA=copy.copy(U_1)
 invLanda2=inverse(Landa2)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda2.setLabel([2,-2])
 invLanda3.setLabel([1,-1])
 invLanda4.setLabel([-4,4])
 GsA=(((GsA*invLanda2)*invLanda3)*invLanda4)
 GsA.permute([51,1,2,85,4],3)
 
 
 U.transpose()
 U_1.transpose()
 U=(Theta *U)*U_1
 U.permute([52,85,84,-8,-9],3)
 invLanda3=inverse(Landa3)
 invLanda8=inverse(Landa8)
 invLanda7=inverse(Landa[6])
 invLanda1=inverse(Landa[0])

 invLanda3.setLabel([-8,8])
 invLanda8.setLabel([-9,9])
 invLanda1.setLabel([-85,85])
 invLanda7.setLabel([-84,84])


 GsB=(U*invLanda3)*invLanda8*(invLanda1*invLanda7)
 #print GsB.printDiagram() 
 GsB.permute([52,-85,-84,8,9],3)
 GsB=max_ten(GsB) 
 
 
 
 

 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])
 GsD.setLabel([0,1,2,3,4])


 blk_qnums = GsA.blockQnum()
 Gamma[0].assign(GsA.bond())
 for qnum in blk_qnums:
  Gamma[0].putBlock(qnum,GsA.getBlock(qnum))
  
 blk_qnums = GsB.blockQnum()
 Gamma[1].assign(GsB.bond())
 for qnum in blk_qnums:
  Gamma[1].putBlock(qnum,GsB.getBlock(qnum))
 
 blk_qnums = GsD.blockQnum()
 Gamma[2].assign(GsD.bond())
 for qnum in blk_qnums:
  Gamma[2].putBlock(qnum,GsD.getBlock(qnum))
 
 #print "hi_inner", GsA[6], GsB[4]
 #print Landa_tem*(1.0/Landa_tem.norm())
 #return GsA, GsB, Landa_tem*(1.0/Landa_tem.norm()) 

def update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd):
 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]
 #print D_dim


 Gamma_a=copy.copy(Gamma[0])
 Gamma_c=copy.copy(Gamma[1])
 Gamma_d=copy.copy(Gamma[2])

 Gamma_a.setLabel([54,1,2,3,4])
 Gamma_c.setLabel([55,-9,8,12,-2])
 Gamma_d.setLabel([56,-12,10,11,-7])

 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([51,52,53,54,55,56])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa4p=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa5p=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])


 Landa1.setLabel([3,-3])
 Landa2.setLabel([-2,2])
 Landa3.setLabel([-1,1])
 Landa4.setLabel([4,-4])
 Landa4p.setLabel([-8,8])
 Landa8.setLabel([-10,10])
 Landa7.setLabel([-7,7])
 Landa5.setLabel([9,-9])
 Landa5p.setLabel([11,-11])
 Landa6.setLabel([12,-12])


 
 Gamma_a=(((Gamma_a*Landa3)*Landa1)*Landa4)
# Gamma_a.permute([54,2,-1,-4,-3],2)
# 
# l_uni,q_uni=lq_parity(Gamma_a)

# q_uni.setLabel([80,81,-1,-4,-3])
# l_uni.setLabel([54,2,80,81])
# l_uni.permute([54,2,80,81],2)


 
 Gamma_d=(((Gamma_d*Landa8)*Landa7)*Landa5p)
# Gamma_d.permute([56,-12,7,-11,-10],2)

# r_uni, qq_uni =lq_parity(Gamma_d)
# 
# r_uni.setLabel([56,-12,82,83])
# r_uni.permute([56,-12,82,83],2)

# qq_uni.setLabel([82,83,7,-11,-10])
 
 
# qqt_uni=copy.copy( qq_uni)
# qqt_uni.transpose()
# A=l_uni*qq_uni
# A.permute([5,6,-7,-8,-9],2)
 #print A.elemCmp(Right),A,Right

 Theta=(Gamma_a*Hamiltonian)*(((((Gamma_c*Landa4p)*Landa6)*Landa2)*Landa5)*Gamma_d)

 Theta.permute([51,-1,-3,-4,9,-8,52,53,-10,-11,7],4)

 if fixbond_itebd is 'off':
  U,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[2].bond(4).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)


 blk_qnums = LA.blockQnum()
 Landa[1].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[1].putBlock(qnum,LA.getBlock(qnum))


 Landa[1].setLabel([0,1])
 Landa[1].permute([1,0],1)
 Landa[1].setLabel([0,1])


 U.setLabel([51,-1,-3,-4,84])

 GsA=copy.copy(U)
 GsA.permute([51,-1,84,-3,-4],3)


 invLanda3=inverse(Landa3)
 invLanda1=inverse(Landa1)
 invLanda4=inverse(Landa4)
 invLanda3.setLabel([1,-1])
 invLanda1.setLabel([-3,3])
 invLanda4.setLabel([-4,4])
 
 GsA=(((GsA*invLanda3)*invLanda1)*invLanda4)
 GsA.permute([51,1,84,3,4],3)





 Theta.permute([53,-10,-11,7,51,-1,-3,-4,9,-8,52],4)


 if fixbond_itebd is 'off':
  U_1,V,LA=setTruncation_long(Theta,D_dim) 
 elif fixbond_itebd is 'on':
  D=[]
  #print Gamma_a.printDiagram() 
  q_D=Gamma[0].bond(3).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
 #blk_qnums = q_D
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)
  D_dim=0
  for i in xrange(len(D)):
   D_dim+=D[i]
   
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U_1=uni10.UniTensor([Theta.bond(0), Theta.bond(1),Theta.bond(2),Theta.bond(3), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(4), Theta.bond(5),Theta.bond(6),Theta.bond(7),Theta.bond(8),Theta.bond(9),Theta.bond(10)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
#   print qnum
   #print U.printDiagram()
   svds[qnum] = Theta.getBlock(qnum).svd()
   U_1.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])     )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)


 
 blk_qnums = LA.blockQnum()
 Landa[5].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[5].putBlock(qnum,LA.getBlock(qnum))

 Landa[5].setLabel([0,1])
 Landa[5].permute([1,0],1)
 Landa[5].setLabel([0,1])



 U_1.setLabel([53,-10,-11,7,85])
 GsD=copy.copy(U_1)
 GsD.permute([53,85,-10,-11,7],3)
 invLanda8=inverse(Landa8)
 invLanda5=inverse(Landa5)
 invLanda7=inverse(Landa7)
 invLanda8.setLabel([10,-10])
 invLanda5.setLabel([-11,11])
 invLanda7.setLabel([7,-7])
 GsD=(((GsD*invLanda7)*invLanda5)*invLanda8)
 GsD.permute([53,85,10,11,-7],3)
 
 
 
 
 U.transpose()
 U_1.transpose()
 U=(Theta *U)*U_1
 U.permute([52,9,-8,85,84],3)
 invLanda4=inverse(Landa4)
 invLanda5=inverse(Landa5)
 invLanda4.setLabel([8,-8])
 invLanda5.setLabel([-9,9])

 invLanda2=inverse(Landa[1])
 invLanda6=inverse(Landa[5])
 invLanda2.setLabel([84,-84])
 invLanda6.setLabel([85,-85])
 #print GsC.printDiagram(), invLanda6.printDiagram(), invLanda2.printDiagram()

 GsC=(U*invLanda4)*invLanda5*(invLanda2*invLanda6)
 GsC.permute([52,-9,8,-85,-84],3)
 GsC=max_ten(GsC) 
 

 
 
 GsA.setLabel([0,1,2,3,4])
 GsD.setLabel([0,1,2,3,4])
 GsC.setLabel([0,1,2,3,4])


 blk_qnums = GsA.blockQnum()
 Gamma[0].assign(GsA.bond())
 for qnum in blk_qnums:
  Gamma[0].putBlock(qnum,GsA.getBlock(qnum))
  
 blk_qnums = GsC.blockQnum()
 Gamma[1].assign(GsC.bond())
 for qnum in blk_qnums:
  Gamma[1].putBlock(qnum,GsC.getBlock(qnum))
 
 blk_qnums = GsD.blockQnum()
 Gamma[2].assign(GsD.bond())
 for qnum in blk_qnums:
  Gamma[2].putBlock(qnum,GsD.getBlock(qnum))
 









def update_ulink_eff(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd):

 if fixbond_itebd is 'on':
  D=[]
  q_D=Gamma[0].bond(2).Qlist()
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  #blk_qnums = q_D
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   D.append(dim)


 
# print q_D
# print  D








 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]


 Gamma_a=copy.copy(Gamma[0])
 Gamma_b=copy.copy(Gamma[1])
 Gamma_a.setLabel([0,1,2,3,4])
 Gamma_b.setLabel([5,6,7,8,9])
 Hamiltonian=copy.copy(U)
 Hamiltonian.setLabel([10,11,0,5])
 Landa1=copy.copy(Landa[0])
 Landa2=copy.copy(Landa[1])
 Landa3=copy.copy(Landa[2])
 Landa4=copy.copy(Landa[3])
 Landa5=copy.copy(Landa[4])
 Landa6=copy.copy(Landa[5])
 Landa7=copy.copy(Landa[6])
 Landa8=copy.copy(Landa[7])
 Landa4p=copy.copy(Landa[3])


 Landa1.setLabel([3,-3])
 Landa2.setLabel([9,2])
 Landa3.setLabel([-1,1])
 Landa4.setLabel([4,-4])
 Landa5.setLabel([-6,6])
 Landa4p.setLabel([-7,7])
 Landa6.setLabel([8,-8])




 Gamma_a=Gamma_a*(Landa1*Landa3*Landa4)
 Gamma_a.permute([0,2,-4,-1,-3],2)
 #print "Left", Left.printDiagram() 
 
 l_uni, q_uni=lq_parity(Gamma_a)

 q_uni.setLabel([20,200,-4,-1,-3])

 l_uni.setLabel([0,2,20,200])
 l_uni.permute([0,2,20,200],2)

 
 Gamma_b=Gamma_b*(Landa5*Landa4p*Landa6)
 Gamma_b.permute([-7,-8,-6,5,9],3)

 qq_uni, r_uni=qr_parity(Gamma_b)

 r_uni.setLabel([40,400,5,9])
 r_uni.permute([40,400,5,9],3)
 qq_uni.setLabel([-7,-8,-6,40,400])


 
 
 
 Theta=(r_uni*Landa2*Hamiltonian)*l_uni
 Theta.permute([11,40,400,10,20,200],3)

 if fixbond_itebd is 'off':
  U,V,LA=setTruncation(Theta,D_dim) 
 elif fixbond_itebd is 'on':
#########################################################################
  count=0
  bdi = uni10.Bond(uni10.BD_IN, q_D)
  bdo = uni10.Bond(uni10.BD_OUT, q_D)
  LA=uni10.UniTensor([bdi, bdo])
  U=uni10.UniTensor([Theta.bond(0), Theta.bond(1), Theta.bond(2), bdo])
  V=uni10.UniTensor([bdi, Theta.bond(3), Theta.bond(4), Theta.bond(5)])
  svds = {}
  blk_qnums = Theta.blockQnum()
  degs = bdi.degeneracy()
  for qnum, dim in degs.iteritems():
   svds[qnum] = Theta.getBlock(qnum).svd()
   U.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), D[count]))
   V.putBlock(qnum, svds[qnum][2].resize(D[count], svds[qnum][2].col()))
   LA.putBlock(qnum, svds[qnum][1].resize(D[count], D[count])   )
   count+=1
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
##############################################################################

 #LA.setLabel([2,9])
 #LA.permute([9,2],1)
 blk_qnums = LA.blockQnum()
 Landa[1].assign(LA.bond()) 
 for qnum in blk_qnums:
  Landa[1].putBlock(qnum,LA.getBlock(qnum))
 
 
 U.setLabel([11,40,400,9])
 U.permute([40,400,11,9],2)

 V.setLabel([2,10,20,200],)
 V.permute([10,2,20,200],2)
 
 GsA=V*q_uni
 GsB=U*qq_uni


 
 invLanda1=inverse(Landa1)
 invLanda3=inverse(Landa3)
 invLanda4=inverse(Landa4)
 invLanda1.setLabel([-3,3])
 invLanda3.setLabel([1,-1])
 invLanda4.setLabel([-4,4])
 GsA=GsA*(invLanda1*invLanda3*invLanda4)

 invLanda4=inverse(Landa4)
 invLanda6=inverse(Landa6)
 invLanda5=inverse(Landa5)


 invLanda4.setLabel([7,-7])
 invLanda5.setLabel([6,-6])
 invLanda6.setLabel([-8,8])
 GsB=GsB*(invLanda5*invLanda4*invLanda6)
 
 
 GsA.permute([10,1,2,3,4],3)
 GsB.permute([11,6,7,8,9],3)
 
 GsA.setLabel([0,1,2,3,4])
 GsB.setLabel([0,1,2,3,4])

 blk_qnums = GsA.blockQnum()
 Gamma[0].assign(GsA.bond()) 
 for qnum in blk_qnums:
  Gamma[0].putBlock(qnum,GsA.getBlock(qnum))
  
 blk_qnums = GsB.blockQnum()
 Gamma[1].assign(GsB.bond()) 
 for qnum in blk_qnums:
  Gamma[1].putBlock(qnum,GsB.getBlock(qnum)) 





