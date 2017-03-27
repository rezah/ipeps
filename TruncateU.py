import pyUni10 as uni10
import sys
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random
import copy
import time

def setTruncation(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),theta.bond(2), bdo_mid])
    GB.assign([bdi_mid, theta.bond(3), theta.bond(4),theta.bond(5)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
    return GA, GB, LA


def svd_parity(theta):

    LA=uni10.UniTensor([theta.bond(0), theta.bond(1)])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1)])
    GB=uni10.UniTensor([theta.bond(0), theta.bond(1)])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
    for qnum in blk_qnums:
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0])
        GB.putBlock(qnum, svd[2])
        LA.putBlock(qnum, svd[1])
    return GA, LA,GB



def svd_parity1(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])
    GB=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])


    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        if (GA.getBlock(qnum).row() == svds[qnum][0].row()) and (GA.getBlock(qnum).col() == svds[qnum][0].col()):
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
        else:
         GA.putBlock(qnum, svds[qnum][0].resize(GA.getBlock(qnum).row(),GA.getBlock(qnum).col()))
         GB.putBlock(qnum, svds[qnum][2].resize(GB.getBlock(qnum).row(),GB.getBlock(qnum).col()))
         LA.putBlock(qnum, svds[qnum][1].resize(LA.getBlock(qnum).row(),LA.getBlock(qnum).col()))


#    print LA
    return GA, LA, GB

def svd_parity2(theta):

    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())
    bd2=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())
    bd3=uni10.Bond(uni10.BD_OUT,theta.bond(2).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),bd1,bd2,bd3])
    LA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),bd1,bd2,bd3])
    GB=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])


    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        print GA.getBlock(qnum).col(), GA.getBlock(qnum).row(), svds[qnum][0].col(), svds[qnum][0].row()
        GA.putBlock(qnum, svds[qnum][0])
        
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

#    print LA
    return GA, LA, GB







def   Sqrt(Landa):
  Landa_cp=copy.copy(Landa)
  blk_qnums=Landa.blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.getBlock(qnum).col())
   Landa_cpm=Landa_cp.getBlock(qnum)
   Landam=Landa_cp.getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     if Landam[i*D+j] > 1.0e-12:
      Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
     else:
      Landa_cpm[i*D+j]=0
   Landa_cp.putBlock(qnum,Landa_cpm)
  return Landa_cp 

def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1)
  invLt = uni10.Matrix(D, D1)
  invLt=Landa2.getBlock(qnum)
  for i in xrange(D):
    for j in xrange(D1):
     invL2[i*D1+j] = 0 if ((invLt[i*D1+j].real) < 1.0e-12) else (1.00 / (invLt[i*D1+j].real))
  invLanda2.putBlock(qnum,invL2)
 return invLanda2


def lq_parity(theta):
    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())
    bd2=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())    
    bd3=uni10.Bond(uni10.BD_OUT,theta.bond(2).Qlist())    
    
    
    LA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),bd1,bd2, bd3])
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA













def sv_merge(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
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
                    if (ori_svs[cur1] > 1.0e-12):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > 1.0e-12):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > 1.0e-12):
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
       elif (sv_mat[0] > 1.0e-12):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs


