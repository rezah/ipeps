import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time

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




def setTruncation(theta, chi):

    LA=uni10.UniTensor([theta.bond(0), theta.bond(2)])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1), theta.bond(2)])
    GB=uni10.UniTensor([theta.bond(0), theta.bond(2), theta.bond(3)])
    svds = {}
    blk_qnums = theta.blockQnum()

#    print '\n', blk_qnums
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        #print_inline(svds[qnum][1])
        dim_svd.append(int(svds[qnum][1].col()))
        #print svds[qnum][1]
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
        #print svs
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    #print dims#,blk_qnums;
    
    dims_val=0
    for i in xrange(len(dims)):
     dims_val+=dims[i]
     
    #print dims_val
    if dims_val < chi:
      dims=Renew_dim(dims,dims_val,chi,dim_svd) 

    dims_val=0
    for i in xrange(len(dims)):
     dims_val+=dims[i]
    #print dims, dims_val

    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([GA.bond(0), GA.bond(1), bdo_mid])
    GB.assign([bdi_mid, GB.bond(1), GB.bond(2)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
#    sv_mat = uni10.Matrix(bdi_mid.dim(), bdo_mid.dim(), svs, True)
#    norm = sv_mat.norm()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )

#    print LA
    return GA, GB, LA







def setTruncation1(theta, chi,q_chi):
    bdi = uni10.Bond(uni10.BD_IN, q_chi)
    bdo = uni10.Bond(uni10.BD_OUT, q_chi)
    LA=uni10.UniTensor([bdi, bdo])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1), bdo])
    GB=uni10.UniTensor([bdi, theta.bond(2), theta.bond(3)])
    svds = {}
    blk_qnums = theta.blockQnum()
    degs = bdi.degeneracy()

    dim_svd=[]
    for qnum, dim in degs.iteritems():
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0].resize(svds[qnum][0].row(), dim))
        GB.putBlock(qnum, svds[qnum][2].resize(dim, svds[qnum][2].col()))
        LA.putBlock(qnum, svds[qnum][1].resize(dim, dim)  )

    return GA, GB, LA

















def print_inline(s):
 print '\n'
 for i in xrange(int(s.col())):
  print s[i]
 print '\n'


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
    #print svs, bidxs,'\n','\n'
    return svs, bidxs


