import pyUni10 as uni10
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import math
import copy
import time
import Move
import MoveCorboz
import MoveFull
import basicA
import basicB
import basicC
import numpy as np
from numpy import linalg as LA


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


def max_ten(a):

 if ( MaxAbs(a) < 0.50e-1) or (MaxAbs(a) > 0.50e+1)   :
  a=a*(1.00/MaxAbs(a));
 else: a=a;
 
 return a


def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(0.5)*uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0])
  return Mat 

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(0.5)*uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);
  return Mat 

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(0.5)*uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);
  return Mat 


def matIden():
    spin_t=0.5
    dimT = int(2*spin_t + 1)
    Mat=uni10.Matrix(dimT, dimT,[1,0,0,1])
    return Mat


#cx_mat X1(3,3);  X1.zeros();         X1(0,1)=1; X1(1,0)=1; X1(1,2)=1,X1(2,1)=1;
#cx_mat Z1(3,3);  Z1.zeros();         Z1(0,0)=1; Z1(1,1)=0; Z1(2,2)=-1;
#cx_mat Y1(3,3);  Y1.zeros();         Y1(0,1)=-1; Y1(1,0)=1; Y1(1,2)=-1,Y1(2,1)=1;


#def matSx():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=(1.0/(2.0**(0.5)))*uni10.Matrix(dimT, dimT,[0, 1.0, 0 ,1.0,0, 1.0,0,1.0,0])
#    return Mat 
#    
#def matSy():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=(1.0/(2.0**(0.5)))*uni10.Matrix(dimT, dimT,[0,-1.0,0,1.0,0,-1.0,0,1.0,0])
#    return Mat 
#   

#def matSz():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=uni10.Matrix(dimT, dimT,[1,0,0,0,0,0,0,0,-1])
#    return Mat

#def matIden():
#    spin_t=1
#    dimT = int(2*spin_t + 1)
#    Mat=uni10.Matrix(dimT, dimT,[1,0,0,0,1,0,0,0,1])
#    return Mat



def Init_env(Env):
 c1=copy.copy(Env[0])
 c2=copy.copy(Env[1])
 c3=copy.copy(Env[2]) 
 c4=copy.copy(Env[3]) 
 Ta1=copy.copy(Env[4])
 Ta2=copy.copy(Env[5])
 Ta3=copy.copy(Env[6])
 Ta4=copy.copy(Env[7])
 Tb1=copy.copy(Env[8])
 Tb2=copy.copy(Env[9])
 Tb3=copy.copy(Env[10])
 Tb4=copy.copy(Env[11])
 return  c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4




def CorrelationH(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,distance_final,fileCorr,fileCorrLength):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 if Model is "Heisenberg":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)(-1.0)*uni10.otimes(sy,sy)
  H=uni10.otimes(sz,iden)+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1=uni10.otimes(iden,sz)+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 if Model is "Heisenberg_Z2":
  #print d_phys
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  iden = matIden()
  szt=uni10.otimes(sz,iden)
  H.setRawElem(szt)
  szt1=uni10.otimes(iden,sz)
  H1.setRawElem(szt1)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])
 if Model is "Heisenberg_U1":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH_tem=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy)
  H_tem=uni10.otimes(sz,iden)#+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1_tem=uni10.otimes(iden,sz)#+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setRawElem(HH_tem)
  H.setRawElem(H_tem)
  H1.setRawElem(H1_tem)
  Iden=copy.copy(H)
  Iden.identity()
  #print HH_tem, HH, H_tem, H, H1_tem, H1
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 vec_left=make_vleft(Tb4,Ta4,c1,c4)
 vec_right=make_vright(Ta2,Tb2,c2,c3)
 ap=make_ap_openindex(a_u)
 bp=make_ap_openindex(b_u)
 cp=make_ap_openindex(c_u)
 dp=make_ap_openindex(d_u)
 dis_val_list=[]
 Corr_val_list=[]
 dis_val_list1=[]
 Corr_val_list1=[]
 dis_val_list2=[]
 Corr_val_list2=[]
 dis_val_list3=[]
 Corr_val_list3=[]


# vec_right_copy=copy.copy(vec_right)
# Corr_length=ED_right(c2, Ta2, Tb2, c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3,vec_right_copy)
# print "Corr_length",Corr_length
# fileCorrLength.write(str(Corr_length)  + "\n")
# fileCorrLength.flush()



 print "\n"
#######################a-a#####################################################

 vec_left_1=Make_first_vecleft_a(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)
 vec_right_1=Make_first_vecright_a(vec_right, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list.append(2)
 Corr_val_list.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list.append(dis_val)
  Corr_val_list.append(Corr_val) 
###################################################################################

 print "\n"
#######################b-b#####################################################

 vec_left_1=Make_first_vecleft_b(vec_left, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)
 vec_right_1=Make_first_vecright_b(vec_right, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list1.append(2)
 Corr_val_list1.append(Corr_val)

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1.append(dis_val)
  Corr_val_list1.append(Corr_val) 
###################################################################################


 print "\n"
#######################c-c#####################################################

 vec_left_1=Make_first_vecleft_c(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)
 vec_right_1=Make_first_vecright_c(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list2.append(2)
 Corr_val_list2.append(Corr_val)

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2.append(dis_val)
  Corr_val_list2.append(Corr_val) 
###################################################################################

 print "\n"
#######################d-d#####################################################

 vec_left_1=Make_first_vecleft_d(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)
 vec_right_1=Make_first_vecright_d(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
 dis_val_list3.append(2)
 Corr_val_list3.append(Corr_val)

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3.append(dis_val)
  Corr_val_list3.append(Corr_val) 
###################################################################################

############################################################################################################


 dis_val_listo=[]
 Corr_val_listo=[]
 dis_val_list1o=[]
 Corr_val_list1o=[]
 dis_val_list2o=[]
 Corr_val_list2o=[]
 dis_val_list3o=[]
 Corr_val_list3o=[]

 print "\n"

######################a-bo#####################################################

 vec_left_1=Make_first_vecleft_a(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)
 vec_right_1=Make_first_vecright_b(vec_right, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,3)
 dis_val_listo.append(3)
 Corr_val_listo.append(Corr_val) 

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_listo.append(dis_val)
  Corr_val_listo.append(Corr_val) 
##################################################################################


 print "\n"

######################b-ao#####################################################

 vec_left_1=Make_first_vecleft_b(vec_left, Ta1, Ta3,Tb1,Tb3,a,bp,c,d)
 vec_right_1=Make_first_vecright_a(vec_right, Ta1, Ta3,Tb1,Tb3,ap,b,c,d)

 #Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,1)
 #dis_val_list1.append(2)
 #Corr_val_list1.append(Corr_val)

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1o.append(dis_val)
  Corr_val_list1o.append(Corr_val) 
##################################################################################


 print "\n"

#######################c-do#####################################################

 vec_left_1=Make_first_vecleft_c(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)
 vec_right_1=Make_first_vecright_d(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,3)
 dis_val_list2o.append(3)
 Corr_val_list2o.append(Corr_val)

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2o.append(dis_val)
  Corr_val_list2o.append(Corr_val) 
###################################################################################

 print "\n"

#######################d-co#####################################################

 vec_left_1=Make_first_vecleft_d(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,c,dp)
 vec_right_1=Make_first_vecright_c(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,cp,d)

# Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,2)
# dis_val_list3.append(2)
# Corr_val_list3.append(Corr_val)

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_left_1=Make_midle_vecleft(vec_left_1, Ta1, Ta3,Tb1,Tb3,a,b,c,d)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_right_1, vec_left_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3o.append(dis_val)
  Corr_val_list3o.append(Corr_val) 
###################################################################################


 print dis_val_list,'\n,\n'
 print Corr_val_list,'\n,\n'
 print Corr_val_list1,'\n,\n'
 print Corr_val_list2,'\n,\n'
 print Corr_val_list3,'\n,\n'

 print dis_val_listo,'\n,\n'
 print Corr_val_listo,'\n,\n'
 print Corr_val_list1o,'\n,\n'
 print Corr_val_list2o,'\n,\n'
 print Corr_val_list3o,'\n,\n'


 Corr_val_list_ave=[ (sum(t)*(1.0/4.0)) for t in zip(Corr_val_list, Corr_val_list1, Corr_val_list2, Corr_val_list3)]
 print Corr_val_list_ave,'\n,\n'


 Corr_val_list_avo=[ (sum(t)*(1.0/4.0)) for t in zip(Corr_val_listo, Corr_val_list1o, Corr_val_list2o, Corr_val_list3o)]
 print Corr_val_list_avo,'\n,\n'


 Corr_val_list_final=Corr_val_list_ave+Corr_val_list_avo
 dis_val_list_final=dis_val_list+dis_val_listo

 print '\n,\n'
 print dis_val_list_final
 print Corr_val_list_final
 print '\n,\n'


 for i in xrange(len(dis_val_list_final)):
  fileCorr.write(str(dis_val_list_final[i])  + " " + str(Corr_val_list_final[i]) + "\n")
  fileCorr.flush()

def CorrelationV(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,distance_final,fileCorr,fileCorrLength):
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)


 if Model is "Heisenberg":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)(-1.0)*uni10.otimes(sy,sy)
  H=uni10.otimes(sz,iden)+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1=uni10.otimes(iden,sz)+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 if Model is "Heisenberg_Z2":
  #print d_phys
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  iden = matIden()
  szt=uni10.otimes(sz,iden)
  H.setRawElem(szt)
  szt1=uni10.otimes(iden,sz)
  H1.setRawElem(szt1)
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])
 if Model is "Heisenberg_U1":
  bdi = uni10.Bond(uni10.BD_IN, d_phys)
  bdo = uni10.Bond(uni10.BD_OUT, d_phys)
  HH = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H = uni10.UniTensor([bdi, bdi, bdo, bdo])
  H1 = uni10.UniTensor([bdi, bdi, bdo, bdo])
  sz = matSz()
  sx = matSx()
  sy = matSy()
  iden = matIden()
  HH_tem=uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy)
  H_tem=uni10.otimes(sz,iden)#+uni10.otimes(sx,iden)#+(-1.0)*uni10.otimes(sy,iden)
  H1_tem=uni10.otimes(iden,sz)#+uni10.otimes(iden,sx)#+(-1.0)*uni10.otimes(iden,sy)
  HH.setRawElem(HH_tem)
  H.setRawElem(H_tem)
  H1.setRawElem(H1_tem)
  Iden=copy.copy(H)
  Iden.identity()
  #print HH_tem, HH, H_tem, H, H1_tem, H1
  HH.setLabel([-10,-20,10,20])
  H.setLabel([-10,-20,10,20])
  H1.setLabel([-10,-20,10,20])
  Iden.setLabel([-10,-20,10,20])

 vec_down=make_down(c4,Ta3, Tb3,c3)
 vec_up=make_up(c1,Tb1, Ta1,c2)
 ap=make_ap_openindex(a_u)
 bp=make_ap_openindex(b_u)
 cp=make_ap_openindex(c_u)
 dp=make_ap_openindex(d_u)
 dis_val_list=[]
 Corr_val_list=[]
 dis_val_list1=[]
 Corr_val_list1=[]
 dis_val_list2=[]
 Corr_val_list2=[]
 dis_val_list3=[]
 Corr_val_list3=[]

# Corr_length=ED_right(c2, Ta2, Tb2, c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3,vec_right)
# print "Corr_length",Corr_length
 vec_up_copy=copy.copy(vec_up)
 Corr_length=ED_up(c1,Tb1, Ta1,c2, a, b, c, d, Tb2, Ta2, Ta4, Tb4,vec_up_copy)
 print "Corr_length",Corr_length
 fileCorrLength.write(str(Corr_length)  + "\n")
 fileCorrLength.flush()


 print "\n"
#######################a-a#####################################################
 vec_down_1=Make_first_vecdown_a(vec_down, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_a(vec_up, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
 dis_val_list.append(2)
 Corr_val_list.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list.append(dis_val)
  Corr_val_list.append(Corr_val) 
###################################################################################

 print "\n"
########################b-b#####################################################

 vec_down_1=Make_first_vecdown_b(vec_down, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_b(vec_up, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
 dis_val_list1.append(2)
 Corr_val_list1.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1.append(dis_val)
  Corr_val_list1.append(Corr_val) 
####################################################################################


 print "\n"
########################c-c#####################################################

 vec_down_1=Make_first_vecdown_c(vec_down, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_c(vec_up, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
 dis_val_list2.append(2)
 Corr_val_list2.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2.append(dis_val)
  Corr_val_list2.append(Corr_val) 
####################################################################################

 print "\n"
########################d-d#####################################################

 vec_down_1=Make_first_vecdown_d(vec_down, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_d(vec_up, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
 dis_val_list3.append(2)
 Corr_val_list3.append(Corr_val) 

 dis_val=2
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3.append(dis_val)
  Corr_val_list3.append(Corr_val) 
####################################################################################

#############################################################################################################


 dis_val_listo=[]
 Corr_val_listo=[]
 dis_val_list1o=[]
 Corr_val_list1o=[]
 dis_val_list2o=[]
 Corr_val_list2o=[]
 dis_val_list3o=[]
 Corr_val_list3o=[]

 print "\n"

#######################a-co#####################################################

 vec_down_1=Make_first_vecdown_a(vec_down, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_c(vec_up, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)

# Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
# dis_val_listo.append(2)
# Corr_val_listo.append(Corr_val) 

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_listo.append(dis_val)
  Corr_val_listo.append(Corr_val) 

###################################################################################


 print "\n"

#######################c-ao#####################################################
 vec_down_1=Make_first_vecdown_c(vec_down, a, b, cp, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_a(vec_up, ap, b, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,3)
 dis_val_list1o.append(3)
 Corr_val_list1o.append(Corr_val) 

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list1o.append(dis_val)
  Corr_val_list1o.append(Corr_val) 
###################################################################################


 print "\n"

########################b-do#####################################################

 vec_down_1=Make_first_vecdown_b(vec_down, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_d(vec_up, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)

# Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,2)
# dis_val_list2.append(2)
# Corr_val_list2.append(Corr_val) 

 dis_val=1
 for i in xrange(distance_final+1):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list2o.append(dis_val)
  Corr_val_list2o.append(Corr_val) 
####################################################################################

 print "\n"

########################d-bo#####################################################

 vec_down_1=Make_first_vecdown_d(vec_down, a, b, c, dp, Tb2, Ta2, Ta4, Tb4)
 vec_up_1=Make_first_vecup_b(vec_up, a, bp, c, d, Tb2, Ta2, Ta4, Tb4)

 Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,3)
 dis_val_list3o.append(3)
 Corr_val_list3o.append(Corr_val) 

 dis_val=3
 for i in xrange(distance_final):
  dis_val+=2
  vec_down_1=Make_midle_vecdown(vec_down_1,a,b,c,d,Tb2, Ta2, Ta4, Tb4)
  Corr_val=Corr_val_function(Iden,HH,H,H1,vec_down_1, vec_up_1,dis_val)
  print dis_val, Corr_val
  dis_val_list3o.append(dis_val)
  Corr_val_list3o.append(Corr_val) 
####################################################################################


 print dis_val_list,'\n,\n'
 print Corr_val_list,'\n,\n'
 print Corr_val_list1,'\n,\n'
 print Corr_val_list2,'\n,\n'
 print Corr_val_list3,'\n,\n'

 print dis_val_listo,'\n,\n'
 print Corr_val_listo,'\n,\n'
 print Corr_val_list1o,'\n,\n'
 print Corr_val_list2o,'\n,\n'
 print Corr_val_list3o,'\n,\n'


 Corr_val_list_ave=[ (sum(t)*(1.0/4.0)) for t in zip(Corr_val_list, Corr_val_list1, Corr_val_list2, Corr_val_list3)]
 print Corr_val_list_ave,'\n,\n'


 Corr_val_list_avo=[ (sum(t)*(1.0/4.0)) for t in zip(Corr_val_listo, Corr_val_list1o, Corr_val_list2o, Corr_val_list3o)]
 print Corr_val_list_avo,'\n,\n'


 Corr_val_list_final=Corr_val_list_ave+Corr_val_list_avo
 dis_val_list_final=dis_val_list+dis_val_listo

 print '\n,\n'
 print dis_val_list_final
 print Corr_val_list_final
 print '\n,\n'


 for i in xrange(len(dis_val_list_final)):
  fileCorr.write(str(dis_val_list_final[i])  + " " + str(Corr_val_list_final[i]) + "\n")
  fileCorr.flush()


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
 #print w,"\n,\n"#, v
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
 #print A
 for i in xrange(D):
  if (i == 0) or ( max_val < abs(A[i]) ):
   max_val = abs(A[i])
   index=i
 return max_val, index, A[index] 


##############################################################################################
def ED_right(c2, Ta2, Tb2, c3, a, b, c, d, Tb1, Ta1, Ta3, Tb3,Vec_uni):

 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D=",  D

 m=10
 W=2
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
  Lambda, index, Lambda_comp =find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'r0', Lambda, Lambda_comp
  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  Vec_FL=Q*eigvec
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


 E1L=copy.copy(E1)
 Vec_FL=Vec_FL*(1.00/Vec_FL.norm())
 m=10
 W=2
 num=0
 E1=0
 p=0
 Vec_F.randomize()
 Vec_F=Vec_F*(1.00/Vec_F.norm())
 while p  <  (W+1):
  #print "norm", p, Vec_F.norm(),Vec_F[0], Vec_F[1], Vec_F[2], Vec_F[3] 
  r=copy.copy(Vec_F)
 
  Vec_FL_trans=copy.copy(Vec_FL)
  Vec_FL_trans.transpose()
  dot_vec=Vec_FL_trans*r
  dot_val=dot_vec.trace()
  r=r+(-1.00*dot_val*Vec_FL)
  
  
  #r = r* (1.00/r.norm()) 
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()
  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3)

   Vec_FL_trans=copy.copy(Vec_FL)
   Vec_FL_trans.transpose()
   dot_vec=Vec_FL_trans*r
   dot_val=dot_vec.trace()
   r=r+(-1.00*dot_val*Vec_FL)



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
  Lambda, index, Lambda_comp=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'r1', Lambda, Lambda_comp
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

 Length=abs(E1/E1L)
 Length_val=-2.0*(1.00/math.log(Length))
 print "Length", Length,Length_val 
 return Length_val


def ED_up(c1,Tb1, Ta1,c2, a, b, c, d, Tb2, Ta2, Ta4, Tb4,Vec_uni):


 Vec_F=Vec_uni.getBlock()
 D=Vec_F.row()
 #print "D=",  D

 m=10
 W=2
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
  Lambda, index, Lambda_comp=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'u0', Lambda,Lambda_comp

  Q=make_Q(q_vec)
  Q.resize(D,m)
  Vec_F=Q*eigvec
  Vec_FL=Q*eigvec
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

 E1L=copy.copy(E1)
 Vec_FL=Vec_FL*(1.00/Vec_FL.norm())
 m=10
 W=2
 num=0
 E1=0
 p=0
 Vec_F.randomize()
 Vec_F=Vec_F*(1.00/Vec_F.norm())
 while p  <  (W+1):
  #print "norm", p, Vec_F.norm(),Vec_F[0], Vec_F[1], Vec_F[2], Vec_F[3] 
  r=copy.copy(Vec_F)
 
  Vec_FL_trans=copy.copy(Vec_FL)
  Vec_FL_trans.transpose()
  dot_vec=Vec_FL_trans*r
  dot_val=dot_vec.trace()
  r=r+(-1.00*dot_val*Vec_FL)
  
  
  #r = r* (1.00/r.norm()) 
  q_vec=[]
  q_vec.append(copy.copy(r))
  h=uni10.Matrix(m,m)
  h.set_zero()
  for j in xrange(m):
   vec_tem=copy.copy(q_vec[j])
   Vec_uni.putBlock(vec_tem)
   r=Multi_u(Vec_uni,a, b, c, d, Tb2, Ta2, Ta4, Tb4)

   Vec_FL_trans=copy.copy(Vec_FL)
   Vec_FL_trans.transpose()
   dot_vec=Vec_FL_trans*r
   dot_val=dot_vec.trace()
   r=r+(-1.00*dot_val*Vec_FL)



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
  Lambda, index,Lambda_comp=find_maxindex(D_eig[0])
  eigvec=return_vec(D_eig[1], index )
  print 'u1', Lambda,Lambda_comp
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



 Length=abs(E1/E1L)
 Length_val=-2.0*(1.00/math.log(Length))
 print "Length", Length,Length_val 
 return Length_val

def Multi_r(Vec_uni,a,b,c,d,Tb1,Ta1,Ta3,Tb3):
 CTM_1 = uni10.Network("Network1/Right.net")
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
 CTM_1 = uni10.Network("Network1/Up.net")
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









############################################################################################




















def  Corr_val_function(Iden,HH,H,H1,vec_right, vec_left,dis_val):

 vec_left.setLabel([10,23, 16, -16 , 9,-9, 2,-10])
 vec_right.setLabel([20,23,16,-16 , 9,-9,2,-20])
 Cor_norm=(vec_left*Iden)*vec_right
 Corr_val=(vec_left*HH)*vec_right
 Corr_val1=(vec_left*H)*vec_right
 Corr_val2=(vec_left*H1)*vec_right

 print  dis_val, Corr_val[0]/Cor_norm[0],Corr_val1[0]/Cor_norm[0],Corr_val2[0]/Cor_norm[0], Cor_norm[0] 

 print dis_val, (Corr_val[0]/Cor_norm[0]),  ((Corr_val1[0]*Corr_val2[0])/(Cor_norm[0]*Cor_norm[0]))  

 val=(Corr_val[0]/Cor_norm[0]) - ((Corr_val1[0]*Corr_val2[0])/(Cor_norm[0]*Cor_norm[0]))  

 return val


def  make_ap_openindex(a_u):
 a_u.setLabel([10,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.transpose()
 a_uc.setLabel([-3,-4,-10,-1,-2])
 result=a_uc*a_u
 result.permute([10,1,-1,2,-2,3,-3,4,-4,-10], 5)
 return result


def make_vleft(Tb4,Ta4,c1,c4):

 Tb4.setLabel([-5,3,-3,4])
 Ta4.setLabel([1,2,-2,-5])
 c1.setLabel([4,6])
 c4.setLabel([1,5])
 vec_left=(Tb4*Ta4)*(c1*c4)
 vec_left.permute([5,2,-2,3,-3,6],0)
 return vec_left

def  make_vright(Ta2,Tb2,c2,c3):

 Ta2.setLabel([-5,3,-3,4])
 Tb2.setLabel([1,2,-2,-5])
 c2.setLabel([6,4])
 c3.setLabel([5,1])
 vec_right=(Ta2*Tb2)*(c2*c3)
 vec_right.permute([5,2,-2,3,-3,6],6)
 return vec_right


def  Make_first_vecleft_a(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d):

 CTM_1 = uni10.Network("Network1/LeftCorra.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_a(vec_right, Ta1, Ta3,Tb1,Tb3,ap,b,c,d):

 CTM_1 = uni10.Network("Network1/RightCorra.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec

def  Make_first_vecleft_b(vec_left, Ta1, Ta3,Tb1,Tb3,a,bp,c,d):

 CTM_1 = uni10.Network("Network1/LeftCorrb.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_b(vec_right, Ta1, Ta3,Tb1,Tb3,a,bp,c,d):

 CTM_1 = uni10.Network("Network1/RightCorrb.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec


def  Make_first_vecleft_c(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,cp,d):

 CTM_1 = uni10.Network("Network1/LeftCorrc.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_c(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,cp,d):

 CTM_1 = uni10.Network("Network1/RightCorrc.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec


def  Make_first_vecleft_d(vec_left, Ta1, Ta3,Tb1,Tb3,a,b,c,dp):

 CTM_1 = uni10.Network("Network1/LeftCorrd.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec


def Make_first_vecright_d(vec_right, Ta1, Ta3,Tb1,Tb3,a,b,c,dp):

 CTM_1 = uni10.Network("Network1/RightCorrd.net")
 CTM_1.putTensor('vec_right',vec_right)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 vec.permute([10,21,14,-14 , 7,-7,0,-10],6)
 vec=max_ten(vec)
 return vec

def  Make_midle_vecleft(vec_left, Ta1, Ta3,Tb1,Tb3,ap,b,c,d):

 CTM_1 = uni10.Network("Network1/LeftCorr1.net")
 CTM_1.putTensor('vec_left',vec_left)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta3',Ta3)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb3',Tb3)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec.permute([10,23, 16, -16 , 9,-9, 2,-10],0)
 vec=max_ten(vec)
 return vec



############################################################################################

def  make_ap_openindex(a_u):
 a_u.setLabel([10,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.transpose()
 a_uc.setLabel([-3,-4,-10,-1,-2])
 result=a_uc*a_u
 result.permute([10,1,-1,2,-2,3,-3,4,-4,-10], 5)
 return result


def make_down(c4,Ta3, Tb3,c3):
 Tb3.setLabel([-5,3,-3,4])
 Ta3.setLabel([1,2,-2,-5])
 c3.setLabel([4,6])
 c4.setLabel([5,1])
 vec_down=(Tb3*Ta3)*(c3*c4)
 vec_down.permute([5, 2, -2 , 3, -3 , 6],0)
 return vec_down

def  make_up(c1,Tb1, Ta1,c2):
 Ta1.setLabel([-5,3,-3,4])
 Tb1.setLabel([1,2,-2,-5])
 c2.setLabel([4,6])
 c1.setLabel([5,1])
 vec_up=(Ta1*Tb1)*(c1*c2)
 vec_up.permute([5,2,-2,3,-3,6],6)
 return vec_up




def  Make_first_vecdown_a(vec_down, ap, b, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/DownCorra.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec




def Make_first_vecup_a(vec_up, ap, b, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/UpCorra.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',ap)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 return vec


def  Make_midle_vecdown(vec_down, a,b,c,d,Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/DownCorr1.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 vec=max_ten(vec)
 return vec



def  Make_first_vecdown_b(vec_down, a, bp, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/DownCorrb.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec




def Make_first_vecup_b(vec_up, a, bp, c, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/UpCorrb.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',bp)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 return vec



def  Make_first_vecdown_c(vec_down, a, b, cp, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/DownCorrc.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec

def Make_first_vecup_c(vec_up, a, b, cp, d, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/UpCorrc.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',cp)
 CTM_1.putTensor('d',d)
 vec=CTM_1.launch()
 return vec




def  Make_first_vecdown_d(vec_down, a, b, c, dp, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/DownCorrd.net")
 CTM_1.putTensor('vec_down',vec_down)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 #print CTM_1
 #vec.permute([10,3, 4, -4 , 5, -5 , 6,-10],8)
 return vec

def Make_first_vecup_d(vec_up, a, b, c, dp, Tb2, Ta2, Ta4, Tb4):

 CTM_1 = uni10.Network("Network1/UpCorrd.net")
 CTM_1.putTensor('vec_up',vec_up)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Ta4',Ta4)
 CTM_1.putTensor('Tb2',Tb2)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 CTM_1.putTensor('c',c)
 CTM_1.putTensor('d',dp)
 vec=CTM_1.launch()
 return vec















#################################################################################################



















