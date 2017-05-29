import pyUni10 as uni10
import copy
import time
import basic
import numpy as np
#import scipy as sp
#import line_profiler

#@profile
def Var_cab(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0,N_grad, Opt_method,plist,MPO_list,Inv_method,N_svd,N_env,method,check_step):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 
 
 if Corner_method is 'CTM':
#  c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1=basic.make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1)
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D,H0,d_phys)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'three',N_env)
  #basic.Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'three',N_env)

 Env=basic.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a,b,c,d=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u)

 a_up=copy.copy(a_u) 
 b_up=copy.copy(b_u) 
 c_up=copy.copy(c_u) 
 d_up=copy.copy(d_u) 


 if method is "SVD_mpo":
  Do_optimization_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, MPO_list, c_u, a_u, b_u, plist, Inv_method, N_grad,N_svd,Gauge,check_step)
  equall_dis_plist(plist)
  c_up,a_up,b_up=recover(c_u, a_u, b_u, plist, MPO_list)
  normal_plist(plist)

 if method is "SVD":
  a_up, b_up, c_up, d_up=Do_optimization_svd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Inv_method,N_svd,Gauge,check_step)

 if method is "Grad":
  a_up, b_up, c_up, d_up=Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Gauge)


# Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
# print "DisFFinal", Dis_val

 for q in xrange(2):
  c_up, d_up=equall_dis_H(c_up, d_up)
  c_up, a_up=equall_dis_V(c_up, a_up)
  a_up, b_up=equall_dis_H(a_up, b_up)
  c_up, d_up=equall_dis_H(c_up, d_up)
  d_up, b_up=equall_dis_V(d_up, b_up)
  a_up, b_up=equall_dis_H(a_up, b_up)
  d_up, b_up=equall_dis_V(d_up, b_up)
  c_up, a_up=equall_dis_V(c_up, a_up)
  Maxa=MaxAbs(a_up)
  Maxb=MaxAbs(b_up)
  Maxc=MaxAbs(c_up)
  Maxd=MaxAbs(d_up)
  print Maxa, Maxb, Maxc, Maxd

 
# Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
# print "DisFFinal", Dis_val
 
 a_up=max_ten(a_up)
 b_up=max_ten(b_up) 
 c_up=max_ten(c_up)
 d_up=max_ten(d_up)

 ap=basic.make_ab(a_up)
 bp=basic.make_ab(b_up)
 cp=basic.make_ab(c_up)
 dp=basic.make_ab(d_up)


 return a_up, b_up, c_up, d_up, ap, bp, cp, dp

##@profile
def Var_abd(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0,N_grad, Opt_method,plist,MPO_list,Inv_method,N_svd,N_env,method,check_step):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 
 if Corner_method is 'CTM':
#  c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1=basic.make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1)
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D,H0,d_phys)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'three1',N_env)
  #basic.Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a_u,b_u,c_u,d_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'three1',N_env)

 Env=basic.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u)

 a_up=copy.copy(a_u) 
 b_up=copy.copy(b_u) 
 c_up=copy.copy(c_u) 
 d_up=copy.copy(d_u) 
 

 
 #randomize_plist(plist)
 if method is "SVD_mpo":
  Do_optimization_MPO1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,N_grad,N_svd,Gauge,check_step)
  equall_dis_plist1(plist)
  a_up,b_up,d_up=recover1( a_u, b_u, d_u, plist, MPO_list)
  normal_plist1(plist)
 

 if method is "SVD":
  a_up, b_up, c_up, d_up=Do_optimization_svd1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Inv_method,N_svd,Gauge,check_step)

 if method is "Grad":
  a_up, b_up, c_up, d_up=Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Gauge)

# Dis_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
# print "DisF", Dis_val

 for q in xrange(2):
  a_up, b_up=equall_dis_H(a_up, b_up)
  c_up, d_up=equall_dis_H(c_up, d_up)
  c_up, a_up=equall_dis_V(c_up, a_up)
  d_up, b_up=equall_dis_V(d_up, b_up)
  c_up, d_up=equall_dis_H(c_up, d_up)
  d_up, b_up=equall_dis_V(d_up, b_up)
  a_up, b_up=equall_dis_H(a_up, b_up)
  c_up, a_up=equall_dis_V(c_up, a_up)
  Maxa=MaxAbs(a_up)
  Maxb=MaxAbs(b_up)
  Maxc=MaxAbs(c_up)
  Maxd=MaxAbs(d_up)
  print Maxa, Maxb, Maxc, Maxd
 
# Dis_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
# print "DisF", Dis_val
 

 a_up=max_ten(a_up)
 b_up=max_ten(b_up) 
 c_up=max_ten(c_up)
 d_up=max_ten(d_up)

 ap=basic.make_ab(a_up)
 bp=basic.make_ab(b_up)
 cp=basic.make_ab(c_up)
 dp=basic.make_ab(d_up)


 return a_up, b_up,c_up,d_up, ap,bp,cp,dp



def produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u):

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,-2,3])
 E1=c1*Tb1
 E1.permute([0,2,-2,3],2)
 E1.setLabel([1,2,-2,3])
 E1.permute([1,2,-2,3],3)
 
 c2.setLabel([1,0])
 Ta1.setLabel([3,2,-2,1])
 E2=c2*Ta1
 E2.permute([0,2,-2,3],2)
 E2.setLabel([5,4,-4,3])
 E2.permute([3,4,-4,5],4)
 
 E3=copy.copy(Ta2)
 E3.setLabel([7,6,-6,5])
 E3.permute([5,6,-6,7],3)
 

 c3.setLabel([8,7])
 Tb2.setLabel([7,15,-15,22])
 E4=(c3*Tb2)
 E4.permute([8,15,-15,22],2)
 E4.setLabel([9,8,-8,7])
 E4.permute([7,8,-8,9],3)

 E5=copy.copy(Tb3)
 E5.setLabel([11,10,-10,9])
 E5.permute([9,10,-10,11],2)
 
 c4.setLabel([11,10])
 Ta3.setLabel([10,12,-12,13])
 E6=c4*Ta3
 E6.permute([11,12,-12,13],1)
 E6.setLabel([13,12,-12,11])
 E6.permute([11,12,-12,13],1)

 E7=copy.copy(Ta4)
 E7.setLabel([13,14,-14,15])
 E7.permute([13,14,-14,15],1)

 E8=copy.copy(Tb4)
 E8.setLabel([15,16,-16,1])
 E8.permute([15,16,-16,1],1)


 b.setLabel([18,-18,20,-20,6,-6,4,-4])
 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])
 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 Norm=(((((E1*E8)*(a))*((E7*E6)*(c))))*(((E2*E3)*(b))))*((E4*E5)*d)
 if Norm[0] < 0: E1=-1.0*E1;

 while Norm[0]<1.0e-4: 
  E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=checking_norm(E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d)
  b.setLabel([18,-18,20,-20,6,-6,4,-4])
  c.setLabel([14,-14,12,-12,19,-19,17,-17])
  a.setLabel([16,-16,17,-17,18,-18,2,-2])
  d.setLabel([19,-19,10,-10,8,-8,20,-20])
  Norm=(((((E1*E8)*(a))*((E7*E6)*(c))))*(((E2*E3)*(b))))*((E4*E5)*d)
  #print Norm[0]


 while Norm[0]>1.0e+5: 
  E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=checking_norm(E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d)
  b.setLabel([18,-18,20,-20,6,-6,4,-4])
  c.setLabel([14,-14,12,-12,19,-19,17,-17])
  a.setLabel([16,-16,17,-17,18,-18,2,-2])
  d.setLabel([19,-19,10,-10,8,-8,20,-20])
  Norm=(((((E1*E8)*(a))*((E7*E6)*(c))))*(((E2*E3)*(b))))*((E4*E5)*d)
  #print Norm[0]

 return E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d 


def checking_norm(E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d):
  
 b.setLabel([18,-18,20,-20,6,-6,4,-4])
 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])
 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 Norm=(((((E1*E8)*(a))*((E7*E6)*(c))))*(((E2*E3)*(b))))*((E4*E5)*d)
 if Norm[0] < 0: E1=-1.0*E1;
 if Norm[0] < 1.0e-4:
  #print "Norm[0] < 1.0e+1 ", Norm[0]
  E1*=1.5
  E8*=1.5
  E2*=1.5
  E3*=1.5
  E4*=1.5
  E6*=1.5
  E7*=1.5
#  a_u*=5.0
#  b_u*=5.0
#  c_u*=5.0
#  d_u*=5.0
#  a=basic.make_ab(a_u)
#  b=basic.make_ab(b_u)
#  c=basic.make_ab(c_u)
#  d=basic.make_ab(d_u)
 elif Norm[0] > 1.e+5:
  #print "Norm[0] > 1.e+8 ", Norm[0]
  E1*=0.5
  E8*=0.5
  E2*=0.5
  E3*=0.5
  E4*=0.5
  E6*=0.5
  E7*=0.5
#  a_u*=0.20
#  b_u*=0.20
#  c_u*=0.20
#  d_u*=0.20
#  a=basic.make_ab(a_u)
#  b=basic.make_ab(b_u)
#  c=basic.make_ab(c_u)
#  d=basic.make_ab(d_u)
 return E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d 


  
 
def energy_cab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, U):


 U.setLabel([-54,-55,-56,54,55,56])

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])

 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,55,-16,-17])
 b_d.setLabel([-6,-4,56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])

 Val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))

###########################################################################################

 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])


 E_val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d)))*U)*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))

 return  E_val[0]/Val[0]


def energy_abd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, U):

 U.setLabel([-55,-56,-57,55,56,57])

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])

 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,55,-16,-17])
 b_d.setLabel([-6,-4,56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])

 Val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))

###########################################################################################

 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])

 E_val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E2*E3)*(b_u*b_d)*U)))*((E4*E5)*(d_u*d_d))

 return  E_val[0]/Val[0]


 
###@profile
def  Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U):

 U.setLabel([-54,-55,-56,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,54,55,56])
 U_2=U*H

 U.setLabel([-54,-55,-56,54,55,56])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-54,-55,-56,54,55,56])

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])

 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])

# Val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d)))*U_2)*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))
# Val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d)))*Iden)*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))

# print Val[0], Val
###########################################################################################

 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
 
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])
 Val1=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp))
 #print Val1[0]

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


 Val2=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d)))*U)*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_up*d_d))
# Val3=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u*d_dp))
 #print Val2[0]
 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=Val1[0]-Val2[0]-Val3[0]
 val_f=Val1[0]-2.0*Val2[0]#-Val3[0]

 if Val1[0]<0: 
  print "Dis, N<0"
  val_f=-Val1[0]-2.0*Val2[0]

 return  val_f

##@profile
def Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U, N_grad, Opt_method, Gauge):


  a_um=copy.copy(a_up) 
  b_um=copy.copy(b_up) 
  c_um=copy.copy(c_up) 
  d_um=copy.copy(d_up) 

  time_val=0
  Es=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
  Ef=0
  E2_val=0
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=1.e+18
  count=0
  D_list=[0]*4
  H_list=[0]*4
  H_a=0; H_b=0;H_c=0;H_d=0;
  for i in xrange(N_grad):
   count+=1
   t0=time.time()

   E1_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   Ef=E1_val
   print 'E1=', E1_val, abs((E_previous-E1_val)/E1_val), count, time_val


   D_a,D_b,D_c,D_d=basic.Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Gauge)
   D_a=(-1.0)*D_a
   D_b=(-1.0)*D_b
   D_c=(-1.0)*D_c
   D_d=(-1.0)*D_d

   if i is 0:
    H_a=D_a
    H_b=D_b
    H_c=D_c
    H_d=D_d
   else:
    Z_a=D_a+(-1.0)*D_list[0]
    Z_b=D_b+(-1.0)*D_list[1]
    Z_c=D_c+(-1.0)*D_list[2]
    Z_d=D_d+(-1.0)*D_list[3]
    A=Z_a*D_a
    B=Z_b*D_b
    C=Z_c*D_c
    D=Z_d*D_d
    A1=D_list[0]*D_list[0]
    A2=D_list[1]*D_list[1]
    A3=D_list[2]*D_list[2]
    A4=D_list[3]*D_list[3]
    Gamma_grad=(A[0]+B[0]+C[0]+D[0]) / (A1[0]+A2[0]+A3[0]+A4[0])
    if Opt_method is 'ST':Gamma_grad=0;
    H_a=D_a+(Gamma_grad)*H_list[0]
    H_b=D_b+(Gamma_grad)*H_list[1]
    H_c=D_c+(Gamma_grad)*H_list[2]
    H_d=D_d+(Gamma_grad)*H_list[3]

#    A=D_a*D_list[0]
#    B=D_b*D_list[1]
#    C=D_c*D_list[2]
#    D=D_d*D_list[3]
#    check=A[0]+B[0]+C[0]+D[0] 
#    print "check", check 

   D_list[0]=copy.copy(D_a)
   D_list[1]=copy.copy(D_b)
   D_list[2]=copy.copy(D_c)
   D_list[3]=copy.copy(D_d)

   H_list[0]=copy.copy(H_a)
   H_list[1]=copy.copy(H_b)
   H_list[2]=copy.copy(H_c)
   H_list[3]=copy.copy(H_d)



   A=D_a*H_a
   B=D_b*H_b
   C=D_c*H_c
   D=D_d*H_d
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-12:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break

   #print E1_val>E_previous  
   if (E1_val>E_previous):
    print "break, not satisfied", E1_val, E_previous
    a_up=copy.copy(a_um) 
    b_up=copy.copy(b_um) 
    c_up=copy.copy(c_um) 
    d_up=copy.copy(d_um)
    break
   else:
    a_um=copy.copy(a_up) 
    b_um=copy.copy(b_up) 
    c_um=copy.copy(c_up) 
    d_um=copy.copy(d_up) 

   E_previous=E1_val
   
   if abs(Norm_Z) < 1.0e-11:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%15)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   #print "Gamma", Gamma
   while Break_loop is 1:
    count+=1
    a_ut=a_up+(2.00)*Gamma*H_a
    b_ut=b_up+(2.00)*Gamma*H_b
    c_ut=c_up+(2.00)*Gamma*H_c
    d_ut=d_up+(2.00)*Gamma*H_d
    E2_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_ut, b_ut, c_ut, d_ut, U)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+12 or  abs(Gamma)>1.0e+12 :
     print "break1", E1_val, abs((0.5)*Norm_Z*Gamma), E2_val, Gamma
     Gamma=1
     break
    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0

   Break_loop=1
   while Break_loop is 1:
    count+=1
    a_ut=a_up+(1.00)*Gamma*H_a
    b_ut=b_up+(1.00)*Gamma*H_b
    c_ut=c_up+(1.00)*Gamma*H_c
    d_ut=d_up+(1.00)*Gamma*H_d
    E2_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_ut, b_ut, c_ut, d_ut, U)
    if abs((0.5)*Norm_Z*Gamma) <1.0e-16 or  (abs((E1_val-E2_val)/E2_val))<1.0e-16 or abs(Gamma)<1.0e-16 :
     print "break2", E1_val, E2_val, Gamma, abs((0.5)*Norm_Z*Gamma), (abs((E1_val-E2_val)/E2_val))
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0


   a_up=a_up+(1.00)*Gamma*H_a
   b_up=b_up+(1.00)*Gamma*H_b
   c_up=c_up+(1.00)*Gamma*H_c
   d_up=d_up+(1.00)*Gamma*H_d
   time_val=time.time() - t0

  return a_up, b_up, c_up, d_up


def Do_optimization_Grad_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,N_grad):
  Es=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u, plist)
  Ef=0
  E2_val=0
  plist_first=[plist[i] for i in xrange(len(plist))]
  plist_tem=[0]*len(plist)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=1.e+18
  count=0
  for i in xrange(N_grad):
   count+=1
   E1_val=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   Ef=E1_val
   #Store_plist(plist)

   print 'E=', E1_val, count
   D_rf=basic.Obtain_grad_four_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   D_r=[ (-1.00)*D_rf[i] for i in xrange(len(D_rf)) ]

   A=D_r[0]*D_r[0]
   B=D_r[1]*D_r[1]
   C=D_r[2]*D_r[2]
   D=D_r[3]*D_r[3]
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-12:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-13:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break


   #print E1_val,  E_previous  
   if (E1_val>E_previous):
    print "break, not satisfied", E1_val, E_previous
    for i in xrange(len(plist_first)):
      plist[i]=plist_first[i]
    break
   else:
    plist_first=[plist[i] for i in xrange(len(plist))]
  
   E_previous=E1_val
   
   if Norm_Z < 1.0e-12:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%10)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   while Break_loop is 1:
    count+=1
    plist_tem[0]=plist[0]+(2.00)*Gamma*D_r[0]
    plist_tem[1]=plist[1]+(2.00)*Gamma*D_r[1]
    plist_tem[2]=plist[2]+(2.00)*Gamma*D_r[2]
    plist_tem[3]=plist[3]+(2.00)*Gamma*D_r[3]
    E2_val=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist_tem)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 or abs(Gamma)<1.0e-17 :
     print "break", E1_val, E2_val, Gamma
     Gamma=0.10
     break

    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
   
   Break_loop=1
   while Break_loop is 1:
    count+=1
    plist_tem[0]=plist[0]+(1.00)*Gamma*D_r[0]
    plist_tem[1]=plist[1]+(1.00)*Gamma*D_r[1]
    plist_tem[2]=plist[2]+(1.00)*Gamma*D_r[2]
    plist_tem[3]=plist[3]+(1.00)*Gamma*D_r[3]
    E2_val=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist_tem)
    if abs((0.5)*Norm_Z*Gamma) <1.0e-15 or  abs(E1_val-E2_val)<1.0e-15 or abs(Gamma)<1.0e-15 :
     print "break", E1_val, E2_val, Gamma
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0


   plist[0]=plist[0]+(1.00)*Gamma*D_r[0]
   plist[1]=plist[1]+(1.00)*Gamma*D_r[1]
   plist[2]=plist[2]+(1.00)*Gamma*D_r[2]
   plist[3]=plist[3]+(1.00)*Gamma*D_r[3]


def Do_optimization_Grad_MPO1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,N_grad):
  Es=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u, plist)
  Ef=0
  E2_val=0
  plist_first=[plist[i] for i in xrange(len(plist))]
  plist_tem=[0]*len(plist)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=1.e+18
  count=0
  for i in xrange(N_grad):
   count+=1
   E1_val=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   Ef=E1_val
   #Store_plist(plist)

   print 'E=', E1_val, count
   D_rf=basic.Obtain_grad_four_MPO1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   D_r=[ (-1.00)*D_rf[i] for i in xrange(len(D_rf)) ]

   A=D_r[0]*D_r[0]
   B=D_r[1]*D_r[1]
   C=D_r[2]*D_r[2]
   D=D_r[3]*D_r[3]
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-12:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-13:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break


   #print E1_val,  E_previous  
   if (E1_val>=E_previous):
    print "break, not satisfied", E1_val, E_previous
    for i in xrange(len(plist_first)):
      plist[i]=plist_first[i]
    break
   else:
    plist_first=[plist[i] for i in xrange(len(plist))]
  
   E_previous=E1_val
   
   if Norm_Z < 1.0e-12:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%10)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   while Break_loop is 1:
    count+=1
    plist_tem[0]=plist[0]+(2.00)*Gamma*D_r[0]
    plist_tem[1]=plist[1]+(2.00)*Gamma*D_r[1]
    plist_tem[2]=plist[2]+(2.00)*Gamma*D_r[2]
    plist_tem[3]=plist[3]+(2.00)*Gamma*D_r[3]
    E2_val=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist_tem)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 or abs(Gamma)<1.0e-17:
     print "break1", E1_val, E2_val, Gamma
     Gamma=0.1
     break

    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
   
   Break_loop=1
   while Break_loop is 1:
    count+=1
    plist_tem[0]=plist[0]+(1.00)*Gamma*D_r[0]
    plist_tem[1]=plist[1]+(1.00)*Gamma*D_r[1]
    plist_tem[2]=plist[2]+(1.00)*Gamma*D_r[2]
    plist_tem[3]=plist[3]+(1.00)*Gamma*D_r[3]
    E2_val=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist_tem)
    if abs((0.5)*Norm_Z*Gamma) <1.0e-15 or  abs(E1_val-E2_val)<1.0e-15 or abs(Gamma)<1.0e-15 :
     print "break2", E1_val, E2_val, Gamma
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0


   plist[0]=plist[0]+(1.00)*Gamma*D_r[0]
   plist[1]=plist[1]+(1.00)*Gamma*D_r[1]
   plist[2]=plist[2]+(1.00)*Gamma*D_r[2]
   plist[3]=plist[3]+(1.00)*Gamma*D_r[3]









###@profile
def  Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U):

 U.setLabel([-55,-56,-57,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,55,56,57])
 U_2=U*H

 U.setLabel([-55,-56,-57,55,56,57])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-55,-56,-57,55,56,57])

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])

 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])

 #Val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d)))*U_2)*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))

###########################################################################################

 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
 
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])
 Val1=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp))

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 Val2=(((((E1*E8)*(a_up*a_d))*(((E2*E3)*(b_up*b_d)))*U)*((E4*E5)*(d_up*d_d)))* ((E7*E6)*(c_up*c_d)))
 #Val3=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u*d_dp))
 #print Val, Val1, Val2, Val3
 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=Val1[0]-Val2[0]-Val3[0]
 val_f=Val1[0]-2.0*Val2[0]#-Val3[0]

 if Val1[0]<0: 
  print "Dis,N<0"
  val_f=-Val1[0]-2.0*Val2[0]#-Val3[0]
 
 #print Val2[0], Val3[0]
 return  val_f
##@profile
def Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Gauge):

  time_val=0
  a_um=copy.copy(a_up) 
  b_um=copy.copy(b_up) 
  c_um=copy.copy(c_up) 
  d_um=copy.copy(d_up) 

  Es=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
  Ef=0
  E2_val=0
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=0
  count=0
  D_list=[0]*4
  H_list=[0]*4
  H_a=0; H_b=0;H_c=0;H_d=0;
  for i in xrange(N_grad):
   t0=time.time()
   count+=1
   E1_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   Ef=E1_val
   print 'E2=', E1_val, abs((E_previous-E1_val)/E1_val), count, time_val
   
   D_a,D_b,D_c,D_d=basic.Obtain_grad_four1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Gauge)
   D_a=(-1.0)*D_a
   D_b=(-1.0)*D_b
   D_c=(-1.0)*D_c
   D_d=(-1.0)*D_d


   if i is 0:
    H_a=D_a
    H_b=D_b
    H_c=D_c
    H_d=D_d
   else:
    Z_a=D_a+(-1.0)*D_list[0]
    Z_b=D_b+(-1.0)*D_list[1]
    Z_c=D_c+(-1.0)*D_list[2]
    Z_d=D_d+(-1.0)*D_list[3]
    A=Z_a*D_a
    B=Z_b*D_b
    C=Z_c*D_c
    D=Z_d*D_d
    A1=D_list[0]*D_list[0]
    A2=D_list[1]*D_list[1]
    A3=D_list[2]*D_list[2]
    A4=D_list[3]*D_list[3]
    Gamma_grad=(A[0]+B[0]+C[0]+D[0]) / (A1[0]+A2[0]+A3[0]+A4[0])
    if Opt_method is 'ST':Gamma_grad=0;
    H_a=D_a+(Gamma_grad)*H_list[0]
    H_b=D_b+(Gamma_grad)*H_list[1]
    H_c=D_c+(Gamma_grad)*H_list[2]
    H_d=D_d+(Gamma_grad)*H_list[3]

#    A=D_a*D_list[0]
#    B=D_b*D_list[1]
#    C=D_c*D_list[2]
#    D=D_d*D_list[3]
#    check=A[0]+B[0]+C[0]+D[0] 
#    print "check", check 

   D_list[0]=copy.copy(D_a)
   D_list[1]=copy.copy(D_b)
   D_list[2]=copy.copy(D_c)
   D_list[3]=copy.copy(D_d)

   H_list[0]=copy.copy(H_a)
   H_list[1]=copy.copy(H_b)
   H_list[2]=copy.copy(H_c)
   H_list[3]=copy.copy(H_d)



   A=D_a*H_a
   B=D_b*H_b
   C=D_c*H_c
   D=D_d*H_d
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-12:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break


   #print E1_val,  E_previous  
   if (E1_val>E_previous):
    print "break, not satisfied", E1_val, E_previous
    a_up=copy.copy(a_um) 
    b_up=copy.copy(b_um) 
    c_up=copy.copy(c_um) 
    d_up=copy.copy(d_um)
    break
   else:
    a_um=copy.copy(a_up) 
    b_um=copy.copy(b_up) 
    c_um=copy.copy(c_up) 
    d_um=copy.copy(d_up) 

      
   E_previous=E1_val
   
   if abs(Norm_Z) < 1.0e-11:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%15)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   #print "Gamma", Gamma
   while Break_loop is 1:
    count+=1
    a_ut=a_up+(2.00)*Gamma*H_a
    b_ut=b_up+(2.00)*Gamma*H_b
    c_ut=c_up+(2.00)*Gamma*H_c
    d_ut=d_up+(2.00)*Gamma*H_d
    E2_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_ut, b_ut, c_ut, d_ut, U)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 :
     print "break1", E1_val, abs((0.5)*Norm_Z*Gamma), E2_val, Gamma
     Gamma=1.0
     break
    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0

   Break_loop=1
   while Break_loop is 1:
    count+=1
    a_ut=a_up+(1.00)*Gamma*H_a
    b_ut=b_up+(1.00)*Gamma*H_b
    c_ut=c_up+(1.00)*Gamma*H_c
    d_ut=d_up+(1.00)*Gamma*H_d
    E2_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_ut, b_ut, c_ut, d_ut, U)
    if abs((0.5)*Norm_Z*Gamma) <1.0e-16 or  (abs((E1_val-E2_val)/E2_val))<1.0e-16 or abs(Gamma)<1.0e-16 :
     print "break2", E1_val, E2_val, Gamma, abs((0.5)*Norm_Z*Gamma), (abs((E1_val-E2_val)/E2_val))
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0


   a_up=a_up+(1.00)*Gamma*H_a
   b_up=b_up+(1.00)*Gamma*H_b
   c_up=c_up+(1.00)*Gamma*H_c
   d_up=d_up+(1.00)*Gamma*H_d
   time_val=time.time() - t0

  return a_up, b_up, c_up, d_up


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





 
def svd_parity(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
    bd3=uni10.Bond(uni10.BD_IN,theta.bond(6).Qlist())
    bd4=uni10.Bond(uni10.BD_IN,theta.bond(7).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
    LA=uni10.UniTensor([bd1,bd2,bd3,bd4,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
    GB=uni10.UniTensor([bd1,bd2,bd3,bd4,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

#    print LA
    return GA, LA, GB

def svd_parity4(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(6).Qlist())
    bd3=uni10.Bond(uni10.BD_IN,theta.bond(7).Qlist())
    bd4=uni10.Bond(uni10.BD_IN,theta.bond(8).Qlist())
    bd5=uni10.Bond(uni10.BD_IN,theta.bond(9).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7),theta.bond(8),theta.bond(9)])
    LA=uni10.UniTensor([bd1,bd2,bd3,bd4,bd5,theta.bond(5),theta.bond(6),theta.bond(7),theta.bond(8),theta.bond(9)])
    GB=uni10.UniTensor([bd1,bd2,bd3,bd4,bd5,theta.bond(5),theta.bond(6),theta.bond(7),theta.bond(8),theta.bond(9)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

#    print LA
    return GA, LA, GB


def svd_parity3(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3)])
    LA=uni10.UniTensor([bd1,theta.bond(3)])
    GB=uni10.UniTensor([bd1,theta.bond(3)])
    svds = {}
    blk_qnums = theta.blockQnum()
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

#    print LA
    return GA, LA, GB

def svd_parity1(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])
    GB=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])
    svds = {}
    blk_qnums = theta.blockQnum()
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

#    print LA
    return GA, LA, GB



###@profile
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
      invL2[i] = 0 if ((invLt[i].real) < 1.0e-9) else (1.00 / (invLt[i].real))

  invLanda2.putBlock(qnum,invL2)
 return invLanda2
 
def equall_dis_H(a_u, b_u):

 A=copy.copy(b_u)
 A.setLabel([54,18,20,6,4])
 A.permute([18,20,6,4,54],1)
 
 l,qq=lq_parity1(A) 
 
 l.setLabel([18,80])
 qq.setLabel([80,20,6,4,54])

 A=copy.copy(a_u)
 A.setLabel([55,16,17,18,2])
 A.permute([2,16,17,55,18],4)
 
 q,r=qr_parity1(A) 
  
 q.setLabel([2,16,17,55,81])
 r.setLabel([81,18])

 Teta=l*r
 Teta.permute([81,80],1)
 U,s,V=svd_parity2(Teta)

 U.setLabel([81,18])
 s.setLabel([18,-18])
 V.setLabel([-18,80])

 s=Sqrt(s)
 
 U=U*s
 V=s*V

 U.permute([81,-18],1)
 U.setLabel([81,18])
 V.permute([18,80],1)
 
 b_up=qq*V
 a_up=q*U
 b_up.permute([54,18,20,6,4],3)
 a_up.permute([55,16,17,18,2],3)
#################################################################################
 return  a_up, b_up 

def equall_dis_V(c_u, a_u):

 A=copy.copy(a_u)
 A.setLabel([55,16,17,18,2])
 A.permute([17,55,16,18,2],1)
 
 l,qq=lq_parity1(A) 
 
 l.setLabel([17,80])
 qq.setLabel([80,55,16,18,2])

 A=copy.copy(c_u)
 A.setLabel([53,14,12,19,17])
 A.permute([53,14,12,19,17],4)
 q,r=qr_parity1(A) 
 
 
 q.setLabel([53,14,12,19,81])
 r.setLabel([81,17])

 Teta=l*r
 Teta.permute([81,80],1)
 U,s,V=svd_parity2(Teta)

 U.setLabel([81,17])
 s.setLabel([17,-17])
 V.setLabel([-17,80])

 s=Sqrt(s)
 
 U=U*s
 V=s*V

 U.permute([81,-17],1)
 U.setLabel([81,17])
 V.permute([17,80],1)
 
 a_up=qq*V
 c_up=q*U
 c_up.permute([53,14,12,19,17],3)
 a_up.permute([55,16,17,18,2],3)

 return c_up, a_up


def equall_dis1(a_up, c_up, d_up):

 A=copy.copy(d_up)
 A.setLabel([54,19,10,8,20])
 A.permute([19,10,8,20,54],1)
 
 l,qq=lq_parity1(A) 
 
 l.setLabel([19,80])
 qq.setLabel([80,10,8,20,54])

 A=copy.copy(c_up)
 A.setLabel([55,14,12,19,17])
 A.permute([55,14,12,17,19],4)
 
 q,r=qr_parity1(A) 
  
 q.setLabel([55,14,12,17,81])
 r.setLabel([81,19])

 Teta=l*r
 Teta.permute([81,80],1)
 U,s,V=svd_parity2(Teta)

 U.setLabel([81,19])
 s.setLabel([19,-19])
 V.setLabel([-19,80])

 s=Sqrt(s)
 
 U=U*s
 V=s*V

 U.permute([81,-19],1)
 U.setLabel([81,19])
 V.permute([19,80],1)
 
 d_up=qq*V
 c_up=q*U
 c_up.permute([55,14,12,19,17],3)
 d_up.permute([54,19,10,8,20],3)
#################################################################################
 A=copy.copy(a_up)
 A.setLabel([55,16,17,18,2])
 A.permute([17,55,16,18,2],1)
 
 l,qq=lq_parity1(A) 
 
 l.setLabel([17,80])
 qq.setLabel([80,55,16,18,2])

 A=copy.copy(c_up)
 A.setLabel([53,14,12,19,17])
 A.permute([53,14,12,19,17],4)
 q,r=qr_parity1(A) 
 
 
 q.setLabel([53,14,12,19,81])
 r.setLabel([81,17])

 Teta=l*r
 Teta.permute([81,80],1)
 U,s,V=svd_parity2(Teta)

 U.setLabel([81,17])
 s.setLabel([17,-17])
 V.setLabel([-17,80])

 s=Sqrt(s)
 
 U=U*s
 V=s*V

 U.permute([81,-17],1)
 U.setLabel([81,17])
 V.permute([17,80],1)
 
 a_up=qq*V
 c_up=q*U
 c_up.permute([53,14,12,19,17],3)
 a_up.permute([55,16,17,18,2],3)
 return a_up, c_up, d_up

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
def lq_parity(theta):
#    bd1=copy.copy(theta.bond(0))
#    bd2=copy.copy(theta.bond(1))
#    bd1.change(uni10.BD_OUT)
#    bd2.change(uni10.BD_OUT)
    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())
    bd2=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())    
    
    LA=uni10.UniTensor([theta.bond(0),theta.bond(1),bd1,bd2])
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4)])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA



def lq_parity1(theta):
    #bd1=copy.copy(theta.bond(0))
    #bd1.change(uni10.BD_OUT)
    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())

    
    LA=uni10.UniTensor([theta.bond(0),bd1])
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4)])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA


def lq_parity2(theta):
    #bd1=copy.copy(theta.bond(0))
    #bd1.change(uni10.BD_OUT)
    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())

    
    LA=uni10.UniTensor([theta.bond(0),bd1])
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3)])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA






    
def qr_parity1(theta):

    #bd1=copy.copy(theta.bond(3))
    #bd1.change(uni10.BD_IN)
    bd1=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4)])
    LA=uni10.UniTensor([bd1, theta.bond(4)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).qr()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])

#    print LA
    return GA, LA


def qr_parity2(theta):

    #bd1=copy.copy(theta.bond(3))
    #bd1.change(uni10.BD_IN)
    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3)])
    LA=uni10.UniTensor([bd1, theta.bond(3)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).qr()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])

#    print LA
    return GA, LA





def svd_parity2(theta):

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
#    print LA
    return GA, LA,GB


def   Sqrt(Landa):
  Landa_cp=copy.copy(Landa)
  blk_qnums=Landa.blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.getBlock(qnum).col())
   Landa_cpm=Landa_cp.getBlock(qnum,True)
   Landam=Landa_cp.getBlock(qnum,True)
   #print Landa_cpm[0], Landa_cpm[1], Landa_cpm[2], Landa_cpm[3]
   for i in xrange(D):
      if Landam[i] > 1.0e-12:
       Landa_cpm[i]=Landam[i]**(1.00/2.00)
      else:
       Landa_cpm[i]=0
   Landa_cp.putBlock(qnum,Landa_cpm)
  return Landa_cp 



def Sqrt_mat(e):
 d=int(e.row())
 for q in xrange(d):
   #print e[q] 
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e  

def randomize_plist(plist):
 for i in xrange(len(plist)):
  plist[i].randomize()




def equall_dis_plist(plist):
 plist0, plist1=equall_dis_qr(plist[0],plist[1])
 plist0.setLabel([62,58,57,17])
 plist1.setLabel([17,64,58,57]) 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64,58,57],1)
 plist[0]=plist0
 plist[1]=plist1
 
 plist0, plist1=equall_dis_qr(plist[2],plist[3])
 plist0.setLabel([66,59,60,18])
 plist1.setLabel([18,68,59,60]) 
 plist0.permute([66,59,60,18],3)
 plist1.permute([18,68,59,60],1)
 plist[2]=plist0
 plist[3]=plist1

def normal_plist(plist):
 if  (MaxAbs(plist[0]) > 1.0e+6) or (MaxAbs(plist[0]) < 1.0e-6):
    print "hi0", MaxAbs(plist[0])
    q, r=qr_parity2(plist[0])
    q.setLabel([62,58,57,17])
    plist[0]=q

 if  (MaxAbs(plist[2]) > 1.0e+6) or (MaxAbs(plist[2]) < 1.0e-6):
    print "hi2", MaxAbs(plist[2])
    q, r=qr_parity2(plist[2])
    q.setLabel([66,59,60,18])
    plist[2]=q

 if  (MaxAbs(plist[1]) > 1.0e+6) or (MaxAbs(plist[1]) < 1.0e-6):
    print "hi1", MaxAbs(plist[1])
    l, q=lq_parity2(plist[1])
    q.setLabel([17,64,58,57])
    plist[1]=q

 if  (MaxAbs(plist[3]) > 1.0e+6) or (MaxAbs(plist[3]) < 1.0e-6):
    print "hi3", MaxAbs(plist[3])
    l, q=lq_parity2(plist[3])
    q.setLabel([18,68,59,60])
    plist[3]=q

 for i in xrange(len(plist)):
   plist[i]=max_ten(plist[i])
 

def normal_plist1(plist):

 if  (MaxAbs(plist[0]) > 1.0e+6) or (MaxAbs(plist[0]) < 1.0e-6):
    print "hi0", MaxAbs(plist[0])
    q, r=qr_parity2(plist[0])
    q.setLabel([62,58,57,18])
    plist[0]=q

 if  (MaxAbs(plist[3]) > 1.0e+6) or (MaxAbs(plist[3]) < 1.0e-6):
    print "hi3", MaxAbs(plist[3])
    plist[3].permute([60,59,68,20],3)
    q, r=qr_parity2(plist[3])
    q.setLabel([60,59,68,20])
    q.permute([60,59,20,68],3)
    plist[3]=q

 if  (MaxAbs(plist[1]) > 1.0e+6) or (MaxAbs(plist[1]) < 1.0e-6):
    print "hi1", MaxAbs(plist[1])
    l, q=lq_parity2(plist[1])
    q.setLabel([18,64,58,57])
    plist[1]=q

 if  (MaxAbs(plist[2]) > 1.0e+6) or (MaxAbs(plist[2]) < 1.0e-6):
    print "hi2", MaxAbs(plist[2])
    plist[2].permute([20,60,59,66],1)
    l, q=lq_parity2(plist[2])
    q.setLabel([20,60,59,66])
    q.permute([66,60,59,20],1)    
    plist[2]=q

 for i in xrange(len(plist)):
   plist[i]=max_ten(plist[i])



def equall_dis_plist1(plist):
 plist0, plist1=equall_dis_qr(plist[0],plist[1])
 plist0.setLabel([62,58,57,18])
 plist1.setLabel([18,64,58,57]) 
 plist0.permute([62,58,57,18],3)
 plist1.permute([18,64,58,57],1)
 plist[0]=plist0
 plist[1]=plist1
 
 
 plist[3].permute([60,59,68,20],3)
 plist[2].permute([20,60,59,66],1)
 
 plist0, plist1=equall_dis_qr(plist[3],plist[2])
 plist0.setLabel([60,59,68,20])
 plist1.setLabel([20,60,59,66]) 
 plist0.permute([60,59,20,68],3)
 plist1.permute([66,60,59,20],1)
 plist[3]=plist0
 plist[2]=plist1





def equall_dis_qr(plist0, plist1 ):
 
 plist0.setLabel([62,58,57,17])
 plist1.setLabel([17,64,58,57])
 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64,58,57],1)
 
 q, r=qr_parity2(plist0)
 q.setLabel([62,58,57,1])
 r.setLabel([1,0])
 

 l, qq=lq_parity2(plist1)
 qq.setLabel([-1,64,58,57])
 l.setLabel([0,-1])

 teta=l*r
 teta.permute([1,-1],1)

 U, s, V =svd_parity2(teta)
 s=Sqrt(s)

 s.setLabel([0,17])
 U.setLabel([1,0])
 U=U*s

 s.setLabel([17,0])
 V.setLabel([0,-1])
 V=s*V

 plist0=q*U
 plist1=V*qq
 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64,58,57],1)
 return plist0, plist1
 


def initialize_plist(a_u, b_u, c_u, MPO_list): 

 plist=[]
 bd1=uni10.Bond(uni10.BD_IN,c_u.bond(4).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(1).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(2).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,c_u.bond(4).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([62,58,57,17])
 Uni_ten.randomize()
 Uni_ten, s, V=svd_parity3(Uni_ten)
 Uni_ten.setLabel([62,58,57,17])
 #Uni_ten.identity()
 plist.append(Uni_ten)
 
 plist.append(copy.copy(plist[0]))
 plist[1].transpose()
 plist[1].setLabel([17,64,58,57])
 
 
 bd1=uni10.Bond(uni10.BD_IN,a_u.bond(3).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[1].bond(3).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[1].bond(4).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,a_u.bond(3).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([66,59,60,18])
 Uni_ten.randomize()
 #Uni_ten.identity()
 Uni_ten, s, V=svd_parity3(Uni_ten)
 Uni_ten.setLabel([66,59,60,18])
 plist.append(Uni_ten)
 plist.append(copy.copy(plist[2]))
 plist[3].transpose()
 plist[3].setLabel([18,68,59,60])

 return plist


def initialize_plist1(a_u, b_u, d_u, MPO_list): 

 plist=[]
 
 bd1=uni10.Bond(uni10.BD_IN,a_u.bond(3).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(1).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(2).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,a_u.bond(3).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([62,58,57,18])
 Uni_ten.identity()
 Uni_ten.randomize()
 Uni_ten, s, V=svd_parity3(Uni_ten)
 Uni_ten.setLabel([62,58,57,18])
 plist.append(Uni_ten)

 plist.append(copy.copy(plist[0]))
 plist[1].transpose()
 plist[1].setLabel([18,64,58,57])

 
 bd1=uni10.Bond(uni10.BD_IN,d_u.bond(4).Qlist())
 bd2=uni10.Bond(uni10.BD_OUT,MPO_list[2].bond(0).Qlist())
 bd3=uni10.Bond(uni10.BD_OUT,MPO_list[2].bond(1).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,d_u.bond(4).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([66,60,59,20])
 Uni_ten.identity()
 Uni_ten.randomize()
 Uni_ten.transpose()
 Uni_ten, s, V=svd_parity3(Uni_ten)
 Uni_ten.transpose()
 Uni_ten.setLabel([66,60,59,20])
 plist.append(Uni_ten)

 plist.append(copy.copy(plist[2]))
 plist[3].transpose()
 plist[3].setLabel([60,59,20,68])

 return plist



##@profile
def Do_optimization_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,N_grad,N_svd,Gauge,check_step):
 
 plist_f=[copy.copy(plist[i])  for i in xrange(len(plist)) ]
 
 Distance_val=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
 print 'Dis0', Distance_val
 
 Res=1
 Res1=Distance_val
 count=0
 check_val=0
 for q in xrange(N_svd[0]):
  t0=time.time()
  optimum_1230(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f,Gauge,check_step)
  #optimum_0123(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f)
  #optimum_3201(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f)
  #optimum_2013(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f)


  Distance_val=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  Res=Res1
  Res1=Distance_val
  print 'Dis-mpo1=', Distance_val, abs(Res1-Res) / abs(Res), q, time.time() - t0
#  print  Res1 , Res,Res1 < Res 

  if Res1 < Res:
    plist_f[0]=plist[0]
    plist_f[1]=plist[1]
    plist_f[2]=plist[2]
    plist_f[3]=plist[3]
  else: 
    plist[0]=plist_f[0]
    plist[1]=plist_f[1]
    plist[2]=plist_f[2]
    plist[3]=plist_f[3]
    check_val=1
    break
  count+=1
  if count > 30: print 'Num_Opt > 30'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < N_svd[1]): 
    print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    check_val=0
    break
  else:
    if (abs(Distance_val) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     print 'break, Dis', Distance_val, abs(Res1-Res)
     break
# Distance_val=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
# print 'DisFMPO', Distance_val
 #if check_val!=0: 
  #Do_optimization_Grad_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,N_grad)


def optimum_0123(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f,Gauge):

  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  if Dis2 <= Dis1:
   plist_f[0]=plist[0]
  else: 
   print "0", Dis1, Dis2, Dis2 <= Dis1 
   plist[0]=plist_f[0]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  if Dis2 <= Dis1:
   plist_f[1]=plist[1]
  else: 
   print "1", Dis1, Dis2, Dis2 <= Dis1
   plist[1]=plist_f[1]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  if Dis2 < Dis1:
   plist_f[2]=plist[2]
  else: 
   print "2", Dis1, Dis2, Dis2 <= Dis1
   plist[2]=plist_f[2]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  if Dis2 <= Dis1:
   plist_f[3]=plist[3]
  else: 
   print "3", Dis1, Dis2, Dis2 <= Dis1
   plist[3]=plist_f[3]
 
def optimum_1230(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f,Gauge,check_step):


  if check_step is 'on':
   Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   #print "1", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[1]=plist[1]
   else: 
    plist[1]=plist_f[1]


   Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   #print "2", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[2]=plist[2]
   else: 
    plist[2]=plist_f[2]



   Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   #print "3", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[3]=plist[3]
   else: 
    plist[3]=plist_f[3]


   Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   #print "0", Dis1, Dis2, Dis2 <= Dis1 
   if Dis2 <= Dis1:
    plist_f[0]=plist[0]
   else: 
    plist[0]=plist_f[0]
  elif check_step is 'off':
   optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)
   optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge)



def optimum_3201(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f):

  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "3", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[3]=plist[3]
  else: 
   plist[3]=plist_f[3]

  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "2", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[2]=plist[2]
  else: 
   plist[2]=plist_f[2]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "1", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[1]=plist[1]
  else: 
   plist[1]=plist_f[1]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "0", Dis1, Dis2, Dis2 <= Dis1 
  if Dis2 <= Dis1:
   plist_f[0]=plist[0]
  else: 
   plist[0]=plist_f[0]

def optimum_2013(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,plist_f):

  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "2", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[2]=plist[2]
  else: 
   plist[2]=plist_f[2]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "0", Dis1, Dis2, Dis2 <= Dis1 
  if Dis2 <= Dis1:
   plist_f[0]=plist[0]
  else: 
   plist[0]=plist_f[0]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "1", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[1]=plist[1]
  else: 
   plist[1]=plist_f[1]


  Dis1=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist,Inv_method)
  Dis2=Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  #print "3", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[3]=plist[3]
  else: 
   plist[3]=plist_f[3]

##@profile
def optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)


###########################################################
 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])



 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.transpose()
 b_d1.transpose()
 a_d1.transpose()

 c_d1.setLabel([-58,-57,-19,-17,-54,-14,-12])
 a_d1.setLabel([-18,-2,-59,-60,-55,-16,-17,-58,-57])
 b_d1.setLabel([-6,-4,-56,-18,-20,-59,-60])




##########################################################
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])

 c_ut=((c_u*(MPO_list[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,62,58,57],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.transpose()
 b_dt.transpose()
 a_dt.transpose()


 c_dt.setLabel([-19,-62,-58,-57,-54,-14,-12])
 a_dt.setLabel([-18,-2,-55,-16,-17])
 b_dt.setLabel([-6,-4,-56,-18,-20])

 A=((((E1*E8)*(a_ut*a_dt))*((E2*E3)*(b_ut*b_dt)))*(((E4*E5)*d)))*((E7*E6)*(c_ut*c_dt))
 A.permute([-62,-58,-57,-17,62,58,57,17],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 
# A=N_Positiv(A)
# A.setLabel([-62,-58,-57,-17,62,58,57,17])



 Ap=(((((E1*E8)*(a_u1*a_dt))) * ((E2*E3)*(b_u1*b_dt)))*(((E4*E5)*d)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-58,-57,-17],4)

 if Gauge is "Fixed":
  Ap1=(((((E1*E8)*(a_ut*a_d1)))*((E2*E3)*(b_ut*b_d1)))*(((E4*E5)*d)))*(((E7*E6)*(c_ut*c_d1)))
  Ap1.permute([62,58,57,17],0)
  Ap1.transpose()
  Ap=Ap+Ap1

 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)
   S=inverse(S)
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])
   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([62,58,57,17,-62,-58,-57,-17])
   A=A_inv*Ap
   A.permute([62,58,57,17],3)
   plist[0]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([62,58,57,17])
   A.permute([62,58,57,17],3)   
   plist[0]=A 
 
 
 
#######################################################



##@profile
def optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)


###########################################################
 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])



 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.transpose()
 b_d1.transpose()
 a_d1.transpose()

 c_d1.setLabel([-58,-57,-19,-17,-54,-14,-12])
 a_d1.setLabel([-18,-2,-59,-60,-55,-16,-17,-58,-57])
 b_d1.setLabel([-6,-4,-56,-18,-20,-59,-60])




##########################################################
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])

 c_ut=((c_u*(plist[0]*MPO_list[0])))
 
 a_ut=(a_u*(plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,64,58,57,18,2],5)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.transpose()
 b_dt.transpose()
 a_dt.transpose()


 c_dt.setLabel([-19,-17,-54,-14,-12])
 a_dt.setLabel([-18,-2,-55,-16,-64,-58,-57])
 b_dt.setLabel([-6,-4,-56,-18,-20])


 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-17,-64,-58,-57,17,64,58,57],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1

# A=N_Positiv(A)
# A.setLabel([-17,-64,-58,-57,17,64,58,57])


 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-17,-64,-58,-57],4)

 if Gauge is "Fixed":
  A1=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
  A1.permute([17,64,58,57],0)
  A1.transpose()
  Ap=Ap+A1
 #Ap.setLabel([-17,-58,-57,-64])

 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)
   S=inverse(S)
   
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])
   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([17,64,58,57,-17,-64,-58,-57])
   A=A_inv*Ap
   A.permute([17,64,58,57],1)
   plist[1]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([17,64,58,57])
   A.permute([17,64,58,57],1)   
   plist[1]=A 







##@profile
def optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)


###########################################################
 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])



 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.transpose()
 b_d1.transpose()
 a_d1.transpose()

 c_d1.setLabel([-58,-57,-19,-17,-54,-14,-12])
 a_d1.setLabel([-18,-2,-59,-60,-55,-16,-17,-58,-57])
 b_d1.setLabel([-6,-4,-56,-18,-20,-59,-60])




##########################################################
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])

 a_ut=(a_u*(plist[1]*MPO_list[1]))
 c_ut=((c_u*(plist[0]*MPO_list[0])))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,66,59,60,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)


 c_dt.transpose()
 b_dt.transpose()
 a_dt.transpose()


 c_dt.setLabel([-19,-17,-54,-14,-12])
 a_dt.setLabel([-66,-59,-60,-2,-55,-16,-17])
 b_dt.setLabel([-6,-4,-56,-18,-20])

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-66,-59,-60,-18,66,59,60,18],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1

# A=N_Positiv(A)
# A.setLabel([-66,-59,-60,-18,66,59,60,18])
    
 Ap=((((E4*E5)*d)*((E2*E3)*(b_u1*b_dt)))*((E7*E6)*(c_u1*c_dt)))*((E1*E8)*(a_u1*a_dt))
 
 Ap.permute([-66,-59,-60,-18],4)

 if Gauge is "Fixed":           
  A1=((((E4*E5)*d)*((E2*E3)*(b_ut*b_d1)))*((E7*E6)*(c_ut*c_d1)))*((E1*E8)*(a_ut*a_d1))
  A1.permute([66,59,60,18],0)
  A1.transpose()
  Ap=Ap+A1

 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)
   S=inverse(S)
   
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])

   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([66,59,60,18,-66,-59,-60,-18])
   A=A_inv*Ap
   A.permute([66,59,60,18],3)
   plist[2]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([66,59,60,18])
   A.permute([66,59,60,18],3)   
   plist[2]=A 









##@profile
def optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,MPO_list,c_u,a_u,b_u,plist,Inv_method,Gauge):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)


###########################################################
 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])



 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.transpose()
 b_d1.transpose()
 a_d1.transpose()

 c_d1.setLabel([-58,-57,-19,-17,-54,-14,-12])
 a_d1.setLabel([-18,-2,-59,-60,-55,-16,-17,-58,-57])
 b_d1.setLabel([-6,-4,-56,-18,-20,-59,-60])




##########################################################
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])

 c_ut=((c_u*(plist[0]*MPO_list[0])))
 
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,68,59,60,20,6,4],5)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)


 c_dt.transpose()
 b_dt.transpose()
 a_dt.transpose()


 c_dt.setLabel([-19,-17,-54,-14,-12])
 a_dt.setLabel([-18,-2,-55,-16,-17])
 b_dt.setLabel([-6,-4,-56,-68,-59,-60,-20])

 A=((((E1*E8)*(a_ut*a_dt)) *((E7*E6)*(c_ut*c_dt))) * ((((E4*E5)*d)))) * ((E2*E3)*(b_ut*b_dt))
 A.permute([-18,-68,-59,-60,18,68,59,60],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1

# A=N_Positiv(A)
# A.setLabel([-18,-68,-59,-60,18,68,59,60])


 Ap=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-68,-59,-60],4)

 if Gauge is "Fixed":
  A1=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)))*((E2*E3)*(b_ut*b_d1))
  A1.permute([18,68,59,60],0)
  A1.transpose()
  Ap=Ap+A1



 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)
   S=inverse(S)
   
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([18,68,59,60,-18,-68,-59,-60])

   A=A_inv*Ap
   #A.transpose()
   A.permute([18,68,59,60],1)
   plist[3]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([18,68,59,60])
   A.permute([18,68,59,60],1)   
   plist[3]=A 





##@profile
def  Dis_ff(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)

###########################################################
 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])



 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)
 
 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.transpose()
 b_d1.transpose()
 a_d1.transpose()


 c_d1.setLabel([-58,-57,-19,-17,-54,-14,-12])
 a_d1.setLabel([-18,-2,-59,-60,-55,-16,-17,-58,-57])
 b_d1.setLabel([-6,-4,-56,-18,-20,-59,-60])

##########################################################
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])


 c_ut=((c_u*(plist[0]*MPO_list[0])))
 #print plist[1].printDiagram(), plist[2].printDiagram(), MPO_list[1].printDiagram(), a_u.printDiagram()  
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.transpose()
 b_dt.transpose()
 a_dt.transpose()


 c_dt.setLabel([-19,-17,-54,-14,-12])
 a_dt.setLabel([-18,-2,-55,-16,-17])
 b_dt.setLabel([-6,-4,-56,-18,-20])


 Val=(((((E1*E8)*(a_ut*a_dt))*((E7*E6)*(c_ut*c_dt))))*(((E2*E3)*(b_ut*b_dt))))*((E4*E5)*d)
 #print 'Val=',Val
 #Val1=(((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c_u1*c_d1))))*((E2*E3)*(b_u1*b_d1)))*(((E4*E5)*d))
 #print 'Val1=',Val1
 Val2=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E2*E3)*(b_ut*b_d1))))*((E4*E5)*d)
 #print 'Val2=',Val2
 #Val3=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E2*E3)*(b_u1*b_dt))))*((E4*E5)*d)
 #print 'Val3=',Val3

 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 val_f=Val[0]-2.00*Val2[0]
 return  val_f

def  recover( c_u, a_u, b_u, plist, MPO_list):

 MPO_list0=[copy.copy(MPO_list[i]) for i in xrange(len(MPO_list))]
 MPO_list0[0].setLabel([-54,58,57,54])
 MPO_list0[1].setLabel([58,57,-55,59,60,55])
 MPO_list0[2].setLabel([60,59,-56,56]) 
 
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])

 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 c_ut.permute([-54,14,12,19,17],3)
  
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 a_ut.permute([-55,16,17,18,2],3)

 b_ut=((b_u*(plist[3]*MPO_list0[2])))
 b_ut.permute([-56,18,20,6,4],3)
 
 return c_ut,a_ut,b_ut

def Store_plist(plist):
 plist[0].save("Store/plist[0]")
 plist[1].save("Store/plist[1]")
 plist[2].save("Store/plist[2]")
 plist[3].save("Store/plist[3]")

def Reload_plist(plist):
 plist[0]=uni10.UniTensor("Store/plist[0]")
 plist[1]=uni10.UniTensor("Store/plist[1]")
 plist[2]=uni10.UniTensor("Store/plist[2]")
 plist[3]=uni10.UniTensor("Store/plist[3]")
 return plist


def Store_plist1(plist):
 plist[0].save("Store/plist1[0]")
 plist[1].save("Store/plist1[1]")
 plist[2].save("Store/plist1[2]")
 plist[3].save("Store/plist1[3]")

def Reload_plist1(plist):
 plist[0]=uni10.UniTensor("Store/plist1[0]")
 plist[1]=uni10.UniTensor("Store/plist1[1]")
 plist[2]=uni10.UniTensor("Store/plist1[2]")
 plist[3]=uni10.UniTensor("Store/plist1[3]")
 return plist



##@profile
def Do_optimization_MPO1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,N_grad,N_svd,Gauge,check_step):
 plist_f=[copy.copy(plist[i])  for i in xrange(len(plist)) ]
 
 Distance_val=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
 print 'Dis0', Distance_val
 
 Res=1
 Res1=Distance_val
 count=0
 check_val=0
 for q in xrange(N_svd[0]):
  t0=time.time()

  optimum_1230s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f,Gauge,check_step)
  #optimum_0123s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f)
#  optimum_3201s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f)
#  optimum_2013s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f)


  Distance_val=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  Res=Res1
  Res1=Distance_val
  print 'Dis-mpo2=', Distance_val, abs(Res1-Res) / abs(Res), q, time.time() - t0
  if Res1 < Res:
    plist_f[0]=plist[0]
    plist_f[1]=plist[1]
    plist_f[2]=plist[2]
    plist_f[3]=plist[3]
  else: 
    plist[0]=plist_f[0]
    plist[1]=plist_f[1]
    plist[2]=plist_f[2]
    plist[3]=plist_f[3]
    check_val=1
    break

  count+=1
  if count > 40: print 'Num_Opt > 40'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < N_svd[1]): 
    print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     print 'break, Dis', Distance_val, abs(Res1-Res)
     break

# Distance_val=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
# print 'DisMPOf', Distance_val


# if check_val !=0:
#  Do_optimization_Grad_MPO1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,N_grad)



def optimum_0123s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f,Gauge):

  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  if Dis2 <= Dis1:
   plist_f[0]=plist[0]
  else: 
   print "0", Dis1, Dis2, Dis2 <= Dis1
   plist[0]=plist_f[0]

  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  if Dis2 <= Dis1:
   plist_f[1]=plist[1]
  else:
   print "1", Dis1, Dis2, Dis2 <= Dis1
   plist[1]=plist_f[1]


  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  if Dis2 <= Dis1:
   plist_f[2]=plist[2]
  else: 
   print "3", Dis1, Dis2, Dis2 <= Dis1
   plist[2]=plist_f[2]


  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  if Dis2 <= Dis1:
   plist_f[3]=plist[3]
  else: 
   print "2", Dis1, Dis2, Dis2 <= Dis1
   plist[3]=plist_f[3]





 
def optimum_1230s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,
MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f,Gauge,check_step):

  if check_step is 'on':
   Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   #print "1", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[1]=plist[1]
   else:
    plist[1]=plist_f[1]

   Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   #print "3", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[2]=plist[2]
   else: 
    plist[2]=plist_f[2]

   Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   #print "2", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[3]=plist[3]
   else: 
    plist[3]=plist_f[3]


   Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
   #print "0", Dis1, Dis2, Dis2 <= Dis1
   if Dis2 <= Dis1:
    plist_f[0]=plist[0]
   else: 
    plist[0]=plist_f[0]
  elif check_step is 'off':
   optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)
   optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge)



def optimum_3201s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f):

  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "2", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[3]=plist[3]
  else: 
   plist[3]=plist_f[3]

  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "3", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[2]=plist[2]
  else: 
   plist[2]=plist_f[2]



  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "0", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[0]=plist[0]
  else: 
   plist[0]=plist_f[0]


  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "1", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[1]=plist[1]
  else:
   plist[1]=plist_f[1]



def optimum_2013s(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,plist_f):

  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "3", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[2]=plist[2]
  else: 
   plist[2]=plist_f[2]



  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "0", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[0]=plist[0]
  else: 
   plist[0]=plist_f[0]


  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "1", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[1]=plist[1]
  else:
   plist[1]=plist_f[1]


  Dis1=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,a_u,b_u,d_u,plist,Inv_method)
  Dis2=Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist)
  #print "2", Dis1, Dis2, Dis2 <= Dis1
  if Dis2 <= Dis1:
   plist_f[3]=plist[3]
  else: 
   plist[3]=plist_f[3]










##@profile
def  Dis_ff1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist):

 c.setLabel([14,-14,12,-12,19,-19,17,-17])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 d_u1=copy.copy(d_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 b_u1.setLabel([55,18,20,6,4])
 d_u1.setLabel([56,19,10,8,20])

 a_u1=(a_u1*(MPO_list[0]))
 b_u1=((b_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 b_u1.permute([-55,18,20,58,57,6,4,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 b_d1=copy.copy(b_u1)
 d_d1=copy.copy(d_u1)

 a_d1.transpose()
 b_d1.transpose()
 d_d1.transpose()
 
 a_d1.setLabel([-58,-57,-18,-2,-54,-16,-17])
 b_d1.setLabel([-6,-4,-59,-60,-55,-18,-20,-58,-57])
 d_d1.setLabel([-8,-20,-56,-19,-10,-59,-60])
 
##########################################################
 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[1]*plist[3]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,18,20,6,4],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-18,-20])
 d_dt.setLabel([-8,-20,-56,-19,-10])
# print a_ut.printDiagram(), a_u.printDiagram(), b_ut.printDiagram(), b_u.printDiagram()
# print d_ut.printDiagram(), d_u.printDiagram()

 
 Val=((((E1*E8)*(a_ut*a_dt))*((E7*E6)*(c))))*(((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b_ut*b_dt)))
 #print 'Val=',Val
# Val1=((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c))))*(((E4*E5)*(d_d1*d_u1))*((E2*E3)*(b_d1*b_u1)))
# print 'Val1=',Val1
 Val2=(((((E1*E8)*(a_ut*a_d1))*((E2*E3)*(b_ut*b_d1))))*(((E4*E5)*(d_ut*d_d1))))*((E7*E6)*(c))
 #print 'Val2=',Val2
# Val3=((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c))))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b_u1*b_dt)))
# print 'Val3=',Val3
 
 val_f=Val[0]-2.0*Val2[0]
 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]

 if Val[0]<0: 
  print "Dis, N<0"
  val_f=-Val[0]-2.0*Val2[0]
 
 return  val_f



##@profile
def optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge):
 c.setLabel([14,-14,12,-12,19,-19,17,-17])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 d_u1=copy.copy(d_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 b_u1.setLabel([55,18,20,6,4])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 b_u1=((b_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 b_u1.permute([-55,18,20,58,57,6,4,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 b_d1=copy.copy(b_u1)
 d_d1=copy.copy(d_u1)

 a_d1.transpose()
 b_d1.transpose()
 d_d1.transpose()
 
 a_d1.setLabel([-58,-57,-18,-2,-54,-16,-17])
 b_d1.setLabel([-6,-4,-59,-60,-55,-18,-20,-58,-57])
 d_d1.setLabel([-8,-20,-56,-19,-10,-59,-60])
 
##########################################################
 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(MPO_list[0])))
 b_ut=(b_u*(plist[1]*plist[3]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))

 a_ut.permute([-54,16,17,62,58,57,2],3)
 b_ut.permute([-55,18,20,6,4],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-62,-58,-57,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-18,-20])
 d_dt.setLabel([-8,-20,-56,-19,-10])

 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b_ut*b_dt)))*(((E7*E6)*(c))))*((E1*E8)*(a_ut*a_dt))
 A.permute([-62,-58,-57,-18,62,58,57,18],4)


 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 
# A=N_Positiv(A)
# A.setLabel([-62,-58,-57,-17,62,58,57,17])



 Ap=((((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b_u1*b_dt)))*(((E7*E6)*(c))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-62,-58,-57,-18],4)

 if Gauge is "Fixed":
  Ap1=((((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b_ut*b_d1)))*(((E7*E6)*(c))))*((E1*E8)*(a_ut*a_d1))
  Ap1.permute([62,58,57,18],0)
  Ap1.transpose()
  Ap=Ap+Ap1


 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   #print S
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)

   S=inverse(S)
   
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([62,58,57,18,-62,-58,-57,-18])

  #######################################################
   A=A_inv*Ap
   A.permute([62,58,57,18],3)
   plist[0]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([62,58,57,18])
   A.permute([62,58,57,18],3)   
   plist[0]=A 





##@profile
def optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge):
 c.setLabel([14,-14,12,-12,19,-19,17,-17])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 d_u1=copy.copy(d_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 b_u1.setLabel([55,18,20,6,4])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 b_u1=((b_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 b_u1.permute([-55,18,20,58,57,6,4,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 b_d1=copy.copy(b_u1)
 d_d1=copy.copy(d_u1)

 a_d1.transpose()
 b_d1.transpose()
 d_d1.transpose()
 

 a_d1.setLabel([-58,-57,-18,-2,-54,-16,-17])
 b_d1.setLabel([-6,-4,-59,-60,-55,-18,-20,-58,-57])
 d_d1.setLabel([-8,-20,-56,-19,-10,-59,-60])
 
 
 
 
##########################################################
 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[3]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,64,58,57,20,6,4],5)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-64,-58,-57,-20])
 d_dt.setLabel([-8,-20,-56,-19,-10])


 A=(((((E4*E5)*(d_ut*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_dt))))*((E2*E3)*(b_ut*b_dt))

 A.permute([-18,-64,-58,-57,18,64,58,57],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1

# A=N_Positiv(A)
# A.setLabel([-17,-64,-58,-57,17,64,58,57])


 Ap=(((((E4*E5)*(d_u1*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_u1*a_dt))))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-64,-58,-57],4)


 if Gauge is "Fixed":
  A1=(((((E4*E5)*(d_ut*d_d1))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_d1))))*((E2*E3)*(b_ut*b_d1))
  A1.permute([18,64,58,57],0)
  A1.transpose()
  Ap=Ap+A1

 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)

   S=inverse(S)
   
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([18,64,58,57,-18,-64,-58,-57])

   A=A_inv*Ap
   #A.transpose()
   A.permute([18,64,58,57],1)
   plist[1]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([18,64,58,57])
   A.permute([18,64,58,57],1)   
   plist[1]=A 





##@profile
def optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge):
 c.setLabel([14,-14,12,-12,19,-19,17,-17])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 d_u1=copy.copy(d_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 b_u1.setLabel([55,18,20,6,4])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 b_u1=((b_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 b_u1.permute([-55,18,20,58,57,6,4,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 b_d1=copy.copy(b_u1)
 d_d1=copy.copy(d_u1)

 a_d1.transpose()
 b_d1.transpose()
 d_d1.transpose()
 

 a_d1.setLabel([-58,-57,-18,-2,-54,-16,-17])
 b_d1.setLabel([-6,-4,-59,-60,-55,-18,-20,-58,-57])
 d_d1.setLabel([-8,-20,-56,-19,-10,-59,-60])
 
 
 
##########################################################
 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[1]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,18,68,59,60,6,4],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-59,-60,-6,-4,-55,-18,-68])
 d_dt.setLabel([-8,-20,-56,-19,-10])

 A=(((((E4*E5)*(d_ut*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_dt))))*((E2*E3)*(b_ut*b_dt))
 A.permute([-60,-59,-20,-68,60,59,20,68],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1

# A=N_Positiv(A)
# A.setLabel([-17,-64,-58,-57,17,64,58,57])


 Ap=(((((E4*E5)*(d_u1*d_dt))))*(((E7*E6)*(c))*((E1*E8)*(a_u1*a_dt))))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-60,-59,-20,-68],4)

 if Gauge is "Fixed":
  A1=(((((E4*E5)*(d_ut*d_d1))))*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_d1))))*((E2*E3)*(b_ut*b_d1))
  A1.permute([60,59,20,68],0)
  A1.transpose()
  Ap=Ap+A1
 #Ap.setLabel([-17,-58,-57,-64])


 if Inv_method is 'SVD':
   U, S, V=svd_parity(A)
   U.transpose()
   V.transpose()
   S=Normal_Singulars(S)

   S=inverse(S)
   
   U.setLabel([8,9,10,11,12,13,14,15])
   S.setLabel([4,5,6,7,8,9,10,11])
   V.setLabel([0,1,2,3,4,5,6,7])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,12,13,14,15],4)
   A_inv.setLabel([60,59,20,68,-60,-59,-20,-68])

   A=A_inv*Ap
   A.permute([60,59,20,68],3)
   plist[3]=A
 elif Inv_method is 'CG':
   A=solve_linear_eq(A,Ap)
   A.setLabel([60,59,20,68])
   A.permute([60,59,20,68],3)   
   plist[3]=A 





##@profile
def optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,MPO_list,a_u,b_u,d_u,plist,Inv_method,Gauge):
 c.setLabel([14,-14,12,-12,19,-19,17,-17])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 d_u1=copy.copy(d_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 b_u1.setLabel([55,18,20,6,4])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 b_u1=((b_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 b_u1.permute([-55,18,20,58,57,6,4,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 b_d1=copy.copy(b_u1)
 d_d1=copy.copy(d_u1)

 a_d1.transpose()
 b_d1.transpose()
 d_d1.transpose()
 

 a_d1.setLabel([-58,-57,-18,-2,-54,-16,-17])
 b_d1.setLabel([-6,-4,-59,-60,-55,-18,-20,-58,-57])
 d_d1.setLabel([-8,-20,-56,-19,-10,-59,-60])
 
##########################################################
 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[1]*plist[3]*MPO_list[1]))
 d_ut=((d_u*(MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,18,20,6,4],3)
 d_ut.permute([-56,19,10,59,60,8,66],5)

 a_dt=copy.copy(a_ut)
 b_dt=copy.copy(b_ut)
 d_dt=copy.copy(d_ut)

 a_dt.transpose()
 b_dt.transpose()
 d_dt.transpose()


 a_dt.setLabel([-18,-2,-54,-16,-17])
 b_dt.setLabel([-6,-4,-55,-18,-20])
 d_dt.setLabel([-8,-66,-56,-19,-10,-59,-60])

 A=((((E2*E3)*(b_ut*b_dt)) )*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_dt))))* (((E4*E5)*(d_ut*d_dt)))
 A.permute([-66,-60,-59,-20,66,60,59,20],4)

 if Gauge is "Fixed":
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1

# A=N_Positiv(A)
# A.setLabel([-18,-68,-59,-60,18,68,59,60])


 Ap=((((E2*E3)*(b_u1*b_dt)) )*(((E7*E6)*(c))*((E1*E8)*(a_u1*a_dt))))* (((E4*E5)*(d_u1*d_dt)))
 Ap.permute([-66,-60,-59,-20],4)

 if Gauge is "Fixed":
  A1=((((E2*E3)*(b_ut*b_d1)) )*(((E7*E6)*(c))*((E1*E8)*(a_ut*a_d1))))* (((E4*E5)*(d_ut*d_d1)))
  A1.permute([66,60,59,20],0)
  A1.transpose()
  Ap=Ap+A1

 if Inv_method is 'SVD':
    U, S, V=svd_parity(A)
    U.transpose()
    V.transpose()
    S=Normal_Singulars(S)
    S=inverse(S)
    
    U.setLabel([8,9,10,11,12,13,14,15])
    S.setLabel([4,5,6,7,8,9,10,11])
    V.setLabel([0,1,2,3,4,5,6,7])


    A_inv=V*S*U
    A_inv.permute([0,1,2,3,12,13,14,15],4)
    A_inv.setLabel([66,60,59,20,-66,-60,-59,-20])

    A=A_inv*Ap
    A.permute([66,60,59,20],1)
    plist[2]=A

 elif Inv_method is 'CG':
    A=solve_linear_eq(A,Ap)
    A.setLabel([66,60,59,20])
    A.permute([66,60,59,20],1)   
    plist[2]=A 



def  recover1( a_u,b_u, d_u, plist, MPO_list):

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56]) 

 a_u.setLabel([54,16,17,62,2])
 b_u.setLabel([55,64,68,6,4])
 d_u.setLabel([56,19,10,8,66])


 a_ut=((a_u*(plist[0]*MPO_list[0])))
 b_ut=(b_u*(plist[1]*plist[3]*MPO_list[1]))
 d_ut=((d_u*(plist[2]*MPO_list[2])))

 a_ut.permute([-54,16,17,18,2],3)
 b_ut.permute([-55,18,20,6,4],3)
 d_ut.permute([-56,19,10,8,20],3)
 
 return a_ut, b_ut,d_ut 
 

#@profile
def Do_optimization_svd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Inv_method,N_svd,Gauge,check_step):

 a_up_first=copy.copy(a_up) 
 b_up_first=copy.copy(b_up)
 c_up_first=copy.copy(c_up)
 d_up_first=copy.copy(d_up)

 a_up_first1=copy.copy(a_up) 
 b_up_first1=copy.copy(b_up)
 c_up_first1=copy.copy(c_up)
 d_up_first1=copy.copy(d_up)


 checking_val=0

 Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
 print Distance_val
 Res=10
 Res1=Distance_val
 count=0
 test=0
 for q in xrange(N_svd[0]):
  t0=time.time()

  if check_step is "on":
   Dis1=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   a_up=opt_a(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      a_up_first=copy.copy(a_up) 
   else: 
      print "Fail0:a", Dis1, Dis2, Dis2 <= Dis1  
      test=1; break;
      a_up=copy.copy(a_up_first) 
   Dis1=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   b_up=opt_b(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      b_up_first=copy.copy(b_up) 
   else: 
      print "Fail0:b", Dis1, Dis2, Dis2 <= Dis1
      test=1; break;  
      b_up=copy.copy(b_up_first) 
   Dis1=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   c_up=opt_c(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      c_up_first=copy.copy(c_up) 
   else: 
      print "Fail0:c", Dis1, Dis2, Dis2 <= Dis1  
      test=1; break;
      c_up=copy.copy(c_up_first) 
   Dis1=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   d_up=opt_d(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      d_up_first=copy.copy(d_up) 
   else: 
      print "Fail0:d", Dis1, Dis2, Dis2 <= Dis1  
      test=1; break;
      d_up=copy.copy(d_up_first) 
  elif check_step is "off":
   a_up=opt_a(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   b_up=opt_b(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   c_up=opt_c(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   d_up=opt_d(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   
   
  Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
  print 'Dis_svd1=', Distance_val, abs(Res1-Res) / abs(Res), q, time.time() - t0
#  print test

  Res=Res1
  Res1=Distance_val
  #print Res, Res1, Res1 < Res 
  #test=q
  if Res1 < Res:
    a_up_first=copy.copy(a_up) 
    b_up_first=copy.copy(b_up)
    c_up_first=copy.copy(c_up)
    d_up_first=copy.copy(d_up)
  else: 
    a_up=copy.copy(a_up_first) 
    b_up=copy.copy(b_up_first)
    c_up=copy.copy(c_up_first)
    d_up=copy.copy(d_up_first)
    #test=q
    break
 
  count+=1
  if count > 100: print 'Num_Opt > 30'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-8) or ((abs(Res1-Res) / abs(Res)) < N_svd[1]): 
    #print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val) < 1.00e-8) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break

# print test
# if test != 0:
##  a_up=copy.copy(a_up_first1)
##  b_up=copy.copy(b_up_first1)
##  c_up=copy.copy(c_up_first1)
##  d_up=copy.copy(d_up_first1)
#  a_up, b_up, c_up, d_up=Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Gauge)
 return a_up, b_up, c_up, d_up


#@profile
def opt_a(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-54,-55,-56,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,54,55,56])
 U_2=U*H

 U.setLabel([-54,-55,-56,54,55,56])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([-54,2,-56,1,55,0])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([1,-55,0,-54,2,-56])

 Iden=Idena*Idenb
 Iden.permute([-55,55],1)
 Iden.identity()
 #print Iden
 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---a---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E1*E8)*Iden)
 #print A.printDiagram()
 A.permute([-55,-16,-17,-18,-2,55,16,17,18,2],5)

 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

# t0=time.time()


 Ap=((((((E2*E3)*(b_u*b_dp)))*((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u*c_dp)*U))*((E1*E8)*(a_u))
 Ap.permute([-55,-16,-17,-18,-2],5)
 
 if Gauge is "Fixed": 
  Ap1=((((((E2*E3)*(b_up*b_d)))*((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_up*c_d)*U))*((E1*E8)*(a_d))
  Ap1.permute([55,16,17,18,2],0)
  Ap1.transpose()
  Ap=Ap+Ap1

 
 if Inv_method is "SVD":
   U, S, V=svd_parity4(A)
   U.transpose()
   V.transpose()
   S=inverse(S)
   
   U.setLabel([10,11,12,13,14,15,16,17,18,19])
   S.setLabel([5,6,7,8,9,10,11,12,13,14])
   V.setLabel([0,1,2,3,4,5,6,7,8,9])

   A_inv=V*S*U
   A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
   A_inv.setLabel([55,16,17,18,2,-55,-16,-17,-18,-2])
   A=A_inv*Ap
   A.permute([55,16,17,18,2],3)

 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([55,16,17,18,2])
   A.permute([55,16,17,18,2],3)
 
 return A


#@profile
def opt_b(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-54,-55,-56,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,54,55,56])
 U_2=U*H

 U.setLabel([-54,-55,-56,54,55,56])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([-54,2,-51,1,0,56])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([1,0,-56,-54,2,-51])

 Iden=Idena*Idenb
 Iden.permute([-56,56],1)
 Iden.identity()
 #print Iden
 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
 
 
 
###########################################---b---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))*(((E2*E3)*(Iden)))
 #print A.printDiagram()
 A.permute([-56,-18,-20,-6,-4,56,18,20,6,4],5)

 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])


# t0=time.time()


 Ap=(((((E1*E8)*(a_u*a_dp)*U)*((E7*E6)*(c_u*c_dp))))*((E4*E5)*(d_u*d_dp)))*(((E2*E3)*(b_u)))
 Ap.permute([-56,-18,-20,-6,-4],5)
 
 if Gauge is "Fixed": 
  Ap1=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*((E4*E5)*(d_up*d_d)))*(((E2*E3)*(b_d)))
  Ap1.permute([56,18,20,6,4],0)
  Ap1.transpose()
  Ap=Ap+Ap1


 if Inv_method is "SVD":
   U, S, V=svd_parity4(A)
   U.transpose()
   V.transpose()
   S=inverse(S)
   
   U.setLabel([10,11,12,13,14,15,16,17,18,19])
   S.setLabel([5,6,7,8,9,10,11,12,13,14])
   V.setLabel([0,1,2,3,4,5,6,7,8,9])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
   A_inv.setLabel([56,18,20,6,4,-56,-18,-20,-6,-4])

   A=A_inv*Ap
  # A.transpose()
  # A=solve_linear_eq(A,Ap)
  # A.setLabel([56,18,20,6,4])
   A.permute([56,18,20,6,4],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([56,18,20,6,4])
   A.permute([56,18,20,6,4],3)
 return A


#@profile
def opt_c(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-54,-55,-56,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,54,55,56])
 U_2=U*H

 U.setLabel([-54,-55,-56,54,55,56])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([-23,2,-51,54,0,56])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([-54,0,56,-23,2,-51])

 Iden=Idena*Idenb
 Iden.permute([-54,54],1)
 Iden.identity()
 #print Iden
 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
 
###########################################---c---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))* ((E7*E6)*(Iden))
 #print A.printDiagram()
 A.permute([-54,-14,-12,-19,-17,54,14,12,19,17],5)
 
 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])



 Ap=(((((E1*E8)*(a_u*a_dp))*((E2*E3)*(b_u*b_dp)))*U)*(((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u))
 Ap.permute([-54,-14,-12,-19,-17],5)
 
 if Gauge is "Fixed": 
  Ap1=(((((E1*E8)*(a_up*a_d))*((E2*E3)*(b_up*b_d)))*U)*(((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_d))
  Ap1.permute([54,14,12,19,17],0)
  Ap1.transpose()
  Ap=Ap+Ap1



 if Inv_method is "SVD":
   U, S, V=svd_parity4(A)
   U.transpose()
   V.transpose()
   S=inverse(S)
   
   U.setLabel([10,11,12,13,14,15,16,17,18,19])
   S.setLabel([5,6,7,8,9,10,11,12,13,14])
   V.setLabel([0,1,2,3,4,5,6,7,8,9])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
   A_inv.setLabel([54,14,12,19,17,-54,-14,-12,-19,-17])

   A=A_inv*Ap
   A.permute([54,14,12,19,17],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([54,14,12,19,17])
   A.permute([54,14,12,19,17],3)

 return A

#@profile
def opt_d(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-54,-55,-56,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,54,55,56])
 U_2=U*H

 U.setLabel([-54,-55,-56,54,55,56])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([-23,2,-51,57,0,56])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([-57,0,56,-23,2,-51])

 Iden=Idena*Idenb
 Iden.permute([-57,57],1)
 Iden.identity()
 #print Iden
 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([-57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,-54,-14,-12])
 d_d.setLabel([-8,-20,57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
 
###########################################---d---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(Iden))
 #print A.printDiagram()
 A.permute([-57,-19,-10,-8,-20,57,19,10,8,20],5)

 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])





 Ap=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u))
 Ap.permute([-57,-19,-10,-8,-20],5)
 

 if Gauge is "Fixed": 
  Ap1=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d)))*U)*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_d))
  Ap1.permute([57,19,10,8,20],0)
  Ap1.transpose()
  Ap=Ap+Ap1



 if Inv_method is "SVD":
    U, S, V=svd_parity4(A)
    U.transpose()
    V.transpose()
    S=inverse(S)
    
    U.setLabel([10,11,12,13,14,15,16,17,18,19])
    S.setLabel([5,6,7,8,9,10,11,12,13,14])
    V.setLabel([0,1,2,3,4,5,6,7,8,9])


    A_inv=V*S*U
    A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
    A_inv.setLabel([57,19,10,8,20,-57,-19,-10,-8,-20])

    A=A_inv*Ap
    A.permute([57,19,10,8,20],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([57,19,10,8,20])
   A.permute([57,19,10,8,20],3)


 return A






#@profile
def Do_optimization_svd1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Inv_method,N_svd,Gauge, check_step):

 a_up_first=copy.copy(a_up) 
 b_up_first=copy.copy(b_up)
 c_up_first=copy.copy(c_up)
 d_up_first=copy.copy(d_up)

 a_up_first1=copy.copy(a_up) 
 b_up_first1=copy.copy(b_up)
 c_up_first1=copy.copy(c_up)
 d_up_first1=copy.copy(d_up)


 checking_val=0

 Distance_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
 print Distance_val
 Res=10
 Res1=Distance_val
 count=0
 test=0
 for q in xrange(N_svd[0]):
  t0=time.time()

  if check_step is "on":
   Dis1=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   a_up=opt_a1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      a_up_first=copy.copy(a_up) 
   else: 
      print "Fail:a", Dis1, Dis2,Dis2 <= Dis1  
      test=1; break;
      a_up=copy.copy(a_up_first) 

   Dis1=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   b_up=opt_b1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      b_up_first=copy.copy(b_up) 
   else: 
      print "Fail:b", Dis1, Dis2, Dis2 <= Dis1  
      test=1; break;
      b_up=copy.copy(b_up_first) 


   Dis1=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   d_up=opt_d1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      d_up_first=copy.copy(d_up) 
   else:
      print "Fail:d", Dis1, Dis2, Dis2 <= Dis1  
      test=1; break;
      d_up=copy.copy(d_up_first) 


   Dis1=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   c_up=opt_c1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   Dis2=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   if Dis2 <= Dis1:
      c_up_first=copy.copy(c_up) 
   else:
      print "Fail:d", Dis1, Dis2, Dis2 <= Dis1  
      test=1; break;
      c_up=copy.copy(c_up_first) 
  elif check_step is "off":
   a_up=opt_a1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   b_up=opt_b1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   d_up=opt_d1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)
   c_up=opt_c1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge)


   Distance_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   print 'Dis_svd2=', Distance_val, abs(Res1-Res) / abs(Res), q, time.time() - t0

   Res=Res1
   Res1=Distance_val
   #print Res, Res1, Res1 < Res 
   #test=q
   if Res1 < Res:
     a_up_first=copy.copy(a_up) 
     b_up_first=copy.copy(b_up)
     c_up_first=copy.copy(c_up)
     d_up_first=copy.copy(d_up)
   else: 
     a_up=copy.copy(a_up_first) 
     b_up=copy.copy(b_up_first)
     c_up=copy.copy(c_up_first)
     d_up=copy.copy(d_up_first)
     #test=q
     break
 
  count+=1
  if count > 100: print 'Num_Opt > 30'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-8) or ((abs(Res1-Res) / abs(Res)) < N_svd[1]): 
    #print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val) < 1.00e-8) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break

# if test != 0:
##  a_up=copy.copy(a_up_first1)
##  b_up=copy.copy(b_up_first1)
##  c_up=copy.copy(c_up_first1)
##  d_up=copy.copy(d_up_first1)
#  a_up, b_up, c_up, d_up=Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method,Gauge)

 return a_up, b_up, c_up, d_up


#@profile
def opt_a1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-55,-56,-57,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,55,56,57])
 U_2=U*H

 U.setLabel([-55,-56,-57,55,56,57])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([10,20,30,55,0,1])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([-55,0,1,10,20,30])

 Iden=Idena*Idenb
 Iden.permute([-55,55],1)
 Iden.identity()

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---a---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E1*E8)*Iden)
 #print A.printDiagram()
 A.permute([-55,-16,-17,-18,-2,55,16,17,18,2],5)
 
 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])

# t0=time.time()


 Ap=((((((E2*E3)*(b_u*b_dp)))*((E4*E5)*(d_u*d_dp)))*U)*((E7*E6)*(c_u*c_dp)))*((E1*E8)*(a_u))
 #print Ap.printDiagram()
 Ap.permute([-55,-16,-17,-18,-2],5)
 

 if Gauge is "Fixed": 
  Ap1=((((((E2*E3)*(b_up*b_d)))*((E4*E5)*(d_up*d_d)))*U)*((E7*E6)*(c_up*c_d)))*((E1*E8)*(a_d))
  Ap1.permute([55,16,17,18,2],0)
  Ap1.transpose()
  Ap=Ap+Ap1



 if Inv_method is "SVD":
    U, S, V=svd_parity4(A)
    U.transpose()
    V.transpose()
    S=inverse(S)
    
    U.setLabel([10,11,12,13,14,15,16,17,18,19])
    S.setLabel([5,6,7,8,9,10,11,12,13,14])
    V.setLabel([0,1,2,3,4,5,6,7,8,9])


    A_inv=V*S*U
    A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
    A_inv.setLabel([55,16,17,18,2,-55,-16,-17,-18,-2])

    A=A_inv*Ap
    A.permute([55,16,17,18,2],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([55,16,17,18,2])
   A.permute([55,16,17,18,2],3)

 return A


#@profile
def opt_b1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-55,-56,-57,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,55,56,57])
 U_2=U*H

 U.setLabel([-55,-56,-57,55,56,57])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([10,20,30,0,56,1])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([0,-56,1,10,20,30])

 Iden=Idena*Idenb
 Iden.permute([-56,56],1)
 Iden.identity()

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---b---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))*(((E2*E3)*(Iden)))
 #print A.printDiagram()
 A.permute([-56,-18,-20,-6,-4,56,18,20,6,4],5)

 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])



# t0=time.time()


 Ap=((((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp))))*((E4*E5)*(d_u*d_dp)))*U)*(((E2*E3)*(b_u)))
 Ap.permute([-56,-18,-20,-6,-4],5)
 
 if Gauge is "Fixed": 
  Ap1=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*((E4*E5)*(d_up*d_d)))*(((E2*E3)*(b_d)))
  Ap1.permute([56,18,20,6,4],0)
  Ap1.transpose()
  Ap=Ap+Ap1


 if Inv_method is "SVD":
   U, S, V=svd_parity4(A)
   U.transpose()
   V.transpose()
   S=inverse(S)
   
   U.setLabel([10,11,12,13,14,15,16,17,18,19])
   S.setLabel([5,6,7,8,9,10,11,12,13,14])
   V.setLabel([0,1,2,3,4,5,6,7,8,9])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
   A_inv.setLabel([56,18,20,6,4,-56,-18,-20,-6,-4])

   A=A_inv*Ap
   A.permute([56,18,20,6,4],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([56,18,20,6,4])
   A.permute([56,18,20,6,4],3)

 return A


#@profile
def opt_d1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-55,-56,-57,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,55,56,57])
 U_2=U*H

 U.setLabel([-55,-56,-57,55,56,57])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([10,20,30,0,1,57])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([0,1,-57,10,20,30])

 Iden=Idena*Idenb
 Iden.permute([-57,57],1)
 Iden.identity()

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---d---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(Iden))
 #print A.printDiagram()
 A.permute([-57,-19,-10,-8,-20,57,19,10,8,20],5)
 
 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


# t0=time.time()


 Ap=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp))))*(((E2*E3)*(b_u*b_dp))*U))*((E4*E5)*(d_u))
 Ap.permute([-57,-19,-10,-8,-20],5)
 
 if Gauge is "Fixed": 
  Ap1=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d))))*(((E2*E3)*(b_up*b_d))*U))*((E4*E5)*(d_d))
  Ap1.permute([57,19,10,8,20],0)
  Ap1.transpose()
  Ap=Ap+Ap1



 if Inv_method is "SVD":
   U, S, V=svd_parity4(A)
   U.transpose()
   V.transpose()
   S=inverse(S)
   
   U.setLabel([10,11,12,13,14,15,16,17,18,19])
   S.setLabel([5,6,7,8,9,10,11,12,13,14])
   V.setLabel([0,1,2,3,4,5,6,7,8,9])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
   A_inv.setLabel([57,19,10,8,20,-57,-19,-10,-8,-20])


   A=A_inv*Ap
   A.permute([57,19,10,8,20],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([57,19,10,8,20])
   A.permute([57,19,10,8,20],3)

 return A

#@profile
def opt_c1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,Inv_method,Gauge):

 U.setLabel([-55,-56,-57,0,1,2])
 H=copy.copy(U)
 H.setLabel([0,1,2,55,56,57])
 U_2=U*H

 U.setLabel([-55,-56,-57,55,56,57])
 Idenb=uni10.UniTensor(U.bond())
 Idenb.identity()
 Idenb.setLabel([10,20,30,0,1,57])

 Idena=uni10.UniTensor(U.bond())
 Idena.identity()
 Idena.setLabel([0,1,-57,10,20,30])

 Iden=Idena*Idenb
 Iden.permute([-57,57],1)
 Iden.setLabel([-54,54])
 Iden.identity()

 a_u.setLabel([55,16,17,18,2])
 b_u.setLabel([56,18,20,6,4])
 c_u.setLabel([-54,14,12,19,17])
 d_u.setLabel([57,19,10,8,20])


 a_d=copy.copy(a_u)
 b_d=copy.copy(b_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 c_d.transpose()
 b_d.transpose()
 a_d.transpose()
 d_d.transpose()
 
 a_d.setLabel([-18,-2,-55,-16,-17])
 b_d.setLabel([-6,-4,-56,-18,-20]) 
 c_d.setLabel([-19,-17,54,-14,-12])
 d_d.setLabel([-8,-20,-57,-19,-10])


 a_up.setLabel([55,16,17,18,2])
 b_up.setLabel([56,18,20,6,4])
 c_up.setLabel([54,14,12,19,17])
 d_up.setLabel([57,19,10,8,20])

 b_dp=copy.copy(b_up)
 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 b_dp.transpose()
 a_dp.transpose()
 c_dp.transpose()
 d_dp.transpose()
###########################################---d---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))* ((E7*E6)*(Iden))
 #print A.printDiagram()
 A.permute([-54,-14,-12,-19,-17,54,14,12,19,17],5)
 
 if Gauge is "Fixed": 
  A1=copy.copy(A)
  A1.transpose()
  A=A+A1
 


 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


# t0=time.time()


 Ap=(((((E1*E8)*(a_u*a_dp))*((E2*E3)*(b_u*b_dp)))*U)*(((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u))
 Ap.permute([-54,-14,-12,-19,-17],5)
 
 if Gauge is "Fixed": 
  Ap1=(((((E1*E8)*(a_up*a_d))*((E2*E3)*(b_up*b_d)))*U)*(((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_d))
  Ap1.permute([54,14,12,19,17],0)
  Ap1.transpose()
  Ap=Ap+Ap1



 if Inv_method is "SVD":
   U, S, V=svd_parity4(A)
   U.transpose()
   V.transpose()
   S=inverse(S)
   
   U.setLabel([10,11,12,13,14,15,16,17,18,19])
   S.setLabel([5,6,7,8,9,10,11,12,13,14])
   V.setLabel([0,1,2,3,4,5,6,7,8,9])


   A_inv=V*S*U
   A_inv.permute([0,1,2,3,4,15,16,17,18,19],5)
   A_inv.setLabel([54,14,12,19,17,-54,-14,-12,-19,-17])


   A=A_inv*Ap
   A.permute([54,14,12,19,17],3)
 elif Inv_method is "CG":
   A=solve_linear_eq(A,Ap)
   A.setLabel([54,14,12,19,17])
   A.permute([54,14,12,19,17],3)

 return A



def Mat_np_to_Uni(Mat_np):
 d0=np.size(Mat_np,0)
 d1=np.size(Mat_np,1)
 Mat_uni=uni10.Matrix(d0,d1)
 for i in xrange(d0):
  for j in xrange(d1):
   Mat_uni[i*d1+j]=Mat_np[i,j]
 return  Mat_uni


 
def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np



###@profile
def solve_linear_eq(A,Ap):
 Ap_h=copy.copy(Ap)
 Ap_h.transpose()
 Result=uni10.UniTensor(Ap.bond())
 blk_qnums = A.blockQnum()
 #print blk_qnums, Ap.printDiagram()
 blk_qnums1 = Ap.blockQnum()
   
 for qnum in blk_qnums:
   if qnum in blk_qnums1:
    A_mat=A.getBlock(qnum)
    Ap_mat=Ap.getBlock(qnum)
    A_np=Mat_uni_to_np(A_mat)
    b_np=Mat_uni_to_np(Ap_mat)
    x_np=np.linalg.lstsq(A_np, b_np)[0] 
    #x_np=sp.linalg.lstsq(A_np, b_np)[0] 
    x=Mat_np_to_Uni(x_np)
    Result.putBlock(qnum, x)
 return Result


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
 if ( Max_val < 0.5e-1) or (Max_val > 0.5e+1)   :

  if Max_val >= 1:
   print ">1",Max_val
   a=a*(1.00/Max_val)
  if Max_val < 1: 
   print "<1",Max_val
   a=a*(1.00/Max_val)

#  if Max_val >= 1:
#   print ">1",Max_val, Max_val**(1./2.)
#   a=a*(1.00/(Max_val**(1./2.)))
#  if Max_val < 1: 
#   print "<1",Max_val
#   a=a*(1.00/Max_val)

 else: a=a;
 return a
 
def Normal_Singulars(S):
# Max_val=MaxAbs(S)
# print "MaxAbs", Max_val
# if Max_val > 1.0e+3 or Max_val < 1.0e-1:
#  S=S*(1.00/Max_val)
# Max_val=MaxAbs(S)
# print "MaxAbs1", Max_val

 return S



