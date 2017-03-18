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
import MoveCorboz
import MoveFull
import basic_FU
import basic


def Var_cab(a_u, b_u, c_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic_FU.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 #t0=time.time()
 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 #print time.time() - t0, "CTM-H, Left"

 Env=basic_FU.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 #t0=time.time()
 E1, E2, E3, E4, E5, E6, E7, E8=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 #print time.time() - t0, "Env, Left"
 #test_env(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d)
 #test_energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u)
 #E=energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,a_u,b_u,c_u)
 #print "Energy", E
 
 
 N_u, l_u, r_u, q_u, qq_u = Qr_lQ_decom(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u)
 c_up, b_up = reproduce_cb(r_u, l_u, q_u, qq_u)
 #test_energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_up,b_up)
 N_u=N_Positiv(N_u)
 
 l_up=copy.copy(l_u)
 r_up=copy.copy(r_u)
 
 
 l_up, r_up=Do_optimization_Full(N_u, l_u, r_u, l_up, r_up, U)

 #Dis_val=Dis_fQR(N_u, l_u, r_u, l_up, r_up, U )
 #print "DisFFF", Dis_val

 c_up, b_up =reproduce_cb(r_up, l_up, q_u, qq_u)

 #Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u,c_up,b_up,a_up)
 #print "Dis_final", Dis_val

 c_up, b_up, a_up=equall_dis(c_up,b_up,a_u) 
 #a_up=copy.copy(a_u)
 #a=basic.make_ab(a_u)
 #Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u,a_u,c_up,b_up,a_up)
 #print "Dis_final", Dis_val


 c_up=basic.max_ten(c_up)
 b_up=basic.max_ten(b_up)
 a_up=basic.max_ten(a_up)

 
 cp=basic.make_ab(c_up)
 bp=basic.make_ab(b_up)
 ap=basic.make_ab(a_up)

 return  a_up, b_up, c_up , ap, bp, cp 



def Var_cad(a_u, b_u, d_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Positive,Corner_method):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic_FU.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 #t0=time.time()
 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 #print time.time() - t0, "CTM-H, Left"

 Env=basic_FU.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 #t0=time.time()
 E1, E2, E3, E4, E5, E6, E7, E8=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 #print time.time() - t0, "Env, Left"
 #test_env(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d)
 #test_energy_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u)
 #E=energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,a_u,b_u,c_u)
 #print "Energy", E
 
 
 N_u, l_u, r_u, q_u, qq_u = Qr_lQ_decom_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u)
 a_up,d_up = reproduce_ad(r_u, l_u, q_u, qq_u)

 #test_energy_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_up,d_up)
 
 
 N_u=N_Positiv_1(N_u)
 
 l_up=copy.copy(l_u)
 r_up=copy.copy(r_u)
 
 #Dis_val=Dis_fQR_1(N_u, l_u, r_u, l_up, r_up, U )
 #print "Dis", Dis_val


 
 
 l_up, r_up=Do_optimization_Full_1(N_u, l_u, r_u, l_up, r_up, U)

 #Dis_val=Dis_fQR_1(N_u, l_u, r_u, l_up, r_up, U )
 #print "DisFFF", Dis_val

 a_up,d_up = reproduce_ad(r_up, l_up, q_u, qq_u)

 d_up, a_up, b_up=equall_dis_1(d_up,a_up,b_u) 
 #b_up=copy.copy(b_u)
 #b=basic.make_ab(b_up)
 #Dis_val=Dis_f_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u,b_u,a_up,d_up,b_up)
 #print "Dis_final", Dis_val


 d_up=basic.max_ten(d_up)
 b_up=basic.max_ten(b_up)
 a_up=basic.max_ten(a_up)

 
 dp=basic.make_ab(d_up)
 bp=basic.make_ab(b_up)
 ap=basic.make_ab(a_up)

 return  a_up, b_up, d_up, ap, bp, dp 


def produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys):

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

 return E1, E2, E3, E4, E5, E6, E7, E8 

def test_env(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d):


 a.setLabel([16,-16,17,-17,18,-18,2,-2])
 b.setLabel([18,-18,20,-20,6,-6,4,-4])
 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 d.setLabel([19,-19,10,-10,8,-8,20,-20])



 Val=(((E1*E8)*a)*((E7*E6)*c))*(((E4*E5)*d)*((E2*E3)*b))
 #Val=((E1*E8)*a)
 
 print "Env_test", Val[0]
def  test_energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])


 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 

 b_u.setLabel([54,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-6,-4,52,-18,-20])


 c_u.setLabel([53,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.transpose()
 c_d.setLabel([-19,-17,51,-14,-12])


 Val=((((E1*E8)*(a))*((E7*E6)*(c_u*c_d)))*Iden)*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 Norm_f=Val
 print 'Env_test=', Val[0]

 Val=((((E1*E8)*(a))*((E7*E6)*(c_u*c_d)))*U)*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 print 'E=',Val[0]/Norm_f[0]

def  energy_cb(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])


 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 

 b_u.setLabel([54,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-6,-4,52,-18,-20])


 c_u.setLabel([53,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.transpose()
 c_d.setLabel([-19,-17,51,-14,-12])


 Val=((((E1*E8)*(a))*((E7*E6)*(c_u*c_d)))*Iden)*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 Norm_f=Val

 Val=((((E1*E8)*(a))*((E7*E6)*(c_u*c_d)))*U)*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 return Val[0]/Norm_f[0]


def  test_energy_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u):

 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 b.setLabel([18,-18,20,-20,6,-6,4,-4])


 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 

 a_u.setLabel([53,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-18,-2,51,-16,-17])


 d_u.setLabel([54,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,52,-19,-10])


 Val=((((E1*E8)*(a_d*a_u))*((E7*E6)*(c)))*Iden)*(((E4*E5)*(d_d*d_u))*((E2*E3)*(b)))
 Norm_f=Val
 print 'Env_test=', Val[0]

 Val=((((E1*E8)*(a_d*a_u))*((E7*E6)*(c)))*U)*(((E4*E5)*(d_d*d_u))*((E2*E3)*(b)))
 print 'E=',Val[0]/Norm_f[0]

def  energy_ad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u):

 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 b.setLabel([18,-18,20,-20,6,-6,4,-4])


 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 

 a_u.setLabel([53,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-18,-2,51,-16,-17])


 d_u.setLabel([54,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,52,-19,-10])


 Val=((((E1*E8)*(a_d*a_u))*((E7*E6)*(c)))*Iden)*(((E4*E5)*(d_d*d_u))*((E2*E3)*(b)))
 Norm_f=Val

 Val=((((E1*E8)*(a_d*a_u))*((E7*E6)*(c)))*U)*(((E4*E5)*(d_d*d_u))*((E2*E3)*(b)))
 return Val[0]/Norm_f[0]



def Qr_lQ_decom(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])

 A=copy.copy(c_u)
 A.setLabel([53,14,12,19,17])
 A.permute([14,12,19,53,17],3)
 
 
 
 q,r=qr_parity(A) 
 
 
 q.setLabel([14,12,19,82,83])
 r.setLabel([82,83,53,17])
 q.permute([14,12,19,82,83],2)
 r.permute([82,83,53,17],3)
 
 q_d=copy.copy(q)
 q_d.transpose()
 q_d.setLabel([-19,-82,-83,-14,-12])
 
 A=copy.copy(b_u)
 A.setLabel([54,18,20,6,4])
 A.permute([54,18,20,6,4],2)
 
 l,qq=lq_parity(A) 
 
 l.setLabel([54,18,80,81])
 qq.setLabel([80,81,20,6,4])
 l.permute([54,18,80,81],2)
 qq.permute([80,81,20,6,4],3)
 
 qq_d=copy.copy(qq)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-80,-81,-20])
 
 N=(((E4*E5)*d)*((E2*E3)*(qq*qq_d)))*(((E1*E8)*(a))*((E7*E6)*(q*q_d)))
 N.permute([80,81,17,-82,-83,-18,-80,-81,-17,82,83,18],6)
 
 
 ###testing####
 r_d=copy.copy(r)
 r_d.transpose()
 l_d=copy.copy(l)
 l_d.transpose()
 
 
 r_d.setLabel([-17,-82,-83,53])
 l_d.setLabel([-80,-81,54,-18])

# Norm=((N*(r*l))*(r_d*l_d))
# print "Norm-QR", Norm[0]
 return N, l, r, q, qq


def Qr_lQ_decom_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u):
 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 b.setLabel([18,-18,20,-20,6,-6,4,-4])





 A=copy.copy(a_u)
 A.setLabel([53,16,17,18,2])
 A.permute([16,17,2,53,18],3)
 
 
 
 q,r=qr_parity(A) 
 
 
 q.setLabel([16,17,2,82,83])
 r.setLabel([82,83,53,18])
 q.permute([16,17,82,83,2],2)
 r.permute([82,83,53,18],3)
 
 q_d=copy.copy(q)
 q_d.transpose()
 q_d.setLabel([-82,-83,-2,-16,-17])
 
 A=copy.copy(d_u)
 A.setLabel([54,19,10,8,20])
 A.permute([19,10,8,54,20],3)
 
 qq , l=qr_parity(A) 
 
 l.setLabel([80,81,54,20])
 qq.setLabel([19,10,8,80,81])
 l.permute([80,81,54,20],3)
 qq.permute([19,10,8,80,81],2)
 
 qq_d=copy.copy(qq)
 qq_d.transpose()
 qq_d.setLabel([-8,-80,-81,-19,-10])
 
 N=(((E4*E5)*(qq*qq_d))*((E2*E3)*(b)))*(((E1*E8)*(q*q_d))*((E7*E6)*(c)))
 N.permute([-80,-81,-20,-82,-83,-18,80,81,20,82,83,18],6)
 

 ###testing####
 r_d=copy.copy(r)
 r_d.transpose()
 l_d=copy.copy(l)
 l_d.transpose()
 ##############
 
 r_d.setLabel([-18,-82,-83,53])
 l_d.setLabel([-20,-80,-81,54])

# Norm=((N*(r*l))*(r_d*l_d))
# print "Norm-QR", Norm[0]
 return N, l, r, q, qq



 
def reproduce_cb(r_u, l_u, q_u, qq_u):

 c_up=q_u*r_u
 c_up.permute([53,14,12,19,17],3)


 b_up=qq_u*l_u
 b_up.permute([54,18,20,6,4],3)


 return c_up, b_up
 
def reproduce_ad(r_u, l_u, q_u, qq_u):

 a_up=q_u*r_u
 a_up.permute([53,16,17,18,2],3)


 d_up=qq_u*l_u
 d_up.permute([54,19,10,8,20],3)


 return a_up, d_up

 
def  Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,c_u,b_u,a_u,c_up,b_up,a_up):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])


 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=U*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([51,52,53,54])
 
 

 a_up.setLabel([54,16,17,18,2])
 a_dp=copy.copy(a_up)
 a_dp.transpose()
 a_dp.setLabel([-18,-2,54,-16,-17])

 a_u.setLabel([54,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-18,-2,54,-16,-17])
  

 b_u.setLabel([54,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-6,-4,52,-18,-20])


 c_u.setLabel([53,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.transpose()
 c_d.setLabel([-19,-17,51,-14,-12])
 
 b_up.setLabel([54,18,20,6,4])
 b_dp=copy.copy(b_up)
 b_dp.transpose()
 b_dp.setLabel([-6,-4,52,-18,-20])


 c_up.setLabel([53,14,12,19,17])
 c_dp=copy.copy(c_up)
 c_dp.transpose()
 c_dp.setLabel([-19,-17,51,-14,-12])
 
 
 Val=(((((E1*E8)*(a))*((E7*E6)*(c_u*c_d)))*H)*(((E2*E3)*(b_u*b_d))))*((E4*E5)*d)
 print 'Val=',Val[0]
 Val1=(((((E1*E8)*(a_dp*a_up))*((E7*E6)*(c_up*c_dp)))*Iden)*((E2*E3)*(b_up*b_dp)))*(((E4*E5)*d))
 print 'Val1=',Val1[0]
 Val2=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*d)
 print 'Val2=',Val2[0]
 Val3=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d)))*U)*(((E2*E3)*(b_up*b_d))))*((E4*E5)*d)
 print 'Val3=',Val3[0]
 
 val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
# val_f=Val[0]-2.00*Val2[0]#-Val3[0]

 return  val_f 

def  Dis_f_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U,a_u,d_u,b_u,a_up,d_up,b_up):

 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 b.setLabel([18,-18,20,-20,6,-6,4,-4])




 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=U*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([51,52,53,54])
 
 

 a_up.setLabel([53,16,17,18,2])
 a_dp=copy.copy(a_up)
 a_dp.transpose()
 a_dp.setLabel([-18,-2,51,-16,-17])

 a_u.setLabel([53,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-18,-2,51,-16,-17])
  

 b_u.setLabel([55,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-6,-4,55,-18,-20])

 b_up.setLabel([55,18,20,6,4])
 b_dp=copy.copy(b_up)
 b_dp.transpose()
 b_dp.setLabel([-6,-4,55,-18,-20])


 d_u.setLabel([54,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,52,-19,-10])
 


 d_up.setLabel([54,19,10,8,20])
 d_dp=copy.copy(d_up)
 d_dp.transpose()
 d_dp.setLabel([-8,-20,52,-19,-10])
 
 
 Val=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c)))*H)*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))
 print 'Val=',Val[0]
 Val1=(((((E1*E8)*(a_dp*a_up))*((E7*E6)*(c)))*Iden)*((E2*E3)*(b_up*b_dp)))*(((E4*E5)*(d_up*d_dp)))
 print 'Val1=',Val1[0]
 Val2=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u*d_dp))
 print 'Val2=',Val2[0]
 Val3=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c)))*U)*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_up*d_d))
 print 'Val3=',Val3[0]
 
 val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
# val_f=Val[0]-2.00*Val2[0]#-Val3[0]

 return  val_f 











def Dis_fQR(N, l, r, lp, rp, U ):
 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=U*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([51,52,53,54])




 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d.setLabel([-17,-82,-83,51])
 l_d.setLabel([-80,-81,52,-18])

 r_dp=copy.copy(rp)
 l_dp=copy.copy(lp)
 l_dp.transpose()
 r_dp.transpose()
 
 
 r_dp.setLabel([-17,-82,-83,51])
 l_dp.setLabel([-80,-81,52,-18])


 Val=((N*(r*l))*(r_d*l_d)*H)
 #print Val[0]
 Val1=((N*(rp*lp))*(r_dp*l_dp)*Iden)
 #print Val1[0]
 Val2=((N*(rp*lp))*(r_d*l_d)*U)
 #print Val2[0]
 Val3=((N*(r*l))*(r_dp*l_dp)*U)
 #print Val3[0]
 Val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]

 return Val_f


def Dis_fQR_1(N, l, r, lp, rp, U ):
 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=U*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([51,52,53,54])

# r.permute([82,83,53,18],3)
# l.permute([80,81,54,20],3)

 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d.setLabel([-18,-82,-83,51])
 l_d.setLabel([-20,-80,-81,52])

 r_dp=copy.copy(rp)
 l_dp=copy.copy(lp)
 l_dp.transpose()
 r_dp.transpose()
 
 
 r_dp.setLabel([-18,-82,-83,51])
 l_dp.setLabel([-20,-80,-81,52])


 Val=((N*(r*l))*(r_d*l_d)*H)
 #print Val
 Val1=((N*(rp*lp))*(r_dp*l_dp)*Iden)
 #print Val1
 Val2=((N*(rp*lp))*(r_d*l_d)*U)
 #print Val2
 Val3=((N*(r*l))*(r_dp*l_dp)*U)
 #print Val3
 Val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]

 return Val_f





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

def Do_optimization_Full(N_u, l_u, r_u, l_up, r_up, U):
 
 Res=1
 Res1=2
 count=0
 #Distance_val=Dis_fQR(N_u, l_u, r_u, l_up, r_up, U)
 for q in xrange(20):
  #print "\n", "\n"
  Distance_val=Dis_fQR(N_u, l_u, r_u, l_up, r_up, U)
  #print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q
  r_up=optimum_0(N_u, l_u, r_u, l_up, r_up, U)
  l_up=optimum_1(N_u, l_u, r_u, l_up, r_up, U)
  Distance_val=Dis_fQR(N_u, l_u, r_u, l_up, r_up, U)

  Res=Res1
  Res1=Distance_val
  count+=1
  if count > 20: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-9) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-9) or (  abs(Res1-Res) < 1.00e-11  ): 
     print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break
 Distance_val=Dis_fQR(N_u, l_u, r_u, l_up, r_up, U)
 #print 'Dis_f', Distance_val, abs(Res1-Res) / abs(Res), q

 return l_up, r_up

def Do_optimization_Full_1(N_u, l_u, r_u, l_up, r_up, U):
 
 Res=1
 Res1=2
 count=0
 #Distance_val=Dis_fQR_1(N_u, l_u, r_u, l_up, r_up, U)
 for q in xrange(20):
  #print "\n", "\n"
  Distance_val=Dis_fQR_1(N_u, l_u, r_u, l_up, r_up, U)
  #print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q
  r_up=optimum_00(N_u, l_u, r_u, l_up, r_up, U)
  l_up=optimum_11(N_u, l_u, r_u, l_up, r_up, U)

  Res=Res1
  Res1=Distance_val
  count+=1
  if count > 20: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break
 Distance_val=Dis_fQR_1(N_u, l_u, r_u, l_up, r_up, U)
 #print 'Dis_f', Distance_val, abs(Res1-Res) / abs(Res), q

 return l_up, r_up




def optimum_0(N, l, r, lp, rp, U):
 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d.setLabel([-17,-82,-83,51])
 l_d.setLabel([-80,-81,52,-18])

 r_dp=copy.copy(rp)
 l_dp=copy.copy(lp)
 l_dp.transpose()
 r_dp.transpose()
 
 r_dp.setLabel([-17,-82,-83,51])
 l_dp.setLabel([-80,-81,52,-18])

 A2=(((lp*l_dp)*N)*Iden)
 A2.permute([-82,-83,51,-17,82,83,53,17],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-82,-83,51,-17,82,83,53,17])
 
 A3=((r)*U*(l*l_dp))*N
 A3.permute([-82,-83,51,-17],0)
 A3p=((r_d)*U*(l_d*lp))*N
 A3p.permute([82,83,53,17],4)
 A3p.transpose()
 #A3=addition_symmetric(A3,A3p)

 A3=A3+A3p
 
 A3.setLabel([-82,-83,51,-17])
 
 U, S, V=svd_parity(A2)
 #print U.printDiagram()
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([82,83,53,17,-82,-83,51,-17])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([82,83,53,17])
 A.permute([82,83,53,17],3)

 rf=copy.copy(A)

 return rf
 
def optimum_1(N, l, r, lp, rp, U):
 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d.setLabel([-17,-82,-83,51])
 l_d.setLabel([-80,-81,52,-18])

 r_dp=copy.copy(rp)
 l_dp=copy.copy(lp)
 l_dp.transpose()
 r_dp.transpose()
 
 r_dp.setLabel([-17,-82,-83,51])
 l_dp.setLabel([-80,-81,52,-18])

 A2=(((rp*r_dp)*N)*Iden)
 A2.permute([52,-18,-80,-81,54,18,80,81],4)
 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 A2.setLabel([52,-18,-80,-81,54,18,80,81])


 A3=((l)*U*(r*r_dp))*N
 A3.permute([52,-18,-80,-81],0)
 A3p=((l_d)*U*(rp*r_d))*N
 A3p.permute([54,18,80,81],4)
 A3p.transpose()
 A3=A3+A3p
 A3.setLabel([52,-18,-80,-81])
 
 
 U, S, V=svd_parity(A2)

 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([54,18,80,81,52,-18,-80,-81])




 A=A3*A2_inv
 A.setLabel([54,18,80,81])
 A.permute([54,18,80,81],2)


 lf=copy.copy(A)

 return lf


def optimum_00(N, l, r, lp, rp, U):
 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d.setLabel([-18,-82,-83,51])
 l_d.setLabel([-20,-80,-81,52])

 r_dp=copy.copy(rp)
 l_dp=copy.copy(lp)
 l_dp.transpose()
 r_dp.transpose()
 
 
 r_dp.setLabel([-18,-82,-83,51])
 l_dp.setLabel([-20,-80,-81,52])

 A2=(((lp*l_dp)*N)*Iden)
 #print A2.printDiagram(), r.printDiagram(),l.printDiagram()
 A2.permute([-82,-83,51,-18,82,83,53,18],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-82,-83,51,-18,82,83,53,18])
 
 A3=((r)*U*(l*l_dp))*N
 A3.permute([-82,-83,51,-18],0)
 A3p=((r_d)*U*(l_d*lp))*N
 A3p.permute([82,83,53,18],4)
 A3p.transpose()
 #A3=addition_symmetric(A3,A3p)

 A3=A3+A3p
 
 A3.setLabel([-82,-83,51,-18])
 
 U, S, V=svd_parity(A2)
 #print U.printDiagram()
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([82,83,53,18,-82,-83,51,-18])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([82,83,53,18])
 A.permute([82,83,53,18],3)

 rf=copy.copy(A)

 return rf
 
def optimum_11(N, l, r, lp, rp, U):
 U.setLabel([51,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54])

 r_d=copy.copy(r)
 l_d=copy.copy(l)
 l_d.transpose()
 r_d.transpose()
 
 
 r_d.setLabel([-18,-82,-83,51])
 l_d.setLabel([-20,-80,-81,52])

 r_dp=copy.copy(rp)
 l_dp=copy.copy(lp)
 l_dp.transpose()
 r_dp.transpose()
 
 
 r_dp.setLabel([-18,-82,-83,51])
 l_dp.setLabel([-20,-80,-81,52])

 A2=(((rp*r_dp)*N)*Iden)
 A2.permute([52,-20,-80,-81,54,20,80,81],4)
 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 A2.setLabel([52,-20,-80,-81,54,20,80,81])


 A3=((l)*U*(r*r_dp))*N
 A3.permute([52,-20,-80,-81],0)
 A3p=((l_d)*U*(rp*r_d))*N
 A3p.permute([54,20,80,81],4)
 A3p.transpose()
 A3=A3+A3p
 A3.setLabel([52,-20,-80,-81])
 
 
 U, S, V=svd_parity(A2)

 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([54,20,80,81,52,-20,-80,-81])




 A=A3*A2_inv
 A.setLabel([54,20,80,81])
 A.permute([80,81,54,20],3)


 lf=copy.copy(A)

 return lf
 
def equall_dis(c_up, b_up, a_u):

 A=copy.copy(b_up)
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



 
 
 return c_up,b_up,a_up



def equall_dis_1(d_up,a_up,b_u):

 A=copy.copy(b_u)
 A.setLabel([54,18,20,6,4])
 A.permute([18,20,6,4,54],1)
 
 l,qq=lq_parity1(A) 
 
 l.setLabel([18,80])
 qq.setLabel([80,20,6,4,54])
 



 A=copy.copy(a_up)
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

 A=copy.copy(b_up)
 A.setLabel([55,18,20,6,4])
 A.permute([20,55,18,6,4],1)
 
 l,qq=lq_parity1(A) 
 
 l.setLabel([20,80])
 qq.setLabel([80,55,18,6,4])
 



 A=copy.copy(d_up)
 A.setLabel([53,19,10,8,20])
 A.permute([53,19,10,8,20],4)
 
 
 
 q,r=qr_parity1(A) 
 
 
 q.setLabel([53,19,10,8,81])
 r.setLabel([81,20])

 Teta=l*r
 Teta.permute([81,80],1)
 U,s,V=svd_parity2(Teta)

 U.setLabel([81,20])
 s.setLabel([20,-20])
 V.setLabel([-20,80])

 s=Sqrt(s)
 
 U=U*s
 V=s*V

 U.permute([81,-20],1)
 U.setLabel([81,20])
 V.permute([20,80],1)
 
 b_up=qq*V
 d_up=q*U
 b_up.permute([55,18,20,6,4],3)
 d_up.permute([53,19,10,8,20],3)



 
 
 return d_up,a_up,b_up





 
















def  Dis_final(E1, E2, E3, E4, E5, E6, E7, E8,c_u,a_u,b_u,c_up,a_up,b_up,d,MPO_list):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 MPO_list[0].setLabel([51,58,57,54])
 MPO_list[1].setLabel([58,57,52,59,60,55])
 MPO_list[2].setLabel([60,59,53,56])
 
 MPO_list0=[copy.copy(MPO_list[i]) for i in xrange(len(MPO_list))]
 MPO_list0[0].setLabel([-54,58,57,54])
 MPO_list0[1].setLabel([58,57,-55,59,60,55])
 MPO_list0[2].setLabel([60,59,-56,56])

 MPO_list1=[copy.copy(MPO_list[i]) for i in xrange(len(MPO_list))]
 MPO_list1[0].setLabel([-54,-58,-57,51])
 MPO_list1[1].setLabel([-58,-57,-55,-59,-60,52])
 MPO_list1[2].setLabel([-60,-59,-56,53])


 a_u1=copy.copy(a_u)
 b_u1=copy.copy(b_u)
 c_u1=copy.copy(c_u)

 a_d1=copy.copy(a_u)
 b_d1=copy.copy(b_u)
 c_d1=copy.copy(c_u)
 
###########################################################
 c_u1.setLabel([54,14,12,19,17])
 a_u1.setLabel([55,16,17,18,2])
 b_u1.setLabel([56,18,20,6,4])


 c_d1.setLabel([51,-14,-12,-19,-17])
 a_d1.setLabel([52,-16,-17,-18,-2])
 b_d1.setLabel([53,-18,-20,-6,-4])


 c_u1=((c_u1*(MPO_list0[0])))
 a_u1=(a_u1*(MPO_list0[1]))
 b_u1=((b_u1*(MPO_list0[2])))

 c_d1=((c_d1*(MPO_list1[0])))
 a_d1=(a_d1*(MPO_list1[1]))
 b_d1=((b_d1*(MPO_list1[2])))

########################################################
######################################################
 a_up.setLabel([-55,16,17,18,2])
 a_dp=copy.copy(a_up)
 a_dp.setLabel([-55,-16,-17,-18,-2])


 b_up.setLabel([-56,18,20,6,4])
 b_dp=copy.copy(b_up)
 b_dp.setLabel([-56,-18,-20,-6,-4])


 c_up.setLabel([-54,14,12,19,17])
 c_dp=copy.copy(c_up)
 c_dp.setLabel([-54,-14,-12,-19,-17])



 Val=((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E4*E5)*d)*((E2*E3)*(b_up*b_dp)))
 #print 'Val=',Val
 Val1=((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c_u1*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_d1)))
 #print 'Val1=',Val1
 Val2=((((E1*E8)*(a_up*a_d1))*((E7*E6)*(c_up*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_up*b_d1)))
 #print 'Val2=',Val2
 Val3=((((E1*E8)*(a_u1*a_dp))*((E7*E6)*(c_u1*c_dp))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dp)))
 #print 'Val3=',Val3
 
 val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=abs(Val[0])-2.0*Val2[0]#-Val3[0]

 return  val_f 





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
def sqrt_general(N2):
  N_init=copy.copy(N2)
  blk_qnums = N2.blockQnum()
  for qnum in blk_qnums:
   M=N2.getBlock(qnum)
   eig=M.eigh()
   
   e=Sqrt_mat(eig[0])
   U_trans=copy.copy(eig[1])
   U_trans.transpose()
   M=U_trans*e*eig[1]
   N_init.putBlock(qnum,M)
  return N_init
def Sqrt_mat(e):
 d=int(e.row())
 
 for q in xrange(d):
   #print e[q] 
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e  
def N_Positiv(N):
 N.setLabel([80,81,17,-82,-83,-18,-80,-81,-17,82,83,18])
 N.permute([-17,-80,-81,-82,-83,-18,17,80,81,82,83,18 ], 6)
 N1=copy.copy(N)
 N1.transpose()
 N=(N+N1)*(1.00/2.00)
 N1=copy.copy(N)
 N1.setLabel( [17,80,81,82,83,18,0,1,2,3,4,5 ] )
 N=N*N1
 N.permute([-17,-80,-81,-82,-83,-18,0,1,2,3,4,5],6)
 N_final=sqrt_general(N)
 N_final.setLabel([-17,-80,-81,-82,-83,-18,17,80,81,82,83,18])
 N_final.permute([80,81,17,-82,-83,-18,-80,-81,-17,82,83,18 ], 6)
 return N_final              
def N_Positiv_1(N):
 N.setLabel([-80,-81,-20,-82,-83,-18,80,81,20,82,83,18])
 N.permute([-80,-81,-20,-82,-83,-18,80,81,20,82,83,18], 6)
 N1=copy.copy(N)
 N1.transpose()
 N=(N+N1)*(1.00/2.00)
 N1=copy.copy(N)
 N1.setLabel( [80,81,20,82,83,18,0,1,2,3,4,5 ] )
 N=N*N1
 N.permute([-80,-81,-20,-82,-83,-18,0,1,2,3,4,5],6)
 N_final=sqrt_general(N)
 N_final.setLabel([-80,-81,-20,-82,-83,-18,80,81,20,82,83,18])
 N_final.permute([-80,-81,-20,-82,-83,-18,80,81,20,82,83,18], 6)
 return N_final               
