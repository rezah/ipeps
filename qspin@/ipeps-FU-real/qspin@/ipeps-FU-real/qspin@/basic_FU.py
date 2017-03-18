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
import basic


def  Do_optimization_Full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 
 Res=1
 Res1=2
 count=0
 Distance_val=Distance_val=Distance(l, r, lp, rp ,N_uni,U)
 for q in xrange(20):
  #print 'Dis', Distance_val[0], abs(Res1-Res) / abs(Res)
  rp, rp_d=r_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
  lp, lp_d=l_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
  Distance_val=Distance(l, r, lp, rp ,N_uni,U)
  Res=Res1
  Res1=Distance_val[0]
  count+=1
  if count > 50: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val[0]) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    #print 'break, Dis', Distance_val[0], (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break

 return rp, rp_d, lp, lp_d



def  Do_optimization_Grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 
  Es=Distance(l, r, lp, rp ,N_uni,U)
  Ef=[0]
  E2=[0.0]
  lp_first=copy.copy(lp)
  rp_first=copy.copy(rp)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=[0]
  count=0
  for i in xrange(60):
   count+=1
   E1=Distance(l, r, lp, rp ,N_uni,U)
   Ef=E1
   #print 'E=', E1[0], count
   #print lp, lp_d
   D_r, D_l=Obtain_grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U)
   D_r=(-1.00)*D_r
   D_l=(-1.00)*D_l
#   A=(D_r.norm()+D_l.norm())*(2.00)
#   D_r=D_r*(1.00/A)
#   D_l=D_l*(1.00/A)
   
   A=D_r*D_r
   B=D_l*D_l
   
   Norm_Z=A[0]+B[0]
   
   #print 'Norm', Norm_Z
   if (E1[0]<E_previous[0]) or (i is 0):
    if (abs(E1[0]) > 1.0e-10):
     if abs((E_previous[0]-E1[0])/E1[0]) < 1.0e-11:
      print 'Differnance Satisfied!', E_previous[0], E1[0], abs((E_previous[0]-E1[0])/E1[0]), i
      break
     else: 
      if abs((E_previous[0]-E1[0])) < 1.0e-11:
       print 'Differnance Satisfied!', E_previous[0], E1[0], abs((E_previous[0]-E1[0])), i
       break
      
   E_previous[0]=E1[0]
   
   if Norm_Z < 1.0e-8:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    count+=1
    rp_tem=rp+(2.00)*Gamma*D_r
    lp_tem=lp+(2.00)*Gamma*D_l
    E2=Distance(l, r, lp_tem, rp_tem ,N_uni,U)
    if E1[0]-E2[0] >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
     
   
   
   Break_loop=1
   while Break_loop is 1:
    count+=1
    rp_tem=rp+(1.00)*Gamma*D_r
    lp_tem=lp+(1.00)*Gamma*D_l
    E2=Distance(l, r, lp_tem, rp_tem ,N_uni,U)
    if abs((0.5)*Norm_Z*Gamma) <1.0e-15 or  abs(E1[0]-E2[0])<1.0e-15 :
     rp=rp+(1.00)*Gamma*D_r
     lp=lp+(1.00)*Gamma*D_l
     rp_d=copy.copy(rp)
     lp_d=copy.copy(lp)
     rp_d.setLabel([-3,-20,-1])
     lp_d.setLabel([-2,-40,-3])
     break
     
    if E1[0]-E2[0] < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   rp=rp+(1.00)*Gamma*D_r
   lp=lp+(1.00)*Gamma*D_l
   rp_d=copy.copy(rp)
   lp_d=copy.copy(lp)
   rp_d.setLabel([-3,-20,-1,-500])
   lp_d.setLabel([-2,-600,-40,-3])


  if(Ef[0] > Es[0]):
   print 'SD method, Fail, f<s', Ef[0], Es[0] 
   lp_first=copy.copy(lp_first)
   rp_first=copy.copy(rp_first)


  lp_d=copy.copy(lp)
  rp_d=copy.copy(rp)

  rp_d.setLabel([-3,-20,-1,-500])
  lp_d.setLabel([-2,-600,-40,-3])


  return rp, rp_d, lp, lp_d


def Obtain_grad(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-20,-40,20,40])



 A2=(((lp*lp_d*rp)*N_uni)*Iden)
 A2.permute([-3,-20,-1,-500],0)
 D_r=copy.copy(A2)
 A2=(((lp*lp_d*rp_d)*N_uni)*Iden)
 A2.permute([3,20,1,500],0)
 D_r=D_r+A2
 A3=((r)*U*(l*lp_d))*N_uni
 A3.permute([-3,-20,-1,-500],0)
 D_r=D_r+(-1.00)*A3
 A3p=((r_d)*U*(l_d*lp))*N_uni
 A3p.permute([3,20,1,500],0)
 D_r=D_r+(-1.00)*A3p
 D_r.setLabel([3,20,1,500])
 D_r.permute([3,20,1,500],2)


 A2=(((rp*rp_d*lp)*N_uni)*Iden)
 A2.permute([-2,-600,-40,-3],0)
 D_l=copy.copy(A2)
 A2=(((rp*rp_d*lp_d)*N_uni)*Iden)
 A2.permute([2,600,40,3],0)
 D_l=A2+D_l

 A3=((l)*U*(r*rp_d))*N_uni
 A3.permute([-2,-600,-40,-3],0)
 D_l=D_l+(-1.00)*A3
 A3p=((l_d)*U*(rp*r_d))*N_uni
 A3p.permute([2,600,40,3],0)
 D_l=D_l+(-1.00)*A3p

 D_l.setLabel([2,600,40,3])
 D_l.permute([2,600,40,3],3)


 return D_r, D_l
def Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):

 D_r=[0]*4
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


##########################################################
 a_u.setLabel([55,16,64,66,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([52,-16,-64,-66,-2])

 b_u.setLabel([56,68,20,6,4])
 b_d=copy.copy(b_u)
 b_d.setLabel([53,-68,-20,-6,-4])

 c_u.setLabel([54,14,12,19,62])
 c_d=copy.copy(c_u)
 c_d.setLabel([51,-14,-12,-19,-62])
##################################  1  #####################################################
 c_ut=((c_u*(MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

 c_dt=((c_d*(MPO_list1[0])))
 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*(((E1*E8)*(a_ut*a_dt))))*((E7*E6)*(c_ut*c_dt))

 A.permute([-62,-58,-57,-17,62,58,57,17],4)

 xt=copy.copy(plistd[0])
 xt.permute([-62,-58,-57,-17],4)
 D_r[0]=xt*A
 D_r[0].permute([62,58,57,17],0)

 x=copy.copy(plist[0])
 x.permute([62,58,57,17],0)
 A.transpose()
 A=x*A
 A.permute([-62,-58,-57,-17],0)
 D_r[0]=D_r[0]+A
##########################################################################################
 
 c_dt=((c_d*(MPO_list1[0])))
 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 Ap=((((E1*E8)*(a_u1*a_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-58,-57,-17],4)
 Ap.transpose()
 D_r[0]=D_r[0]+(-1.00)*Ap

 c_ut=((c_u*(MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))


 Ap=((((E1*E8)*(a_ut*a_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1)))*(((E7*E6)*(c_ut*c_d1)))
 Ap.permute([62,58,57,17],0)
 D_r[0]=D_r[0]+(-1.00)*Ap
 D_r[0].permute([62,58,57,17],3)
############################################################################################

##################################   2   #####################################################
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 a_ut=(a_u*(plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 a_dt=(a_d*(plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))
 A.permute([-17,-64,-58,-57,17,64,58,57],4)

 xt=copy.copy(plistd[1])
 xt.permute([-17,-64,-58,-57],4)
 D_r[1]=xt*A
 D_r[1].permute([17,64,58,57],0)

 x=copy.copy(plist[1])
 x.permute([17,64,58,57],0)
 A.transpose()
 A=x*A
 A.permute([-17,-64,-58,-57],0)
 D_r[1]=D_r[1]+A
##########################################################################################
 

 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 a_dt=(a_d*(plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-17,-64,-58,-57],4)
 Ap.transpose()
 D_r[1]=D_r[1]+(-1.00)*Ap

 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 a_ut=(a_u*(plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))


 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([17,64,58,57],0)
 D_r[1]=D_r[1]+(-1.00)*Ap
 D_r[1].permute([17,64,58,57],1)
############################################################################################


##################################  3  #####################################################
 a_ut=(a_u*(plist[1]*MPO_list0[1]))
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

 a_dt=(a_d*(plistd[1]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))
 A.permute([-66,-59,-60,-18,66,59,60,18],4)

 xt=copy.copy(plistd[2])
 xt.permute([-66,-59,-60,-18],4)
 D_r[2]=xt*A
 D_r[2].permute([66,59,60,18],0)

 x=copy.copy(plist[2])
 x.permute([66,59,60,18],0)
 A.transpose()
 A=x*A
 A.permute([-66,-59,-60,-18],0)
 D_r[2]=D_r[2]+A
##########################################################################################
 

 a_dt=(a_d*(plistd[1]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-66,-59,-60,-18],4)
 Ap.transpose()
 D_r[2]=D_r[2]+(-1.00)*Ap

 a_ut=(a_u*(plist[1]*MPO_list0[1]))
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))


 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([66,59,60,18],0)
 D_r[2]=D_r[2]+(-1.00)*Ap
 D_r[2].permute([66,59,60,18],3)
############################################################################################

##################################  4  #####################################################
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(MPO_list0[2])))

 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(MPO_list1[2])))


 A=((((E1*E8)*(a_ut*a_dt)) *((E7*E6)*(c_ut*c_dt))) * ((((E4*E5)*d)))) * ((E2*E3)*(b_ut*b_dt))

 A.permute([-18,-68,-59,-60,18,68,59,60],4)

 xt=copy.copy(plistd[3])
 xt.permute([-18,-68,-59,-60],4)
 D_r[3]=xt*A
 D_r[3].permute([18,68,59,60],0)

 x=copy.copy(plist[3])
 x.permute([18,68,59,60],0)
 A.transpose()
 A=x*A
 A.permute([-18,-68,-59,-60],0)
 D_r[3]=D_r[3]+A
##########################################################################################
 

 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(MPO_list1[2])))

 Ap=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-68,-59,-60],4)
 Ap.transpose()
 D_r[3]=D_r[3]+(-1.00)*Ap

 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 b_ut=((b_u*(MPO_list0[2])))

 
 Ap=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)))*((E2*E3)*(b_ut*b_d1))
 Ap.permute([18,68,59,60],0)
 D_r[3]=D_r[3]+(-1.00)*Ap
 D_r[3].permute([18,68,59,60],1)
############################################################################################


 return D_r 

def Obtain_grad_four1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,a_u,c_u,d_u,plist, plistd):

 D_r=[0]*4
 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 U.setLabel([51,52,53,54,55,56])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54,55,56])
 Iden1=copy.copy(Iden)
 Iden1.setLabel([51,52,53,-54,-55,-56])



 MPO_list[0].setLabel([51,58,57,54])
 MPO_list[1].setLabel([58,57,52,59,60,55])
 MPO_list[2].setLabel([60,59,53,56])
 
# H=MPO_list[0]*MPO_list[1]*MPO_list[2]
# H.permute([51,52,53,54,55,56],3)
# print "test-MPO", H.elemCmp(U) 
  
 
 
 MPO_list0=[copy.copy(MPO_list[i]) for i in xrange(len(MPO_list))]
 MPO_list0[0].setLabel([-54,58,57,54])
 MPO_list0[1].setLabel([58,57,-55,59,60,55])
 MPO_list0[2].setLabel([60,59,-56,56])

 MPO_list1=[copy.copy(MPO_list[i]) for i in xrange(len(MPO_list))]
 MPO_list1[0].setLabel([-54,-58,-57,51])
 MPO_list1[1].setLabel([-58,-57,-55,-59,-60,52])
 MPO_list1[2].setLabel([-60,-59,-56,53])


 a_u.setLabel([54,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([54,-16,-17,-18,-2])

 c_u.setLabel([55,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.setLabel([55,-14,-12,-19,-17])


 d_u.setLabel([56,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.setLabel([56,-19,-10,-8,-20])
 

 a_u0=copy.copy(a_u)
 d_u0=copy.copy(d_u)
 c_u0=copy.copy(c_u)
 
 a_d1=copy.copy(a_u)
 d_d1=copy.copy(d_u)
 c_d1=copy.copy(c_u)
 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)
###########################################################
 c_d1.setLabel([-55,-14,-12,-19,-17])
 a_d1.setLabel([-54,-16,-17,-18,-2])
 d_d1.setLabel([-56,-19,-10,-8,-20])

 c_u1.setLabel([-55,14,12,19,17])
 a_u1.setLabel([-54,16,17,18,2])
 d_u1.setLabel([-56,19,10,8,20])
##########################################################
 a_u0.setLabel([54,16,68,18,2])
 a_d0=copy.copy(a_u0)
 a_d0.setLabel([51,-16,-68,-18,-2])

 c_u0.setLabel([55,14,12,62,66])
 c_d0=copy.copy(c_u0)
 c_d0.setLabel([52,-14,-12,-62,-66])


 d_u0.setLabel([56,64,10,8,20])
 d_d0=copy.copy(d_u0)
 d_d0.setLabel([53,-64,-10,-8,-20])


##################################  1  #####################################################
 a_ut=((a_u0*(plist[3]*MPO_list0[0])))
 c_ut=(c_u0*(plist[2]*MPO_list0[1]))
 d_ut=((d_u0*(plist[1]*MPO_list0[2])))

 a_dt=((a_d0*(plistd[3]*MPO_list1[0])))
 c_dt=(c_d0*(plistd[2]*MPO_list1[1]))
 d_dt=((d_d0*(plistd[1]*MPO_list1[2])))

 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))*(((E1*E8)*(a_ut*a_dt))))*((E7*E6)*(c_ut*c_dt))


 A.permute([-62,-59,-60,-19,62,59,60,19],4)

 xt=copy.copy(plistd[0])
 xt.permute([-62,-59,-60,-19],4)
 D_r[0]=xt*A
 D_r[0].permute([62,59,60,19],0)

 x=copy.copy(plist[0])
 x.permute([62,59,60,19],0)
 A.transpose()
 A=x*A
 A.permute([-62,-59,-60,-19],0)
 D_r[0]=D_r[0]+A
##########################################################################################
 
 a_dt=((a_d0*(plistd[3]*MPO_list1[0])))
 c_dt=(c_d0*(plistd[2]*MPO_list1[1]))
 d_dt=((d_d0*(plistd[1]*MPO_list1[2])))

 Ap=((((E1*E8)*(a_u1*a_dt))))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-59,-60,-19],4)
 Ap.transpose()
 D_r[0]=D_r[0]+(-1.00)*Ap

 a_ut=((a_u0*(plist[3]*MPO_list0[0])))
 c_ut=(c_u0*(plist[2]*MPO_list0[1]))
 d_ut=((d_u0*(plist[1]*MPO_list0[2])))


 Ap=((((E1*E8)*(a_ut*a_d1))))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b)))*(((E7*E6)*(c_ut*c_d1)))
 Ap.permute([62,59,60,19],0)

 D_r[0]=D_r[0]+(-1.00)*Ap
 D_r[0].permute([62,59,60,19],3)
############################################################################################

##################################   2   #####################################################
 a_ut=((a_u0*(plist[3]*MPO_list0[0])))
 c_ut=(c_u0*(plist[2]*plist[0]*MPO_list0[1]))
 d_ut=((d_u0*(MPO_list0[2])))

 a_dt=((a_d0*(plistd[3]*MPO_list1[0])))
 c_dt=(c_d0*(plistd[2]*plistd[0]*MPO_list1[1]))
 d_dt=((d_d0*(MPO_list1[2])))


 A=((((E1*E8)*(a_ut*a_dt))*((E2*E3)*(b)))*((E7*E6)*(c_ut*c_dt)))*((E4*E5)*(d_ut*d_dt))

 A.permute([-19,-59,-60,-64,19,59,60,64],4)

 xt=copy.copy(plistd[1])
 xt.permute([-19,-59,-60,-64],4)
 D_r[1]=xt*A
 D_r[1].permute([19,59,60,64],0)

 x=copy.copy(plist[1])
 x.permute([19,59,60,64],0)
 A.transpose()
 A=x*A
 A.permute([-19,-59,-60,-64],0)
 D_r[1]=D_r[1]+A
##########################################################################################
 
 a_dt=((a_d0*(plistd[3]*MPO_list1[0])))
 c_dt=(c_d0*(plistd[2]*plistd[0]*MPO_list1[1]))
 d_dt=((d_d0*(MPO_list1[2])))

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E1*E8)*(a_u1*a_dt))*((E2*E3)*(b))))*((E4*E5)*(d_u1*d_dt))
 Ap.permute([-19,-59,-60,-64],4)
 Ap.transpose()
 D_r[1]=D_r[1]+(-1.00)*Ap

 a_ut=((a_u0*(plist[3]*MPO_list0[0])))
 c_ut=(c_u0*(plist[2]*plist[0]*MPO_list0[1]))
 d_ut=((d_u0*(MPO_list0[2])))


 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E1*E8)*(a_ut*a_d1))*((E2*E3)*(b))))*((E4*E5)*(d_ut*d_d1))
 Ap.permute([19,59,60,64],0)

 D_r[1]=D_r[1]+(-1.00)*Ap
 D_r[1].permute([19,59,60,64],1)
############################################################################################


##################################  3  #####################################################
 a_ut=((a_u0*(plist[3]*MPO_list0[0])))
 c_ut=(c_u0*(plist[0]*MPO_list0[1]))
 d_ut=((d_u0*(plist[1]*MPO_list0[2])))

 a_dt=((a_d0*(plistd[3]*MPO_list1[0])))
 c_dt=(c_d0*(plistd[0]*MPO_list1[1]))
 d_dt=((d_d0*(plistd[1]*MPO_list1[2])))

 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))*((E1*E8)*(a_ut*a_dt)))*((E7*E6)*(c_ut*c_dt))

 A.permute([-66,-58,-57,-17,66,58,57,17],4)

 xt=copy.copy(plistd[2])
 xt.permute([-66,-58,-57,-17],4)
 D_r[2]=xt*A
 D_r[2].permute([66,58,57,17],0)

 x=copy.copy(plist[2])
 x.permute([66,58,57,17],0)
 A.transpose()
 A=x*A
 A.permute([-66,-58,-57,-17],0)
 D_r[2]=D_r[2]+A
##########################################################################################
 
 a_dt=((a_d0*(plistd[3]*MPO_list1[0])))
 c_dt=(c_d0*(plistd[0]*MPO_list1[1]))
 d_dt=((d_d0*(plistd[1]*MPO_list1[2])))

 Ap=((((E1*E8)*(a_u1*a_dt)))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b))))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-66,-58,-57,-17],4)
 Ap.transpose()
 D_r[2]=D_r[2]+(-1.00)*Ap

 a_ut=((a_u0*(plist[3]*MPO_list0[0])))
 c_ut=(c_u0*(plist[0]*MPO_list0[1]))
 d_ut=((d_u0*(plist[1]*MPO_list0[2])))

 Ap=((((E1*E8)*(a_ut*a_d1)))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b))))*(((E7*E6)*(c_ut*c_d1)))
 Ap.permute([66,58,57,17],0)
 D_r[2]=D_r[2]+(-1.00)*Ap
 D_r[2].permute([66,58,57,17],3)
############################################################################################

##################################  4  #####################################################
 a_ut=((a_u0*(MPO_list0[0])))
 c_ut=(c_u0*(plist[0]*plist[2]*MPO_list0[1]))
 d_ut=((d_u0*(plist[1]*MPO_list0[2])))

 a_dt=((a_d0*(MPO_list1[0])))
 c_dt=(c_d0*(plistd[0]*plistd[2]*MPO_list1[1]))
 d_dt=((d_d0*(plistd[1]*MPO_list1[2])))

 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-17,-58,-57,-68,17,58,57,68],4)

 xt=copy.copy(plistd[3])
 xt.permute([-17,-58,-57,-68],4)
 D_r[3]=xt*A
 D_r[3].permute([17,58,57,68],0)

 x=copy.copy(plist[3])
 x.permute([17,58,57,68],0)
 A.transpose()
 A=x*A
 A.permute([-17,-58,-57,-68],0)
 D_r[3]=D_r[3]+A
##########################################################################################
 
 a_dt=((a_d0*(MPO_list1[0])))
 c_dt=(c_d0*(plistd[0]*plistd[2]*MPO_list1[1]))
 d_dt=((d_d0*(plistd[1]*MPO_list1[2])))

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b))))*((E1*E8)*(a_u1*a_dt))
 
 Ap.permute([-17,-58,-57,-68],4)
 Ap.transpose()
 D_r[3]=D_r[3]+(-1.00)*Ap

 a_ut=((a_u0*(MPO_list0[0])))
 c_ut=(c_u0*(plist[0]*plist[2]*MPO_list0[1]))
 d_ut=((d_u0*(plist[1]*MPO_list0[2])))

 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([17,58,57,68],0)
 D_r[3]=D_r[3]+(-1.00)*Ap
 D_r[3].permute([17,58,57,68],1)
############################################################################################


 return D_r 





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








def  r_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 #H1.transpose()
 #H1.setLabel([-20,-40,30,50])
 #print H1, U
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-20,-40,20,40])
 A2=(((lp*lp_d)*N_uni)*Iden)
 A2.permute([-3,-20,-1,-500,3,20,1,500],4)
 
 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 #A2=addition_symmetric(A2,A2_trans)
 
 A2.setLabel([-3,-20,-1,-500,3,20,1,500])
 
 A3=((r)*U*(l*lp_d))*N_uni
 A3.permute([-3,-20,-1,-500],0)
 A3p=((r_d)*U*(l_d*lp))*N_uni
 A3p.permute([3,20,1,500],4)
 A3p.transpose()
 #A3=addition_symmetric(A3,A3p)

 A3=A3+A3p
 
 A3.setLabel([-3,-20,-1,-500])
 
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
 A2_inv.setLabel([3,20,1,500,-3,-20,-1,-500])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([3,20,1,500])
 A.permute([3,20,1,500],2)

 rf=copy.copy(A)
 rf_d=copy.copy(rf)
 rf_d.transpose()
 rf_d.setLabel([-1,-500,-3,-20])

 return rf, rf_d

def  l_optimum_full(l, r, l_d, r_d, lp, rp, lp_d, rp_d ,N_uni,U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 U1=copy.copy(U)

 H1.transpose()
 H1.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-20,-40,20,40])

 A2=(((rp*rp_d)*N_uni)*Iden)
 A2.permute([-2,-600,-40,-3,2,600,40,3],4)
 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 A2.setLabel([-2,-600,-40,-3,2,600,40,3])


 A3=((l)*H1*(r*rp_d))*N_uni
 A3.permute([-2,-600,-40,-3],0)
 A3p=((l_d)*U*(rp*r_d))*N_uni
 A3p.permute([2,600,40,3],4)
 A3p.transpose()
 A3=A3+A3p
 A3.setLabel([-2,-600,-40,-3])
 
 
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
 A2_inv.setLabel([2,600,40,3,-2,-600,-40,-3])




 A=A3*A2_inv
 A.setLabel([2,600,40,3])
 A.permute([2,600,40,3],3)


 lf=copy.copy(A)
 lf_d=copy.copy(lf)
 lf_d.transpose()
 lf_d.setLabel([-3,-2,-600,-40])

 return lf, lf_d





def putBlock_Total_trans(E1,E2):
    blk_qnums = E1.blockQnum()
    for qnum in blk_qnums:
        A=E1.getBlock(qnum)
        A.transpose()
        E2.putBlock(qnum,A)
    return E2



def inv(svd):
 bdi=uni10.Bond(uni10.BD_IN, svd.col())
 bdo=uni10.Bond(uni10.BD_OUT, svd.row())
 Landa=uni10.UniTensor([bdi,bdo])
 Landa.putBlock(svd)
 Landa=basic.inverse(Landa)
 return Landa
 
def distance_iden(A2,A2_inv):
 A=A2*A2_inv
 d=int(A.row())
 matrixI=uni10.Matrix(d, d)
 for i in xrange(d):
  for j in xrange(d):
   if(i==j):
    matrixI[i*d+j]=1.0
   else:
    matrixI[i*d+j]=0.0
 sum=0 
 for i in xrange(d):
  for j in xrange(d):
   sum=abs(A[i*d+j]+(-1.00)*matrixI[i*d+j])
 return sum
 
 
def  recover(l, r, q,qq):
 a_u=q*r
 a_u.permute([20,4,5,3,2],3)
 a_u.setLabel([0,1,2,3,4])
 #print a_ut.similar(a_u)#,a_ut.printDiagram() 
 
 b_u=qq*l
 b_u.permute([40,3,4,5,6],3)
 b_u.setLabel([0,1,2,3,4])
 #print b_ut.similar(b_u) 
 return a_u, b_u





def  test_energy_lr(N_uni, l, l_d, r, r_d, q,qq,U,E1, E2, E3, E4, E5,E6,a_u,b_u):
 U.setLabel([-20,-40,20,40])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-20,-40,20,40])

 
 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-2,-3,-20,-13,-1])

 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-5,-6,-40,-2,-4])

 B=(((r*r_d)*U)*(l*l_d))*(N_uni)
 A=(((r*r_d)*Iden)*(l*l_d))*(N_uni)
 print 'E=', B[0]/A[0]
 print 'NormQR=', A[0]
 
 
 a_ut=q*r
 a_ut.permute([20,4,5,3,2],3)
 print a_ut.elemCmp(a_u) 
 
 b_ut=qq*l
 b_ut.permute([40,3,4,5,6],3)
 
 print b_ut.elemCmp(b_u)






def  test_energy(E1, E2, E3, E4, E5,E6, a, b, c1,c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4, a_u,b_u,U):


 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-2,-3,20,-13,-1])
 #print a_u.printDiagram()
 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-5,-6,40,-2,-4])
 A=(((E2*(a_u*a_d))*E1)*E3)*(((((b_u*b_d)*E5)*E6)*E4))
 Norm_f=A
 print 'Norm=', A[0]

 a_u.setLabel([20,13,1,2,3])
 a_d.setLabel([-2,-3,-20,-13,-1])

 b_u.setLabel([40,2,4,5,6])
 b_d.setLabel([-5,-6,-40,-2,-4,])


 U.setLabel([-20,-40,20,40])
 #print U.printDiagram(), a_u.printDiagram(),a_d.printDiagram(), b_d.printDiagram(),b_u.printDiagram()

 #print a_u.printDiagram(), b_u.printDiagram() 
 B=((((E2*(a_u*a_d))*E1)*E3)*U)*(((((b_u*b_d)*E5)*E6)*E4))
 print 'E=', B[0]/A[0]


def test_env(E1, E2, E3, E4, E5,E6, a, b, c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 a.setLabel([13,-13,1,-1,2,-2,3,-3])
 b.setLabel([2,-2,4,-4,5,-5,6,-6])
 A=(((E2*a)*E1)*E3)*((((b*E5)*E6)*E4))
 print 'norm=', A[0]


def Distance(l, r, lp, rp ,N_uni, U):
 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,30,50])
  
 U.setLabel([30,50,20,40])
 H=U*H1
 H.permute([-20,-40,20,40],2)
 #H.setLabel([-20,-40,20,40]) 
 
 H1.setLabel([-20,-40,20,40]) 
 U.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-20,-40,20,40])
 
 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)
 rp_d.transpose()
 lp_d.transpose()
 
 rp_d.setLabel([-1,-500,-3,-20])
 lp_d.setLabel([-3,-2,-600,-40])

 
 l_d=copy.copy(l)
 r_d=copy.copy(r)
 r_d.transpose()
 l_d.transpose()

 r_d.setLabel([-1,-500,-3,-20])
 l_d.setLabel([-3,-2,-600,-40])
 
 
 
 A1=((((r*r_d)*H)*(N_uni))*l*l_d)
 A2=((((rp*rp_d)*Iden)*(N_uni))*lp*lp_d)
 A3=((((r*rp_d)*U)*(N_uni))*l*lp_d)
 A4=((((rp*r_d)*H1)*(N_uni))*lp*l_d)
 A=A1+A2+(-1.00)*A3+(-1.00)*A4
 return A

def initialize_lrprime(l, r, l_d, r_d, N_uni):
 lp=copy.copy(l)
 rp=copy.copy(r)
 lp_d=copy.copy(l_d)
 rp_d=copy.copy(r_d)

# lp_d.transpose()
# rp_d.transpose()
# 
# rp_d.setLabel([-1,-500,-3,-20])
# lp_d.setLabel([-3,-2,-600,-40])
  
 return  lp, rp, lp_d, rp_d

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

def putBlock_Total(E1,E2):
    blk_qnums = E1.blockQnum()
    for qnum in blk_qnums:
        A=E1.getBlock(qnum)
        E2.putBlock(qnum,A)
    return E2


def make_XtranX(N_final):
  blk_qnums = N_final.blockQnum()
  X_tran_u=uni10.UniTensor(N_final.bond())
  X_u=uni10.UniTensor(N_final.bond())

  for qnum in blk_qnums:
   eig=N_final.getBlock(qnum)
   eig=eig.eigh() 
   e=Positiv(eig[0])
   for q in xrange(int(eig[0].row())):
      e[q]=(e[q]**(1.00/2.00)) 
   X=e*eig[1]
   eig[1].transpose()
   X_tran=eig[1]*e
   X_tran_u.putBlock(qnum,X_tran)
   X_u.putBlock(qnum,X)
   
#  X_u_tem=copy.copy(X_tran_u)
#  X_u_tem.transpose()
#  X_u=putBlock_Total(X_u_tem,X_u)
  return X_u, X_tran_u


def reproduceX(X_u,X_tran_u,qq_u,q_u):
 bdiB=X_tran_u.bond(0)
 bdoB=X_u.bond(1)
 
 X_u.setLabel([0,1,2,3,4,5,6,7])

 X_u.permute([4,5,0,1,2,3,6,7],2)
# bdout=X_u.bond(0)
# bdout1=X_u.bond(1)
# 
# bdout.change(uni10.BD_OUT)
# bdout1.change(uni10.BD_OUT)

 bdout=uni10.Bond(uni10.BD_OUT,X_u.bond(0).Qlist())
 bdout1=uni10.Bond(uni10.BD_OUT,X_u.bond(1).Qlist())
 
 l1_uni=uni10.UniTensor([X_u.bond(0),X_u.bond(1),bdout,bdout1])
 l1_uni_inv=uni10.UniTensor([X_u.bond(0),X_u.bond(1),bdout,bdout1])
 
 blk_qnums = X_u.blockQnum()
 for qnum in blk_qnums:
  lq1=X_u.getBlock(qnum).lq()
  l1=lq1[0]
  l1_uni.putBlock(qnum,l1)
  l1_inv=l1.inverse()
  l1_uni_inv.putBlock(qnum,l1_inv)

 
 X_u.permute([0,1,2,3,4,5,6,7],6)

# bdin=X_u.bond(6)
# bdin1=X_u.bond(7)
# 
# bdin.change(uni10.BD_IN)
# bdin1.change(uni10.BD_IN)
 bdin=uni10.Bond(uni10.BD_IN,X_u.bond(6).Qlist())
 bdin1=uni10.Bond(uni10.BD_IN,X_u.bond(7).Qlist())


 r1_uni_inv=uni10.UniTensor([bdin,bdin1,X_u.bond(6),X_u.bond(7)])
 r1_uni=uni10.UniTensor([bdin,bdin1,X_u.bond(6),X_u.bond(7)])

 blk_qnums = X_u.blockQnum()
 for qnum in blk_qnums:
  qr1=X_u.getBlock(qnum).qr()
  r1=qr1[1]
  r1_uni.putBlock(qnum,r1)
  r1_inv=r1.inverse()
  r1_uni_inv.putBlock(qnum,r1_inv)
 
 
 r1_uni_inv.setLabel([2,600,-2,-600])

 qq_u=qq_u*r1_uni_inv
 qq_u.permute([-2,-600,4,5,6],2)
 qq_u.setLabel([2,600,4,5,6])

 
 
 l1_uni_inv.setLabel([-1,-500,1,500])
 q_u=q_u*l1_uni_inv
 q_u.permute([2,4,5,-1,-500],3)
 q_u.setLabel([2,4,5,1,500])
 
 X_u.permute([0,1,2,3,4,5,6,7],4)
 X_u.setLabel([-10,-11,-12,-13,1,500,2,600])
 
 X_u=X_u*l1_uni_inv*r1_uni_inv
 X_u.permute([-10,-11,-12,-13,-1,-500,-2,-600],4)
 X_u.setLabel([-10,-11,-12,-13,1,500,2,600])
 return X_u, q_u, qq_u,r1_uni,l1_uni
def addition_symmetric(N,N1):
 N2=copy.copy(N)
 N_final=copy.copy(N)

 blk_qnums = N2.blockQnum()
 for qnum in blk_qnums:
  M2=N2.getBlock(qnum)
  M1=N1.getBlock(qnum)
  M=(M1+M2)*(0.5)
  N_final.putBlock(qnum,M)

 return N_final

def initialize_Positiv_lrprime(l, r, l_d, r_d, N_uni, U1, D, d_phys, q_u, qq_u,a_u,b_u,Positive):
 U1.setLabel([-20,-40,20,40])
 N=copy.copy(N_uni)
 N.setLabel([1,500,-1,-500,2,600,-2,-600])
 #print N.printDiagram()
 N.permute([-1,-500,-2,-600,1,500,2,600], 4)
 N1=copy.copy(N)
 N1.transpose()
 #print N.printDiagram(),N1.printDiagram(),N,N1
 N=(N+N1)*(1.00/2.00)
 #N=addition_symmetric(N,N1)
 if Positive is 'Restrict':
  N1=copy.copy(N)
  N.setLabel([-1,-500,-2,-600,1,500,2,600] )
  N1.setLabel([1,500,2,600,10,11,12,13] )
  N2=N*N1
  N2.permute([-1,-500,-2,-600,10,11,12,13],4)
  N_final=sqrt_general(N2)
  
 X_u, X_tran_u=make_XtranX(N_final)

#################################################


 X_u, q_u, qq_u,r1_uni,l1_uni=reproduceX(X_u,X_tran_u,qq_u,q_u)

######################################################
 X_u.setLabel([-10,-11,-12,-13,1,500,2,600])

 X_tran_u=copy.copy(X_u)
 X_tran_u.transpose()
 X_tran_u.setLabel([-1,-500,-2,-600,-10,-11,-12,-13])
 
# X_u_tem=copy.copy(X_u)
# X_u_tem.transpose()
# X_tran_u=putBlock_Total(X_u_tem,X_tran_u)
# X_tran_u.setLabel([-1,-500,-2,-600,-10,-11,-12,-13])

 N=X_tran_u*X_u

 N.permute([1,500,-1,-500,2,600,-2,-600],4)
 #print N.printDiagram()
 l1_uni.setLabel([1,500,-1,-500])
 r=r*l1_uni
 r.permute([3,20,-1,-500],2)
 r.setLabel([3,20,1,500])
 r_d=copy.copy(r)
 r_d.transpose()
 r_d.setLabel([-1,-500,-3,-20])
 
 
 r1_uni.setLabel([-2,-600,2,600])
 l=r1_uni*l
 l.permute([-2,-600,40,3],3)
 l.setLabel([2,600,40,3])
 l_d=copy.copy(l)
 l_d.transpose()
 l_d.setLabel([-3,-2,-600,-40])

 lp=copy.copy(l)
 rp=copy.copy(r)
 lp_d=copy.copy(l_d)
 rp_d=copy.copy(r_d)

 lp_d.transpose()
 lp_d.transpose()

 rp_d.setLabel([-1,-500,-3,-20])
 lp_d.setLabel([-3,-2,-600,-40])

#######################################################
# 
# D_dim=0
# for i in xrange(len(D)):
#  D_dim+=D[i]
# #print "D_dim", D_dim

# U1.setLabel([-20,-40,20,40])
# lp=copy.copy(l)
# rp=copy.copy(r)
# Teta=U1*lp*rp
# Teta.permute([1,500,-20,2,600,-40],3)
# U,s,V=svd_parity1(Teta,D_dim)
# 
# #print "s", s
# s=Sqrt(s)
# #print "sqrt(s)", s
# 
# U.setLabel([0,1,2,3])
# s.setLabel([3,6])
# V.setLabel([6,9,10,11])
# U=U*s
# V=s*V

# U.permute([0,1,2,6],3)
# V.permute([3,9,10,11],1)
# 
# U.setLabel([1,500,20,3])
# rp=copy.copy(U)
# rp.permute([3,20,1,500],2)
# 
# V.setLabel([3,2,600,40])  
# lp=copy.copy(V)
# lp.permute([2,600,40,3],3)

# lp_d=copy.copy(lp)
# rp_d=copy.copy(rp)

# lp_d.transpose()
# rp_d.transpose()
# 
# rp_d.setLabel([-1,-500,-3,-20])
# lp_d.setLabel([-3,-2,-600,-40])
 ###################################
 a_u=q_u*r
 a_u.permute([20,4,5,3,2],3)
 
 b_u=qq_u*l
 b_u.permute([40,3,4,5,6],3)
 
 return  lp, rp, lp_d, rp_d, N, l, r, l_d, r_d, q_u, qq_u, a_u, b_u 



def Positiv(e):
 d=int(e.row())
 if (e[0] > 0) and (e[d-1] > 0):
  return e

 if (e[0] < 0) and (e[d-1] < 0):
  return (-1.00)*e

 if abs(e[0]) >  abs(e[d-1]) :
   print 'e_small > e_large',  e, e[0],  e[d-1]
   e=(-1.00)*e


 for q in xrange(d):
  if e[q] < 0: e[q]=0;
 return e










def  distance_two(A,B ):
 d=int(A.row()*A.col())
 sum=0.0
 for q in xrange(d):
  sum+=abs(A[q]-B[q])
 return sum








def lq_parity1(theta):
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


def Equall_Dist(l, r,D,d_phys):

 lp=copy.copy(l)
 rp=copy.copy(r)
 rp.permute([3,20,1,500],1)
 #lq_r=rp.getBlock().lq()
 l1, q1=lq_parity1(rp)
 
 l1.setLabel([3,-1])
 q1.setLabel([-1,20,1,500])
 
 lp.permute([2,600,40,3],3)
 #qr_l=lp.getBlock().qr()
 qq1, r1=qr_parity1(lp)
 
 r1.setLabel([-2,3])
 qq1.setLabel([2,600,40,-2])
 
 
 

 Teta=l1*r1
 Teta.permute([-2,-1],1)
 U,s,V=svd_parity2(Teta)

 U.setLabel([0,1])
 s.setLabel([1,2])
 V.setLabel([2,3])

 s=Sqrt(s)
 U=U*s
 V=s*V

 U.permute([0,2],1)
 V.permute([1,3],1)
 U.setLabel([-2,3])
 V.setLabel([3,-1])



 lp=U*qq1
 rp=q1*V
 
 rp.permute([3,20,1,500],2)
 lp.permute([2,600,40,3],3)
 
  
 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)
 lp_d.transpose()
 rp_d.transpose()
 
 rp_d.setLabel([-1,-500,-3,-20])
 lp_d.setLabel([-3,-2,-600,-40])
 
 return  lp, rp, lp_d, rp_d


def svd_parity1(theta,chi):

    LA=uni10.UniTensor([theta.bond(0), theta.bond(3)])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1), theta.bond(2),theta.bond(3)])
    GB=uni10.UniTensor([theta.bond(0), theta.bond(3), theta.bond(4),theta.bond(5)])
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
#    #print dims, dims_val

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
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim) * (1.00/norm)  )
#    print LA
    return GA, LA,GB




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

def initialize_SVD_lrprime(l, r, l_d, r_d, N_uni,U1,D,d_phys):
 
 D_dim=0
 for i in xrange(len(D)):
  D_dim+=D[i]
 print "D_dim", D_dim

 U1.setLabel([-20,-40,20,40])
 lp=copy.copy(l)
 rp=copy.copy(r)
 Teta=U1*lp*rp
 Teta.permute([1,500,-20,2,600,-40],3)
 U,s,V=svd_parity1(Teta,D_dim)
 
 #print "s", s
 s=Sqrt(s)
 #print "sqrt(s)", s
 
 U.setLabel([0,1,2,3])
 s.setLabel([3,6])
 V.setLabel([6,9,10,11])
 U=U*s
 V=s*V

 U.permute([0,1,2,6],3)
 V.permute([3,9,10,11],1)
 
 U.setLabel([1,500,20,3])
 rp=copy.copy(U)
 rp.permute([3,20,1,500],2)
 
 V.setLabel([3,2,600,40])  
 lp=copy.copy(V)
 lp.permute([2,600,40,3],3)

 lp_d=copy.copy(lp)
 rp_d=copy.copy(rp)
 lp_d.transpose()
 rp_d.transpose()
 rp_d.setLabel([-1,-500,-3,-20])
 lp_d.setLabel([-3,-2,-600,-40])
 
 return  lp, rp, lp_d, rp_d


 
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

 
def Qr_lQ_decom(a_u,b_u, E1, E2, E3, E4, E5,E6,D,d_phys):
 a_uc=copy.copy(a_u)
 a_uc.setLabel([20,13,1,2,3])
 a_uc.permute([3,13,1,20,2],3)
 
 
 
 q_uni,r_uni=qr_parity(a_uc) 
 
 
 q_uni.setLabel([3,13,1,100,500])
 r_uni.setLabel([100,500,20,2])



 b_uc=copy.copy(b_u)
 b_uc.setLabel([40,2,4,5,6])
 b_uc.permute([40,2,4,5,6],2)



 l_uni,qq_uni=lq_parity(b_uc)

 l_uni.setLabel([40,2,200,600])
 qq_uni.setLabel([200,600,4,5,6])


 r_uni_d=copy.copy(r_uni)
 r_uni_d.transpose()
 r_uni_d.setLabel([-20,-2,-100,-500])

 l_uni_d=copy.copy(l_uni)
 l_uni_d.transpose()
 l_uni_d.setLabel([-200,-600,-40,-2])

 q_uni_d=copy.copy(q_uni)
 q_uni_d.transpose()
 q_uni_d.setLabel([-100,-500,-3,-13,-1])
 qq_uni_d=copy.copy(qq_uni)
 qq_uni_d.transpose()
 qq_uni_d.setLabel([-4,-5,-6,-200,-600])
 N=((((E1*q_uni)*E2)*q_uni_d)*E3)*((((E6*qq_uni)*E5)*qq_uni_d)*E4)
 N.permute([100,500,-100,-500,200,600,-200,-600],4)
 #print "hi", N.printDiagram()

 N.setLabel([1,500,-1,-500,2,600,-2,-600])
 #r_uni.setLabel([100,20,2])
 r_uni.setLabel([1,500,20,3])
 r_uni.permute([3,20,1,500],2)
 
 r_uni_d=copy.copy(r_uni)
 r_uni_d.transpose()
 r_uni_d.setLabel([-1,-500,-3,-20,])

 #l_uni.setLabel([40,2,200])
 l_uni.setLabel([40,3,2,600])
 l_uni.permute([2,600,40,3],3)
 
 l_uni_d=copy.copy(l_uni)
 l_uni_d.transpose()
 l_uni_d.setLabel([-3,-2,-600,-40])

 #q_uni.setLabel([3,13,1,100])
 q_uni.setLabel([2,4,5,1,500])
 #qq_uni.setLabel([200,4,5,6])
 qq_uni.setLabel([2,600,4,5,6])
 
 
 return   N, l_uni, l_uni_d, r_uni, r_uni_d, q_uni,qq_uni
 
 
def  proper_bond(E1, E2, E3, E4, E5,E6,D,d_phys):
 bd=uni10.Bond(E1.bond()[1].type(), D)
 A=uni10.UniTensor([E1.bond()[0],bd,bd,E1.bond()[2]])
 A.putBlock(E1.getBlock())
 E1=copy.copy(A)
 E1.setLabel([7,3,-3,8])
 
 bd=uni10.Bond(E2.bond()[1].type(), D)
 A=uni10.UniTensor([E2.bond()[0],bd,bd,E2.bond()[2]])
 A.putBlock(E2.getBlock())
 E2=copy.copy(A)
 E2.setLabel([12,13,-13,7])


 bd=uni10.Bond(E3.bond()[1].type(), D)
 A=uni10.UniTensor([E3.bond()[0],bd,bd,E3.bond()[2]])
 A.putBlock(E3.getBlock())
 E3=copy.copy(A)
 E3.setLabel([11,1,-1,12])


 bd=uni10.Bond(E4.bond()[1].type(), D)
 A=uni10.UniTensor([E4.bond()[0],bd,bd,E4.bond()[2]])
 A.putBlock(E4.getBlock())
 E4=copy.copy(A)
 E4.setLabel([11,4,-4,10])


 bd=uni10.Bond(E5.bond()[1].type(), D)
 A=uni10.UniTensor([E5.bond()[0],bd,bd,E5.bond()[2]])
 A.putBlock(E5.getBlock())
 E5=copy.copy(A)
 E5.setLabel([10,5,-5,9])

 bd=uni10.Bond(E6.bond()[1].type(), D)
 A=uni10.UniTensor([E6.bond()[0],bd,bd,E6.bond()[2]])
 A.putBlock(E6.getBlock())
 E6=copy.copy(A)
 E6.setLabel([9,6,-6,8])

 return E1, E2, E3, E4, E5,E6
 

def produce_Env_Hac(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys):

 c1.setLabel([0,1])
 Tb4.setLabel([3,2,-2,0])
 E2=c1*Tb4
 E2.permute([3,2,-2,1],3)
 E2.setLabel([10,2,-2,12])
 E2.setLabel([8,6,-6,9])
 E2.permute([8,6,-6,9],3)
 
 c4.setLabel([0,1])
 Ta3.setLabel([1,2,-2,3])
 E4=c4*Ta3
 E4.permute([0,2,-2,3],2)
 E4.setLabel([9,6,-6,8])
 E4.setLabel([7,13,-13,12])
 E4.permute([12,13,-13,7],1)
 
 E1=Tb1
 E1.setLabel([12,1,-1,13])
 E1.setLabel([9,5,-5,10])
 E1.permute([10,5,-5,9],3)

 E3=Ta4
 E3.setLabel([9,5,-5,10])
 E3.setLabel([7,3,-3,8])
 E3.permute([7,3,-3,8],3)
 
 
 c2.setLabel([5,6])
 Ta1.setLabel([21,20,-20,5])
 Ta2.setLabel([22,19,-19,6])
 b.setLabel([4,-4,23,-23,19,-19,20,-20])
 E6=(((c2*Ta2)*Ta1)*b)
 #E6.combineBond([23,22])
 E6.permute([21,4,-4,23,-23,22],4)
 E6.setLabel([13,3,-3,11,-11,22])
 E6.setLabel([10,4,-4,11,-11,22])
 E6.permute([11,-11,22,4,-4,10],3)
 
 
 
 c3.setLabel([8,7])
 Tb2.setLabel([7,15,-15,22])
 Tb3.setLabel([13,14,-14,8])
 d.setLabel([16,-16,14,-14,15,-15,23,-23])
 E5=(((c3*Tb3)*Tb2)*d)
 #E5.combineBond([23,22])
 E5.permute([13,16,-16,23,-23,22],4)
 E5.setLabel([8,7,-7,11,-11,22])
 E5.setLabel([12,1,-1,11,-11,22])
 E5.permute([11,-11,22,1,-1,12],0)
 return E1, E2, E3, E4, E5,E6









 
 
def produce_Env_Hab(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys):

 c1.setLabel([0,1])
 Tb1.setLabel([1,2,-2,3])
 E1=c1*Tb1
 E1.permute([0,2,-2,3],2)
 E1.setLabel([7,3,-3,8])
 E1.permute([7,3,-3,8],3)
 
 c2.setLabel([1,0])
 Ta1.setLabel([3,2,-2,1])
 E6=c2*Ta1
 E6.permute([0,2,-2,3],2)
 E6.setLabel([9,6,-6,8])
 E6.permute([8,6,-6,9],4)
 
 E5=Ta2
 E5.setLabel([10,5,-5,9])
 E5.permute([10,5,-5,9],3)
 
 E2=Tb4
 E2.setLabel([12,13,-13,7])
 E2.permute([12,13,-13,7],1)

 
 c3.setLabel([8,7])
 Tb2.setLabel([7,15,-15,22])
 Tb3.setLabel([13,14,-14,8])
 d.setLabel([16,-16,14,-14,15,-15,23,-23])
 E4=(((c3*Tb3)*Tb2)*d)
 E4.permute([13,16,-16,23,-23,22],3)
 E4.setLabel([11,16,-16,4,-4,10])
 E4.permute([11,16,-16,4,-4,10],3)


 
 c4.setLabel([11,10])
 Ta4.setLabel([11,17,-17,18])
 Ta3.setLabel([10,12,-12,13])
 c.setLabel([17,-17,12,-12,16,-16,19,-19])
 E3=(((c4*Ta4)*Ta3)*c)
 E3.permute([13,16,-16,19,-19,18],0)
 E3.setLabel([11,16,-16,1,-1,12])

 return E1, E2, E3, E4, E5,E6
 
 
 
 
def final_test_distance(ap_u, bp_u, a_u, b_u,E1, E2, E3, E4, E5,E6,U,N_uni):


 U.setLabel([-20,-40,20,40])
 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([-20,-40,30,50])
  
 U.setLabel([30,50,20,40])
 H=U*H1
 H.permute([-20,-40,20,40],2)
 #H.setLabel([-20,-40,20,40]) 
 
 H1.setLabel([-20,-40,20,40]) 
 U.setLabel([-20,-40,20,40])
 
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-20,-40,20,40])
 
 
 a_u.setLabel([20,13,1,2,3])
 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-2,-3,-20,-13,-1])

 b_u.setLabel([40,2,4,5,6])
 b_d=copy.copy(b_u)
 b_d.transpose()
 b_d.setLabel([-5,-6,-40,-2,-4])


 ap_u.setLabel([20,13,1,2,3])
 ap_d=copy.copy(ap_u)
 ap_d.transpose()
 ap_d.setLabel([-2,-3,-20,-13,-1])

 bp_u.setLabel([40,2,4,5,6])
 bp_d=copy.copy(bp_u)
 bp_d.transpose()
 bp_d.setLabel([-5,-6,-40,-2,-4])
 #N=(((((E2*)*E1)*E3)*Iden)*((((E5)*E6)*E4)))

# a_u=q*r
# a_u.permute([20,4,5,3,2],3)
# a_u.setLabel([0,1,2,3,4])
# #print a_ut.similar(a_u)#,a_ut.printDiagram() 
# 
# b_u=qq*l
# b_u.permute([40,3,4,5,6],3)
# b_u.setLabel([0,1,2,3,4])
# #print b_ut.similar(b_u) 
 
 A1=(((((E2*(a_u*a_d))*E1)*E3)*H)*(((((b_u*b_d)*E5)*E6)*E4))) 
 A2=((((E2*(ap_u*ap_d))*E1)*E3)*Iden)*(((((bp_u*bp_d)*E5)*E6)*E4))
 A3=((((E2*(a_u*ap_d))*E1)*E3)*U)*(((((b_u*bp_d)*E5)*E6)*E4))
 A4=((((E2*(ap_u*a_d))*E1)*E3)*H1)*(((((bp_u*b_d)*E5)*E6)*E4))

 A=A1+A2+(-1.00)*A3+(-1.00)*A4
 return A
 
 
 
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

def Sqrt_mat(e):
 d=int(e.row())
 
 for q in xrange(d):
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e

def Inv_mat(e):
 d=int(e.row())
 
 for q in xrange(d):
   if e[q] > 1.0e-11:  
    e[q]=((e[q])**(-1.00))
   else:  
    e[q]=0.0 
 return e
 
 
 
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

def reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env):
  
  Env[0]=copy.copy(c1)
  Env[1]=copy.copy(c2)
  Env[2]=copy.copy(c3) 
  Env[3]=copy.copy(c4) 
  Env[4]=copy.copy(Ta1)
  Env[5]=copy.copy(Ta2)
  Env[6]=copy.copy(Ta3)
  Env[7]=copy.copy(Ta4)
  Env[8]=copy.copy(Tb1)
  Env[9]=copy.copy(Tb2)
  Env[10]=copy.copy(Tb3)
  Env[11]=copy.copy(Tb4)
  
def Rand_env_slight(Env):
 
 for i in xrange(len(Env)):
  Env_tem=copy.copy(Env[i]) 
  Env_tem.randomize()
  Env[i]=Env[i]+0.01*Env_tem

def Rand_env_total(Env):
 Env1=copy.copy(Env)
 for i in xrange(len(Env1)):
  Env1[i]=copy.copy(Env[i])
  Env1[i].randomize()
 
 return  Env1
  

