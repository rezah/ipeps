import pyUni10 as uni10
import copy
import time
import basic
#import line_profiler

def Var_cab(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0,N_grad, Opt_method):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 
 #t0=time.time()
 if Corner_method is 'CTM':
#  c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1=basic.make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1)
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D,H0,d_phys)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
  #basic.Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
 #print time.time() - t0, "CTM-H, Left"

 Env=basic.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 #t0=time.time()
 E1, E2, E3, E4, E5, E6, E7, E8=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 #print time.time() - t0, "Env, Left"
 a_up=copy.copy(a_u) 
 b_up=copy.copy(b_u) 
 c_up=copy.copy(c_u) 
 d_up=copy.copy(d_u) 


# Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
# print "DisF", Dis_val

 a_up, b_up, c_up, d_up=Inite_opt(E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,a, b, c, d)


 
 a_up, b_up, c_up, d_up=Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method)



 c_up, b_up, a_up=equall_dis(c_up,b_up,a_up) 
 a_up, c_up, d_up=equall_dis1(a_up, c_up, d_up)
 
 Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
 print "Dis", Dis_val
 

 a_up=basic.max_ten(a_up)
 b_up=basic.max_ten(b_up) 
 c_up=basic.max_ten(c_up)
 d_up=basic.max_ten(d_up)

 ap=basic.make_ab(a_up)
 bp=basic.make_ab(b_up)
 cp=basic.make_ab(c_up)
 dp=basic.make_ab(d_up)


 return a_up, b_up, c_up, d_up, ap, bp, cp, dp


def Var_abd(a_u, b_u,c_u,d_u,a,b,c,d,Env,D,U,d_phys,chi,Gauge,Corner_method,H0,N_grad, Opt_method):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)
 
 #t0=time.time()
 if Corner_method is 'CTM':
#  c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1=basic.make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1)
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D,H0,d_phys)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
  #basic.Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
 #print time.time() - t0, "CTM-H, Left"

 Env=basic.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 #t0=time.time()
 E1, E2, E3, E4, E5, E6, E7, E8=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 #print time.time() - t0, "Env, Left"
 a_up=copy.copy(a_u) 
 b_up=copy.copy(b_u) 
 c_up=copy.copy(c_u) 
 d_up=copy.copy(d_u) 
 

 
# Dis_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
# print "Dis", Dis_val
 
 #plist=Reload_plist(plist)
 
 
 a_up, b_up, c_up, d_up=Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method)

 c_up, b_up, a_up=equall_dis(c_up,b_up,a_up) 
 a_up, c_up, d_up=equall_dis1(a_up, c_up, d_up)
 
 Dis_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
 print "Dis", Dis_val
 

 a_up=basic.max_ten(a_up)
 b_up=basic.max_ten(b_up) 
 c_up=basic.max_ten(c_up)
 d_up=basic.max_ten(d_up)

 ap=basic.make_ab(a_up)
 bp=basic.make_ab(b_up)
 cp=basic.make_ab(c_up)
 dp=basic.make_ab(d_up)


 return a_up, b_up,c_up,d_up, ap,bp,cp,dp



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


 b.setLabel([18,-18,20,-20,6,-6,4,-4])
 c.setLabel([14,-14,12,-12,19,-19,17,-17])
 a.setLabel([16,-16,17,-17,18,-18,2,-2])
 d.setLabel([19,-19,10,-10,8,-8,20,-20])
 Norm=(((((E1*E8)*(a))*((E7*E6)*(c))))*(((E2*E3)*(b))))*((E4*E5)*d)
 #print "Norm", Norm, Norm[0]
 if Norm[0] < 0: print "Norm<0"; E1=-1.0*E1;

 return E1, E2, E3, E4, E5, E6, E7, E8 


  
 
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


 E_val=(((((E1*E8)*(a_u*a_d)*U)*((E7*E6)*(c_u*c_d))))*(((E2*E3)*(b_u*b_d))))*((E4*E5)*(d_u*d_d))

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


 
#@profile 
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


 Val2=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_up*d_d))
# Val3=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u*d_dp))
 #print Val2[0]
 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=Val1[0]-Val2[0]-Val3[0]
 val_f=Val1[0]-2.0*Val2[0]#-Val3[0]

 if Val1[0]<0: print "gggggggggggg"
 #print Val2[0], Val3[0]
 return  val_f

#@profile
def Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U, N_grad, Opt_method):

  a_uf=copy.copy(a_up) 
  b_uf=copy.copy(b_up) 
  c_uf=copy.copy(c_up) 
  d_uf=copy.copy(d_up) 

  Es=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
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
   count+=1

   E1_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   Ef=E1_val
   print 'E=', E1_val, count

#   if i%4==0 and i >0:
#    a_up, b_up, c_up, d_up=Inite_opt(E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U, a, b, c, d)

   D_a,D_b,D_c,D_d=basic.Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
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
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 :
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

  if( Ef > Es):
   print 'SD method, Fail, f<s', Ef, Es 
   a_up=copy.copy(a_uf) 
   b_up=copy.copy(b_uf) 
   c_up=copy.copy(c_uf) 
   d_up=copy.copy(d_uf) 

  return a_up, b_up, c_up, d_up


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


 Val2=(((((E1*E8)*(a_up*a_d)*U)*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_up*d_d)))* ((E7*E6)*(c_up*c_d)))
 #Val3=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp)))*U)*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u*d_dp))
 #print Val, Val1, Val2, Val3
 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=Val1[0]-Val2[0]-Val3[0]
 val_f=Val1[0]-2.0*Val2[0]#-Val3[0]

 if Val1[0]<0: print "gggggggggggg"
 #print Val2[0], Val3[0]
 return  val_f

def Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U,N_grad, Opt_method):


  a_uf=copy.copy(a_up) 
  b_uf=copy.copy(b_up) 
  c_uf=copy.copy(c_up) 
  d_uf=copy.copy(d_up) 

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
   count+=1
   E1_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
   Ef=E1_val
   print 'E=', E1_val, count
   
   D_a,D_b,D_c,D_d=basic.Obtain_grad_four1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
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

  if( Ef > Es):
   print 'SD method, Fail, f<s', Ef, Es 
   a_up=copy.copy(a_uf) 
   b_up=copy.copy(b_uf) 
   c_up=copy.copy(c_uf) 
   d_up=copy.copy(d_uf) 

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

 return c_up, b_up, a_up


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
def Sqrt_mat(e):
 d=int(e.row())
 for q in xrange(d):
   #print e[q] 
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e  

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

#def N_Positiv(N):
# N.permute([17,-80,-81,-82,-83,84,85,-17,80,81,82,83,-84,-85],7)
# N.permute([-17,-80,-81,-82,-83,-84,-85,17,80,81,82,83,84,85], 7)
# N1=copy.copy(N)
# N1.transpose()
# N=(N+N1)*(1.00/2.00)
# N1=copy.copy(N)
# N1.setLabel([17,80,81,82,83,84,85,0,1,2,3,4,5,6])
# N=N*N1
# N.permute([-17,-80,-81,-82,-83,-84,-85,0,1,2,3,4,5,6],7)
# N_final=sqrt_general(N)
# N_final.setLabel([-17,-80,-81,-82,-83,-84,-85,17,80,81,82,83,84,85])
# N_final.permute([17,-80,-81,-82,-83,84,85,-17,80,81,82,83,-84,-85],7)
# return N_final              

#@profile
def Inite_opt(E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U, a, b, c, d):

  r_u, r1_u, r2_u, q_u, qq_u, qqq_u = Qr_lQ_decomf(a_u,b_u,c_u)
  r_u, r1_u, r2_ug, q_u, qq_u, qqq_ug = Qr_lQ_decomg(a_u,b_u,c_u)
  N_up, r_up, r1_up, r2_up, q_up, qq_up, qqq_up = Qr_lQ_decomP(E1, E2, E3, E4, E5, E6, E7, E8, d_up, U,a_up,b_up,c_up)

  N_iter=20
  E_0=1.0
  E_1=2.0
#  Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U)
#  print "Dis_first", Dis_val
  
  
 ##################################################################################################
  for i in xrange(N_iter): 

   N_up, r2_up, qqq_up = Qr_lQ_decomVV(E1, E2, E3, E4, E5, E6, E7, E8, d_up, U, a_up, q_up, qq_up)


   N_uper, N_down = make_env_N(E1, E2, E3, E4, E5, E6, E7, E8, q_u, qq_u, qqq_u, q_up, qq_up, qqq_up,d_u, d_up)

   E_0=E_1
   E_1=Dis_fQR(N_up,N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )
   print "E=", E_1, abs((E_0-E_1) / E_0)

   if (( abs((E_0-E_1) / E_0) ) < 1.0e-6): 
    print 'break', E_0, E_1, abs((E_0-E_1) / E_0)
    break


#   print "M1",Dis_fQR(N_up,N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )


   r_up, r1_up, r2_up=Do_optimization_Full(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )

   a_up=reproduce_a(r2_up, qqq_up)

#   print "E=", Dis_fQR(N_up,N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U ) 
################################################################################################
   N_up, r2_up, qqq_up = Qr_lQ_decomV(E1, E2, E3, E4, E5, E6, E7, E8, d_up, U,a_up,q_up,qq_up)


   N_uper, N_down = make_env_Ng(E1, E2, E3, E4, E5, E6, E7, E8, q_u, qq_u, qqq_ug, q_up, qq_up, qqq_up,d_u, d_up)

#   print "E=", Dis_fQRg(N_up,N_uper, N_down, r_u, r1_u, r2_ug,r_up, r1_up, r2_up, U ) 

   r_up, r1_up, r2_up=Do_optimization_Fullg(N_up, N_uper, N_down, r_u, r1_u, r2_ug,r_up, r1_up, r2_up, U )

   a_up=reproduce_a(r2_up, qqq_up)
#   print "M2", Dis_fQRg(N_up,N_uper, N_down, r_u, r1_u, r2_ug,r_up, r1_up, r2_up, U )



  b_up, c_up =reproduce_bc(r_up, r1_up,q_up, qq_up)

  return a_up, b_up, c_up, d_up


def Qr_lQ_decomP(E1, E2, E3, E4, E5, E6, E7, E8, d_u, U,a_u,b_u,c_u):
 d_u.setLabel([57,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,57,-19,-10])
 
 
############################################################################
 A=copy.copy(c_u)
 A.setLabel([52,14,12,19,17])
 A.permute([14,12,19,52,17],3)
 
 q,r=qr_parity(A) 
 
 q.setLabel([14,12,19,80,81])
 r.setLabel([80,81,52,17])
 q.permute([14,12,19,80,81],2)
 r.permute([80,81,52,17],3)
 
 q_d=copy.copy(q)
 q_d.transpose()
 q_d.setLabel([-19,-80,-81,-14,-12])
####################################################################################
 A=copy.copy(a_u)
 A.setLabel([53,16,17,18,2])
 A.permute([16,17,2,53,18],3)
 
 qqq,r2=qr_parity(A) 
 
 qqq.setLabel([16,17,2,82,83])
 r2.setLabel([82,83,53,18])
 qqq.permute([16,17,82,83,2],2)
 r2.permute([82,83,53,18],3)
 
 qqq_d=copy.copy(qqq)
 qqq_d.transpose()
 qqq_d.setLabel([-82,-83,-2,-16,-17])
########################################## 
 A=copy.copy(b_u)
 A.setLabel([54,18,20,6,4])
 A.permute([54,18,20,6,4],2)
 
 r1, qq=lq_parity(A) 
 
 r1.setLabel([54,18,84,85])
 qq.setLabel([84,85,20,6,4])
 r1.permute([54,18,84,85],2)
 qq.permute([84,85,20,6,4],3)
 
 qq_d=copy.copy(qq)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-84,-85,-20])
################################################


 N=((((E4*E5)*(d_u*d_d))*((E2*E3)*(qq*qq_d)))*((E7*E6)*(q*q_d)))*(((E1*E8)*(qqq*qqq_d)))

 N.permute([17,-80,-81,-82,-83,84,85,-17,80,81,82,83,-84,-85],7)
 ###testing####
# r_d=copy.copy(r)
# r1_d=copy.copy(r1)
# r2_d=copy.copy(r2)

# r_d.transpose()
# r1_d.transpose()
# r2_d.transpose()

# r_d.setLabel([-17,-80,-81,52])
# r1_d.setLabel([-84,-85,54,-18])
# r2_d.setLabel([-18,-82,-83,53])
# 
# 
# Norm=(((N*(r*r_d))*(r1*r1_d)))*(r2*r2_d)
# print "Norm-QR", Norm[0], Norm.printDiagram()

 return N, r, r1, r2, q, qq, qqq


def Qr_lQ_decomVV(E1, E2, E3, E4, E5, E6, E7, E8, d_u, U, a_u, q, qq):
 d_u.setLabel([57,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,57,-19,-10])
############################################################################ 
 q_d=copy.copy(q)
 q_d.transpose()
 q_d.setLabel([-19,-80,-81,-14,-12])
####################################################################################
 A=copy.copy(a_u)
 A.setLabel([53,16,17,18,2])
 A.permute([16,17,2,53,18],3)
 
 qqq,r2=qr_parity(A) 
 
 qqq.setLabel([16,17,2,82,83])
 r2.setLabel([82,83,53,18])
 qqq.permute([16,17,82,83,2],2)
 r2.permute([82,83,53,18],3)
 
 qqq_d=copy.copy(qqq)
 qqq_d.transpose()
 qqq_d.setLabel([-82,-83,-2,-16,-17])
##########################################  
 qq_d=copy.copy(qq)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-84,-85,-20])
################################################


 N=((((E4*E5)*(d_u*d_d))*((E2*E3)*(qq*qq_d)))*((E7*E6)*(q*q_d)))*(((E1*E8)*(qqq*qqq_d)))
 N.permute([17,-80,-81,-82,-83,84,85,-17,80,81,82,83,-84,-85],7)

 return N, r2, qqq












     
def Qr_lQ_decomV(E1, E2, E3, E4, E5, E6, E7, E8, d_u, U,a_u,q,qq):
 d_u.setLabel([57,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,57,-19,-10])
############################################################################
 q_d=copy.copy(q)
 q_d.transpose()
 q_d.setLabel([-19,-80,-81,-14,-12])
####################################################################################
 A=copy.copy(a_u)
 A.setLabel([53,16,17,18,2])
 A.permute([53,17,16,18,2],2)
 
 r2, qqq=lq_parity(A) 
 
 qqq.setLabel([82,83,16,18,2])
 r2.setLabel([53,17,82,83])
 qqq.permute([16,82,83,18,2],3)
 r2.permute([53,17,82,83],2)
 
 qqq_d=copy.copy(qqq)
 qqq_d.transpose()
 qqq_d.setLabel([-18,-2,-16,-82,-83])
########################################## 
 qq_d=copy.copy(qq)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-84,-85,-20])
################################################


 
 N=((((E4*E5)*(d_u*d_d))*((E2*E3)*(qq*qq_d)))*((E7*E6)*(q*q_d)))*(((E1*E8)*(qqq*qqq_d)))

 N.permute([18,-80,-81,82,83,84,85,-18,80,81,-82,-83,-84,-85],7)
 ###testing####
# r_d=copy.copy(r)
# r1_d=copy.copy(r1)
# r2_d=copy.copy(r2)

# r_d.transpose()
# r1_d.transpose()
# r2_d.transpose()

# r_d.setLabel([-17,-80,-81,52])
# r1_d.setLabel([-84,-85,54,-18])
# r2_d.setLabel([-82,-83,53,-17])
# 
# 
# Norm=(((N*(r*r_d))*(r1*r1_d)))*(r2*r2_d)
# print "Norm-QR", Norm[0], Norm.printDiagram()

 return N, r2, qqq



def Qr_lQ_decomf(a_u,b_u,c_u):
############################################################################
 A=copy.copy(c_u)
 A.setLabel([52,14,12,19,17])
 A.permute([14,12,19,52,17],3)
 
 q,r=qr_parity(A) 
 
 q.setLabel([14,12,19,80,81])
 r.setLabel([80,81,52,17])
 q.permute([14,12,19,80,81],2)
 r.permute([80,81,52,17],3)
 
####################################################################################
 A=copy.copy(a_u)
 A.setLabel([53,16,17,18,2])
 A.permute([16,17,2,53,18],3)
 
 qqq,r2=qr_parity(A) 
 
 qqq.setLabel([16,17,2,82,83])
 r2.setLabel([82,83,53,18])
 qqq.permute([16,17,82,83,2],2)
 r2.permute([82,83,53,18],3)
 
########################################## 
 A=copy.copy(b_u)
 A.setLabel([54,18,20,6,4])
 A.permute([54,18,20,6,4],2)
 
 r1, qq=lq_parity(A) 
 
 r1.setLabel([54,18,84,85])
 qq.setLabel([84,85,20,6,4])
 r1.permute([54,18,84,85],2)
 qq.permute([84,85,20,6,4],3)
 
 qq_d=copy.copy(qq)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-84,-85,-20])
################################################
 return r, r1, r2, q, qq, qqq





def Qr_lQ_decomg(a_u,b_u,c_u):
############################################################################
 A=copy.copy(c_u)
 A.setLabel([52,14,12,19,17])
 A.permute([14,12,19,52,17],3)
 
 q,r=qr_parity(A) 
 
 
 q.setLabel([14,12,19,80,81])
 r.setLabel([80,81,52,17])
 q.permute([14,12,19,80,81],2)
 r.permute([80,81,52,17],3)
 
####################################################################################
 A=copy.copy(a_u)
 A.setLabel([53,16,17,18,2])
 A.permute([53,17,16,18,2],2)
 
 r2, qqq=lq_parity(A) 
 
 qqq.setLabel([82,83,16,18,2])
 r2.setLabel([53,17,82,83])
 qqq.permute([16,82,83,18,2],3)
 r2.permute([53,17,82,83],2)
 
########################################## 
 A=copy.copy(b_u)
 A.setLabel([54,18,20,6,4])
 A.permute([54,18,20,6,4],2)
 
 r1, qq=lq_parity(A) 
 
 r1.setLabel([54,18,84,85])
 qq.setLabel([84,85,20,6,4])
 r1.permute([54,18,84,85],2)
 qq.permute([84,85,20,6,4],3)
 
################################################
 return r, r1, r2, q, qq, qqq



#@profile
def make_env_N(E1, E2, E3, E4, E5, E6, E7, E8, q_u, qq_u, qqq_u, q_up, qq_up, qqq_up,d_u, d_up):
 d_u.setLabel([57,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,57,-19,-10])



 d_up.setLabel([57,19,10,8,20])
 d_dp=copy.copy(d_up)
 d_dp.transpose()
 d_dp.setLabel([-8,-20,57,-19,-10])

 q_d=copy.copy(q_u)
 q_d.transpose()
 q_d.setLabel([-19,-80,-81,-14,-12])

 qq_d=copy.copy(qq_u)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-84,-85,-20])


 qqq_d=copy.copy(qqq_u)
 qqq_d.transpose()
 qqq_d.setLabel([-82,-83,-2,-16,-17])


 q_dp=copy.copy(q_up)
 q_dp.transpose()
 q_dp.setLabel([-19,-80,-81,-14,-12])

 qq_dp=copy.copy(qq_up)
 qq_dp.transpose()
 qq_dp.setLabel([-6,-4,-84,-85,-20])


 qqq_dp=copy.copy(qqq_up)
 qqq_dp.transpose()
 qqq_dp.setLabel([-82,-83,-2,-16,-17])

 N_upper=((((E4*E5)*(d_up*d_d))*((E2*E3)*(qq_up*qq_d)))*((E7*E6)*(q_up*q_d)))*(((E1*E8)*(qqq_up*qqq_d)))

 N_down=((((E4*E5)*(d_u*d_dp))*((E2*E3)*(qq_u*qq_dp)))*((E7*E6)*(q_u*q_dp)))*(((E1*E8)*(qqq_u*qqq_dp)))


 return N_upper, N_down
#@profile
def make_env_Ng(E1, E2, E3, E4, E5, E6, E7, E8, q_u, qq_u, qqq_u, q_up, qq_up, qqq_up,d_u, d_up):
 d_u.setLabel([57,19,10,8,20])
 d_d=copy.copy(d_u)
 d_d.transpose()
 d_d.setLabel([-8,-20,57,-19,-10])



 d_up.setLabel([57,19,10,8,20])
 d_dp=copy.copy(d_up)
 d_dp.transpose()
 d_dp.setLabel([-8,-20,57,-19,-10])

 q_d=copy.copy(q_u)
 q_d.transpose()
 q_d.setLabel([-19,-80,-81,-14,-12])

 qq_d=copy.copy(qq_u)
 qq_d.transpose()
 qq_d.setLabel([-6,-4,-84,-85,-20])


 qqq_d=copy.copy(qqq_u)
 qqq_d.transpose()
 qqq_d.setLabel([-18,-2,-16,-82,-83])

 q_dp=copy.copy(q_up)
 q_dp.transpose()
 q_dp.setLabel([-19,-80,-81,-14,-12])

 qq_dp=copy.copy(qq_up)
 qq_dp.transpose()
 qq_dp.setLabel([-6,-4,-84,-85,-20])


 qqq_dp=copy.copy(qqq_up)
 qqq_dp.transpose()
 qqq_dp.setLabel([-18,-2,-16,-82,-83])

 N_upper=((((E4*E5)*(d_up*d_d))*((E2*E3)*(qq_up*qq_d)))*((E7*E6)*(q_up*q_d)))*(((E1*E8)*(qqq_up*qqq_d)))

 N_down=((((E4*E5)*(d_u*d_dp))*((E2*E3)*(qq_u*qq_dp)))*((E7*E6)*(q_u*q_dp)))*(((E1*E8)*(qqq_u*qqq_dp)))

 return N_upper, N_down


def Do_optimization_Full(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):

 r_up_first=copy.copy(r_up) 
 r1_up_first=copy.copy(r1_up)
 r2_up_first=copy.copy(r2_up)
 checking_val=0

 Distance_val=Dis_fQR(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )
# print 'Dis0', Distance_val 
 
 Res=10
 Res1=Distance_val
 count=0
 for q in xrange(3):
  r_up=optimum_0(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U)
  r1_up=optimum_1(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U)
  r2_up=optimum_2(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U)

  Distance_val=Dis_fQR(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )
#  print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q

  Res=Res1
  Res1=Distance_val
#  print Res, Res1, Res1 < Res 

  if Res1 < Res:
    r_up_first=copy.copy(r_up) 
    r1_up_first=copy.copy(r1_up)
    r2_up_first=copy.copy(r2_up)
  else: 
    r_up=copy.copy(r_up_first) 
    r1_up=copy.copy(r1_up_first)
    r2_up=copy.copy(r2_up_first)
    break
 
  count+=1
  if count > 100: print 'Num_Opt > 30'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-8) or ((abs(Res1-Res) / abs(Res)) < 1.00e-9): 
    #print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val) < 1.00e-8) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break


 return r_up, r1_up, r2_up




def Do_optimization_Fullg(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):

 r_up_first=copy.copy(r_up) 
 r1_up_first=copy.copy(r1_up)
 r2_up_first=copy.copy(r2_up)
 checking_val=0

 Distance_val=Dis_fQRg(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )
 #print 'Dis0', Distance_val
 
 Res=10
 Res1=Distance_val
 count=0
 for q in xrange(3):
  #print "\n", "\n"
  r_up=optimum_0g(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U)
  r1_up=optimum_1g(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U)
  r2_up=optimum_2g(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U)

  Distance_val=Dis_fQRg(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )
#  print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q

  Res=Res1
  Res1=Distance_val
#  print Res, Res1, Res1 < Res 
  if Res1 < Res:
    r_up_first=copy.copy(r_up) 
    r1_up_first=copy.copy(r1_up)
    r2_up_first=copy.copy(r2_up)
  else: 
    r_up=copy.copy(r_up_first) 
    r1_up=copy.copy(r1_up_first)
    r2_up=copy.copy(r2_up_first)
    break


  if (q>1) and (Res1 > Res) and ((abs(Res1-Res) / abs(Res)) > 1.00e-4): checking_val=1;

  count+=1
  if count > 100: print 'Num_Opt > 30'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-8) or ((abs(Res1-Res) / abs(Res)) < 1.00e-9): 
    #print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-8) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break
 #Dis_val=Dis_fQRg(N_up, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U )
 #print 'DisF', Distance_val, abs(Res1-Res) / abs(Res), q


# if checking_val != 0:
#  r_up=copy.copy(r_up_first) 
#  r1_up=copy.copy(r1_up_first)
#  r2_up=copy.copy(r2_up_first)

 return r_up, r1_up, r2_up



 
def Dis_fQR(N_u,N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U ):
 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([52,53,54,1,2,3])
 H=U*H1
 H.permute([-52,-53,-54,1,2,3],3)
 H.setLabel([-52,-53,-54,52,53,54])


 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-18,-82,-83,-53])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-18,-82,-83,-53])

 #print N.printDiagram(), H.printDiagram(), r.printDiagram() 
 Val=((r_d*(r2_d*r1_d)))*(N_u*(r_u*(r2_u*r1_u)*H))
# print Val[0]
 Val1=((r_dp*(r2_dp*r1_dp)))*(N_u*(r_up*(r2_up*r1_up)*Iden))
# print Val1[0]
 Val2=((r_dp*(r2_dp*r1_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))
# print Val2[0]
 Val3=((r_d*(r2_d*r1_d)))*(N_uper*((r_up*(r2_up*r1_up))*U))
# print Val3[0]
# Val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #Val_f=Val1[0]-Val2[0]-Val3[0]
# Val_f=Val1[0]-2.0*Val2[0]#-Val3[0]
 Val_f=Val1[0]-2.0*Val3[0]#-Val3[0]

 return Val_f
 


def Dis_fQRg(N_u,N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U ):
 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 H1=copy.copy(U)
 H1.transpose()
 H1.setLabel([52,53,54,1,2,3])
 H=U*H1
 H.permute([-52,-53,-54,1,2,3],3)
 H.setLabel([-52,-53,-54,52,53,54])


 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-82,-83,-53,-17])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-82,-83,-53,-17])

 #print N.printDiagram(), H.printDiagram(), r.printDiagram() 
 Val=((r_d*(r2_d*r1_d)))*(N_u*(r_u*(r2_u*r1_u)*H))
# print Val[0]
 Val1=((r_dp*(r2_dp*r1_dp)))*(N_u*(r_up*(r2_up*r1_up)*Iden))
# print  Val1[0]
 Val2=((r_dp*(r2_dp*r1_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))
# print Val2[0]
 Val3=((r_d*(r2_d*r1_d)))*(N_uper*((r_up*(r2_up*r1_up))*U))
# print Val3[0]
# Val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #Val_f=Val1[0]-Val2[0]-Val3[0]
# Val_f=Val1[0]-2.0*Val2[0]#-Val3[0]
 Val_f=Val1[0]-2.0*Val3[0]#-Val3[0]

 return Val_f


 
 
def optimum_0(N_u, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):
 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-18,-82,-83,-53])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-18,-82,-83,-53])

 A2=(((r2_dp*r1_dp)))*(N_u*((r2_up*r1_up)*Iden))        
 A2.permute([-80,-81,-52,-17,80,81,52,17],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-80,-81,-52,-17,80,81,52,17])

 A3=(((r2_dp*r1_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))        
 A3.permute([-80,-81,-52,-17],0)
 A3p=((r_d*(r2_d*r1_d)))*(N_uper*(((r2_up*r1_up))*U))
 A3p.permute([80,81,52,17],4)
 A3p.transpose()

 A3=A3+A3p
 
 A3.setLabel([-80,-81,-52,-17])
 
 U, S, V=svd_parity(A2)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([80,81,52,17,-80,-81,-52,-17])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([80,81,52,17])
 A.permute([80,81,52,17],3)

 rf=copy.copy(A)

 return rf
def optimum_1(N_u, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):

 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-18,-82,-83,-53])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-18,-82,-83,-53])

 A2=((r_dp*(r2_dp)))*(N_u*(r_up*(r2_up)*Iden))        
 A2.permute([-54,-18,-84,-85,54,18,84,85],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-54,-18,-84,-85,54,18,84,85])

 A3=((r_dp*(r2_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))        
 A3.permute([-54,-18,-84,-85],0)
 A3p=((r_d*(r2_d*r1_d)))*(N_uper*((r_up*(r2_up))*U))
 A3p.permute([54,18,84,85],4)
 A3p.transpose()

 A3=A3+A3p
 
 A3.setLabel([-54,-18,-84,-85])
 
 U, S, V=svd_parity(A2)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([54,18,84,85,-54,-18,-84,-85])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([54,18,84,85])
 A.permute([54,18,84,85],2)

 rf=copy.copy(A)

 return rf


def optimum_2(N_u, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):

 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-18,-82,-83,-53])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-18,-82,-83,-53])

 A2=((r_dp*(r1_dp)))*(N_u*(r_up*(r1_up)*Iden))        
 A2.permute([-82,-83,-53,-18,82,83,53,18],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-82,-83,-53,-18,82,83,53,18])

 A3=((r_dp*(r1_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))        
 A3.permute([-82,-83,-53,-18],0)
 A3p=((r_d*(r2_d*r1_d)))*(N_uper*((r_up*(r1_up))*U))
 A3p.permute([82,83,53,18],4)
 A3p.transpose()

 A3=A3+A3p
 
 A3.setLabel([-82,-83,-53,-18])
 
 U, S, V=svd_parity(A2)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([82,83,53,18,-82,-83,-53,-18])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([82,83,53,18])
 A.permute([82,83,53,18],3)

 rf=copy.copy(A)

 return rf



def optimum_0g(N_u, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):
 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-82,-83,-53,-17])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-82,-83,-53,-17])

 A2=(((r2_dp*r1_dp)))*(N_u*((r2_up*r1_up)*Iden))        
 A2.permute([-80,-81,-52,-17,80,81,52,17],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-80,-81,-52,-17,80,81,52,17])

 A3=(((r2_dp*r1_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))        
 A3.permute([-80,-81,-52,-17],0)
 A3p=((r_d*(r2_d*r1_d)))*(N_uper*(((r2_up*r1_up))*U))
 A3p.permute([80,81,52,17],4)
 A3p.transpose()

 A3=A3+A3p
 
 A3.setLabel([-80,-81,-52,-17])
 
 U, S, V=svd_parity(A2)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([80,81,52,17,-80,-81,-52,-17])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([80,81,52,17])
 A.permute([80,81,52,17],3)

 rf=copy.copy(A)

 return rf
def optimum_1g(N_u, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):

 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-82,-83,-53,-17])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-82,-83,-53,-17])

 A2=((r_dp*(r2_dp)))*(N_u*(r_up*(r2_up)*Iden))        
 A2.permute([-54,-18,-84,-85,54,18,84,85],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-54,-18,-84,-85,54,18,84,85])

 A3=((r_dp*(r2_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))        
 A3.permute([-54,-18,-84,-85],0)
 A3p=((r_d*(r2_d*r1_d)))*(N_uper*((r_up*(r2_up))*U))
 A3p.permute([54,18,84,85],4)
 A3p.transpose()

 A3=A3+A3p
 
 A3.setLabel([-54,-18,-84,-85])
 
 U, S, V=svd_parity(A2)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([54,18,84,85,-54,-18,-84,-85])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([54,18,84,85])
 A.permute([54,18,84,85],2)

 rf=copy.copy(A)

 return rf


def optimum_2g(N_u, N_uper, N_down, r_u, r1_u, r2_u,r_up, r1_up, r2_up, U):

 U.setLabel([-52,-53,-54,52,53,54])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([-52,-53,-54,52,53,54])

 r_d=copy.copy(r_u)
 r1_d=copy.copy(r1_u)
 r2_d=copy.copy(r2_u)

 r_d.transpose()
 r1_d.transpose()
 r2_d.transpose()

 r_d.setLabel([-17,-80,-81,-52])
 r1_d.setLabel([-84,-85,-54,-18])
 r2_d.setLabel([-82,-83,-53,-17])

 r_dp=copy.copy(r_up)
 r1_dp=copy.copy(r1_up)
 r2_dp=copy.copy(r2_up)

 r_dp.transpose()
 r1_dp.transpose()
 r2_dp.transpose()

 r_dp.setLabel([-17,-80,-81,-52])
 r1_dp.setLabel([-84,-85,-54,-18])
 r2_dp.setLabel([-82,-83,-53,-17])

 A2=((r_dp*(r1_dp)))*(N_u*(r_up*(r1_up)*Iden))        
 A2.permute([-53,-17,-82,-83,53,17,82,83],4)

 A2_trans=copy.copy(A2)
 A2_trans.transpose()
 A2=A2+A2_trans
 
 A2.setLabel([-53,-17,-82,-83,53,17,82,83])

 A3=((r_dp*(r1_dp)))*(N_down*((r_u*(r2_u*r1_u))*U))        
 A3.permute([-53,-17,-82,-83],0)
 A3p=((r_d*(r2_d*r1_d)))*(N_uper*((r_up*(r1_up))*U))
 A3p.permute([53,17,82,83],4)
 A3p.transpose()

 A3=A3+A3p
 
 A3.setLabel([-53,-17,-82,-83])
 
 U, S, V=svd_parity(A2)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A2_inv=V*S*U
 A2_inv.permute([0,1,2,3,12,13,14,15],4)
 A2_inv.setLabel([53,17,82,83,-53,-17,-82,-83])
 
 #distance_iden_val=distance_iden(A2_mat,A2_inv)
 #print 'distance1=', distance_iden_val
 #print A2.getBlock()*A2_inv


 A=A2_inv*A3
 A.setLabel([53,17,82,83])
 A.permute([53,17,82,83],2)

 rf=copy.copy(A)

 return rf



def reproduce_bc(r_up, r1_up,q_up, qq_up):


 b_up=qq_up*r1_up
 b_up.permute([54,18,20,6,4],3)


 c_up=q_up*r_up
 c_up.permute([52,14,12,19,17],3)


 return b_up, c_up







def reproduce_a(r2_up, qqq_up):

 a_up=qqq_up*r2_up
 a_up.permute([53,16,17,18,2],3)
 return  a_up 

#def equall_dis(a_up, b_up):

# A=copy.copy(b_up)
# A.setLabel([54,18,20,6,4])
# A.permute([18,20,6,4,54],1)
# 
# l,qq=lq_parity1(A) 
# 
# l.setLabel([18,80])
# qq.setLabel([80,20,6,4,54])

# A=copy.copy(a_up)
# A.setLabel([55,16,17,18,2])
# A.permute([2,16,17,55,18],4)
# 
# q,r=qr_parity1(A) 
# 
# q.setLabel([2,16,17,55,81])
# r.setLabel([81,18])

# Teta=l*r
# Teta.permute([81,80],1)
# U,s,V=svd_parity2(Teta)

# U.setLabel([81,18])
# s.setLabel([18,-18])
# V.setLabel([-18,80])

# s=Sqrt(s)
# 
# U=U*s
# V=s*V

# U.permute([81,-18],1)
# U.setLabel([81,18])
# V.permute([18,80],1)
# 
# b_up=qq*V
# a_up=q*U
# b_up.permute([54,18,20,6,4],3)
# a_up.permute([55,16,17,18,2],3)
# 
# return a_up,b_up






