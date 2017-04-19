import pyUni10 as uni10
import copy
import time
import basic


def Var_cab(c_u, a_u, b_u,a,b,c,d,Env,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic.Init_env(Env)

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 #t0=time.time()
 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMFull(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 #print time.time() - t0, "CTM-H, Left"

 Env=basic.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)


 #t0=time.time()
 E1, E2, E3, E4, E5, E6, E7, E8=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 #print time.time() - t0, "Env, Left"

 
 Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
 print "Dis", Dis_val
 
 #plist=Reload_plist(plist)
 
 
 #Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
 Do_optimization_Full(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)

 Store_plist(plist)
 
 c_up,a_up,b_up=recover(c_u,a_u,b_u,plist, MPO_list)
 
 c_up, b_up, a_up=equall_dis(c_up,b_up,a_up) 
 
 Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
 print "Dis_f", Dis_val
 
 Dis_val=Dis_final(E1, E2, E3, E4, E5, E6, E7, E8,c_u,a_u,b_u,c_up,a_up,b_up,d,MPO_list)
 print "Dis_final", Dis_val
 
 c_up=basic.max_ten(c_up)
 a_up=basic.max_ten(a_up)
 b_up=basic.max_ten(b_up)

 cp=basic.make_ab(c_up)
 ap=basic.make_ab(a_up)
 bp=basic.make_ab(b_up)

 return c_up, a_up, b_up, cp, ap, bp


def Var_acd(a_u, c_u, d_u,a,b,c,d,Env,D,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=basic.Init_env(Env)

 #Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 #t0=time.time()
 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=basic.corner_transfer_matrix_twosite_CTMRG(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
 if Corner_method is'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,Truncation=basic.corner_transfer_matrix_twosite_CTMFull(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D, Truncation)
 #print time.time() - t0, "CTM-H, Left"

 Env=basic.reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)



 #t0=time.time()
 E1, E2, E3, E4, E5, E6, E7, E8=produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)
 
 Dis_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
 print "Dis", Dis_val
 
 Do_optimization_Full1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
 #Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)

 Store_plist1(plist)
 
 
 a_up,c_up,d_up=recover1( a_u, c_u, d_u, plist, MPO_list)


 a_up, c_up, d_up=equall_dis1(a_up, c_up, d_up) 

 
 Dis_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
 print "Disf", Dis_val
 
 Dis_val=Dis_final1(E1, E2, E3, E4, E5, E6, E7, E8,a_u,c_u,d_u,a_up,c_up,d_up,b,MPO_list)
 print "Dis_final", Dis_val
 
 
 c_up=basic.max_ten(c_up)
 a_up=basic.max_ten(a_up)
 d_up=basic.max_ten(d_up)
 
 
 ap=basic.make_ab(a_up)
 dp=basic.make_ab(d_up)
 cp=basic.make_ab(c_up)

 return a_up, c_up, d_up, ap, cp, dp




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


  
def initialize_plist(a_u, b_u, c_u, MPO_list): 

 plist=[]
 bd1=uni10.Bond(uni10.BD_IN,c_u.bond(4).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(1).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(2).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,c_u.bond(4).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([62,58,57,17])
 Uni_ten.randomize()
 Uni_ten.identity()
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
 Uni_ten.identity()
 plist.append(Uni_ten)
 
 plist.append(copy.copy(plist[2]))
 plist[3].transpose()
 plist[3].setLabel([18,68,59,60])

 return plist
def initialize_plist1(a_u, c_u, d_u, MPO_list): 

 plist=[]
 
 bd1=uni10.Bond(uni10.BD_IN,c_u.bond(3).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[1].bond(3).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[1].bond(4).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,c_u.bond(3).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([62,59,60,19])
 Uni_ten.identity()
 plist.append(Uni_ten)

 plist.append(copy.copy(plist[0]))
 plist[1].transpose()
 plist[1].setLabel([19,64,59,60])

 
 bd1=uni10.Bond(uni10.BD_IN,c_u.bond(4).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,c_u.bond(4).Qlist())
 bd2=uni10.Bond(uni10.BD_OUT,MPO_list[0].bond(1).Qlist())
 bd3=uni10.Bond(uni10.BD_OUT,MPO_list[0].bond(2).Qlist())
 Uni_ten=uni10.UniTensor([bd1,bd4,bd2,bd3])
 Uni_ten.setLabel([66,17,58,57])
 Uni_ten.identity()
 plist.append(Uni_ten)


 plist.append(copy.copy(plist[2]))
 plist[3].transpose()
 plist[3].setLabel([17,58,57,68])

 return plist

 
def energy_cab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


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

 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)
 
 c_d=copy.copy(c_u)
 b_d=copy.copy(b_u)
 a_d=copy.copy(a_u)

 c_d.setLabel([-54,-14,-12,-19,-17])
 a_d.setLabel([-55,-16,-17,-18,-2])
 b_d.setLabel([-56,-18,-20,-6,-4])

 c_u.setLabel([-54,14,12,19,17])
 a_u.setLabel([-55,16,17,18,2])
 b_u.setLabel([-56,18,20,6,4])

 Norm=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E2*E3)*(b_u*b_d))))*((E4*E5)*d)

 Val=(((((E1*E8)*(a_u1*a_d))*((E7*E6)*(c_u1*c_d))))*(((E2*E3)*(b_u1*b_d))))*((E4*E5)*d)

 print  "Energy1", Val, Norm
 return  Val[0]/Norm[0]

def energy_acd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u):

 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])
 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
  
##########################################################
 a_u.setLabel([-54,16,17,18,2])
 c_u.setLabel([-55,14,12,19,17])
 d_u.setLabel([-56,19,10,8,20])


 a_d=copy.copy(a_u)
 c_d=copy.copy(c_u)
 d_d=copy.copy(d_u)

 a_d.setLabel([-54,-16,-17,-18,-2])
 c_d.setLabel([-55,-14,-12,-19,-17])
 d_d.setLabel([-56,-19,-10,-8,-20])

 Norm=(((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E2*E3)*(b))))*((E4*E5)*(d_u*d_d))

 Val=(((((E1*E8)*(a_u1*a_d))*((E7*E6)*(c_u1*c_d))))*(((E2*E3)*(b))))*((E4*E5)*(d_u1*d_d))

 print  "Energy2", Val, Norm
 return  Val[0]/Norm[0]


 
 
def  Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):

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

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])

##########################################################
 a_u.setLabel([55,16,64,66,2])
 b_u.setLabel([56,68,20,6,4])
 c_u.setLabel([54,14,12,19,62])


 c_ut=((c_u*(plist[0]*MPO_list[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-17,-18,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])


 Val=(((((E1*E8)*(a_ut*a_dt))*((E7*E6)*(c_ut*c_dt))))*(((E2*E3)*(b_ut*b_dt))))*((E4*E5)*d)
 #print 'Val=',Val
 #Val1=(((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c_u1*c_d1))))*((E2*E3)*(b_u1*b_d1)))*(((E4*E5)*d))
 #print 'Val1=',Val1
 Val2=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E2*E3)*(b_ut*b_d1))))*((E4*E5)*d)
 #print 'Val2=',Val2
 #Val3=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E2*E3)*(b_u1*b_dt))))*((E4*E5)*d)
 #print 'Val3=',Val3

 #val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 val_f=abs(Val[0])-2.00*Val2[0]
 return  val_f




def  Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):

 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 c_d1=copy.copy(c_u1)
 d_d1=copy.copy(d_u1)

 a_d1.setLabel([-54,-16,-17,-58,-57,-18,-2])
 c_d1.setLabel([-55,-14,-12,-58,-57,-19,-17,-59,-60])
 d_d1.setLabel([-56,-19,-10,-59,-60,-8,-20])
 
 
 
 
##########################################################
 a_u.setLabel([54,16,68,18,2])
 c_u.setLabel([55,14,12,62,66])
 d_u.setLabel([56,64,10,8,20])


 a_ut=((a_u*(plist[3]*MPO_list[0])))
 c_ut=(c_u*(plist[0]*plist[2]*MPO_list[1]))
 d_ut=((d_u*(plist[1]*MPO_list[2])))


 a_ut.permute([-54,16,17,18,2],3)
 c_ut.permute([-55,14,12,19,17],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 c_dt=copy.copy(c_ut)
 d_dt=copy.copy(d_ut)

 a_dt.setLabel([-54,-16,-17,-18,-2])
 c_dt.setLabel([-55,-14,-12,-19,-17])
 d_dt.setLabel([-56,-19,-10,-8,-20])


 Val=((((E1*E8)*(a_ut*a_dt))*((E7*E6)*(c_ut*c_dt))))*(((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))
 #print 'Val=',Val
 Val1=((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c_u1*c_d1))))*(((E4*E5)*(d_d1*d_u1))*((E2*E3)*(b)))
 #print 'Val1=',Val1
 Val2=((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b)))
 #print 'Val2=',Val2
 Val3=((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b)))
 #print 'Val3=',Val3
 
 val_f=abs(Val[0])-2.0*Val2[0]
 val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 
 
 return  val_f



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

def Do_optimization_Full(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):
 
 Res=1
 Res1=2
 count=0
 Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
 print 'Dis0', Distance_val
 for q in xrange(20):
  optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist)
  optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist)
  optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist)
  optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,  MPO_list,c_u,a_u,b_u,plist)


  Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
  Res=Res1
  Res1=Distance_val
  print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q
  
  count+=1
  if count > 50: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break
 Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
 print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q



def Do_optimization_Full1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):
 
 Res=1
 Res1=2
 count=0
 Distance_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
 for q in xrange(20):
  print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q
  optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
  optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
  optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
  optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)


  Distance_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
  Res=Res1
  Res1=Distance_val
  count+=1
  if count > 50: print 'Num_Opt > 50'; break;
  if abs(Res) > 1.00e-10:
   if (abs(Distance_val) < 1.00e-7) or ((abs(Res1-Res) / abs(Res)) < 8.00e-9): 
    print 'break, Dis', Distance_val, (abs(Res1-Res) / abs(Res)), count
    break
  else:
    if (abs(Distance_val[0]) < 1.00e-7) or (  abs(Res1-Res) < 1.00e-11  ): 
     #print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break



def optimum_00(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):
 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 c_d1=copy.copy(c_u1)
 d_d1=copy.copy(d_u1)

 a_d1.setLabel([-54,-16,-17,-58,-57,-18,-2])
 c_d1.setLabel([-55,-14,-12,-58,-57,-19,-17,-59,-60])
 d_d1.setLabel([-56,-19,-10,-59,-60,-8,-20])
 
 
 
 
##########################################################
 a_u.setLabel([54,16,68,18,2])
 c_u.setLabel([55,14,12,62,66])
 d_u.setLabel([56,64,10,8,20])


 a_ut=((a_u*(plist[3]*MPO_list[0])))
 c_ut=(c_u*(plist[2]*MPO_list[1]))
 d_ut=((d_u*(plist[1]*MPO_list[2])))
 
 a_ut.permute([-54,16,17,18,2],3)
 c_ut.permute([-55,14,12,62,59,60,17],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 c_dt=copy.copy(c_ut)
 d_dt=copy.copy(d_ut)

 a_dt.setLabel([-54,-16,-17,-18,-2])
 c_dt.setLabel([-55,-14,-12,-62,-59,-60,-17])
 d_dt.setLabel([-56,-19,-10,-8,-20])
 


 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))*(((E1*E8)*(a_ut*a_dt))))*((E7*E6)*(c_ut*c_dt))


 A.permute([-62,-59,-60,-19,62,59,60,19],4)

 A1=copy.copy(A)
 A1.transpose()
 A=A+A1
 
 Ap=((((E1*E8)*(a_u1*a_dt))))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-59,-60,-19],4)


 A1=((((E1*E8)*(a_ut*a_d1))))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b)))*(((E7*E6)*(c_ut*c_d1)))
 A1.permute([62,59,60,19],0)
 A1.transpose()
 Ap=Ap+A1

 U, S, V=svd_parity(A)
 #print S
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A_inv=V*S*U
 A_inv.permute([0,1,2,3,12,13,14,15],4)
 A_inv.setLabel([62,59,60,19,-62,-59,-60,-19])

#######################################################
 A=A_inv*Ap
 A.permute([62,59,60,19],3)
 plist[0]=A



def optimum_11(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):
 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 c_d1=copy.copy(c_u1)
 d_d1=copy.copy(d_u1)

 a_d1.setLabel([-54,-16,-17,-58,-57,-18,-2])
 c_d1.setLabel([-55,-14,-12,-58,-57,-19,-17,-59,-60])
 d_d1.setLabel([-56,-19,-10,-59,-60,-8,-20])
 
 
 
 
##########################################################
 a_u.setLabel([54,16,68,18,2])
 c_u.setLabel([55,14,12,62,66])
 d_u.setLabel([56,64,10,8,20])


 a_ut=((a_u*(plist[3]*MPO_list[0])))
 c_ut=(c_u*(plist[2]*plist[0]*MPO_list[1]))
 d_ut=((d_u*(MPO_list[2])))

 a_ut.permute([-54,16,17,18,2],3)
 c_ut.permute([-55,14,12,19,17],3)
 d_ut.permute([-56,64,59,60,10,8,20],4)

 a_dt=copy.copy(a_ut)
 c_dt=copy.copy(c_ut)
 d_dt=copy.copy(d_ut)

 a_dt.setLabel([-54,-16,-17,-18,-2])
 c_dt.setLabel([-55,-14,-12,-19,-17])
 d_dt.setLabel([-56,-64,-59,-60,-10,-8,-20])


 A=((((E1*E8)*(a_ut*a_dt))*((E2*E3)*(b)))*((E7*E6)*(c_ut*c_dt)))*((E4*E5)*(d_ut*d_dt))

 A.permute([-19,-64,-59,-60,19,64,59,60],4)

 A1=copy.copy(A)
 A1.transpose()
 A=A+A1



 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E1*E8)*(a_u1*a_dt))*((E2*E3)*(b))))*((E4*E5)*(d_u1*d_dt))
 Ap.permute([-19,-64,-59,-60],4)



 A1=(((((E7*E6)*(c_ut*c_d1))))*(((E1*E8)*(a_ut*a_d1))*((E2*E3)*(b))))*((E4*E5)*(d_ut*d_d1))
 A1.permute([19,64,59,60],0)
 A1.transpose()
 Ap=Ap+A1
 #Ap.setLabel([-17,-58,-57,-64])

 U, S, V=svd_parity(A)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A_inv=V*S*U
 A_inv.permute([0,1,2,3,12,13,14,15],4)
 A_inv.setLabel([19,64,59,60,-19,-64,-59,-60])

 A=A_inv*Ap
 A.permute([19,64,59,60],1)
 plist[1]=A

def optimum_22(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):
 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])
 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 c_d1=copy.copy(c_u1)
 d_d1=copy.copy(d_u1)

 a_d1.setLabel([-54,-16,-17,-58,-57,-18,-2])
 c_d1.setLabel([-55,-14,-12,-58,-57,-19,-17,-59,-60])
 d_d1.setLabel([-56,-19,-10,-59,-60,-8,-20])
 
##########################################################
 a_u.setLabel([54,16,68,18,2])
 c_u.setLabel([55,14,12,62,66])
 d_u.setLabel([56,64,10,8,20])


 a_ut=((a_u*(plist[3]*MPO_list[0])))
 c_ut=(c_u*(plist[0]*MPO_list[1]))
 d_ut=((d_u*(plist[1]*MPO_list[2])))

 a_ut.permute([-54,16,17,18,2],3)
 c_ut.permute([-55,14,12,58,57,19,66],5)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 c_dt=copy.copy(c_ut)
 d_dt=copy.copy(d_ut)

 a_dt.setLabel([-54,-16,-17,-18,-2])
 c_dt.setLabel([-55,-14,-12,-58,-57,-19,-66])
 d_dt.setLabel([-56,-19,-10,-8,-20])


 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))*((E1*E8)*(a_ut*a_dt)))*((E7*E6)*(c_ut*c_dt))

 A.permute([-66,-58,-57,-17,66,58,57,17],4)
 A1=copy.copy(A)
 A1.transpose()
 A=A+A1



 Ap=((((E1*E8)*(a_u1*a_dt)))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b))))*(((E7*E6)*(c_u1*c_dt)))
 
 Ap.permute([-66,-58,-57,-17],4)


 A1=((((E1*E8)*(a_ut*a_d1)))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b))))*(((E7*E6)*(c_ut*c_d1)))
 A1.permute([66,58,57,17],0)
 A1.transpose()
 Ap=Ap+A1


 U, S, V=svd_parity(A)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A_inv=V*S*U
 A_inv.permute([0,1,2,3,12,13,14,15],4)
 A_inv.setLabel([66,58,57,17,-66,-58,-57,-17])

 A=A_inv*Ap
 A.permute([66,17,58,57],1)
 plist[2]=A


def optimum_33(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):
 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])


 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 c_d1=copy.copy(c_u1)
 d_d1=copy.copy(d_u1)

 a_d1.setLabel([-54,-16,-17,-58,-57,-18,-2])
 c_d1.setLabel([-55,-14,-12,-58,-57,-19,-17,-59,-60])
 d_d1.setLabel([-56,-19,-10,-59,-60,-8,-20])
 
 
 
 
##########################################################
 a_u.setLabel([54,16,68,18,2])
 c_u.setLabel([55,14,12,62,66])
 d_u.setLabel([56,64,10,8,20])


 a_ut=((a_u*(MPO_list[0])))
 c_ut=(c_u*(plist[0]*plist[2]*MPO_list[1]))
 d_ut=((d_u*(plist[1]*MPO_list[2])))
 
 a_ut.permute([-54,16,68,58,57,18,2],3)
 c_ut.permute([-55,14,12,19,17],3)
 d_ut.permute([-56,19,10,8,20],3)

 a_dt=copy.copy(a_ut)
 c_dt=copy.copy(c_ut)
 d_dt=copy.copy(d_ut)

 a_dt.setLabel([-54,-16,-68,-58,-57,-18,-2])
 c_dt.setLabel([-55,-14,-12,-19,-17])
 d_dt.setLabel([-56,-19,-10,-8,-20])
 

 A=((((E4*E5)*(d_ut*d_dt))*((E2*E3)*(b)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-17,-58,-57,-68,17,58,57,68],4)
 A1=copy.copy(A)
 A1.transpose()
 A=A+A1


 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*(d_u1*d_dt))*((E2*E3)*(b))))*((E1*E8)*(a_u1*a_dt))
 
 Ap.permute([-17,-58,-57,-68],4)


 A1=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*(d_ut*d_d1))*((E2*E3)*(b))))*((E1*E8)*(a_ut*a_d1))
 A1.permute([17,58,57,68],0)
 A1.transpose()
 Ap=Ap+A1



 U, S, V=svd_parity(A)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A_inv=V*S*U
 A_inv.permute([0,1,2,3,12,13,14,15],4)
 A_inv.setLabel([17,58,57,68,-17,-58,-57,-68])

 A=A_inv*Ap
 
 A.permute([17,58,57,68],3)

 plist[3]=A




def optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])

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

 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])
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

 c_ut=((c_u*(MPO_list[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,62,58,57],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-62,-58,-57])
 a_dt.setLabel([-55,-16,-17,-18,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])

 A=((((E1*E8)*(a_ut*a_dt))*((E2*E3)*(b_ut*b_dt)))*(((E4*E5)*d)))*((E7*E6)*(c_ut*c_dt))
 A.permute([-62,-58,-57,-17,62,58,57,17],4)

 A1=copy.copy(A)
 A1.transpose()
 A=A+A1
 
# A=N_Positiv(A)
# A.setLabel([-62,-58,-57,-17,62,58,57,17])



 Ap=(((((E1*E8)*(a_u1*a_dt))) * ((E2*E3)*(b_u1*b_dt)))*(((E4*E5)*d)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-58,-57,-17],4)


 Ap1=(((((E1*E8)*(a_ut*a_d1)))*((E2*E3)*(b_ut*b_d1)))*(((E4*E5)*d)))*(((E7*E6)*(c_ut*c_d1)))
 Ap1.permute([62,58,57,17],0)
 Ap1.transpose()
 Ap=Ap+Ap1

 U, S, V=svd_parity(A)
 #print S
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A_inv=V*S*U
 A_inv.permute([0,1,2,3,12,13,14,15],4)
 A_inv.setLabel([62,58,57,17,-62,-58,-57,-17])

###################################################
# A_inv1=copy.copy(A_inv)
# A_inv1.setLabel([0,1,2,3,-62,-58,-57,-17])
# Res=A_inv1*A
# Res.permute([0,1,2,3,62,58,57,17],4)
# #print A, Res
# M=A.getBlock()
# #M1=M.inverse()
# svd=M.svd()
# s=basic_FU.Inv_mat(svd[1])
# svd[0].transpose()
# svd[2].transpose()
# M_inv=svd[2]*s*svd[0]
# #print M_inv*M
#######################################################
 A=A_inv*Ap
 A.permute([62,58,57,17],3)
 plist[0]=A



def optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])



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


 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])


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

 c_ut=((c_u*(plist[0]*MPO_list[0])))
 a_ut=(a_u*(plist[2]*MPO_list[1]))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,64,58,57,18,2],5)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-64,-58,-57,-18,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])


 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-17,-64,-58,-57,17,64,58,57],4)

 A1=copy.copy(A)
 A1.transpose()
 A=A+A1

# A=N_Positiv(A)
# A.setLabel([-17,-64,-58,-57,17,64,58,57])


 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-17,-64,-58,-57],4)


 A1=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 A1.permute([17,64,58,57],0)
 A1.transpose()
 Ap=Ap+A1
 #Ap.setLabel([-17,-58,-57,-64])

 U, S, V=svd_parity(A)
 U.transpose()
 V.transpose()
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


def optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):
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

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])


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

 a_ut=(a_u*(plist[1]*MPO_list0[1]))
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,66,59,60,2],3)
 b_ut.permute([-56,18,20,6,4],3)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-17,-66,-59,-60,-2])
 b_dt.setLabel([-56,-18,-20,-6,-4])

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-66,-59,-60,-18,66,59,60,18],4)
 A1=copy.copy(A)
 A1.transpose()
 A=A+A1

# A=N_Positiv(A)
# A.setLabel([-66,-59,-60,-18,66,59,60,18])

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 
 Ap.permute([-66,-59,-60,-18],4)

 A1=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 A1.permute([66,59,60,18],0)
 A1.transpose()
 Ap=Ap+A1


 U, S, V=svd_parity(A)
 U.transpose()
 V.transpose()
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


def optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,MPO_list,c_u,a_u,b_u,plist):
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

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])


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

 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(MPO_list0[2])))

 c_ut.permute([-54,14,12,19,17],3)
 a_ut.permute([-55,16,17,18,2],3)
 b_ut.permute([-56,68,59,60,20,6,4],5)

 c_dt=copy.copy(c_ut)
 b_dt=copy.copy(b_ut)
 a_dt=copy.copy(a_ut)

 c_dt.setLabel([-54,-14,-12,-19,-17])
 a_dt.setLabel([-55,-16,-17,-18,-2])
 b_dt.setLabel([-56,-68,-59,-60,-20,-6,-4])

 A=((((E1*E8)*(a_ut*a_dt)) *((E7*E6)*(c_ut*c_dt))) * ((((E4*E5)*d)))) * ((E2*E3)*(b_ut*b_dt))


 A.permute([-18,-68,-59,-60,18,68,59,60],4)
 A1=copy.copy(A)
 A1.transpose()
 A=A+A1

# A=N_Positiv(A)
# A.setLabel([-18,-68,-59,-60,18,68,59,60])


 Ap=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-68,-59,-60],4)


 A1=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)))*((E2*E3)*(b_ut*b_d1))
 
 A1.permute([18,68,59,60],0)
 A1.transpose()
 Ap=Ap+A1



 U, S, V=svd_parity(A)
 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([8,9,10,11,12,13,14,15])
 S.setLabel([4,5,6,7,8,9,10,11])
 V.setLabel([0,1,2,3,4,5,6,7])


 A_inv=V*S*U
 A_inv.permute([0,1,2,3,12,13,14,15],4)
 A_inv.setLabel([18,68,59,60,-18,-68,-59,-60])

 A=A_inv*Ap
 
 A.permute([18,68,59,60],1)
 plist[3]=A
 
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
 
def  recover1( a_u,c_u, d_u, plist, MPO_list):

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56]) 
 a_u.setLabel([54,16,68,18,2])
 c_u.setLabel([55,14,12,62,66])
 d_u.setLabel([56,64,10,8,20])

 a_ut=((a_u*(plist[3]*MPO_list[0])))
 c_ut=(c_u*(plist[0]*plist[2]*MPO_list[1]))
 d_ut=((d_u*(plist[1]*MPO_list[2])))

 a_ut.permute([-54,16,17,18,2],3)
 c_ut.permute([-55,14,12,19,17],3)
 d_ut.permute([-56,19,10,8,20],3)
 
 return a_ut, c_ut,d_ut 
 
 
 
def  Dis_final(E1, E2, E3, E4, E5, E6, E7, E8,c_u,a_u,b_u,c_up,a_up,b_up,d,MPO_list):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])



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

 c_u1=((c_u1*(MPO_list[0])))
 a_u1=(a_u1*(MPO_list[1]))
 b_u1=((b_u1*(MPO_list[2])))

 c_u1.permute([-54,14,12,58,57,19,17],3)
 a_u1.permute([-55,16,17,58,57,18,2,59,60],5)
 b_u1.permute([-56,18,20,59,60,6,4],5)

 c_d1=copy.copy(c_u1)
 b_d1=copy.copy(b_u1)
 a_d1=copy.copy(a_u1)

 c_d1.setLabel([-54,-14,-12,-58,-57,-19,-17])
 a_d1.setLabel([-55,-16,-17,-58,-57,-18,-2,-59,-60])
 b_d1.setLabel([-56,-18,-20,-59,-60,-6,-4])

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

def  Dis_final1(E1, E2, E3, E4, E5, E6, E7, E8,a_u,c_u,d_u,a_up,c_up,d_up,b,MPO_list):

 b.setLabel([18,-18,20,-20,6,-6,4,-4])

 MPO_list[0].setLabel([-54,58,57,54])
 MPO_list[1].setLabel([58,57,-55,59,60,55])
 MPO_list[2].setLabel([60,59,-56,56])
 
 a_u1=copy.copy(a_u)
 d_u1=copy.copy(d_u)
 c_u1=copy.copy(c_u)

###########################################################
 a_u1.setLabel([54,16,17,18,2])
 c_u1.setLabel([55,14,12,19,17])
 d_u1.setLabel([56,19,10,8,20])

  
 a_u1=(a_u1*(MPO_list[0]))
 c_u1=((c_u1*(MPO_list[1])))
 d_u1=((d_u1*(MPO_list[2]))) 
 
 a_u1.permute([-54,16,17,58,57,18,2],3)
 c_u1.permute([-55,14,12,58,57,19,17,59,60],5)
 d_u1.permute([-56,19,10,59,60,8,20],5)
 
 a_d1=copy.copy(a_u1)
 c_d1=copy.copy(c_u1)
 d_d1=copy.copy(d_u1)

 a_d1.setLabel([-54,-16,-17,-58,-57,-18,-2])
 c_d1.setLabel([-55,-14,-12,-58,-57,-19,-17,-59,-60])
 d_d1.setLabel([-56,-19,-10,-59,-60,-8,-20])
 
 
##########################################################
 a_up.setLabel([-54,16,17,18,2])
 c_up.setLabel([-55,14,12,19,17])
 d_up.setLabel([-56,19,10,8,20])



 a_up.permute([-54,16,17,18,2],3)
 c_up.permute([-55,14,12,19,17],3)
 d_up.permute([-56,19,10,8,20],3)

 a_dp=copy.copy(a_up)
 c_dp=copy.copy(c_up)
 d_dp=copy.copy(d_up)

 a_dp.setLabel([-54,-16,-17,-18,-2])
 c_dp.setLabel([-55,-14,-12,-19,-17])
 d_dp.setLabel([-56,-19,-10,-8,-20])

 Val=((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E4*E5)*(d_up*d_dp))*((E2*E3)*(b)))

 Val1=((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c_u1*c_d1))))*(((E4*E5)*(d_u1*d_d1))*((E2*E3)*(b)))

 Val2=((((E1*E8)*(a_up*a_d1))*((E7*E6)*(c_up*c_d1))))*(((E4*E5)*(d_up*d_d1))*((E2*E3)*(b)))

 Val3=((((E1*E8)*(a_u1*a_dp))*((E7*E6)*(c_u1*c_dp))))*(((E4*E5)*(d_u1*d_dp))*((E2*E3)*(b)))

 
 val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=abs(Val[0])-2.0*Val2[0]#-Val3[0]

 return  val_f 



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

def N_Positiv(N):
 N.setLabel([-1,-2,-3,-4,1,2,3,4])
 N1=copy.copy(N)
 N1.transpose()
 N=(N+N1)*(1.00/2.00)
 N1=copy.copy(N)
 N1.setLabel( [1,2,3,4,5,6,7,8] )
 N=N*N1
 N.permute([-1,-2,-3,-4, 5,6,7,8 ],4)
 N_final=sqrt_general(N)
 N_final.setLabel([-1,-2,-3,-4, 5,6,7,8])
 return N_final              

def Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):
  Es=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u, plist)
  Ef=0
  E2_val=0
  plist_first=[plist[i] for i in xrange(len(plist))]
  plist_tem=[0]*len(plist)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=0
  count=0
  for i in xrange(40):
   count+=1
   E1_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   Ef=E1_val
   Store_plist(plist)

   print 'E=', E1_val, count
   D_rf=basic.Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist)
   D_r=[ (-1.00)*D_rf[i] for i in xrange(len(D_rf)) ]

   A=D_r[0]*D_r[0]
   B=D_r[1]*D_r[1]
   C=D_r[2]*D_r[2]
   D=D_r[3]*D_r[3]
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-15:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break
      
   E_previous=E1_val
   
   if Norm_Z < 1.0e-16:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    count+=1
    plist_tem[0]=plist[0]+(2.00)*Gamma*D_r[0]
    plist_tem[1]=plist[1]+(2.00)*Gamma*D_r[1]
    plist_tem[2]=plist[2]+(2.00)*Gamma*D_r[2]
    plist_tem[3]=plist[3]+(2.00)*Gamma*D_r[3]
    E2_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist_tem)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 :
     print "break", E1_val, E2_val, Gamma
     Gamma=0
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
    E2_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist_tem)
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

  if( Ef > Es):
   print 'SD method, Fail, f<s', Ef, Es 
   for i in xrange(len(plist_first)):
    plist[i]=plist_first[i]
    
     
def Do_optimization_Grad1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist):
  Es=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
  Ef=0
  E2_val=0
  plist_first=[plist[i] for i in xrange(len(plist))]
  plist_tem=[0]*len(plist)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=0
  count=0
  for i in xrange(40):
   count+=1
   E1_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
   Ef=E1_val
   Store_plist(plist)

   print 'E=', E1_val, count
   D_rf=basic.Obtain_grad_four1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist)
   D_r=[ (-1.00)*D_rf[i] for i in xrange(len(D_rf)) ]

   A=D_r[0]*D_r[0]
   B=D_r[1]*D_r[1]
   C=D_r[2]*D_r[2]
   D=D_r[3]*D_r[3]
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-15:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break
      
   E_previous=E1_val
   
   if Norm_Z < 1.0e-16:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    count+=1
    plist_tem[0]=plist[0]+(2.00)*Gamma*D_r[0]
    plist_tem[1]=plist[1]+(2.00)*Gamma*D_r[1]
    plist_tem[2]=plist[2]+(2.00)*Gamma*D_r[2]
    plist_tem[3]=plist[3]+(2.00)*Gamma*D_r[3]
    E2_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist_tem)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 :
     print "break", E1_val, E2_val, Gamma
     Gamma=0
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
    E2_val=Dis_f1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,a_u,c_u,d_u,plist_tem)
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

  if( Ef > Es):
   print 'SD method, Fail, f<s', Ef, Es 
   for i in xrange(len(plist_first)):
    plist[i]=plist_first[i] 

