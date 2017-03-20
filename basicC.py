import pyUni10 as uni10
import copy
import time
import basic


def Var_cab(c_u, a_u, b_u,a,b,c,d,Env,D,U,MPO_list,d_phys,chi,Gauge,Positive,Corner_method,plist, plistd):



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
 test_env(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d)

 test_energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u)
 #E=energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,a_u,b_u,c_u)
 #print "Energy", E
 
 Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
 #print "Dis", Dis_val
 
 plist,plistd=Reload_plist(plist,plistd)
 
 #Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
 
 Do_optimization_Full(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)


 Store_plist(plist,plistd)


 #print plist[3],plist[1],plist[2],plist[0]
 
 cp_u,ap_u,bp_u=recover(c_u,a_u,b_u,plist, plistd, MPO_list)
 
 Dis_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
 print "Dis", Dis_val
 
 Dis_val=Dis_final(E1, E2, E3, E4, E5, E6, E7, E8,c_u,a_u,b_u,cp_u,ap_u,bp_u,d,MPO_list)
 print "Dis_final", Dis_val
 
 


 cp_u=basic.max_ten(cp_u)
 ap_u=basic.max_ten(ap_u)
 bp_u=basic.max_ten(bp_u)

 cp=basic.make_ab(cp_u)
 ap=basic.make_ab(ap_u)
 bp=basic.make_ab(bp_u)



 return cp_u, ap_u, bp_u, cp, ap, bp

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
def  test_energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])


 U.setLabel([51,52,53,54,55,56])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54,55,56])

 MPO_list[0].setLabel([51,58,57,54])
 MPO_list[1].setLabel([58,57,52,59,60,55])
 MPO_list[2].setLabel([60,59,53,56])
 H=MPO_list[0]*MPO_list[1]*MPO_list[2]
 H.permute([51,52,53,54,55,56],3)
 print "test-MPO", H.elemCmp(U) 
 
 
 a_u.setLabel([55,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([52,-16,-17,-18,-2])


 b_u.setLabel([56,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.setLabel([53,-18,-20,-6,-4])


 c_u.setLabel([54,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.setLabel([51,-14,-12,-19,-17])


 Val=((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d)))*Iden)*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 Norm_f=Val
 print 'Env_test=', Val[0]

 Val=((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d)))*U)*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 print 'E=',Val[0]/Norm_f[0]
 
 Val=((((E1*E8)*(a_u*a_d*MPO_list[1]))*((E7*E6)*(c_u*c_d*MPO_list[0]))))*(((E4*E5)*d)*((E2*E3)*(b_u*b_d*MPO_list[2])))
 print 'E=',Val[0]/Norm_f[0]
 
def  energy(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,a_u,b_u,c_u):
 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 U.setLabel([51,52,53,54,55,56])

 MPO_list[0].setLabel([51,58,57,54])
 MPO_list[1].setLabel([58,57,52,59,60,55])
 MPO_list[2].setLabel([60,59,53,56])
 H=MPO_list[0]*MPO_list[1]*MPO_list[2]
 H.permute([51,52,53,54,55,56],3)
 #print "test-MPO", H.elemCmp(U), H 
 
 
 a_u.setLabel([55,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([55,-16,-17,-18,-2])


 b_u.setLabel([56,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.setLabel([56,-18,-20,-6,-4])


 c_u.setLabel([54,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.setLabel([54,-14,-12,-19,-17])
 

 Val=((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))
 Norm_f=Val

 a_u.setLabel([55,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([52,-16,-17,-18,-2])


 b_u.setLabel([56,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.setLabel([53,-18,-20,-6,-4])


 c_u.setLabel([54,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.setLabel([51,-14,-12,-19,-17])
 
 Val=((((E1*E8)*(a_u*a_d*MPO_list[1]))*((E7*E6)*(c_u*c_d*MPO_list[0]))))*(((E4*E5)*d)*((E2*E3)*(b_u*b_d*MPO_list[2])))

 a_u.setLabel([55,16,17,18,2])
 a_d=copy.copy(a_u)
 a_d.setLabel([52,-16,-17,-18,-2])


 b_u.setLabel([56,18,20,6,4])
 b_d=copy.copy(b_u)
 b_d.setLabel([53,-18,-20,-6,-4])


 c_u.setLabel([54,14,12,19,17])
 c_d=copy.copy(c_u)
 c_d.setLabel([51,-14,-12,-19,-17])
 
 Val2=((((E1*E8)*(a_u*a_d))*((E7*E6)*(c_u*c_d))))*(((E4*E5)*d)*((E2*E3)*(b_u*b_d)))*U



 A_val=Val[0]*(1.00/Norm_f[0]) 
 print Norm_f[0] ,Val2[0]*(1.00/Norm_f[0]), Val2
 return A_val
  
def initialize_plist(a_u, b_u, c_u, MPO_list): 

 plist=[]
 plistd=[]
 
 bd1=uni10.Bond(uni10.BD_IN,c_u.bond(4).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(1).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[0].bond(2).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,c_u.bond(4).Qlist())
 #bd4=uni10.Bond(uni10.BD_OUT,3*2*2)
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([62,58,57,17])
 Uni_ten.randomize()
 Uni_ten.identity()
 plist.append(Uni_ten)
 Uni_ten=copy.copy(Uni_ten)
 Uni_ten.setLabel([-62,-58,-57,-17])
 plistd.append(Uni_ten)


 
 bd1=uni10.Bond(uni10.BD_IN,a_u.bond(2).Qlist())
 bd2=uni10.Bond(uni10.BD_OUT,MPO_list[0].bond(1).Qlist())
 bd3=uni10.Bond(uni10.BD_OUT,MPO_list[0].bond(2).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,a_u.bond(2).Qlist())
 #bd1=uni10.Bond(uni10.BD_IN,12)
 Uni_ten=uni10.UniTensor([bd1,bd4,bd2,bd3])
 Uni_ten.setLabel([17,64,58,57])
 Uni_ten.randomize()
 Uni_ten.identity()
 plist.append(Uni_ten)
 Uni_ten=copy.copy(Uni_ten)
 Uni_ten.setLabel([-17,-64,-58,-57])
 plistd.append(Uni_ten)
 
 
 
 
 bd1=uni10.Bond(uni10.BD_IN,a_u.bond(3).Qlist())
 bd2=uni10.Bond(uni10.BD_IN,MPO_list[1].bond(3).Qlist())
 bd3=uni10.Bond(uni10.BD_IN,MPO_list[1].bond(4).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,a_u.bond(3).Qlist())
 #bd4=uni10.Bond(uni10.BD_OUT,3*2*2)
 Uni_ten=uni10.UniTensor([bd1,bd2,bd3,bd4])
 Uni_ten.setLabel([66,59,60,18])
 Uni_ten.randomize()
 Uni_ten.identity()
 plist.append(Uni_ten)
 Uni_ten=copy.copy(Uni_ten)
 Uni_ten.setLabel([-66,-59,-60,-18])
 plistd.append(Uni_ten)
 
 

 bd1=uni10.Bond(uni10.BD_IN,b_u.bond(2).Qlist())
 bd2=uni10.Bond(uni10.BD_OUT,MPO_list[1].bond(3).Qlist())
 bd3=uni10.Bond(uni10.BD_OUT,MPO_list[1].bond(4).Qlist())
 bd4=uni10.Bond(uni10.BD_OUT,b_u.bond(2).Qlist())
 #bd1=uni10.Bond(uni10.BD_IN,12)
 Uni_ten=uni10.UniTensor([bd1,bd4,bd2,bd3])
 Uni_ten.setLabel([18,68,59,60])
 Uni_ten.randomize()
 Uni_ten.identity()
 plist.append(Uni_ten)
 Uni_ten=copy.copy(Uni_ten)
 Uni_ten.setLabel([-18,-68,-59,-60])
 plistd.append(Uni_ten)
 return plist, plistd
 
 
 
def  Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):

 d.setLabel([19,-19,10,-10,8,-8,20,-20])

 U.setLabel([51,52,53,54,55,56])
 Iden=uni10.UniTensor(U.bond())
 Iden.identity()
 Iden.setLabel([51,52,53,54,55,56])
 Iden1=copy.copy(Iden)
 Iden1.setLabel([51,52,53,-54,-55,-56])



 MPO_list[0].setLabel([51,58,57,54])
 MPO_list[1].setLabel([58,57,52,59,60,55])
 MPO_list[2].setLabel([60,59,53,56])
 
 H=MPO_list[0]*MPO_list[1]*MPO_list[2]
 H.permute([51,52,53,54,55,56],3)
 #print "test-MPO", H.elemCmp(U) , U,H
  
 
 
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

 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

# print plistd[0],plistd[1], plistd[2], plistd[3]
# print plist[0],plist[1], plist[2], plist[3]

# print "test_a", a_ut.elemCmp(a_dt), a_ut.elemCmp(a_u1), a_u1.norm(),a_ut.norm()  
# print "test_c", c_ut.elemCmp(c_dt), c_ut.elemCmp(c_u1), c_u1.norm(),c_ut.norm()  
# print "test_b", b_ut.elemCmp(b_dt),b_u1.norm(),b_ut.norm()

 Val=(((((E1*E8)*(a_ut*a_dt))*((E7*E6)*(c_ut*c_dt))))*(((E2*E3)*(b_ut*b_dt))))*((E4*E5)*d)
 #print 'Val=',Val
 Val1=(((((E1*E8)*(a_u1*a_d1))*((E7*E6)*(c_u1*c_d1))))*((E2*E3)*(b_u1*b_d1)))*(((E4*E5)*d))
 #print 'Val1=',Val1
 Val2=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E2*E3)*(b_ut*b_d1))))*((E4*E5)*d)
 #print 'Val2=',Val2
 Val3=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E2*E3)*(b_u1*b_dt))))*((E4*E5)*d)
 #print 'Val3=',Val3

 val_f=Val[0]+Val1[0]-Val2[0]-Val3[0]
 #val_f=Val[0]-2.00*Val2[0]#-Val3[0]

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

def Do_optimization_Full(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):
 
 Res=1
 Res1=2
 count=0
 Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
 for q in xrange(20):
  print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q
  optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
  optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
  optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
  optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)


  Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
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
     print 'break, Dis', Distance_val[0], abs(Res1-Res)
     break
 Distance_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
 print 'Dis', Distance_val, abs(Res1-Res) / abs(Res), q

 return plist, plistd

def   make_plistd(plist,plistd):

 
 plistd[0]=copy.copy(plist[0])
 plistd[0].setLabel([-62,-58,-57,-17])
 
 plistd[1]=copy.copy(plist[1])
 plistd[1].setLabel([-17,-64,-58,-57])
 
 plistd[2]=copy.copy(plist[2])
 plistd[2].setLabel([-66,-59,-60,-18])

 plistd[3]=copy.copy(plist[3])
 plistd[3].setLabel([-18,-68,-59,-60])

 return plistd


def Do_optimization_Grad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):
  Es=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u, plist, plistd)
  Ef=0
  E2_val=0
  plist_first=[plist[i] for i in xrange(len(plist))]
  plistd_first=[plistd[i] for i in xrange(len(plistd))]
  plist_tem=[0]*len(plist)
  plistd_tem=[0]*len(plistd)
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=0
  count=0
  for i in xrange(200):
   count+=1
   E1_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
   Ef=E1_val
   Store_plist(plist,plistd)

   print 'E=', E1_val, count
   D_rf=basic_FU.Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
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
    make_plistd(plist_tem,plistd_tem)
    E2_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist_tem, plistd_tem)
    #print "Gamma", Gamma, E2_val, E1_val  
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+15 or  abs(Gamma)>1.0e+15 :
     print "break", E1_val, E2_val
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
    make_plistd(plist_tem,plistd_tem)
    E2_val=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist_tem, plistd_tem)
    #print "Gamma", Gamma, E2_val, E1_val
    v=Dis_f(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd)
    #print "v", v, plist[0], plist_tem[0]
    if abs((0.5)*Norm_Z*Gamma) <1.0e-15 or  abs(E1_val-E2_val)<1.0e-15 or abs(Gamma)<1.0e-15 :
     print "break", E1_val, E2_val
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0


   plist[0]=plist[0]+(1.00)*Gamma*D_r[0]
   plist[1]=plist[1]+(1.00)*Gamma*D_r[1]
   plist[2]=plist[2]+(1.00)*Gamma*D_r[2]
   plist[3]=plist[3]+(1.00)*Gamma*D_r[3]
   make_plistd(plist,plistd)

  if( Ef > Es):
   print 'SD method, Fail, f<s', Ef, Es 
   for i in xrange(len(plist_first)):
    plist[i]=plist_first[i] 
    plistd[i]=plistd_first[i]








def optimum_0(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):
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

 c_ut=((c_u*(MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))


 c_dt=((c_d*(MPO_list1[0])))
 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))


 A=((((E1*E8)*(a_ut*a_dt))*((E2*E3)*(b_ut*b_dt)))*(((E4*E5)*d)))*((E7*E6)*(c_ut*c_dt))


 A.permute([-62,-58,-57,-17,62,58,57,17],4)

 A1=copy.copy(A)
 A1.transpose()
 A=A+A1
 

 c_dt=((c_d*(MPO_list1[0])))
 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 Ap=(((((E1*E8)*(a_u1*a_dt))) * ((E2*E3)*(b_u1*b_dt)))*(((E4*E5)*d)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-58,-57,-17],4)


 c_ut=((c_u*(MPO_list0[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))


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
 A_d=copy.copy(A)
 A_d.setLabel([-62,-58,-57,-17])

 plistd[0]=copy.copy(A_d)



def optimum_1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):
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
 a_ut=(a_u*(plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 a_dt=(a_d*(plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))


 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-17,-64,-58,-57,17,64,58,57],4)

 A1=copy.copy(A)
 A1.transpose()
 A=A+A1


 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 a_dt=(a_d*(plistd[2]*MPO_list1[1]))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-17,-64,-58,-57],4)

 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 a_ut=(a_u*(plist[2]*MPO_list0[1]))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

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
 A_d=copy.copy(A)
 A_d.setLabel([-17,-64,-58,-57])

 plistd[1]=copy.copy(A_d)


def optimum_2(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):
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

 a_dt=(a_d*(plistd[1]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*((E7*E6)*(c_ut*c_dt)))*((E1*E8)*(a_ut*a_dt))

 A.permute([-66,-59,-60,-18,66,59,60,18],4)
 A1=copy.copy(A)
 A1.transpose()
 A=A+A1


 a_dt=(a_d*(plistd[1]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(plistd[3]*MPO_list1[2])))

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 
 Ap.permute([-66,-59,-60,-18],4)

 a_ut=(a_u*(plist[1]*MPO_list0[1]))
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 b_ut=((b_u*(plist[3]*MPO_list0[2])))

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
 A_d=copy.copy(A)
 A_d.setLabel([-66,-59,-60,-18])

 plistd[2]=copy.copy(A_d)


def optimum_3(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, U, MPO_list,c_u,a_u,b_u,plist, plistd):
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

 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(MPO_list1[2])))

 #print  a_ut.printDiagram(), a_dt.printDiagram(), a_ut[20],a_dt[20], a_ut.elemCmp(a_dt)
 #print  c_ut.printDiagram(), c_dt.printDiagram(), c_ut[30],c_dt[30], c_ut.elemCmp(c_dt)

 A=((((E1*E8)*(a_ut*a_dt)) *((E7*E6)*(c_ut*c_dt))) * ((((E4*E5)*d)))) * ((E2*E3)*(b_ut*b_dt))


 A.permute([-18,-68,-59,-60,18,68,59,60],4)
 A1=copy.copy(A)
 A1.transpose()
 A=A+A1

 a_dt=(a_d*(plistd[1]*plistd[2]*MPO_list1[1]))
 c_dt=((c_d*(plistd[0]*MPO_list1[0])))
 b_dt=((b_d*(MPO_list1[2])))

 Ap=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-68,-59,-60],4)


 a_ut=(a_u*(plist[1]*plist[2]*MPO_list0[1]))
 c_ut=((c_u*(plist[0]*MPO_list0[0])))
 b_ut=((b_u*(MPO_list0[2])))

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
 A_d=copy.copy(A)
 A_d.setLabel([-18,-68,-59,-60])

 plistd[3]=copy.copy(A_d)

 
def  recover( c_u, a_u, b_u, plist, plistd, MPO_list):

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




def Store_plist(plist,plistd):
 plist[0].save("Store/plist[0]")
 plist[1].save("Store/plist[1]")
 plist[2].save("Store/plist[2]")
 plist[3].save("Store/plist[3]")
 plistd[0].save("Store/plistd[0]")
 plistd[1].save("Store/plistd[1]")
 plistd[2].save("Store/plistd[2]")
 plistd[3].save("Store/plistd[3]")

def Reload_plist(plist,plistd):

 plist[0]=uni10.UniTensor("Store/plist[0]")
 plist[1]=uni10.UniTensor("Store/plist[1]")
 plist[2]=uni10.UniTensor("Store/plist[2]")
 plist[3]=uni10.UniTensor("Store/plist[3]")

 plistd[0]=uni10.UniTensor("Store/plistd[0]")
 plistd[1]=uni10.UniTensor("Store/plistd[1]")
 plistd[2]=uni10.UniTensor("Store/plistd[2]")
 plistd[3]=uni10.UniTensor("Store/plistd[3]")

 return plist,plistd
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



