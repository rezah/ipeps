import pyUni10 as uni10
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import copy
import time
import basicitebd
import basic


def itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,
Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iterF, h,Model,q_D):
  
   
 if Model is "Heisenberg":
   H0=basic.Heisenberg0(h[0],h[1])
   H00=basic.Heisenberg00(h[0],h[1])
   H1=basic.Heisenberg1(h[2])
   H2=basic.threebody(h,d_phys)   
 if Model is "Heisenberg_Z2":
   H0=basic.Heisenberg0_Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg00_Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)
   H2=basic.threebody_U1(h,d_phys)   
 if Model is "Heisenberg_U1Z2":
   H0=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)
  
 U = uni10.UniTensor(H0.bond(), "U");
 U0 = uni10.UniTensor(H00.bond(), "U");

 for i in xrange(1,600):

   delta=1.00/pow(2,i) 

   if delta>1.0e-1:
    N_iter=N_iterF
    #delta=1.0e-1
   if delta<1.0e-1 and delta>1.0e-3:
    N_iter=N_iterF
   if delta<1.0e-3  and delta>1.0e-5:
    N_iter=N_iterF
   if delta<1.0e-5:
    break


   print 'delta =', delta
   print "N_iterF=", N_iterF
   
   for q in xrange(N_iterF):


##########################################################################################
#    H0=basic.Heisenberg_U1(h,d_phys)
#    H0=basic.Heisenberg(h)
    U = uni10.UniTensor(H0.bond(), "U");
    blk_qnums = H0.blockQnum()
    for qnum in blk_qnums:
        U.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))


    U0 = uni10.UniTensor(H00.bond(), "U");
    blk_qnums = H00.blockQnum()
    for qnum in blk_qnums:
     U0.putBlock(qnum, uni10.takeExp(-delta, H00.getBlock(qnum)))






#################################################################################
    Gamma=[Gamma_a,Gamma_b]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
  #rlink

    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)
  #ulink
    Gamma=[Gamma_a, Gamma_c]
    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
    basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D)
    #print Landa_2.printDiagram(),Gamma_a.printDiagram(),Gamma_c.printDiagram()


  #rlink
    Gamma=[Gamma_b, Gamma_a]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_5,Landa_6,Landa_2,Landa_4]
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)
  #ulink
    Gamma=[ Gamma_b,Gamma_d]
    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
    basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D)


    Gamma=[Gamma_c,Gamma_d]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
  #rlink
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)
  #ulink
    Gamma=[ Gamma_c,Gamma_a]
    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_7,Landa_8]
    basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D)


  #rlink
    Gamma=[Gamma_d, Gamma_c]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_3,Landa_1,Landa_4,Landa_2]
    basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D)
 #ulink
    Gamma=[ Gamma_d,Gamma_b]
    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_2,Landa_4]
    basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D)
################################################################################################



#    blk_qnums = H2.blockQnum()
#    U = uni10.UniTensor(H2.bond(), "U");
#    for qnum in blk_qnums:
#        U.putBlock(qnum, uni10.takeExp((-1.00/2.00)*delta, H2.getBlock(qnum)))




#################################  abcd  ##############################################################
#    Gamma=[Gamma_a,Gamma_b,Gamma_c]
#    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
#  #rlink
#    basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_c,Gamma_d,Gamma_b]
#    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
#  #rlink
#    basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_a,Gamma_b,Gamma_d]
#    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
#  #rlink
#    basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_a,Gamma_c,Gamma_d]
#    Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
#  #rlink
#    basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

######################################################################################################
#    
#################################  badc   #####################################################    
#    Gamma=[Gamma_b,Gamma_a,Gamma_d]
#    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
#  #rlink
#    basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)
#    
#    
#    Gamma=[Gamma_a,Gamma_b,Gamma_d]
#    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
#  #rlink
#    basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_b,Gamma_a,Gamma_c]
#    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
#  #rlink
#    basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D)


#    Gamma=[Gamma_c,Gamma_a,Gamma_b]
#    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
#  #rlink
#    basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

################################################################################################    
#    
#    
#    
##################################cdab########################################################    
#    Gamma=[Gamma_c,Gamma_d,Gamma_a]
#    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
#  #rlink
#    basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_b,Gamma_a,Gamma_c]
#    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
#  #rlink
#    basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_c,Gamma_d,Gamma_b]
#    Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
#  #rlink
#    basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D)

#    Gamma=[Gamma_d,Gamma_b,Gamma_a]
#    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
#  #rlink
#    basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)
##########################################################################################
#    
####################################     dcba    #################################################    
#    Gamma=[Gamma_d,Gamma_c,Gamma_b]
#    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
#  #rlink
#    basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)
#    
#    Gamma=[Gamma_d,Gamma_c,Gamma_a]
#    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
#  #rlink
#    basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)
#    
#    Gamma=[Gamma_d,Gamma_c,Gamma_a]
#    Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
#  #rlink
#    basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D)
#    
#    
#    Gamma=[Gamma_b,Gamma_d,Gamma_c]
#    Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
#  #rlink
#    basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D)


 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8  



