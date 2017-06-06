import pyUni10 as uni10
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import copy
import time
import basicitebd
import basic


def itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,
Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,d_phys,D,N_iterF, h,Model,q_D,fixbond_itebd,start_itebd, division_itebd,itebd_method):
  
 E_0=1.0
 E_1=2.0
   
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
   H2=basic.threebody_U1_help(h,d_phys)   
 if Model is "Heisenberg_U1Z2":
   H0=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=basic.Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=basic.Heisenberg1_U1(h[2],d_phys)
  
 U = uni10.UniTensor(H0.bond(), "U");
 U0 = uni10.UniTensor(H00.bond(), "U");

 for i in xrange(1,600):
 
   delta=start_itebd/pow(division_itebd,i) 

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





    if itebd_method is 'short':
 ##############################################################################
     Gamma=[Gamma_a,Gamma_b]
     Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
   #rlink

     basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
   #ulink
     Gamma=[Gamma_a, Gamma_c]
     Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
     basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D,fixbond_itebd)
     #print Landa_2.printDiagram()


   #rlink
     Gamma=[Gamma_b, Gamma_a]
     Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_5,Landa_6,Landa_2,Landa_4]
     basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
   #ulink
     Gamma=[ Gamma_b,Gamma_d]
     Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
     basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D,fixbond_itebd)


     Gamma=[Gamma_c,Gamma_d]
     Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
   #rlink
     basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
   #ulink
     Gamma=[ Gamma_c,Gamma_a]
     Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_7,Landa_8]
     basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D,fixbond_itebd)


   #rlink
     Gamma=[Gamma_d, Gamma_c]
     Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_3,Landa_1,Landa_4,Landa_2]
     basicitebd.update_rlink_eff(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
  #ulink
     Gamma=[ Gamma_d,Gamma_b]
     Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_2,Landa_4]
     basicitebd.update_ulink_eff(Gamma,Landa,U0,D,d_phys,q_D,fixbond_itebd)
 #############################################################################################
    
    
    A1=Landa_1.trace().real
    A2=Landa_2.trace().real
    A3=Landa_3.trace().real
    A4=Landa_4.trace().real
    A5=Landa_5.trace().real
    A6=Landa_6.trace().real
    A7=Landa_7.trace().real
    A8=Landa_8.trace().real
    #print A1, Landa_1
    E_0=E_1
    E_1=A1+A2+A3+A4+A5+A6+A7+A8
    print E_0, E_1, abs((E_0-E_1) / E_0), q
    if (( abs((E_0-E_1) / E_0) ) < 1.00e-10) and (q > 2): 
     print 'break', E_0, E_1, abs((E_0-E_1) / E_0), q
     E_0=10.0
     E_1=20.0
     break;


    blk_qnums = H2.blockQnum()
    U = uni10.UniTensor(H2.bond(), "U");
    for qnum in blk_qnums:
        U.putBlock(qnum, uni10.takeExp((-1.00/2.00)*delta, H2.getBlock(qnum)))
    #print H2

    #print Landa_1, '\n', Landa_2, '\n',Landa_3 ,'\n',Landa_4 ,'\n',Landa_5, '\n',Landa_6, '\n',Landa_7, 
    if itebd_method is 'long':

 #################################  abcd  ##############################################################
     #print q
     Gamma=[Gamma_a,Gamma_b,Gamma_c]
     Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
   #rlink
     basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_c,Gamma_d,Gamma_b]
     Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
 #  #rlink
     basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_a,Gamma_b,Gamma_d]
     Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
   #rlink
     basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_a,Gamma_c,Gamma_d]
     Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
   #rlink
     basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

 ######################################################################################################
 #    
 #################################  badc   #####################################################    
     Gamma=[Gamma_b,Gamma_a,Gamma_d]
     Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
   #rlink
     basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
     
     
     Gamma=[Gamma_a,Gamma_b,Gamma_d]
     Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
   #rlink
     basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_b,Gamma_a,Gamma_c]
     Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
    #rlink
     basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)


     Gamma=[Gamma_c,Gamma_a,Gamma_b]
     Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
   #rlink
     basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

 ################################################################################################    
 #    
 #    
 #    
 ##################################cdab########################################################    
     Gamma=[Gamma_c,Gamma_d,Gamma_a]
     Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
   #rlink
     basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_b,Gamma_a,Gamma_c]
     Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
   #rlink
     basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_c,Gamma_d,Gamma_b]
     Landa=[Landa_6,Landa_4,Landa_5,Landa_2,Landa_3,Landa_1,Landa_8,Landa_7]
 #  #rlink
     basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)

     Gamma=[Gamma_d,Gamma_b,Gamma_a]
     Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
   #rlink
     basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
 ##########################################################################################

 ####################################     dcba    #################################################    
     Gamma=[Gamma_d,Gamma_c,Gamma_b]
     Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
   #rlink
     basicitebd.update_rlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
     
     Gamma=[Gamma_d,Gamma_c,Gamma_a]
     Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
   #rlink
     basicitebd.update_rdlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
     
     Gamma=[Gamma_d,Gamma_c,Gamma_a]
     Landa=[Landa_5,Landa_8,Landa_6,Landa_7,Landa_1,Landa_3,Landa_4,Landa_2]
   #rlink
     basicitebd.update_ulink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)
     
     
     Gamma=[Gamma_b,Gamma_d,Gamma_c]
     Landa=[Landa_3,Landa_7,Landa_1,Landa_8,Landa_6,Landa_5,Landa_2,Landa_4]
   #rlink
     basicitebd.update_udlink_eff_long(Gamma,Landa,U,D,d_phys,q_D,fixbond_itebd)


 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8  



