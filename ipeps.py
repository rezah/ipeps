import pyUni10 as uni10
#import matplotlib
#matplotlib.use('pdf')
#import matplotlib.pyplot as plt
#import pylab
import copy
import basic
import itebd
import Fullupdate
###################### Initialize parameters ###########################
#Model="Heisenberg"
#Model="Heisenberg_Z2"
Model="Heisenberg_U1"
#Model="Heisenberg_U1Z2"

#D=[2]
#chi=[20]
##d_phys=[8]
#d_phys=[2]

#D=[2,2]
#chi=[10,10]
#d_phys=[1,1]

D=[1,2,1]
chi=[5,6,10,6,5]
d_phys=[1,1]

#D=[2,2,2,2,2,2]
#chi=[10,10,10,10,10,10]
#d_phys=[1,1]


Inv_method='CG'       #SVD, CG
Opt_method="CG"        # CG,ST
Gauge='Fixed'
Corner_method='CTMFull'   #CTM, CTMRG, CTMFull

Acc_E=1.00e-7
N_iteritebd=100
N_iterFull=50
N_grad=20
N_svd=[30,1.00e-8]
N_env=[30,1.00e-7]

J1_list=[h*0.0100 for h in range(270,400)]
J1_list=[1.0]
J2_list=[0.00]
h_list=[1.0]


###################################################################


q_D, q_chi, q_phys=basic.full_make_bond()


Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.produce_GammaLanda(q_D, q_phys)


Gamma=[Gamma_a,Gamma_b,Gamma_c,Gamma_d]
Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
basic.Initialize_function(Gamma,Landa)
a_u,b_u,c_u,d_u,a,b,c,d=basic.produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d)


Env=basic.produce_env_init(q_chi,q_D)
Env1=b asic.Rand_env_total(Env)
Env2=basic.Rand_env_total(Env)
Env3=basic.Rand_env_total(Env)


for h , J1, J2 in zip( h_list, J1_list, J2_list):
 print h, J1, J2 
 h=[h, J1, J2]
#########################################################################################
 
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.Reload_itebd()
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,q_phys,D,N_iteritebd,h,Model,q_D)
 #basic.Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8)

# print Landa_1, '\n', Landa_2, '\n',Landa_3 ,'\n',Landa_4 ,'\n',Landa_5, '\n',Landa_6, '\n',Landa_7, '\n', Landa_8 

# Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
# a_u,b_u,c_u,d_u,a,b,c,d=basic.produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d)

# a_u,b_u,c_u,d_u,a,b,c,d=basic.increase_norm(a_u,b_u,c_u,d_u,a,b,c,d, 0.0025)

 #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model,N_env)
 #print 'E_toal=', E_value
 #break
#########################################################################################

############################################################################
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.total_random(a_u,b_u,c_u,d_u,a,b,c,d)
 a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full_previous(a_u, b_u, c_u, d_u)
 basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.increase_norm(a_u,b_u,c_u,d_u,a,b,c,d, 0.5)
 #print a_u
 a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,q_phys,D,h,Env,Env1,Env2,Env3,Gauge,Corner_method,N_iterFull,Acc_E,Model,N_grad, Opt_method,Inv_method,N_svd,N_env)

 #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model,N_env)
# print "E_final", E_value
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)



