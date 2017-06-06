import pyUni10 as uni10
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
#d_phys=[8]
#d_phys=[2]
#d_phys=[3]

#D=[2,2]
#chi=[10,10]
#d_phys=[1,1]

D=[1,2,1]
chi=[2,2,2,5,5,10,5,5,2,2,2]
d_phys=[1,1]
#d_phys=[1,1,1]

#D=[2,2,2,2,2,2]
#chi=[10,10,10,10,10,10]
#d_phys=[1,1]

method="SVD_mpo"            #SVD, Grad, SVD_mpo
Inv_method='SVD'         #SVD, CG
Grad_method="CG"        # CG,ST
Gauge='Non-Fixed'
Corner_method='CTMFull'   #CTMRG, CTMFull
check_step='off'            #on, off
fixbond_itebd='off'            #on, off
itebd_method='long'            #long, short
start_itebd=0.50
division_itebd=5

N_grad=120
N_env=[30,1.00e-8]
N_svd=[30,1.00e-8]
N_iteritebd=50
N_iterFull=50
Acc_E=1.00e-7
distance_val=15

J1_list=[h*0.0100 for h in range(270,400)]
J1_list=[1.0]
J2_list=[0.5000]
h_list=[1.00]

###################################################################

q_D, q_chi, q_phys=basic.full_make_bond(Model, D, chi, d_phys)


Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.produce_GammaLanda(q_D, q_phys)

#Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.produce_GammaLanda_manual(q_D, q_phys)



Gamma=[Gamma_a,Gamma_b,Gamma_c,Gamma_d]
Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
basic.Initialize_function(Gamma,Landa)
a_u,b_u,c_u,d_u,a,b,c,d=basic.produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d)


###
#a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
##a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full_previous(a_u, b_u, c_u, d_u)
##print a_u
##a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
#Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
#Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8=basic.produce_gamma_abcd(a_u,b_u,c_u,d_u,Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8)
##

Env=basic.produce_env_init(q_chi,q_D)
Env1=basic.Rand_env_total(Env)
Env2=basic.Rand_env_total(Env)
Env3=basic.Rand_env_total(Env)

fileCorrH = open("Data/CorrelationH.txt", "w")
fileCorrV = open("Data/CorrelationV.txt", "w")
fileCorrLength = open("Data/CorrelationLength.txt", "w")

for h , J1, J2 in zip( h_list, J1_list, J2_list):
 print h, J1, J2 
 h=[h, J1, J2]
#########################################################################################
 
 Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.Reload_itebd()
 #Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,q_phys,D,N_iteritebd,h,Model,q_D,fixbond_itebd,start_itebd, division_itebd,itebd_method)
 basic.Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8)

 print Landa_1, '\n', Landa_2, '\n',Landa_3 ,'\n',Landa_4 ,'\n',Landa_5, '\n',Landa_6, '\n',Landa_7, '\n', Landa_8 

 Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
 a_u,b_u,c_u,d_u,a,b,c,d=basic.produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d)

 a_u,b_u,c_u,d_u,a,b,c,d=basic.increase_norm(a_u,b_u,c_u,d_u,a,b,c,d, 4.00)
 a_u,b_u,c_u,d_u,a,b,c,d=basic.make_equall_dis(a_u,b_u,c_u,d_u,a,b,c,d)

# print a_u, b_u, c_u, d_u


# basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model,N_env)
 print 'E_toal=', E_value
 basic.Store_EnvEnv(Env,Env1,Env2,Env3)
## #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 M_value=basic.M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model)
 print 'M_s=', M_value
# basic.CorrelationH(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,q_phys,chi,Corner_method,Model,distance_val,fileCorrH,fileCorrLength)
# basic.CorrelationV(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,q_phys,chi,Corner_method,Model,distance_val,fileCorrV,fileCorrLength)
 #break
#########################################################################################

############################################################################
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.total_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.total_random1(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full_previous(a_u, b_u, c_u, d_u)
 #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
 #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 #a_u,b_u,c_u,d_u,a,b,c,d=basic.increase_norm(a_u,b_u,c_u,d_u,a,b,c,d, 0.5)
 #print a_u
 a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,q_phys,D,h,Env,Env1,Env2,Env3,Gauge,Corner_method,N_iterFull,Acc_E,Model,N_grad, Grad_method,Inv_method,N_svd,N_env,method,check_step)

 #E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,q_phys,chi,Corner_method,Model,N_env)
 #print "E_final", E_value
 #basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)

 #basic.Store_EnvEnv(Env,Env1,Env2,Env3)
 #basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
 #basic.CorrelationH(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,q_phys,chi,Corner_method,Model,distance_val,fileCorrH,fileCorrLength)
# basic.CorrelationV(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,q_phys,chi,Corner_method,Model,distance_val,fileCorrV,fileCorrLength)


