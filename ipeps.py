import pyUni10 as uni10
import copy
import basic
import basic_H
import itebd
import Fullupdate
###################### Initialize parameters ###########################
#Model="Heisenberg"
#Model="Heisenberg_Z2"
Model="Heisenberg_U1"
#Model="Heisenberg_U1Z2"

# D=[4]
# chi=[20]
# #d_phys=[8]
# #d_phys=[2]
# d_phys=[3]

#D=[2,2]
#chi=[10,10]
#d_phys=[1,1]

D=[1,2,1]
chi=[1,1,2,10,2,1,1]
#d_phys=[1,1]
d_phys=[1,1,1]

#D=[2,2,2,2,2,2]
#chi=[10,10,10,10,10,10]
#d_phys=[1,1]

method="SVD"            #SVD, Grad, SVD_mpo, SVD_QR, SVD_mix,SVD_mg
Inv_method='SVD'         #SVD, CG
Grad_method="CG"        # CG,ST
Gauge='Fixed'
Corner_method='CTMFull'   #CTMRG, CTMFull
check_step='off'            #on, off
fixbond_itebd='off'            #on, off
itebd_method='long'            #long, short
start_itebd=1.0
division_itebd=5
iteration_per_step=1

N_grad=150
N_env=[20,1.00e-8]
N_svd=[10,1.00e-9,20]
N_iteritebd=20
N_iterFull=100
Acc_E=1.00e-7
distance_val=2
Arnoldi_bond=5
J1_x=1.0
J1_y=1.0
J2=0.50
h_z=1.0
h_x=1.0
h_y=1.0
K_1=0.0
K_2=0.0

###################################################################

q_D, q_chi, q_phys=basic.full_make_bond(Model, D, chi, d_phys)
Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.produce_GammaLanda(q_D, q_phys)
#Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.produce_GammaLanda_manual(q_D, q_phys)

Gamma=[Gamma_a,Gamma_b,Gamma_c,Gamma_d]
Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
basic.Initialize_function(Gamma,Landa)
a_u,b_u,c_u,d_u,a,b,c,d=basic.produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d)

Env=basic.produce_env_init(q_chi,q_D)
Env1=basic.Rand_env_total(Env)
Env2=basic.Rand_env_total(Env)
Env3=basic.Rand_env_total(Env)

###
a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
##a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full_previous(a_u, b_u, c_u, d_u)
##print a_u
##a_u,b_u,c_u,d_u,a,b,c,d=basic.slighty_random(a_u,b_u,c_u,d_u,a,b,c,d)
Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8=basic.produce_gamma_abcd(a_u,b_u,c_u,d_u,Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8)
####


fileCorrH = open("Data/CorrelationH.txt", "w")
fileCorrV = open("Data/CorrelationV.txt", "w")
fileCorrLength = open("Data/CorrelationLength.txt", "w")
fileSpectrum = open("Data/Spectrum.txt", "w")
fileSpectrum2 = open("Data/Spectrum2.txt", "w")
fileSpectrum4 = open("Data/Spectrum4.txt", "w")

fileSpectrum6 = open("Data/Spectrum6.txt", "w")
fileSpectrum8 = open("Data/Spectrum8.txt", "w")
fileSpectrum10 = open("Data/Spectrum10.txt", "w")

fileSpectrum0 = open("Data/Spectrum0.txt", "w")
fileSpectrum1 = open("Data/Spectrum1.txt", "w")
fileSpectrum3 = open("Data/Spectrum3.txt", "w")

fileCorrHH = open("Data/CorrelationHH.txt", "w")
fileCorrVV = open("Data/CorrelationVV.txt", "w")
fileCorrLength = open("Data/CorrelationLength.txt", "w")

fileCorrHH1 = open("Data/CorrelationHH1.txt", "w")
fileCorrVV1 = open("Data/CorrelationVV1.txt", "w")
fileCorrLength1 = open("Data/CorrelationLength1.txt", "w")


print J1_x, J1_y, J2,  h_z,h_x, h_y, K_1, K_2 
h_coupling=[J1_x, J1_y, J2,  h_z,h_x, h_y, K_1, K_2]
#########################################################################################
Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=basic.Reload_itebd()
Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8=itebd.itebd_eff(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8,chi,q_phys,D,N_iteritebd,h_coupling,Model,q_D,fixbond_itebd,start_itebd, division_itebd,itebd_method)
basic.Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8)

print Landa_1, '\n', Landa_2, '\n',Landa_3 ,'\n',Landa_4 ,'\n',Landa_5, '\n',Landa_6, '\n',Landa_7, '\n', Landa_8 

Landa=[Landa_1,Landa_2,Landa_3,Landa_4,Landa_5,Landa_6,Landa_7,Landa_8]
a_u,b_u,c_u,d_u,a,b,c,d=basic.produce_abcd_gamma(Landa, Gamma_a,Gamma_b,Gamma_c,Gamma_d)

a_u,b_u,c_u,d_u,a,b,c,d=basic.increase_norm(a_u,b_u,c_u,d_u,a,b,c,d, 6.00)
a_u,b_u,c_u,d_u,a,b,c,d=basic.make_equall_dis(a_u,b_u,c_u,d_u,a,b,c,d)

#a_u,b_u,c_u,d_u,a,b,c,d=basic.Reload_Full()
#print a_u.printDiagram(), b_u.printDiagram(), c_u.printDiagram(), d_u.printDiagram()
#basic.Reload_EnvEnv(Env,Env1,Env2,Env3)
E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h_coupling,q_phys,chi,Corner_method,Model,N_env)
print 'E_toal=', E_value
basic.Store_EnvEnv(Env,Env1,Env2,Env3)
M_value=basic.M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h_coupling,q_phys,chi,Corner_method,Model)
print 'M_s=', M_value, '\n'
basic.Translational_sym(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h_coupling,q_phys,chi,Corner_method,Model)
###########################################################################################

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


a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3=Fullupdate.Full_Update(a_u,b_u,c_u,d_u,a,b,c,d,chi,q_phys,D,h_coupling,Env,Env1,Env2,Env3,Gauge,Corner_method,N_iterFull,Acc_E,Model,N_grad, Grad_method,Inv_method,N_svd,N_env,method,check_step,iteration_per_step)

#E_value=basic.E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h_coupling,q_phys,chi,Corner_method,Model,N_env)
#print "E_final", E_value
#basic.Store_Full(a_u,b_u,c_u,d_u,a,b,c,d)
#basic.Store_EnvEnv(Env,Env1,Env2,Env3)
#basic.Reload_EnvEnv(Env,Env1,Env2,Env3)

#############################################
#basic.CorrelationH(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h_coupling,q_phys,chi,Corner_method,Model,distance_val,fileCorrHH,fileCorrLength)
#basic.CorrelationV(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h_coupling,q_phys,chi,Corner_method,Model,distance_val,fileCorrVV,fileCorrLength)
#basic.CorrelationH(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h_coupling,q_phys,chi,Corner_method,Model,distance_val,fileCorrHH1,fileCorrLength)
#basic.CorrelationV(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h_coupling,q_phys,chi,Corner_method,Model,distance_val,fileCorrVV1,fileCorrLength)
#############################################
#basic_H.CorrelationH(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h_coupling,q_phys,chi,Corner_method,Model,distance_val,fileCorrH,fileCorrLength)
#basic_H.CorrelationV(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h_coupling,q_phys,chi,Corner_method,Model,distance_val,fileCorrV,fileCorrLength)
#################################################
#basic_H.Spectrum(a_u,b_u,c_u,d_u,a,b,c,d,Env,fileSpectrum)
#basic_H.Spectrum_Arnoldi_right0(a_u,b_u,c_u,d_u,a,b,c,d,Env,fileSpectrum,fileSpectrum2,fileSpectrum4,Arnoldi_bond)
#basic_H.Spectrum_Arnoldi_right1(a_u,b_u,c_u,d_u,a,b,c,d,Env,fileSpectrum0,fileSpectrum1,fileSpectrum3,Arnoldi_bond)
#basic_H.Spectrum_Arnoldi_full(a,b,c,d,fileSpectrum6,fileSpectrum8,fileSpectrum10,Arnoldi_bond)




