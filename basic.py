import pyUni10 as uni10
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import copy
import time
import Move
import MoveCorboz
import MoveFull
import basicA
import basicB
import basicC


def produce_GammaLanda(q_D, q_phys):
 bdi = uni10.Bond(uni10.BD_IN, q_D)
 bdo = uni10.Bond(uni10.BD_OUT, q_D)
 bdi_pys = uni10.Bond(uni10.BD_IN, q_phys)

 Gamma_a=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_a")
 Gamma_b=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_b")
 Gamma_c=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_c")
 Gamma_d=uni10.UniTensor([bdi_pys,bdi,bdi,bdo,bdo], "Gamma_d")

 Landa_1=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa_2=uni10.UniTensor([bdi,bdo],"Landa_2")
 Landa_3=uni10.UniTensor([bdi,bdo],"Landa_3")
 Landa_4=uni10.UniTensor([bdi,bdo],"Landa_4")
 Landa_5=uni10.UniTensor([bdi,bdo],"Landa_5")
 Landa_6=uni10.UniTensor([bdi,bdo],"Landa_6")
 Landa_7=uni10.UniTensor([bdi,bdo],"Landa_7")
 Landa_8=uni10.UniTensor([bdi,bdo],"Landa_8")
 
 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8

def produce_env_init(q_chi,q_D):
 c1, c2,c3,c4=makec1(q_chi,q_D)
 Ta1, Tb1=makeTab(q_chi,q_D)
 Ta2, Tb2=makeTab(q_chi,q_D)
 Ta3, Tb3=makeTab1(q_chi,q_D)
 Ta4, Tb4=makeTab1(q_chi,q_D)
 Env=[c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4]
 
 return Env




def Initialize_function(Gamma,Landa):
 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
 #Gamma[0].randomize()
 for i in xrange(len(Gamma)):
  Gamma[i].randomize()
  Gamma[i]=Gamma[i]*(1.00/MaxAbs(Gamma[i]))
  
 for i in xrange(len(Landa)):
  #Landa[i].identity()
  blk_qnums = Landa[i].blockQnum()
  for qnum in blk_qnums:
    M=Landa[i].getBlock(qnum)
    if qnum == q0_even:
     #M[0]=1.100
     M.randomize()
     #M.identity()
     #print "M0", M
    else: 
     #M[0]=0.01
     M.randomize()
     #M.identity()

    Landa[i].putBlock(qnum,M)
  Landa[i]=Landa[i]*(1.00/MaxAbs(Landa[i]))

 
def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0]);

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);


def Heisenberg0(h, J1):
    spin = 0.5
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    ham =J1*0.25*(h*uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy))
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
    H.putBlock(ham)

#    szt0=uni10.otimes(iden,sz)
#    sxt0=uni10.otimes(iden,sx)
#    syt0=uni10.otimes(iden,sy)

#    ident=uni10.otimes(iden,iden)
#    szt=uni10.otimes(ident,sz)
#    sxt=uni10.otimes(ident,sx)
#    syt=uni10.otimes(ident,sy)
#    
#    sztt=uni10.otimes(sz,ident)
#    sxtt=uni10.otimes(sx,ident)
#    sytt=uni10.otimes(sy,ident)

#    szttt=uni10.otimes(szt0,iden)
#    sxttt=uni10.otimes(sxt0,iden)
#    syttt=uni10.otimes(syt0,iden)

#    ham =J1*(h*uni10.otimes(szt,sztt)+uni10.otimes(sxt,sxtt)+(-1.0)*uni10.otimes(syt,sytt))
#    ham =ham + J1*(h*uni10.otimes(szt,szttt)+uni10.otimes(sxt,sxttt)+(-1.0)*uni10.otimes(syt,syttt))

#    szt=uni10.otimes(sz,sz)
#    sxt=uni10.otimes(sx,sx)
#    syt=uni10.otimes(sy,sy)

#    sztt=uni10.otimes(iden,szt)
#    sxtt=uni10.otimes(iden,sxt)
#    sytt=uni10.otimes(iden,syt)

#    szttt=uni10.otimes(szt,iden)
#    sxttt=uni10.otimes(sxt,iden)
#    syttt=uni10.otimes(syt,iden)
#    identt=uni10.otimes(ident,iden)


#    ham =ham + 0.25*J1*(h*uni10.otimes(szttt,identt)+uni10.otimes(sxttt,identt)+(-1.0)*uni10.otimes(syttt,identt))

#    ham =ham + 0.25*J1*(h*uni10.otimes(sztt,identt)+uni10.otimes(sxtt,identt)+(-1.0)*uni10.otimes(sytt,identt))


#    ham =ham + 0.25*J1*(h*uni10.otimes(identt,szttt)+uni10.otimes(identt,sxttt)+(-1.0)*uni10.otimes(identt,syttt))

#    ham =ham + 0.25*J1*(h*uni10.otimes(identt,sztt)+uni10.otimes(identt,sxtt)+(-1.0)*uni10.otimes(identt,sytt))


#    dim = int(8)
#    bdi = uni10.Bond(uni10.BD_IN, dim);
#    bdo = uni10.Bond(uni10.BD_OUT, dim);
#    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
#    H.putBlock(ham)
#    #print H.bond()
    #print H
    return H

def Heisenberg00(h, J1):
    spin = 0.5
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    ham =J1*0.25*(h*uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy))
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
    H.putBlock(ham)


#    szt0=uni10.otimes(iden,sz)
#    sxt0=uni10.otimes(iden,sx)
#    syt0=uni10.otimes(iden,sy)


#    ident=uni10.otimes(iden,iden)
#    szt=uni10.otimes(ident,sz)
#    sxt=uni10.otimes(ident,sx)
#    syt=uni10.otimes(ident,sy)
#    
#    sztt=uni10.otimes(sz,ident)
#    sxtt=uni10.otimes(sx,ident)
#    sytt=uni10.otimes(sy,ident)

#    szttt=uni10.otimes(szt0,iden)
#    sxttt=uni10.otimes(sxt0,iden)
#    syttt=uni10.otimes(syt0,iden)


#    ham =J1*(h*uni10.otimes(sztt,szttt)+uni10.otimes(sxtt,sxttt)+(-1.0)*uni10.otimes(sytt,syttt))
#    ham =ham + J1*(h*uni10.otimes(sztt,szt)+uni10.otimes(sxtt,sxt)+(-1.0)*uni10.otimes(sytt,syt))


#    szt=uni10.otimes(sz,sz)
#    sxt=uni10.otimes(sx,sx)
#    syt=uni10.otimes(sy,sy)

#    sztt=uni10.otimes(iden,szt)
#    sxtt=uni10.otimes(iden,sxt)
#    sytt=uni10.otimes(iden,syt)

#    szttt=uni10.otimes(szt,iden)
#    sxttt=uni10.otimes(sxt,iden)
#    syttt=uni10.otimes(syt,iden)
#    identt=uni10.otimes(ident,iden)


#    ham =ham + 0.25*J1*(h*uni10.otimes(szttt,identt)+uni10.otimes(sxttt,identt)+(-1.0)*uni10.otimes(syttt,identt))

#    ham =ham + 0.25*J1*(h*uni10.otimes(sztt,identt)+uni10.otimes(sxtt,identt)+(-1.0)*uni10.otimes(sytt,identt))


#    ham =ham + 0.25*J1*(h*uni10.otimes(identt,szttt)+uni10.otimes(identt,sxttt)+(-1.0)*uni10.otimes(identt,syttt))

#    ham =ham + 0.25*J1*(h*uni10.otimes(identt,sztt)+uni10.otimes(identt,sxtt)+(-1.0)*uni10.otimes(identt,sytt))


#    dim = int(8)
#    bdi = uni10.Bond(uni10.BD_IN, dim);
#    bdo = uni10.Bond(uni10.BD_OUT, dim);
#    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
#    H.putBlock(ham)
    #print H
    return H




def Heisenberg1(J2):
    spin = 0.5
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    ham =J2*(uni10.otimes(sz,sz)+uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy))
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
    H.putBlock(ham)
    return H


def threebody(h,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi,bdi, bdo, bdo,bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    szt=uni10.otimes(sz,sz)
    sxt=uni10.otimes(sx,sx)
    syt=(-1.0)*uni10.otimes(sy,sy)
    ident=uni10.otimes(iden,iden)

    sztt=uni10.otimes(iden,sz)
    sxtt=uni10.otimes(iden,sx)
    sytt=(-1.0)*uni10.otimes(iden,sy)


    ham =(0.25*h[1])*(h[0]*uni10.otimes(szt,iden)+uni10.otimes(sxt,iden)+uni10.otimes(syt,iden))
    ham =ham + 0.25*h[1]*(h[0]*uni10.otimes(iden,szt)+uni10.otimes(iden,sxt)+uni10.otimes(iden,syt))
    ham =ham + 0.50*h[2]*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

    H.putBlock(ham)
#    H.randomize()
#    H.randomize()
#    H1=copy.copy(H)
#    H1.transpose()
#    H=H+H1
#    print H
#    H.save("Store/H")
#    H=uni10.UniTensor("Store/H")
#    ham=H.getBlock()
#    D=ham.eigh()
#    print "eig", D[0][0], D[0]
#    H.save("Store/H")
#    print H
    return H

def threebody_U1(h,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)

    H = uni10.UniTensor([bdi, bdi,bdi, bdo, bdo,bdo], "Heisenberg")
    sx = matSx()
    sy = matSy()
    sz = matSz()
    iden = uni10.Matrix(2,2, [1, 0, 0, 1])
    szt=uni10.otimes(sz,sz)
    sxt=uni10.otimes(sx,sx)
    syt=(-1.0)*uni10.otimes(sy,sy)
    ident=uni10.otimes(iden,iden)

    sztt=uni10.otimes(iden,sz)
    sxtt=uni10.otimes(iden,sx)
    sytt=(-1.0)*uni10.otimes(iden,sy)


    ham =(0.25*h[1])*(h[0]*uni10.otimes(szt,iden)+uni10.otimes(sxt,iden)+uni10.otimes(syt,iden))
    ham =ham + 0.25*h[1]*(h[0]*uni10.otimes(iden,szt)+uni10.otimes(iden,sxt)+uni10.otimes(iden,syt))
    ham =ham + 0.25*2.00*h[2]*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

    #ham = 0.25*2.00*h[2]*(uni10.otimes(sz,sztt)+uni10.otimes(sx,sxtt)+uni10.otimes(sy,sytt))

    H.setRawElem(ham)
    #print ham, H

    return H



def Heisenberg0_Z2(h,J1,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    blk_qnums=H.blockQnum()
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    M[1]=0.0
    M[2]=0.0
    M[3]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=(-J1)*h
    M[1]=2.0*J1
    M[2]=2.0*J1
    M[3]=(-J1)*h
    H.putBlock(blk_qnums[1],M)
##########################################################################

    return H

def Heisenberg00_Z2(h,J1,d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")

    blk_qnums=H.blockQnum()
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    M[1]=0.0
    M[2]=0.0
    M[3]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=(-J1)*h
    M[1]=2.0*J1
    M[2]=2.0*J1
    M[3]=(-J1)*h
    H.putBlock(blk_qnums[1],M)
##########################################################################

    return H





def Heisenberg0_U1(h,J1, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    #H.randomize()
    #print transverseIsing(h).getBlock()
    #H.setRawElem(transverseIsing(h).getBlock().getElem());
    #H.setRawElem(Heisenberg().getBlock());
    blk_qnums=H.blockQnum()
    blk_qnums[0]
    
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=J1*(-h)
    M[1]=J1*2.0
    M[2]=J1*2.0
    M[3]=J1*(-h)
    H.putBlock(blk_qnums[1],M)


    M=H.getBlock(blk_qnums[2])
    M[0]=J1*h
    H.putBlock(blk_qnums[2],M)
    H=0.25*H
    #print H
    return H



def Heisenberg0_U1Z2(h,J1, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    #H.randomize()
    #print transverseIsing(h).getBlock()
    #H.setRawElem(transverseIsing(h).getBlock().getElem());
    #H.setRawElem(Heisenberg().getBlock());
    blk_qnums=H.blockQnum()
    blk_qnums[0]
    
    M=H.getBlock(blk_qnums[0])
    M[0]=J1*h
    H.putBlock(blk_qnums[0],M)

    M=H.getBlock(blk_qnums[1])
    M[0]=J1*(-h)
    M[1]=J1*2.0
    M[2]=J1*2.0
    M[3]=J1*(-h)
    H.putBlock(blk_qnums[1],M)


    M=H.getBlock(blk_qnums[2])
    M[0]=J1*h
    H.putBlock(blk_qnums[2],M)
    H=H*0.250
    #print H
    return H



def Heisenberg1_U1(J2, d_phys):
    bdi = uni10.Bond(uni10.BD_IN, d_phys)
    bdo = uni10.Bond(uni10.BD_OUT, d_phys)
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    blk_qnums=H.blockQnum()
    blk_qnums[0]
    M=H.getBlock(blk_qnums[0])
    M[0]=J2
    H.putBlock(blk_qnums[0],M)
    M=H.getBlock(blk_qnums[1])
    M[0]=J2*1.00
    M[1]=J2*2.0
    M[2]=J2*2.0
    M[3]=J2*1.00
    H.putBlock(blk_qnums[1],M)
    M=H.getBlock(blk_qnums[2])
    M[0]=J2
    H.putBlock(blk_qnums[2],M)
    #print Heisenberg(h)
    return H






def makeTab(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_IN, D)
 #bdi1.combine(copy.copy(bdi1))
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem0.randomize()
 Tem0*=(1.00/MaxAbs(Tem0))
 Tem1=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem1.randomize()
 Tem1*=(1.00/MaxAbs(Tem1))
 return Tem0, Tem1

def makeTab1(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdi1 = uni10.Bond(uni10.BD_OUT, D)
 #bdi1.combine(copy.copy(bdi1))
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 Tem0=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem0.randomize()
 Tem0*=(1.00/Tem0.norm())
 Tem1=uni10.UniTensor([bdi, bdi1,bdi1, bdo])
 Tem1.randomize()
 Tem1*=(1.00/Tem1.norm())
 return Tem0, Tem1

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

def makec1(chi,D):
 bdi = uni10.Bond(uni10.BD_IN, chi)
 bdo = uni10.Bond(uni10.BD_OUT, chi)
 c1=uni10.UniTensor([bdi, bdo])
 c2=uni10.UniTensor([bdi, bdi])
 c3=uni10.UniTensor([bdi, bdo])
 c4=uni10.UniTensor([bdo, bdo])
 c1.randomize()
 c2.randomize()
 c3.randomize()
 c4.randomize()
 c1*=(1.00/MaxAbs(c1))
 c2*=(1.00/MaxAbs(c2))
 c3*=(1.00/MaxAbs(c3))
 c4*=(1.00/MaxAbs(c4))
 return c1,c2,c3,c4
 
def makeab(Landa,Gamma):
 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa)) ]
 Landa_sq=sqrt(Landa_cp)
 a_u=copy.copy(Gamma)
 Landa_sq[0].setLabel([-1,1])
 Landa_sq[1].setLabel([-2,2])
 Landa_sq[2].setLabel([3,-3])
 Landa_sq[3].setLabel([4,-4])


 a_u.setLabel([0,1,2,3,4])
 
 a_u=(((((a_u*Landa_sq[0])*Landa_sq[1])*Landa_sq[2])*Landa_sq[3]))
 a_u.permute([0,-1,-2,-3,-4],3)

 a_d=copy.copy(a_u)
 a_d.transpose()
 a_d.setLabel([-3,-4,0,-1,-2])

 a_u.setLabel([0,1,2,3,4])

 a=a_u*a_d
 a.permute([1,-1,2,-2,3,-3,4,-4],4)

 return a_u, a


 
def sqrt(Landa):
 Landa_cp=[ copy.copy(Landa[i]) for i in xrange(len(Landa))   ]
 for q in xrange(len(Landa_cp)): 
  blk_qnums=Landa_cp[q].blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp[q].getBlock(qnum).col())
   Landa_cpm=Landa_cp[q].getBlock(qnum)
   Landam=Landa[q].getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
   Landa_cp[q].putBlock(qnum,Landa_cpm)
 return Landa_cp

 
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


def make_ab(a_u):
 a_u.setLabel([0,1,2,3,4])
 a_uc=copy.copy(a_u)
 a_uc.transpose()
 a_uc.setLabel([-3,-4,0,-1,-2])
 result=a_uc*a_u
 result.permute([1,-1,2,-2,3,-3,4,-4], 4)
 return result




def corner_transfer_matrix_twosite(a,b,c,d,chi,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D):
 z1=copy.copy(a)
 z1.identity()
 z2=copy.copy(a)
 z2.randomize()
 z1=z2#+(1.0e-1)*z2
 
 Accuracy=1.00e-7
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 #print  '\n', '\n', 'CTM'
 
 while Loop_iter is 0: 

  c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4=Move.add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D)
  
  norm=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  norm1=Move.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
  E0=E1
  if (abs(norm) > 1.00e-10):
   E1=abs(norm1)/abs(norm)
   if (abs((E0-E1)/E0) < Accuracy):Loop_iter=1;
  else:
   E1=abs(norm1)
   if (abs((E0-E1)) < Accuracy) : print 'Warning: norm~0', E1; Loop_iter=1;
  count+=1
  if (count > 15 ): print 'break! CTM'; break;
  print E1, abs((E0-E1)/E1),norm, count
  #print E1, Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()
 
 
 #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4)
 #print 'CTM', norm
 return c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4



def corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,direct_env):
 z1=copy.copy(a)
 z1.identity()
 z2=copy.copy(a)
 z2.randomize()
 z1=z2#+(1.0e-1)*z2
 
 Accuracy=1.00e-7
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 #print  '\n', '\n', 'CTM' 
 while Loop_iter is 0: 
  t0=time.time()
  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveCorboz.add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D)
  
  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveCorboz.add_left1(c1,c2,c3,c4,Tb1,Ta2,Tb3,Ta4,Ta1,Tb2,Ta3,Tb4,b,a,d,c,chi,D) 

  MoveCorboz.permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveCorboz.add_left1(c1,c4,c3,c2,Ta4,Ta3,Ta2,Ta1,Tb4,Tb3,Tb2,Tb1,a,c,b,d,chi,D)
  

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveCorboz.add_left1(c1,c4,c3,c2,Tb4,Ta3,Tb2,Ta1,Ta4,Tb3,Ta2,Tb1,c,a,d,b,chi,D)

  MoveCorboz.permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)


  
  norm=MoveCorboz.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  norm1=MoveCorboz.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)


#  norm=1.0
#  if direct_env is 'h':
#   norm1=MoveCorboz.Env_energy_h(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
#  elif direct_env is 'v': 
#   norm1=MoveCorboz.Env_energy_v(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)


  
  E0=E1
  if (abs(norm) > 1.00e-10):
   E1=abs(norm1)/abs(norm)
   if (abs((E0-E1)/E0) < Accuracy):Loop_iter=1;
  else:
   E1=abs(norm1)
   if (abs((E0-E1)) < Accuracy) : print 'Warning: norm~0', E1; Loop_iter=1;
  count+=1
  if (count > 20 ): print 'break! CTMRG'; break;
  print E1, abs((E0-E1)/E1),norm, count,  time.time() - t0, "CTMRG"
  #print E1, Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()
  
 #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4)
 #print 'CTM', norm
 return c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4


def corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,direct_env):
 z1=copy.copy(a)
 z1.identity()
 z2=copy.copy(a)
 z2.randomize()
 z1=z2#+(1.0e-1)*z2
 #z1=200*z2
 
 Accuracy=1.00e-7
 E0=20.00
 E1=10.00
 Loop_iter=0
 count=0
 while Loop_iter is 0: 
  t0=time.time()
  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveFull.add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D)

  c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3=MoveFull.add_left1(c1,c2,c3,c4,Tb1,Ta2,Tb3,Ta4,Ta1,Tb2,Ta3,Tb4,b,a,d,c,chi,D) 

  MoveFull.permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveFull.add_left1(c1,c4,c3,c2,Ta4,Ta3,Ta2,Ta1,Tb4,Tb3,Tb2,Tb1,a,c,b,d,chi,D)
  

  c1, Ta1, Tb1, c2, c4, Ta3, Tb3, c3=MoveFull.add_left1(c1,c4,c3,c2,Tb4,Ta3,Tb2,Ta1,Ta4,Tb3,Ta2,Tb1,c,a,d,b,chi,D)

  MoveFull.permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

  
  norm=MoveFull.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  norm1=MoveFull.magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,z1,b,c,d)
  norm=norm[0]
  norm1=norm1[0]
   
#  norm=1.0
#  if direct_env is 'h':
#   norm1=MoveFull.Env_energy_h(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)
#  elif direct_env is 'v':
#   norm1=MoveFull.Env_energy_v(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys)



  E0=E1
  if (abs(norm) > 1.00e-10):
   E1=abs(norm1)/abs(norm)
   if (abs((E0-E1)/E0) < Accuracy):Loop_iter=1;
  else:
   E1=abs(norm1)
   if (abs((E0-E1)) < Accuracy) : print 'Warning: norm~0', E1; Loop_iter=1;
  count+=1
  if (count > 20 ): print 'break! CTMFull'; break;
  print E1, abs((E0-E1)/E1),norm, count, time.time() - t0,"CTMFull"
  #print E1, Truncation[0], abs((E0-E1)/E1)
  #print a.norm(), b.norm(), c.norm(), d.norm()

 #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 
 #print 'CTM', norm
 return c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4
 



def E_total_conv(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model):

############################################################################
# E_val1=Energy_cb(c_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val1", E_val1

# E_val2=Energy_ad(a_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val2", E_val2

# E_val3=Energy_cb(d_u,a_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val3", E_val3

# E_val4=Energy_ad(b_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val4", E_val4
# 
# E_val5=Energy_cb(a_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val5", E_val5

# E_val6=Energy_ad(c_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val6", E_val6

# E_val7=Energy_cb(b_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val7", E_val7

# E_val8=Energy_ad(d_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val8", E_val8
############################################################################
 E_1=Energy_cab(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
 print "E_1", E_1

 E_2=Energy_cab(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
 print "E_2", E_2

 E_3=Energy_cab(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
 print "E_3", E_3

 E_4=Energy_cab(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
 print "E_4", E_4
 
 E_5=Energy_abd(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
 print "E_5", E_5

 E_6=Energy_abd(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
 print "E_6", E_6

 E_7=Energy_abd(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
 print "E_7", E_7

 E_8=Energy_abd(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
 print "E_8", E_8

 #return E_5
# return (E_5+E_6+E_7+E_8)/4.00

# return (E_1+E_2+E_3+E_4)/4.00 
 return (E_1+E_2+E_3+E_4+E_5+E_6+E_7+E_8)/8.00
# return (E_ab+E_ba)/2.0
# return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc) / 4.00)#+((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc)/4.00)


def E_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model):

############################################################################
# E_val1=Energy_cb(c_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val1", E_val1

# E_val2=Energy_ad(a_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val2", E_val2

# E_val3=Energy_cb(d_u,a_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val3", E_val3

# E_val4=Energy_ad(b_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val4", E_val4
# 
# E_val5=Energy_cb(a_u,d_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val5", E_val5

# E_val6=Energy_ad(c_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val6", E_val6

# E_val7=Energy_cb(b_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val7", E_val7

# E_val8=Energy_ad(d_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
# #print "E_val8", E_val8
############################################################################

 E_1=Energy_cab(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
 print "E_1", E_1

 E_2=Energy_cab(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
 print "E_2", E_2

 E_3=Energy_cab(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
 print "E_3", E_3

 E_4=Energy_cab(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
 print "E_4", E_4
 
 E_5=Energy_abd(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model)
 print "E_5", E_5

 E_6=Energy_abd(b_u,a_u,d_u,c_u,b,a,d,c,Env1,D,h,d_phys,chi,Corner_method,Model)
 print "E_6", E_6

 E_7=Energy_abd(c_u,d_u,a_u,b_u,c,d,a,b,Env2,D,h,d_phys,chi,Corner_method,Model)
 print "E_7", E_7

 E_8=Energy_abd(d_u,c_u,b_u,a_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model)
 print "E_8", E_8

 #return E_5
# return (E_5+E_6+E_7+E_8)/4.00

# return (E_1+E_2+E_3+E_4)/4.00 
 return (E_1+E_2+E_3+E_4+E_5+E_6+E_7+E_8)/8.00
# return (E_ab+E_ba)/2.0
# return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc) / 4.00)#+((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_val1+E_val2+E_val3+E_val4) / 4.00) + ((E_val5+E_val6+E_val7+E_val8) / 4.00)
 #return ((E_ca+E_ac+E_db+E_bd) / 4.00) + ((E_ab+E_ba+E_cd+E_dc)/4.00)



def Energy_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)
 
# Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

# if Corner_method is 'CTM':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
# if Corner_method is 'CTMRG':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)
# if Corner_method is 'CTMFull':
#  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D)

# reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 
 if Model is "Heisenberg":
   H0=Heisenberg0(h[0],h[1])
   H00=Heisenberg00(h[0],h[1])
   H1=Heisenberg1(h[2])
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h[0],h[1],d_phys)
   H00=Heisenberg00_Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h[0],h[1],d_phys)
   H00=Heisenberg0_U1(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)




 E_ca=basicB.Energy_ca(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H00,c_u,a_u)

# E_ca=basicB.Energy_ca_positive(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H00,c_u,a_u)
# print "E_ca", E_ca
 return E_ca




def Energy_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model):

 if Model is "Heisenberg":
   H0=Heisenberg0(h[0],h[1])
   H1=Heisenberg1(h[2])
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 #c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Reload_Env()

 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D)
 if Corner_method is 'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
 if Corner_method is 'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
  #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 

 reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4, Tb4,D,d_phys)

 E_ab=basicB.Energy_ab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0, a_u, b_u)

# E_ab=basicB.Energy_ab_positive(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0, a_u, b_u)
# print "E_ab", E_ab
 return E_ab


def Energy_cab(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model):
 
 if Model is "Heisenberg":
   H0=Heisenberg0(h[0],h[1])
   H00=Heisenberg00(h[0],h[1])
   H1=Heisenberg1(h[2])
   H2=threebody(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h[0],h[1],d_phys)
   H00=Heisenberg00_Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h[0],h[1],d_phys)
   H00=Heisenberg0_U1(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)
   H2=threebody_U1(h,d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)
 
 
 #H0.randomize()
 
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 #c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Reload_Env()


 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D)
 if Corner_method is 'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
 if Corner_method is 'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
  #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 

  
  
 reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)
 E1, E2, E3, E4, E5, E6, E7, E8=basicC.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)


 E=basicC.energy_cab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, H2)

 return E

def Energy_abd(a_u,b_u,c_u,d_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model):
 
 if Model is "Heisenberg":
   H0=Heisenberg0(h[0],h[1])
   H00=Heisenberg00(h[0],h[1])
   H1=Heisenberg1(h[2])
   H2=threebody(h,d_phys)
 if Model is "Heisenberg_Z2":
   H0=Heisenberg0_Z2(h[0],h[1],d_phys)
   H00=Heisenberg00_Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_Z2(h[2],d_phys)
 if Model is "Heisenberg_U1":
   H0=Heisenberg0_U1(h[0],h[1],d_phys)
   H00=Heisenberg0_U1(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)
   H2=threebody_U1(h,d_phys)
 if Model is "Heisenberg_U1Z2":
   H0=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H00=Heisenberg0_U1Z2(h[0],h[1],d_phys)
   H1=Heisenberg1_U1(h[2],d_phys)
 
 
 #H0.randomize()
 
 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 #c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Reload_Env()


 Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4)

 if Corner_method is 'CTM':
  c1, c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4=corner_transfer_matrix_twosite(a,b,c,d,chi,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4,Tb4,D)
 if Corner_method is 'CTMRG':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMRG(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
 if Corner_method is 'CTMFull':
  c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4=corner_transfer_matrix_twosite_CTMFull(a_u,b_u,a,b,c,d,chi,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,H0,d_phys,'h')
  #Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4) 

  
  
 reconstruct_env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4,Env)
 E1, E2, E3, E4, E5, E6, E7, E8=basicC.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)


 E=basicC.energy_abd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, H2)

 return E

 

def M_total(a_u,b_u,c_u,d_u,a,b,c,d,Env,Env1,Env2,Env3,D,h,d_phys,chi,Corner_method,Model):

 iden = uni10.Matrix(2,2, [1, 0, 0, 1])
 sx = matSx()
 sy = matSy()
 sz = matSz()

 ident=uni10.otimes(iden,iden)
 identt=uni10.otimes(ident,iden)

 szttt=uni10.otimes(ident,sz)
 sxttt=uni10.otimes(ident,sx)
 syttt=uni10.otimes(ident,sy)
 
 szt=uni10.otimes(sz,ident)
 sxt=uni10.otimes(sx,ident)
 syt=uni10.otimes(sy,ident)

 szt0=uni10.otimes(iden,sz)
 sxt0=uni10.otimes(iden,sx)
 syt0=uni10.otimes(iden,sy)

 sztt=uni10.otimes(szt0,iden)
 sxtt=uni10.otimes(sxt0,iden)
 sytt=uni10.otimes(syt0,iden)


 H0 =uni10.otimes(szt,identt)+uni10.otimes(sxt,identt)+uni10.otimes(syt,identt)

 H0 =H0 + -1.*(uni10.otimes(sztt,identt)+uni10.otimes(sxtt,identt)+uni10.otimes(sytt,identt))

 H0 =H0 + uni10.otimes(szttt,identt)+uni10.otimes(sxttt,identt)+uni10.otimes(syttt,identt)

 H0 =H0 + -1.*(uni10.otimes(identt,szt)+uni10.otimes(identt,sxt)+uni10.otimes(identt,syt))

 H0 =H0 + uni10.otimes(identt,sztt)+uni10.otimes(identt,sxtt)+uni10.otimes(identt,sytt)

 H0 =H0 + -1.*(uni10.otimes(identt,szttt)+uni10.otimes(identt,sxttt)+uni10.otimes(identt,syttt))

 H0*=(1./6.)

 dim = int(8)
 bdi = uni10.Bond(uni10.BD_IN, dim);
 bdo = uni10.Bond(uni10.BD_OUT, dim);
 H =  uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg");
 H.putBlock(H0)


 H0=copy.copy(H)
 H00=copy.copy(H)

##############################################################################
 E_ab=M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H0)
 print E_ab 
 E_ba=M_h(b_u,a_u,b,a,d,c,Env2,D,h,d_phys,chi,Corner_method,Model,H0)
 print E_ba
 E_cd=M_h(c_u,d_u,c,d,a,b,Env1,D,h,d_phys,chi,Corner_method,Model,H0)
 print E_cd
 E_dc=M_h(d_u,c_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,H0)
 print E_dc
 E_bd=M_v(b_u,d_u,d,c,b,a,Env3,D,h,d_phys,chi,Corner_method,Model,H00)
 print E_bd
 E_ca=M_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H00)
 print E_ca
 E_db=M_v(d_u,b_u,b,a,d,c,Env2,D,h,d_phys,chi,Corner_method,Model,H00)
 print E_db
 E_ac=M_v(a_u,c_u,c,d,a,b,Env1,D,h,d_phys,chi,Corner_method,Model,H00)
 print E_ac 
 return ((abs(E_ca)+abs(E_ac)+abs(E_db)+abs(E_bd)) / 8.00) + ((abs(E_ab)+abs(E_ba)+abs(E_cd)+abs(E_dc)) / 8.00)#+((E_val1+E_val2
###############################################################################



def M_v(c_u,a_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H00):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)
 
 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E_ca=basicB.Energy_ca(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H00,c_u,a_u)
 
 return E_ca




def M_h(a_u,b_u,a,b,c,d,Env,D,h,d_phys,chi,Corner_method,Model,H0):

 c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4=Init_env(Env)

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E_ab=basicB.Energy_ab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0,a_u,b_u)
 return E_ab































def Store_itebd(Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8):
 Gamma_a.save("Store/Gamma_a")
 Gamma_b.save("Store/Gamma_b")
 Gamma_c.save("Store/Gamma_c")
 Gamma_d.save("Store/Gamma_d")
 Landa_1.save("Store/Landa_1")
 Landa_2.save("Store/Landa_2")
 Landa_3.save("Store/Landa_3")
 Landa_4.save("Store/Landa_4")
 Landa_5.save("Store/Landa_5")
 Landa_6.save("Store/Landa_6")
 Landa_7.save("Store/Landa_7")
 Landa_8.save("Store/Landa_8")

def Reload_itebd():
 Gamma_a=uni10.UniTensor("Store/Gamma_a")
 Gamma_b=uni10.UniTensor("Store/Gamma_b")
 Gamma_c=uni10.UniTensor("Store/Gamma_c")
 Gamma_d=uni10.UniTensor("Store/Gamma_d")
 Landa_1=uni10.UniTensor("Store/Landa_1")
 Landa_2=uni10.UniTensor("Store/Landa_2")
 Landa_3=uni10.UniTensor("Store/Landa_3")
 Landa_4=uni10.UniTensor("Store/Landa_4")
 Landa_5=uni10.UniTensor("Store/Landa_5")
 Landa_6=uni10.UniTensor("Store/Landa_6")
 Landa_7=uni10.UniTensor("Store/Landa_7")
 Landa_8=uni10.UniTensor("Store/Landa_8")
 return Gamma_a,Gamma_b,Gamma_c,Gamma_d,Landa_1,Landa_2,Landa_3,Landa_4,Landa_5, Landa_6, Landa_7,Landa_8

def Store_Full(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.save("Store/a_u")
 b_u.save("Store/b_u")
 c_u.save("Store/c_u")
 d_u.save("Store/d_u")
 a.save("Store/a")
 b.save("Store/b")
 c.save("Store/c")
 d.save("Store/d")

def Store_Fullp(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.save("Store/ap_u")
 b_u.save("Store/bp_u")
 c_u.save("Store/cp_u")
 d_u.save("Store/dp_u")
 a.save("Store/ap")
 b.save("Store/bp")
 c.save("Store/cp")
 d.save("Store/dp")

def Reload_Fullp():
 ap_u=uni10.UniTensor("Store/ap_u")
 bp_u=uni10.UniTensor("Store/bp_u")
 cp_u=uni10.UniTensor("Store/cp_u")
 dp_u=uni10.UniTensor("Store/dp_u")
 ap=uni10.UniTensor("Store/ap")
 bp=uni10.UniTensor("Store/bp")
 cp=uni10.UniTensor("Store/cp")
 dp=uni10.UniTensor("Store/dp")
 return ap_u,bp_u,cp_u,dp_u,ap,bp,cp,dp


def Reload_Full():
 a_u=uni10.UniTensor("Store/a_u")
 b_u=uni10.UniTensor("Store/b_u")
 c_u=uni10.UniTensor("Store/c_u")
 d_u=uni10.UniTensor("Store/d_u")
 a=uni10.UniTensor("Store/a")
 b=uni10.UniTensor("Store/b")
 c=uni10.UniTensor("Store/c")
 d=uni10.UniTensor("Store/d")
 return a_u,b_u,c_u,d_u,a,b,c,d

def reconstruct_ab(a_u, a_uf):

 blk_qnums = a_u.blockQnum()
 blk_qnumsf = a_uf.blockQnum()

 for qnum in blk_qnums:
  if qnum in blk_qnumsf:
   mat_t=a_u.getBlock(qnum)
   dx=int(mat_t.row())
   dy=int(mat_t.col())
   mat_f=a_uf.getBlock(qnum)
   a_u.putBlock(qnum,mat_f.resize(dx,dy) )
 return a_u


def Reload_Full_previous(a_u, b_u, c_u, d_u):
 a_uf=uni10.UniTensor("Store/a_u")
 b_uf=uni10.UniTensor("Store/b_u")
 c_uf=uni10.UniTensor("Store/c_u")
 d_uf=uni10.UniTensor("Store/d_u")
 a_u=reconstruct_ab(a_u, a_uf)
 b_u=reconstruct_ab(b_u, b_uf)
 c_u=reconstruct_ab(c_u, c_uf)
 d_u=reconstruct_ab(d_u, d_uf)

 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)

 return a_u,b_u,c_u,d_u,a,b,c,d














def Rand_env_total(Env):
 Env1=copy.copy(Env)
 for i in xrange(len(Env1)):
  Env1[i]=copy.copy(Env[i])
  Env1[i].randomize()
 return  Env1


def qr_parity(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
    
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([bd1,bd2, theta.bond(4),theta.bond(5)])

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

    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())
    bd2=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([theta.bond(0),theta.bond(1),bd1,bd2])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

#    print LA
    return  LA, GA


def Decomposition(U):

 H=copy.copy(U)
 H.setLabel([0,1,2,3,4,5])
 H.permute([0,3,4,5,2,1],2)
 
 l,q=lq_parity(H)
 l.setLabel([0,3,-1,-2])

 q.setLabel([-1,-2,4,5,2,1])
 q.permute([-1,-2,4,1,5,2],4)
 qq, r=qr_parity(q)
 qq.setLabel([-1,-2,4,1,-3,-4])
 r.setLabel([-3,-4,5,2]) 
 
 l.permute([0,-1,-2,3],1)
 qq.permute([-1,-2,1,-4,-3,4],3)
 r.permute([-3,-4,2,5],3)
 
 MPO_list=[]
 MPO_list.append(l)
 MPO_list.append(qq)
 MPO_list.append(r)

 #H1=MPO_list[0]*MPO_list[1]*MPO_list[2]

 #H1.permute([0,1,2,3,4,5],3)
 #M=H1.getBlock()
 #D=M.eigh()
 #print D[0]
 #print "Test", H1.elemCmp(U), MPO_list[0],MPO_list[1],MPO_list[2]
 #print "hiihi",H1, D[0][0]
 return MPO_list

def slighty_random(a_u,b_u,c_u,d_u,a,b,c,d):
 rand=copy.copy(a_u)
 rand.randomize()
 a_u=a_u+(0.2)*rand

 rand=copy.copy(b_u)
 rand.randomize()
 b_u=b_u+(0.2)*rand
 
 rand=copy.copy(c_u)
 rand.randomize()
 c_u=c_u+(0.2)*rand
 
 rand=copy.copy(d_u)
 rand.randomize()
 d_u=d_u+(0.2)*rand


 a_u=a_u*(1.00/MaxAbs(a_u)) 
 b_u=b_u*(1.00/MaxAbs(b_u)) 
 c_u=c_u*(1.00/MaxAbs(c_u)) 
 d_u=d_u*(1.00/MaxAbs(d_u)) 
 
 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)
 
 return a_u,b_u,c_u,d_u,a,b,c,d

def total_random(a_u,b_u,c_u,d_u,a,b,c,d):
 a_u.randomize()
 b_u.randomize()
 c_u.randomize()
 d_u.randomize()


 a_u=a_u*(1.00/MaxAbs(a_u)) 
 b_u=b_u*(1.00/MaxAbs(b_u)) 
 c_u=c_u*(1.00/MaxAbs(c_u)) 
 d_u=d_u*(1.00/MaxAbs(d_u)) 
 
 a=make_ab(a_u)
 b=make_ab(b_u)
 c=make_ab(c_u)
 d=make_ab(d_u)

# a=a*(1.00/MaxAbs(a)) 
# b=b*(1.00/MaxAbs(b)) 
# c=c*(1.00/MaxAbs(c)) 
# d=d*(1.00/MaxAbs(d)) 
 
 return a_u,b_u,c_u,d_u,a,b,c,d


def max_ten(a):

 if ( MaxAbs(a) < 0.50e-1) or (MaxAbs(a) > 0.50e+1)   :
  a=a*(1.00/MaxAbs(a));
 else: a=a;
 
 return a

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


def positve_Env(Env,a,b,c,d):
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

 a1=copy.copy(a)
 a1.setLabel([1,-1,2,-2,3,-3,4,-4])
 a1.partialTrace(1,-1)
 a1.partialTrace(4,-4)
 a1.combineBond([2,-2])
 a1.combineBond([3,-3])
 a1.permute([2,3],1)
 Env[0]=a1


 b1=copy.copy(b)
 b1.setLabel([1,-1,2,-2,3,-3,4,-4])
 b1.partialTrace(3,-3)
 b1.partialTrace(4,-4)
 b1.combineBond([1,-1])
 b1.combineBond([2,-2])
 b1.permute([1,2],2)
 Env[1]=b1


 c1=copy.copy(c)
 c1.setLabel([1,-1,2,-2,3,-3,4,-4])
 c1.partialTrace(1,-1)
 c1.partialTrace(2,-2)
 c1.combineBond([3,-3])
 c1.combineBond([4,-4])
 c1.permute([3,4],0)
 Env[2]=c1


 d1=copy.copy(d)
 d1.setLabel([1,-1,2,-2,3,-3,4,-4])
 d1.partialTrace(2,-2)
 d1.partialTrace(3,-3)
 d1.combineBond([1,-1])
 d1.combineBond([4,-4])
 d1.permute([1,4],1)
 Env[3]=d1


 Ta1=copy.copy(Env[4])
 Ta2=copy.copy(Env[5])
 Ta3=copy.copy(Env[6])
 Ta4=copy.copy(Env[7])
 Tb1=copy.copy(Env[8])
 Tb2=copy.copy(Env[9])
 Tb3=copy.copy(Env[10])
 Tb4=copy.copy(Env[11])

 a1=copy.copy(a)
 a1.setLabel([1,-1,2,-2,3,-3,4,-4])
 a1.partialTrace(4,-4)
 a1.combineBond([1,-1])
 a1.combineBond([3,-3])
 a1.permute([1,2,-2,3],3)
 Env[8]=a1
 Env[9]=copy.copy(a1)
 Env[4]=copy.copy(a1)
 Env[5]=copy.copy(a1)

 a1=copy.copy(a)
 a1.setLabel([1,-1,2,-2,3,-3,4,-4])
 a1.partialTrace(1,-1)
 a1.combineBond([2,-2])
 a1.combineBond([4,-4])
 a1.permute([2,3,-3,4],1)
 Env[6]=a1
 Env[7]=a1
 Env[10]=a1
 Env[11]=a1
 return Env


def Store_Env(c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4):
 c1.save("Store/c1")
 c2.save("Store/c2")
 c3.save("Store/c3")
 c4.save("Store/c4")
 Ta1.save("Store/Ta1")
 Ta2.save("Store/Ta2")
 Ta3.save("Store/Ta3")
 Ta4.save("Store/Ta4")
 Tb1.save("Store/Tb1")
 Tb2.save("Store/Tb2")
 Tb3.save("Store/Tb3")
 Tb4.save("Store/Tb4")

def Store_EnvEnv(Env,Env1,Env2,Env3):

 Env[0].save("Store/c11")
 Env[1].save("Store/c21")
 Env[2].save("Store/c31")
 Env[3].save("Store/c41")
 Env[4].save("Store/Ta11")
 Env[5].save("Store/Ta21")
 Env[6].save("Store/Ta31")
 Env[7].save("Store/Ta41")
 Env[8].save("Store/Tb11")
 Env[9].save("Store/Tb21")
 Env[10].save("Store/Tb31")
 Env[11].save("Store/Tb41")

 Env1[0].save("Store/c12")
 Env1[1].save("Store/c22")
 Env1[2].save("Store/c32")
 Env1[3].save("Store/c42")
 Env1[4].save("Store/Ta12")
 Env1[5].save("Store/Ta22")
 Env1[6].save("Store/Ta32")
 Env1[7].save("Store/Ta42")
 Env1[8].save("Store/Tb12")
 Env1[9].save("Store/Tb22")
 Env1[10].save("Store/Tb32")
 Env1[11].save("Store/Tb42")

 Env2[0].save("Store/c13")
 Env2[1].save("Store/c23")
 Env2[2].save("Store/c33")
 Env2[3].save("Store/c43")
 Env2[4].save("Store/Ta13")
 Env2[5].save("Store/Ta23")
 Env2[6].save("Store/Ta33")
 Env2[7].save("Store/Ta43")
 Env2[8].save("Store/Tb13")
 Env2[9].save("Store/Tb23")
 Env2[10].save("Store/Tb33")
 Env2[11].save("Store/Tb43")

 Env3[0].save("Store/c14")
 Env3[1].save("Store/c24")
 Env3[2].save("Store/c34")
 Env3[3].save("Store/c44")
 Env3[4].save("Store/Ta14")
 Env3[5].save("Store/Ta24")
 Env3[6].save("Store/Ta34")
 Env3[7].save("Store/Ta44")
 Env3[8].save("Store/Tb14")
 Env3[9].save("Store/Tb24")
 Env3[10].save("Store/Tb34")
 Env3[11].save("Store/Tb44")

def Reload_EnvEnv(Env,Env1,Env2,Env3):

 Env[0]=uni10.UniTensor("Store/c11")
 Env[1]=uni10.UniTensor("Store/c21")
 Env[2]=uni10.UniTensor("Store/c31")
 Env[3]=uni10.UniTensor("Store/c41")
 Env[4]=uni10.UniTensor("Store/Ta11")
 Env[5]=uni10.UniTensor("Store/Ta21")
 Env[6]=uni10.UniTensor("Store/Ta31")
 Env[7]=uni10.UniTensor("Store/Ta41")
 Env[8]=uni10.UniTensor("Store/Tb11")
 Env[9]=uni10.UniTensor("Store/Tb21")
 Env[10]=uni10.UniTensor("Store/Tb31")
 Env[11]=uni10.UniTensor("Store/Tb41")

 Env1[0]=uni10.UniTensor("Store/c12")
 Env1[1]=uni10.UniTensor("Store/c22")
 Env1[2]=uni10.UniTensor("Store/c32")
 Env1[3]=uni10.UniTensor("Store/c42")
 Env1[4]=uni10.UniTensor("Store/Ta12")
 Env1[5]=uni10.UniTensor("Store/Ta22")
 Env1[6]=uni10.UniTensor("Store/Ta32")
 Env1[7]=uni10.UniTensor("Store/Ta42")
 Env1[8]=uni10.UniTensor("Store/Tb12")
 Env1[9]=uni10.UniTensor("Store/Tb22")
 Env1[10]=uni10.UniTensor("Store/Tb32")
 Env1[11]=uni10.UniTensor("Store/Tb42")

 Env2[0]=uni10.UniTensor("Store/c13")
 Env2[1]=uni10.UniTensor("Store/c23")
 Env2[2]=uni10.UniTensor("Store/c33")
 Env2[3]=uni10.UniTensor("Store/c43")
 Env2[4]=uni10.UniTensor("Store/Ta13")
 Env2[5]=uni10.UniTensor("Store/Ta23")
 Env2[6]=uni10.UniTensor("Store/Ta33")
 Env2[7]=uni10.UniTensor("Store/Ta43")
 Env2[8]=uni10.UniTensor("Store/Tb13")
 Env2[9]=uni10.UniTensor("Store/Tb23")
 Env2[10]=uni10.UniTensor("Store/Tb33")
 Env2[11]=uni10.UniTensor("Store/Tb43")

 Env3[0]=uni10.UniTensor("Store/c14")
 Env3[1]=uni10.UniTensor("Store/c24")
 Env3[2]=uni10.UniTensor("Store/c34")
 Env3[3]=uni10.UniTensor("Store/c44")
 Env3[4]=uni10.UniTensor("Store/Ta14")
 Env3[5]=uni10.UniTensor("Store/Ta24")
 Env3[6]=uni10.UniTensor("Store/Ta34")
 Env3[7]=uni10.UniTensor("Store/Ta44")
 Env3[8]=uni10.UniTensor("Store/Tb14")
 Env3[9]=uni10.UniTensor("Store/Tb24")
 Env3[10]=uni10.UniTensor("Store/Tb34")
 Env3[11]=uni10.UniTensor("Store/Tb44")


def Reload_Env():
 c1=uni10.UniTensor("Store/c1")
 c2=uni10.UniTensor("Store/c2")
 c3=uni10.UniTensor("Store/c3")
 c4=uni10.UniTensor("Store/c4")

 Ta1=uni10.UniTensor("Store/Ta1")
 Ta2=uni10.UniTensor("Store/Ta2")
 Ta3=uni10.UniTensor("Store/Ta3")
 Ta4=uni10.UniTensor("Store/Ta4")

 Tb1=uni10.UniTensor("Store/Tb1")
 Tb2=uni10.UniTensor("Store/Tb2")
 Tb3=uni10.UniTensor("Store/Tb3")
 Tb4=uni10.UniTensor("Store/Tb4")
 return c1,c2,c3,c4, Ta1, Ta2, Ta3, Ta4, Tb1, Tb2, Tb3, Tb4





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

def make_equall_bond(c1, c2,c3,c4, Tb3, Ta3, Ta1, Tb1):

 q_list=[]
 q_list1=[]
 q_list2=[]

 degs = c2.bond(0).degeneracy()
 for qnum, dim in degs.iteritems():
  if qnum.U1() < 0:
   for i in xrange(dim):
    q_list.append(qnum)
    q_list1.append(-qnum)
  if qnum.U1() == 0:
   for i in xrange(dim):
    q_list2.append(qnum)

 q_list1.reverse()
 q_list=q_list+q_list2+q_list1
 bd=uni10.Bond(uni10.BD_IN,q_list)


 blk_qnums = c2.blockQnum()
 bd_list=[bd,c2.bond(1)]
 c2b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c2b.getBlock(qnum).row())
  dy=int(c2b.getBlock(qnum).col())
  c2b.putBlock(qnum,c2.getBlock(qnum).resize(dx,dy) )





 bd=uni10.Bond(uni10.BD_OUT,c2b.bond(0).Qlist())
 blk_qnums = Ta1.blockQnum()
 bd_list=[Ta1.bond(0),Ta1.bond(1),Ta1.bond(2), bd]
 Ta1b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Ta1b.getBlock(qnum).row())
  dy=int(Ta1b.getBlock(qnum).col())
  Ta1b.putBlock(qnum,Ta1.getBlock(qnum).resize(dx,dy) )




 q_list=[]
 q_list1=[]
 q_list2=[]

 
 degs = c1.bond(1).degeneracy()
 for qnum, dim in degs.iteritems():
  if qnum.U1() < 0:
   for i in xrange(dim):
    q_list.append(qnum)
    q_list1.append(-qnum)
  if qnum.U1() == 0:
   for i in xrange(dim):
    q_list2.append(qnum)

 q_list1.reverse()
 q_list=q_list+q_list2+q_list1
 bd=uni10.Bond(uni10.BD_OUT,q_list)
 
 
 blk_qnums = c1.blockQnum()
 bd_list=[c1.bond(0),bd]
 c1b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c1b.getBlock(qnum).row())
  dy=int(c1b.getBlock(qnum).col())
  c1b.putBlock(qnum,c1.getBlock(qnum).resize(dx,dy) )

 
 
 bd=uni10.Bond(uni10.BD_IN,c1b.bond(1).Qlist())
 blk_qnums = Tb1.blockQnum()
 bd_list=[bd,Tb1.bond(1),Tb1.bond(2), Tb1.bond(3)]
 Tb1b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Tb1b.getBlock(qnum).row())
  dy=int(Tb1b.getBlock(qnum).col())
  Tb1b.putBlock(qnum,Tb1.getBlock(qnum).resize(dx,dy) )
  
 
 blk_qnums = c3.blockQnum()
 bd_list=[c2b.bond(0),c3.bond(1)]
 c3b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c3b.getBlock(qnum).row())
  dy=int(c3b.getBlock(qnum).col())
  c3b.putBlock(qnum,c3.getBlock(qnum).resize(dx,dy) )
 
 blk_qnums = c4.blockQnum()
 bd_list=[c4.bond(0),c1b.bond(1)]
 c4b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(c4b.getBlock(qnum).row())
  dy=int(c4b.getBlock(qnum).col())
  c4b.putBlock(qnum,c4.getBlock(qnum).resize(dx,dy) )



 bd=uni10.Bond(uni10.BD_OUT,c3b.bond(0).Qlist())
 blk_qnums = Tb3.blockQnum()
 bd_list=[Tb3.bond(0),Tb3.bond(1),Tb3.bond(2), bd]
 Tb3b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Tb3b.getBlock(qnum).row())
  dy=int(Tb3b.getBlock(qnum).col())
  Tb3b.putBlock(qnum,Tb3.getBlock(qnum).resize(dx,dy) )

 bd=uni10.Bond(uni10.BD_IN,c4b.bond(1).Qlist())
 blk_qnums = Ta3.blockQnum()
 bd_list=[bd,Ta3.bond(1),Ta3.bond(2), Ta3.bond(3)]
 Ta3b=uni10.UniTensor(bd_list)
 for qnum in blk_qnums:
  dx=int(Ta3b.getBlock(qnum).row())
  dy=int(Ta3b.getBlock(qnum).col())
  Ta3b.putBlock(qnum,Ta3.getBlock(qnum).resize(dx,dy) )



 return c1b, c2b, c3b, c4b, Tb3b, Ta3b, Ta1b, Tb1b

def rebond_corner(a,b,c,d,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):

 bd=uni10.Bond(uni10.BD_OUT,a.bond(0).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,a.bond(1).Qlist())

 tempo=copy.copy(Tb4)
 bd_list=[Tb4.bond(0),bd,bd1,Tb4.bond(3)]
 Tb4.assign(bd_list)
 blk_qnums = Tb4.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb4.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb4 dimension"      
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb4.putBlock(qnum,sv_mat.resize(dx,dy) )
 
 #print "Tb4",tempo.printDiagram(),Tb4.printDiagram(),a.printDiagram()
 #print "Tb4", tempo[10],Tb4[10],tempo.elemCmp(Tb4)#,Tb4.printDiagram(),
################################################################

 bd=uni10.Bond(uni10.BD_OUT,c.bond(0).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,c.bond(1).Qlist())
 
 tempo=copy.copy(Ta4)
 bd_list=[Ta4.bond(0),bd,bd1,Ta4.bond(3)]
 Ta4.assign(bd_list)
 blk_qnums = Ta4.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta4.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta4 dimension"      
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta4.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Ta4", Ta4.printDiagram(),Ta4[4]
 #print "Ta4", tempo.elemCmp(Ta4),Tb4.printDiagram(),

##################################################################

 bd=uni10.Bond(uni10.BD_IN,b.bond(4).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,b.bond(5).Qlist())

 tempo=copy.copy(Ta2)
 bd_list=[Ta2.bond(0),bd,bd1,Ta2.bond(3)]
 Ta2.assign(bd_list)
 blk_qnums = Ta2.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta2.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta2 dimension"      
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta2.putBlock(qnum,sv_mat.resize(dx,dy) )
# print "Ta2", Ta2.printDiagram(),Ta2[3]
# print "Ta2", tempo.elemCmp(Ta2),Tb4.printDiagram(),



 bd=uni10.Bond(uni10.BD_IN,d.bond(4).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,d.bond(5).Qlist())


 tempo=copy.copy(Tb2)
 bd_list=[Tb2.bond(0),bd,bd1,Tb2.bond(3)]
 Tb2.assign(bd_list)
 blk_qnums = Tb2.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb2.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb2 dimension"    
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb2.putBlock(qnum,sv_mat.resize(dx,dy) )
# print "Tb2", Tb2.printDiagram(),Tb2[4]
# print "Tb2", tempo.elemCmp(Tb2),Tb4.printDiagram(),

################################################################3
 bd=uni10.Bond(uni10.BD_IN,a.bond(6).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,a.bond(7).Qlist())

 tempo=copy.copy(Tb1)
 bd_list=[Tb1.bond(0),bd,bd1,Tb1.bond(3)]
 Tb1.assign(bd_list)
 blk_qnums = Tb1.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb1.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb1 dimension"  
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb1.putBlock(qnum,sv_mat.resize(dx,dy) )
# print "Tb1", Tb1.printDiagram(),Tb1[4]
# print "Tb1", tempo.elemCmp(Tb1),Tb4.printDiagram(),


 bd=uni10.Bond(uni10.BD_IN,b.bond(6).Qlist())
 bd1=uni10.Bond(uni10.BD_IN,b.bond(7).Qlist())

 tempo=copy.copy(Ta1)
 bd_list=[Ta1.bond(0),bd,bd1,Ta1.bond(3)]
 Ta1.assign(bd_list)
 blk_qnums = Ta1.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta1.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta1 dimension"  
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta1.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Ta1",tempo.printDiagram(),Ta1.printDiagram(),b.printDiagram()
 #print "Ta1", tempo[10],Ta1[10],tempo.elemCmp(Ta1)#,Tb4.printDiagram(),

######################################################
 bd=uni10.Bond(uni10.BD_OUT,c.bond(2).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,c.bond(3).Qlist())

 tempo=copy.copy(Ta3)
 bd_list=[Ta3.bond(0),bd,bd1,Ta3.bond(3)]
 Ta3.assign(bd_list)
 blk_qnums = Ta3.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Ta3.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Ta3 dimension"  
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Ta3.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Ta3", Ta3.printDiagram(),Ta3[4]
 #print "Ta3", tempo.elemCmp(Ta3),Tb4.printDiagram(),
 
 
 
 bd=uni10.Bond(uni10.BD_OUT,d.bond(2).Qlist())
 bd1=uni10.Bond(uni10.BD_OUT,d.bond(3).Qlist())

 tempo=copy.copy(Tb3)
 bd_list=[Tb3.bond(0),bd,bd1,Tb3.bond(3)]
 Tb3.assign(bd_list)
 blk_qnums = Tb3.blockQnum()
 blk_qnums1 = tempo.blockQnum()

 for qnum in blk_qnums:
  mat_t=Tb3.getBlock(qnum)
  dx=int(mat_t.row())
  dy=int(mat_t.col())
  if qnum in blk_qnums1: 
   sv_mat=tempo.getBlock(qnum)
  else:
   print "changing Tb3 dimension"
   sv_mat=uni10.Matrix(dx,dy)
   sv_mat.randomize()
  Tb3.putBlock(qnum,sv_mat.resize(dx,dy) )
 #print "Tb3", Tb3.printDiagram(),Tb3[4]
 #print "Tb3", tempo.elemCmp(Tb3),Tb4.printDiagram(),

################################################ 

 return Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4

def Short_TrotterSteps(N_iterF):
 List_delN=[]

 for i in xrange(6, 1, -1):
  Delta_N=(i*(0.50/10),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(5, 1, -1):
  Delta_N=(i*(1.0/100),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(5, 1, -1):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -1):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)


########################################
# Delta_N=((100.0/1000),N_iterF)
# List_delN.append(Delta_N)

# Delta_N=((9.0/1000),N_iterF)
# List_delN.append(Delta_N)

# Delta_N=((1.0/100),N_iterF)
# List_delN.append(Delta_N)

# Delta_N=((0.90/100),N_iterF)
# List_delN.append(Delta_N)

# Delta_N=((0.80/100),N_iterF)
# List_delN.append(Delta_N)

# Delta_N=((0.70/100),N_iterF)
# List_delN.append(Delta_N)

# Delta_N=((0.40/100),N_iterF)
# List_delN.append(Delta_N)


 return List_delN
 
def Short_TrotterSteps1(N_iterF):
 List_delN=[]

 for i in xrange(5, 0, -2):
  Delta_N=(i*(1.0/100),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -2):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -2):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)

 return List_delN 
 
def Long_TrotterSteps(N_iterF):
 List_delN=[]


 for i in xrange(10, 1, -1):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -1):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)

 return List_delN 
 
 
def Long_TrotterSteps1(N_iterF):
 List_delN=[]


 for i in xrange(10, 1, -2):
  Delta_N=(i*(1.0/1000),N_iterF)
  List_delN.append(Delta_N)

 for i in xrange(10, 0, -2):
  Delta_N=(i*(1.0/10000),N_iterF)
  List_delN.append(Delta_N)

 return List_delN 
 
#@profile
def Obtain_grad_four(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U):

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
###########################################---a---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E1*E8)*(a_dp))
 A.permute([55,16,17,18,2],3)
 D_a=A
 A=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E1*E8)*(a_up))
 A.permute([-18,-2,55,-16,-17],2)
 A.transpose()
 D_a=D_a+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

# t0=time.time()

 A=((((((E2*E3)*(b_up*b_d)))*((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_up*c_d)*U))*((E1*E8)*(a_d))
 A.permute([55,16,17,18,2],3)
 D_a=D_a+(-1.0)*A

# print "1", time.time() - t0


 A=((((((E2*E3)*(b_u*b_dp)))*((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u*c_dp)*U))*((E1*E8)*(a_u))
 A.permute([-18,-2,-55,-16,-17],2)
 A.transpose()
 D_a=D_a+(-1.0)*A
 D_a.transpose()
 D_a.permute([55,16,17,18,2],3)
 
 
#####################################----b----#############################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))*(((E2*E3)*(b_dp)))
 A.permute([56,18,20,6,4],3)
 D_b=A
 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))*(((E2*E3)*(b_up)))
 A.permute([-6,-4,56,-18,-20],2)
 A.transpose()
 D_b=D_b+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

# t0=time.time()

 A=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*((E4*E5)*(d_up*d_d)))*(((E2*E3)*(b_d)))
 A.permute([56,18,20,6,4],3)
 D_b=D_b+(-1.0)*A

# print "2", time.time() - t0


 A=(((((E1*E8)*(a_u*a_dp)*U)*((E7*E6)*(c_u*c_dp))))*((E4*E5)*(d_u*d_dp)))*(((E2*E3)*(b_u)))
 A.permute([-6,-4,-56,-18,-20],2)
 A.transpose()
 D_b=D_b+(-1.0)*A
 D_b.transpose()
 D_b.permute([56,18,20,6,4],3)

##################################---c---#################################################3


 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))* ((E7*E6)*(c_dp))
 A.permute([54,14,12,19,17],3)
 D_c=A
 A=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))* ((E7*E6)*(c_up))
 A.permute([-19,-17,54,-14,-12],2)
 A.transpose()
 D_c=D_c+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

# t0=time.time()

 A=(((((E1*E8)*(a_up*a_d))*((E2*E3)*(b_up*b_d)))*U)*(((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_d))
 A.permute([54,14,12,19,17],3)
 D_c=D_c+(-1.0)*A

# print "3", time.time() - t0


 A=(((((E1*E8)*(a_u*a_dp))*((E2*E3)*(b_u*b_dp)))*U)*(((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u))
 A.permute([-19,-17,-54,-14,-12],2)
 A.transpose()
 D_c=D_c+(-1.0)*A
 D_c.transpose()
 D_c.permute([54,14,12,19,17],3)
 
 
###################################----d----######################################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_dp))
 A.permute([57,19,10,8,20],3)
 D_d=A
 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up))
 A.permute([-8,-20,57,-19,-10],2)
 A.transpose()
 D_d=D_d+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,-54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

# t0=time.time()

 A=(((((E1*E8)*(a_up*a_d)*U)*((E7*E6)*(c_up*c_d))))*(((E2*E3)*(b_up*b_d))))*((E4*E5)*(d_d))
 A.permute([57,19,10,8,20],3)
 D_d=D_d+(-1.0)*A

# print "4", time.time() - t0


 A=(((((E1*E8)*(a_u*a_dp)*U)*((E7*E6)*(c_u*c_dp))))*(((E2*E3)*(b_u*b_dp))))*((E4*E5)*(d_u))
 A.permute([-8,-20,57,-19,-10],2)
 A.transpose()
 D_d=D_d+(-1.0)*A
 D_d.transpose()
 D_d.permute([57,19,10,8,20],3)
 
##################################################################################################
 return D_a, D_b, D_c, D_d 
 
 
def Obtain_grad_four1(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c, d, a_u, b_u, c_u, d_u, a_up, b_up, c_up, d_up, U):

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
###########################################---a---##############################################
 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E1*E8)*(a_dp))
 A.permute([55,16,17,18,2],3)
 D_a=A
 A=(((((E4*E5)*(d_up*d_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E1*E8)*(a_up))
 A.permute([-18,-2,55,-16,-17],2)
 A.transpose()
 D_a=D_a+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=((((((E2*E3)*(b_up*b_d))*U)*((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_up*c_d)))*((E1*E8)*(a_d))
 A.permute([55,16,17,18,2],3)
 D_a=D_a+(-1.0)*A


 A=((((((E2*E3)*(b_u*b_dp))*U)*((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u*c_dp)))*((E1*E8)*(a_u))
 A.permute([-18,-2,-55,-16,-17],2)
 A.transpose()
 D_a=D_a+(-1.0)*A
 D_a.transpose()
 D_a.permute([55,16,17,18,2],3)
 
 
#####################################----b----#############################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))*(((E2*E3)*(b_dp)))
 A.permute([56,18,20,6,4],3)
 D_b=A
 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*((E4*E5)*(d_up*d_dp)))*(((E2*E3)*(b_up)))
 A.permute([-6,-4,56,-18,-20],2)
 A.transpose()
 D_b=D_b+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d))))*((E4*E5)*(d_up*d_d))*U)*(((E2*E3)*(b_d)))
 A.permute([56,18,20,6,4],3)
 D_b=D_b+(-1.0)*A


 A=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp))))*((E4*E5)*(d_u*d_dp))*U)*(((E2*E3)*(b_u)))
 A.permute([-6,-4,-56,-18,-20],2)
 A.transpose()
 D_b=D_b+(-1.0)*A
 D_b.transpose()
 D_b.permute([56,18,20,6,4],3)

##################################---c---#################################################3


 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))* ((E7*E6)*(c_dp))
 A.permute([54,14,12,19,17],3)
 D_c=A
 A=(((((E1*E8)*(a_up*a_dp))*((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up*d_dp)))* ((E7*E6)*(c_up))
 A.permute([-19,-17,54,-14,-12],2)
 A.transpose()
 D_c=D_c+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E2*E3)*(b_up*b_d)))*U)*(((E4*E5)*(d_up*d_d))))*((E7*E6)*(c_d))
 A.permute([54,14,12,19,17],3)
 D_c=D_c+(-1.0)*A


 A=(((((E1*E8)*(a_u*a_dp))*((E2*E3)*(b_u*b_dp)))*U)*(((E4*E5)*(d_u*d_dp))))*((E7*E6)*(c_u))
 A.permute([-19,-17,54,-14,-12],2)
 A.transpose()
 D_c=D_c+(-1.0)*A
 D_c.transpose()
 D_c.permute([54,14,12,19,17],3)
 
 
###################################----d----######################################################

 a_dp.setLabel([-18,-2,55,-16,-17])
 b_dp.setLabel([-6,-4,56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,57,-19,-10])

 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_dp))
 A.permute([57,19,10,8,20],3)
 D_d=A
 A=(((((E1*E8)*(a_up*a_dp))*((E7*E6)*(c_up*c_dp))))*(((E2*E3)*(b_up*b_dp))))*((E4*E5)*(d_up))
 A.permute([-8,-20,57,-19,-10],2)
 A.transpose()
 D_d=D_d+A

 a_dp.setLabel([-18,-2,-55,-16,-17])
 b_dp.setLabel([-6,-4,-56,-18,-20]) 
 c_dp.setLabel([-19,-17,54,-14,-12])
 d_dp.setLabel([-8,-20,-57,-19,-10])


 A=(((((E1*E8)*(a_up*a_d))*((E7*E6)*(c_up*c_d))))*(((E2*E3)*(b_up*b_d))*U))*((E4*E5)*(d_d))
 A.permute([57,19,10,8,20],3)
 D_d=D_d+(-1.0)*A


 A=(((((E1*E8)*(a_u*a_dp))*((E7*E6)*(c_u*c_dp))))*(((E2*E3)*(b_u*b_dp))*U))*((E4*E5)*(d_u))
 A.permute([-8,-20,-57,-19,-10],2)
 A.transpose()
 D_d=D_d+(-1.0)*A
 D_d.transpose()
 D_d.permute([57,19,10,8,20],3)
 
##################################################################################################
 return D_a, D_b, D_c, D_d 
def Obtain_grad_four_MPO(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, MPO_list,c_u,a_u,b_u,plist):

 D_r=[0]*4
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
##################################  1  #####################################################
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
 A=((((E4*E5)*d)*((E2*E3)*(b_ut*b_dt)))*(((E1*E8)*(a_ut*a_dt))))*((E7*E6)*(c_ut*c_dt))

 A.permute([-62,-58,-57,-17,62,58,57,17],4)

 xt=copy.copy(plist[0])
 xt.setLabel([-62,-58,-57,-17])
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
 
 
 Ap=((((E1*E8)*(a_u1*a_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt)))*(((E7*E6)*(c_u1*c_dt)))
 Ap.permute([-62,-58,-57,-17],4)
 Ap.transpose()
 D_r[0]=D_r[0]+(-1.00)*Ap

 Ap=((((E1*E8)*(a_ut*a_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1)))*(((E7*E6)*(c_ut*c_d1)))
 Ap.permute([62,58,57,17],0)
 D_r[0]=D_r[0]+(-1.00)*Ap
 D_r[0].permute([62,58,57,17],3)
############################################################################################

##################################   2   #############################################
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

 xt=copy.copy(plist[1])
 xt.setLabel([-17,-64,-58,-57])
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
 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-17,-64,-58,-57],4)
 Ap.transpose()
 D_r[1]=D_r[1]+(-1.00)*Ap

 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([17,64,58,57],0)
 D_r[1]=D_r[1]+(-1.00)*Ap
 D_r[1].permute([17,64,58,57],1)
############################################################################################


##################################  3  #####################################################
 a_ut=(a_u*(plist[1]*MPO_list[1]))
 c_ut=((c_u*(plist[0]*MPO_list[0])))
 b_ut=((b_u*(plist[3]*MPO_list[2])))

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

 xt=copy.copy(plist[2])
 xt.setLabel([-66,-59,-60,-18]) 
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

 Ap=(((((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)*((E2*E3)*(b_u1*b_dt))))*((E1*E8)*(a_u1*a_dt))
 Ap.permute([-66,-59,-60,-18],4)
 Ap.transpose()
 D_r[2]=D_r[2]+(-1.00)*Ap

 Ap=(((((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)*((E2*E3)*(b_ut*b_d1))))*((E1*E8)*(a_ut*a_d1))
 Ap.permute([66,59,60,18],0)
 D_r[2]=D_r[2]+(-1.00)*Ap
 D_r[2].permute([66,59,60,18],3)
############################################################################################

##################################  4  #####################################################
 c_ut=((c_u*(plist[0]*MPO_list[0])))
 a_ut=(a_u*(plist[1]*plist[2]*MPO_list[1]))
 b_ut=((b_u*(MPO_list[2])))

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

 xt=copy.copy(plist[3])
 xt.setLabel([-18,-68,-59,-60]) 
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
 


 Ap=(((((E1*E8)*(a_u1*a_dt))*((E7*E6)*(c_u1*c_dt))))*(((E4*E5)*d)))*((E2*E3)*(b_u1*b_dt))
 Ap.permute([-18,-68,-59,-60],4)
 Ap.transpose()
 D_r[3]=D_r[3]+(-1.00)*Ap


 
 Ap=(((((E1*E8)*(a_ut*a_d1))*((E7*E6)*(c_ut*c_d1))))*(((E4*E5)*d)))*((E2*E3)*(b_ut*b_d1))
 Ap.permute([18,68,59,60],0)
 D_r[3]=D_r[3]+(-1.00)*Ap
 D_r[3].permute([18,68,59,60],1)
############################################################################################
 return D_r 

