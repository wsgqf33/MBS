import matplotlib.pyplot as plt
import InterestrateModel as IRM
import PrepaymentModel as PPM
import numpy as np
import pandas as pd
import QuantLib as ql
import timeit
import random
import math
from Calibration import *

from os.path import join, dirname, abspath,exists
import pickle


#1. MBS Parameter**********************************************************
# r0=0.42
r0 = 0.078
r0=r0/100
fee=0.005
T = 30
WAC = 0.04
PV0 = 1
dT = 1/360
nsims =2000
t=ql.Date(12,12, 2016)
# 2.Parameter Calibration***************************************************

address=join(dirname(dirname(abspath(__file__))), 'MBS', 'G2PPPara.pkl')
address2 = join(dirname(dirname(abspath(__file__))), 'MBS', 'MCtest.pkl')
if  exists(address):
    print("G2PP Parameter Exist...")
    [a, sigma, b, eta, rho]= pickle.load(open(address, "rb"))
else:
    print("G2PP Parameter Miss, Re Calibration")
    G2ppCal=Calibration()
    [a, sigma, b, eta , rho]=G2ppCal.G2pp(t)
    para=[a, sigma, b, eta , rho]
    pickle.dump(para, open(address, "wb"))
    print("G2PP Parameter Completed")





# 3 Multiplicative Prepayment Model******************************************************************************
# [a, sigma, b,eta,rho]=[0.04190596403020456, 0.004764227424705504, 0.050107680884392575, 0.003037434233114634, 0.03442977268843137]
start_time1 = timeit.default_timer()
g2pp=IRM.G2pp(r0=r0,x0=0,y0=0, T=T+10, dT=dT,a=a,b=b, sigma=sigma,eta=eta,rho=rho,nsims=nsims)
r_sim2=g2pp.EulerDiscretize_R()

numerixPPM = PPM.Prepayment(PV0, WAC, r_sim2, fee, T, "Multi")
numerixPPM.setInterestRateModel(g2pp)
MBSPrice = numerixPPM.MBSPrice()

del numerixPPM
del g2pp

elapsed1 = timeit.default_timer() - start_time1

print('a) HullWhite_Multiplicative_Prepayment_MBS = %f      ****[%f sec]'%(MBSPrice,elapsed1))



### PCA100 Prepayment Model *****************************************************************************************
start_time1 = timeit.default_timer()

g2pp=IRM.G2pp(r0=r0,x0=0,y0=0, T=T+10, dT=dT,a=a,b=b, sigma=sigma,eta=eta,rho=rho,nsims=nsims)
r_sim2 = g2pp.EulerDiscretize_R()

PCA100PPM = PPM.Prepayment(PV0, WAC, r_sim2, fee, T, "PSA")
PCA100PPM.setInterestRateModel(g2pp)
MBSPrice = PCA100PPM.MBSPrice()
del g2pp
del PCA100PPM
elapsed1 = timeit.default_timer() - start_time1

print('b) HullWhite_100PCA_Prepayment_MBS = %f      ****[%f sec]'%(MBSPrice,elapsed1))

### PCA100 Prepayment Model
start_time1 = timeit.default_timer()

g2pp=IRM.G2pp(r0=r0,x0=0,y0=0, T=T+10, dT=dT,a=a,b=b, sigma=sigma,eta=eta,rho=rho,nsims=nsims)
r_sim2 = g2pp.EulerDiscretize_R()

PCA100PPM = PPM.Prepayment(PV0, WAC, r_sim2,fee, T, "PSA1.5")
PCA100PPM.setInterestRateModel(g2pp)
MBSPrice = PCA100PPM.MBSPrice()
del g2pp
del PCA100PPM
elapsed1 = timeit.default_timer() - start_time1

print('c) HullWhite_150PCA_Prepayment_MBS = %f      ****[%f sec]'%(MBSPrice,elapsed1))

### PCA100 Prepayment Model
start_time1 = timeit.default_timer()

g2pp=IRM.G2pp(r0=r0,x0=0,y0=0, T=T+10, dT=dT,a=a,b=b, sigma=sigma,eta=eta,rho=rho,nsims=nsims)
r_sim2 = g2pp.EulerDiscretize_R()

PCA100PPM = PPM.Prepayment(PV0, WAC, r_sim2,fee, T, "PSA2")
PCA100PPM.setInterestRateModel(g2pp)
MBSPrice = PCA100PPM.MBSPrice()
del g2pp
del PCA100PPM
elapsed1 = timeit.default_timer() - start_time1

print('d) HullWhite_200PCA_Prepayment_MBS = %f      ****[%f sec]'%(MBSPrice,elapsed1))
import matplotlib.pyplot as plt
for i in range(10):
    plt.plot(r_sim2[i, :])
plt.suptitle('2 Factor Hull White Interest Rate')
plt.xlabel('Days')
plt.ylabel('2 Factor Hull White Rate')
plt.show()

#******************************************** MC TEST*******************************************
# if exists(address2):
#     print("Already Done MC test")
#     aa =pickle.load(open(address2, "rb"))
#     print(aa.keys())
#
# else:
#     MCresult = {}
#     path = []
#     test = [100,250,500,1000,2000]
#     for npath in test:
#         nsims = npath
#         price = []
#         time = []
#
#         for i in range(20):
#             start_time1 = timeit.default_timer()
#             g2pp = IRM.G2pp(r0=r0, x0=0, y0=0, T=T + 10, dT=dT, a=a, b=b, sigma=sigma, eta=eta, rho=rho, nsims=nsims)
#             r_sim2 = g2pp.EulerDiscretize_R()
#             numerixPPM = PPM.Prepayment(PV0, WAC, r_sim2, fee, T, "Multi")
#             numerixPPM.setInterestRateModel(g2pp)
#             MBSPrice = numerixPPM.MBSPrice()
#
#             elapsed1 = timeit.default_timer() - start_time1
#             time.append(elapsed1)
#             price.append(MBSPrice)
#             print(i, "finish path", nsims, "......")
#             del numerixPPM
#             del g2pp
#         print(nsims, " paths test is done")
#         pt = [price, time]
#
#         MCresult.update({str(nsims): pt})
#
#     pickle.dump(MCresult, open(address2, "wb"))
#     print("MCResults test", MCresult.keys())
#     print(MCresult["100"])
#     print(MCresult["250"])
#     print(MCresult["500"])
#     print(MCresult["1000"])
#     print(MCresult["2000"])



# ### Simple Prepayment Model
# start_time1 = timeit.default_timer()
# g2pp=IRM.G2pp(r0=r0,x0=0,y0=0, T=T+10, dT=dT,a=a,b=b, sigma=sigma,eta=eta,rho=rho,nsims=nsims)
# r_sim2 = g2pp.EulerDiscretize_R()
#
# SimplePPM = PPM.Prepayment(PV0, WAC, r_sim2, fee, T, "Simple")
# SimplePPM.setInterestRateModel(g2pp)
# MBSPrice = SimplePPM.MBSPrice()
#
# del g2pp
# del SimplePPM
# elapsed1 = timeit.default_timer() - start_time1
#
# print('b) HullWhite_Simple_Prepayment_MBS = %f      ****[%f sec]'%(MBSPrice,elapsed1))
#

