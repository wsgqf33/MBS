import numpy as np
import math
import InterestrateModel as IRM
from scipy.optimize import fsolve

import Calibration
# from Calibration import *
from MarketData import *
from DiscountCurve import *
from CDS import *
from HazardRateSim import *

from os.path import join, dirname, abspath,exists
import pickle


class Prepayment:
    def __init__(self, PV0, WAC, rt, fee, T, Model):
        self.Model=Model
        self.__PV0 = PV0
        self.__r = WAC/12
        self.__N = T*12
        self.__f=fee

        self.__rt  = rt
        self.__ndiv = rt.shape[1]
        self.__nsims = rt.shape[0]
        self.__dT = 1/360
        self.__mr = rt

        # self.__default=Calibration.Default(self.__nsims,self.__N)
        # self.__default.setENV()
        # self._Sur = self.__default.SurvivalFun()

        self._Sur = self.getSurvival()

        # print(self._Sur[1][0],self._Sur[1][1],self._Sur[1][2],self._Sur[1][3])
        # print(self._Sur[2][0], self._Sur[2][1], self._Sur[2][2], self._Sur[2][3])

        self.__PV = np.zeros([self.__nsims,self.__N+1])
        self.__PV[:,0] = PV0
        self.__c = np.zeros([self.__nsims, self.__N+1])
        self.__IP = np.zeros([self.__nsims, self.__N+1])
        self.__TPP = np.zeros([self.__nsims, self.__N+1])
        self._SMM=np.zeros([self.__nsims, self.__N+1])

    def setInterestRateModel(self,IR):
        self.__IR = IR


    def CPR(self, t):
        if self.Model == "Multi":
            return self.Ref_Incentive(t) * self.Burnout(t) * self.Seasoning(t) * self.Seasonality(t)
        elif self.Model == "PSA":
            if t <= 30:
                return 0.002 * t
            else:
                return 0.06
        elif self.Model == "PSA1.5":
            if t <= 30:
                return 0.003 * t
            else:
                return 0.09
        elif self.Model == "PSA2":
            if t <= 30:
                return 0.004 * t
            else:
                return 0.12
            return 0

    def SMM(self,t,flag):
        if flag=="Multi" or flag =="PSA" or flag =="PSA1.5" or flag =="PSA2":
            self._SMM[:,t]=1 - np.power(1 - self.CPR(t), 1 / 12)
            return self._SMM[:, t]
        # print(type(self.__mr[:,t]))

        if t==0:
            self._SMM[:, t] = 0
        else:
            fee=self.__f
            logic1 = np.where((self.Mor_R(t) < self.__r * 12 - fee) & (self.Mor_R(t - 1) > self.Mor_R(t)))
            logic2 = np.where((self.Mor_R(t - 1) < self.__r * 12 - fee) & (self.Mor_R(t - 1) > self.Mor_R(t)))

            self._SMM[np.where(self.Mor_R(t) > self.__r * 12 - fee)] = 0
            self._SMM[np.where(self.Mor_R(t) > self.Mor_R(t - 1)), t] = 0

            self._SMM[logic1, t] = 0.05 * (self.__r * 12 - fee - self.Mor_R(t))
            self._SMM[logic2, t] = self._SMM[logic2, t - 1] + 0.04 * (self.__r * 12 - fee - self.Mor_R(t))


            #
            #
            # logic1 = np.where((self.__mr[:, t] < self.__r*12 - fee) & (self.__mr[:, t - 1] > self.__mr[:, t]))
            # logic2 = np.where((self.__mr[:, t - 1] < self.__r*12 - fee) & (self.__mr[:, t - 1] > self.__mr[:, t]))
            #
            # self._SMM[np.where(self.__mr[:, t] > self.__r*12 - fee)]=0
            # self._SMM[np.where(self.__mr[:, t] > self.__mr[:, t - 1]), t] = 0
            #
            # self._SMM[logic1, t] = 0.05*(self.__r*12-fee-self.__mr[logic1, t])
            # self._SMM[logic2, t] = self._SMM[logic2, t-1]+0.04*(self.__r*12-fee-self.__mr[logic2, t])
            #
            # #



            # self._SMM[np.where(self.__mr[:,t-1]>self.__mr[:,t]), t] = 0.05 * (self.__r - 0.275 - self.__mr[:, t])
            # self._SMM[np.where(self.__mr[:,t-1]<self.__r-0.275), t] = 0.05*(self.__r-0.275-self.__mr[:,t])
            # self._SMM[np.where(self.__mr[:,t-1]>self.__mr[:,t]), t] = 0.05 * (self.__r - 0.275 - self.__mr[:, t])
        # if t==0 or self.__r-0.275<self.__mr[:,t] or self.__mr[:,t]>self.__mr[:,t-1]:
        #     self._SMM[:,t]= 0
        # elif self.__mr[:,t]<self.__r-0.275 and self.__mr[:,t-1]>self.__mr[:,t]:
        #     self._SMM[:,t]=0.05*(self.__r-0.275-self.__mr[:,t])
        # else:
        #     self._SMM[:,t]=self._SMM[:,t-1]+0.04*(self.__r-0.275-self.__mr[:,t])

        return self._SMM[:, t]



    # Schedual Principal Payment
    def SPRIN(self, t):  ## t in months
        disc= 1/np.power(1+self.__r,self.__N - (t -1))
        return self.__PV[:,t-1]*self.__r*(1/(1- disc) -1)

    # Principal Prepayment
    def PP(self,t):
        # SMM=1 - np.power(1 - self.CPR(t), 1 / 12)
        return (self.__PV[:,t-1] - self.SPRIN(t)) * self.SMM(t,self.Model)

    # Interest Payment
    def IP(self,t):
        return self.__PV[:,t-1]*self.__r

    def setPV(self,t):
        self.__PV[:,t] = self.__PV[:,t-1] - (self.SPRIN(t) + self.PP(t))


    def r_10yr(self,t):## t in months
        #ind = (int)(((t - 1)) / self.__N) * self.__ndiv  # t is in years. t-1 corresponds to one month before
        ind = (t-1)*30
        return -np.log(self.__IR.BondPrice_ClosedForm(self.__rt[:,ind],t/12,t/12+10))/10


    def r_2yr(self,t):## t in months
        ind = (t-1)*30
        return -np.log(self.__IR.BondPrice_ClosedForm(self.__rt[:,ind],t/12,t/12+2))/10

    def Ref_Incentive(self,t): ## t in months
        R = 12*self.__r
        return 0.28 + 0.14 * np.arctan(-8.57 + 430 * (R - self.r_10yr(t)))

    def Burnout(self,t): ## t in months
        return 0.3 + 0.7*self.__PV[:,t-1]/self.__PV0

    def Seasoning(self,t): ## t in months
        return min(1,t/30)

    def Seasonality(self,t): ## t in months
        SY = [0.98, 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22,1.23]  # SY[0] - corresponds to Dec
        return SY[t % 12]


    def setc(self,t):

        self.__TPP[:,t] = self.SPRIN(t) + self.PP(t)
        self.__IP[:,t] = self.IP(t)
        #self.__c[:,t] =self.SP(t) + self.PP(t) + self.IP(t)
        # self.__c[:, t] = self.__TPP[:,t] + self.__IP[:,t]
        self.__c[:, t] = np.multiply(self.__TPP[:,t] + self.__IP[:,t]*(self.__r/(self.__r+self.__f)),self._Sur[:,t])
        # self.__c[:, t] = self.__TPP[:, t] + self.__IP[:, t]


    def generateC(self):
        for i in range(self.__N):
            self.setPV(i+1)
            self.setc(i+1)

    def setdt(self,x):
        '''
        rtx = self.__rt[:,0:(self.__N+1)] + x
        #self.__d= np.exp(-np.true_divide(rtx.cumsum(1), np.arange(1, rtx.shape[1] + 1)))
        self.__d = np.exp(-rtx.cumsum(1)*self.__dT)
        '''
        self.__d = np.zeros([self.__nsims, self.__N+1])
        rtx = self.__rt +x
        nskip = 30
        for i in range(self.__N):
            #self.__d[:, i+1] = np.exp(-np.mean(self.__rt[:, 0:(nskip * (i+1))], 1)*(i+1)/nskip)
            self.__d[:, i+1] = np.exp(-np.sum(rtx[:, 0:(nskip * (i+1))], 1)/360)



    def MBSPrice(self):
        self.setdt(0)
        self.generateC()
        P= np.sum(self.__c*self.__d,1)
        #print(np.mean(P))
        return np.mean(P)

    def MBS_IOPrice(self):
        self.setdt(0)
        self.generateC()
        P = np.sum(self.__IP * self.__d, 1)
        # print(np.mean(P))
        return np.mean(P)

    def MBS_POPrice(self):
        self.setdt(0)
        self.generateC()
        P = np.sum(self.__TPP * self.__d, 1)
        # print(np.mean(P))
        return np.mean(P)



    def solveSpread(self,x,MBSPrice):
        '''
        dx = np.zeros([self.__nsims,self.__N + 1])
        nskip = (int)(365 / 12)
        rtx = self.__rt + x
        for j in range(self.__N):
            dx[:,j + 1] = np.exp(-np.mean(rtx[:,0:(nskip * (j + 1))]) * (j + 1) / 12)
        '''
        self.setdt(x)
        P = np.mean(np.sum(self.__c * self.__d,1))
        return P - MBSPrice

    def OASSpread(self,MBSPrice):
        self.generateC()
        x = fsolve(self.solveSpread,-0.013,MBSPrice)
        return x

    def getDurationConvexity(self, ybps, MBSPrice):
        x = self.OASSpread(MBSPrice)
        P1 = self.solveSpread(x+ybps/10000,MBSPrice) + MBSPrice
        P2 = self.solveSpread(x-ybps/10000,MBSPrice) + MBSPrice
        duration = 10000 * (P2 - P1) / (2 * ybps * MBSPrice)
        convexity = math.pow(10, 8) * (P2 + P1 - 2 * MBSPrice) / (2 * ybps * ybps * MBSPrice)
        return duration,convexity

    def getSurvival(self):
        ### Hazard Rate Calibration
        address=join(dirname(dirname(abspath(__file__))), 'MBS', 'CIRPara.pkl')

        if  exists(address):
            print("Hazard CIR Parameter Exist...")
            CIRpara = pickle.load(open(address, "rb"))
        else:
            print("Hazard CIR Parameter Miss, Re Calibration...")
            # CDS spreads from Bloomberg (tenors,spreads）
            spreads = {'Date': '12/12/2016', '1': '9.99', '2': '16.92', '3': '21.605', '4': '25.885', '5': '27.45',
                       '7': '32.03'}
            z = MarketData(spreads)
            t = ql.Date(12, 12, 2016)
            CIR = Calibration.CDSCalibration(t, DiscountCurve=OISDisountCurve,
                                 MarketData=z,
                                 CDS=CIRCreditDefaultSwap,
                                 Process="CIR",
                                 Guess=[0.37, 0.032, 0.03, 0.002]
                                 )
            CIRpara = CIR.Calibrate()
            pickle.dump(CIRpara, open(address, "wb"))


        ## Hazard Rate Simulation
        Q = HazardRateSim(CIRpara, self.__nsims, self.__N/12, self.__dT)
        print("Survival Matrix Completed")
        # for i in range(10):
        #     plt.plot(self.__rt[i, :])
        # # plt.suptitle('MBS Survival Rate')
        # plt.xlabel('Days')
        # plt.ylabel('2 Factor Hull White Rate')
        # plt.savefig('IR.eps')
        # plt.show()

        return Q

    def Mor_R(self,t):
        return 0.024+0.2*self.r_2yr(t)+0.6*self.r_10yr(t)