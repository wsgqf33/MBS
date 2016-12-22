import math
import numpy as np
from scipy.stats import ncx2
import random

import matplotlib as plt
import matplotlib.gridspec as gridspec
from numpy import *
import math
import matplotlib.pyplot as plt



class InterestRateModel:
    def __init__(self, r0,T, dT, sigma,nsims):
        self._r0 =r0
        self._sigma = sigma
        self.T = T
        self.dT = dT
        self._ndiv = (int)(T/dT)
        self._nsims = nsims

    @property
    def T(self):
        return self._T


    @T.setter
    def T(self, T):
        if T < 0:
            raise AttributeError('T cannot be negative')
        else:
            self._T = T


    @property
    def dT(self):
        return self._dT


    @dT.setter
    def dT(self, dT):
        if dT < 0 or dT > self._T:
            raise AttributeError('dT cannot be negative or grater than t')
        else:
            self._dT = dT

class Vasicek(InterestRateModel):

    def __init__(self, r0, r_bar, T, dT, sigma,nsims,kappa):
        InterestRateModel.__init__(self, r0,T, dT, sigma,nsims)
        self.__kappa = kappa
        self.__rbar = r_bar

    def EulerDiscretize_R(self):
        self.__r = np.zeros([self._nsims,(self._ndiv+1)])
        self.__r[:,0] = self._r0
        num = 1
        while(num <= self._ndiv):
            z = np.random.normal(0,1,self._nsims)
            self.__r[:,num]= self.__r[:,(num-1)] +  self.__kappa*(self.__rbar - self.__r[:,(num-1)])*self._dT + self._sigma*math.sqrt(self._dT)*z
            num = num+1

        return self.__r


    def getA_B(self,t, T):
        B = 1/self.__kappa*(1- np.exp(-self.__kappa*(T - t)))
        A = np.exp((self.__rbar - (self._sigma**2)/(2*(self.__kappa**2)))*(B - (T-t)) - ((self._sigma*B/2)**2)/self.__kappa)
        return A,B

    def getsigma_p(self,t0,t,T):
        d = np.sqrt((1- np.exp(-2*self.__kappa*(T-t0)))/(2*self.__kappa))
        return (self._sigma*(1- np.exp(-self.__kappa*(T-t)) )*d)

    def BondPrice_ClosedForm(self,r,t,T):
        A,B = self.getA_B(t,T)
        return A*np.exp(-B*r)

class CIR(InterestRateModel):
    def __init__(self, r0, r_bar, T, dT, sigma, nsims, kappa):
        InterestRateModel.__init__(self, r0, T, dT, sigma, nsims)
        self.__kappa = kappa
        self.__rbar = r_bar

    def EulerDiscretize_R(self):
        np.random.seed(100)
        self.__r = np.zeros([self._nsims, (self._ndiv + 1)])
        self.__r[:, 0] = self._r0
        num = 1
        while (num <= self._ndiv):
            z = np.random.normal(0, 1, self._nsims)
            self.__r[:, num] = self.__r[:, (num - 1)] + self.__kappa * (
            self.__rbar - self.__r[:, (num - 1)]) * self._dT + self._sigma * np.sqrt(np.abs(self.__r[:,(num-1)])*self._dT) * z
            num = num + 1

        return self.__r

    def getA_B(self, t, T):
        h1 = np.sqrt(self.__kappa**2 + 2*(self._sigma**2))
        h2 = (self.__kappa+ h1)/2
        h3 = 2*self.__kappa*self.__rbar/(self._sigma**2)
        A = np.power(h1*np.exp(h2*(T-t))/(h2*(np.exp(h1*(T-t))-1)+h1),h3)
        B = (np.exp(h1*(T-t))-1)/(h2*(np.exp(h1*(T-t))-1)+h1)
        return A, B

    def getsigma_p(self, t0, t, T):
        d = np.sqrt((1 - np.exp(-2 * self.__kappa * (T - t0))) / (2 * self.__kappa))
        return (self._sigma * (1 - np.exp(-self.__kappa * (T - t))) * d)

    def BondPrice_ClosedForm(self, r, t, T):
        A, B = self.getA_B(t, T)
        return A * np.exp(-B * r)

    def getChiSq1(self,t,T,S,K):
        theta = np.sqrt(self.__kappa**2 + 2*(self._sigma**2))
        phi = 2*theta/(self._sigma**2*(np.exp(theta*(T-t))-1))
        psi = (self.__kappa + theta)/(self._sigma**2)
        A,B = self.getA_B(T,S)
        r_star = np.log(A/K)/ B
        x = 2*r_star*(phi+ psi + B)
        p = 4*self.__kappa*self.__rbar/(self._sigma**2)
        q= 2*(phi**2)*self._r0*(np.exp(theta*(T-t)))/(phi+psi+B)
        return ncx2.cdf(x,p,q)

    def getChiSq2(self, t, T,S, K):
        theta = np.sqrt(self.__kappa ** 2 + 2 * (self._sigma ** 2))
        phi = 2 * theta / (self._sigma ** 2 * (np.exp(theta * (T - t)) - 1))
        psi = (self.__kappa + theta) / (self._sigma ** 2)
        A, B = self.getA_B(T, S)
        r_star = np.log(A / K) / B
        x = 2 * r_star*(phi + psi)
        p = 4 * self.__kappa * self.__rbar / (self._sigma ** 2)
        q = 2 * (phi ** 2) * self._r0 * (np.exp(theta * (T - t))) / (phi + psi )
        return ncx2.cdf(x, p, q)


class G2pp(InterestRateModel):

    def __init__(self, r0, x0, y0, T, dT, a, b, sigma, eta, rho, nsims):
        InterestRateModel.__init__(self, r0,T, dT, sigma,nsims)
        # self.__phi = phi
        self.__rho = rho

        self.__x0 = x0
        self.__a = a

        self.__y0 = y0
        self.__b = b
        self.__eta = eta

    def EulerDiscretize_R(self):
        self.__r = np.zeros([self._nsims, (self._ndiv + 1)])
        self.__x = np.zeros([self._nsims, (self._ndiv + 1)])
        self.__y = np.zeros([self._nsims, (self._ndiv + 1)])
        self.__phi = np.zeros([self._nsims, (self._ndiv + 1)])


        self.__r[:,0] = self._r0
        self.__x[:,0] = self.__x0
        self.__y[:,0] = self.__y0
        self.__phi[:, 0] = self._r0

        num = 1
        while (num <= self._ndiv):
            z1 = np.random.normal(0, 1, self._nsims)
            z = np.random.normal(0, 1, self._nsims)
            z2 = self.__rho*z1 + np.sqrt(1-self.__rho**2)*z
            self.__x[:, num] = self.__x[:,(num-1)] - self.__a*self.__x[:,(num-1)]*self._dT + self._sigma*z1*np.sqrt(self._dT)
            self.__y[:, num] = self.__y[:,(num-1)] - self.__b*self.__y[:,(num-1)]*self._dT + self.__eta*z2*np.sqrt(self._dT)
            self.__phi[:, num] =self.phi(num)

            num = num + 1

        self.__r = self.__x + self.__y + self._r0
        # return self.__r,self.__x,self.__y
        return self.__r

    def V(self,t,T):
        dt = T-t
        exp1 = ((self._sigma / self.__a) ** 2) * (
        dt + 2 * np.exp(-self.__a * dt) / self.__a - np.exp(-2 * self.__a * dt) / (2 * self.__a) - 3 / (2 * self.__a))
        exp2 = ((self.__eta / self.__b) ** 2) * (
        dt + 2 * np.exp(-self.__b * dt) / self.__b - np.exp(-2 * self.__b * dt) / (2 * self.__b) - 3 / (2 * self.__b))
        exp3 = (2*self.__rho*self._sigma*self.__eta/(self.__a*self.__b)) * (
            dt + (np.exp(-self.__a * dt) - 1) / self.__a +  (np.exp(-self.__b * dt) - 1) / self.__b - (np.exp(-(self.__a+ self.__b) * dt) - 1)/(self.__a+ self.__b))
        return (exp1 + exp2 + exp3)

    def BigSigma(self,T,S):
        exp1 = self._sigma ** 2 / (2*self.__a ** 3) * ((1 - np.exp(-self.__a * (S - T))) ** 2) * (
        1 - np.exp(-2 * self.__a * (T)))

        exp2 = self.__eta ** 2 / (2*self.__b ** 3) * ((1 - np.exp(-self.__b * (S - T))) ** 2) * (
        1 - np.exp(-2 * self.__b * (T)))

        exp3 = 2*self.__rho*self._sigma*self.__eta / (self.__a * self.__b * (self.__a + self.__b)) * (1 - np.exp(-self.__a * (S - T))) * (1 - np.exp(-self.__b * (S - T))) * (
        1 - np.exp(-(self.__a + self.__b) * (S - T)))
        return np.sqrt((exp1 + exp2 + exp3))

    def phi(self,T):
        return self._r0+self._sigma**2/(2*self.__a**2)*(1-math.exp(-self.__a*T))**2+self.__eta**2/(2*self.__b**2)*(1-math.exp(-self.__b*T))**2+self.__rho*self._sigma*self.__eta/(self.__a*self.__b)*(1-math.exp(-self.__b*T))*(1-math.exp(-self.__a*T))


    def BondPrice_ClosedForm(self,r,t, T):
        exp1 = -self.phi(T)*(T-t)
        exp2 = - (1 - np.exp(-self.__a * (T-t))) * self.__x0 / self.__a
        exp3 = - (1 - np.exp(-self.__b * (T-t))) * self.__y0 / self.__b
        exp4 = self.V(t,T)/2
        return np.exp(exp1 + exp2 + exp3 + exp4)
        #
        # exp1 = -self.__phi * (T)
        # exp2 = - (1 - np.exp(-self.__a * T)) * self.__x0 / self.__a
        # exp3 = - (1 - np.exp(-self.__b * T)) * self.__y0 / self.__b
        # exp4 = self.V(0, T) / 2
        # return np.exp(exp1 + exp2 + exp3 + exp4)
    def PlotIR(self):

        self.__ndiv = self.__r.shape[1]
        self.__nsims = self.__r.shape[0]
        print("paths:",self.__nsims,"div:",self.__ndiv)
        print(self.__r)

        fig = plt.figure()

        ax = fig.add_subplot(111)
        ax.set_title('colorMap')
        plt.imshow(np.transpose(self.__r))
        plt.colorbar(orientation='vertical')
        plt.show()
