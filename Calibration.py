import QuantLib as ql
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import xlrd
import numpy as np
from scipy import optimize
from collections import namedtuple
import math

import QuantLib as ql
from CDS import *
from DiscountCurve import *
from MarketData import *


class Calibration(object):
    def __init__(self):
        self.IRdata=self.GetIRdata()
        self.Sdata=self.GetSwaptiondata()

    def SetForwardRate(self,rate):
        self.forwardrate=rate

    def GetIRdata(self):
        wb = xlrd.open_workbook('tmp.xlsx')
        table = wb.sheets()[0]
        list = table.col_values(1)
        list = list[2:]
        return list

    def GetSwaptiondata(self):
        CalibrationData = namedtuple("CalibrationData","start, length, volatility")
        data = [CalibrationData(5, 5, 0.009073),
                CalibrationData(10, 10,0.008322 ),
                CalibrationData(15, 15, 0.006432),
                CalibrationData(20, 20, 0.005982),
                CalibrationData(1, 1,0.006001 ),
                CalibrationData(1, 3, 0.007631),
                CalibrationData(25, 4, 0.006769),
                CalibrationData(5, 25, 0.008167),
                CalibrationData(25, 5, 0.00661),
                CalibrationData(20, 10, 0.00596),
                CalibrationData(10, 20, 0.007111),
                CalibrationData(15, 10, 0.007315)]
        return data

    def CIR(self):
        gamma = np.mean(self.IRdata)
        Rt = [x - gamma for x in self.IRdata]
        n = len(Rt) - 1

        def Residual(phi):
            RSS = []
            for i in range(len(Rt) - 1):
                RSS.append((Rt[i + 1] - phi * Rt[i]) ** 2 / (Rt[i + 1] + gamma))
            return sum(RSS)

        a = optimize.brent(Residual)
        sigma_ = np.sqrt(Residual(a) / n)
        kappa = -np.log(a)
        sigma = np.sqrt(2 * kappa * sigma_ ** 2 / (1 - np.exp(-2 * kappa)))
        return [gamma, sigma, kappa,self.IRdata[len(self.IRdata)-1]]

    def G2pp(self,date):
        today = date
        settlement = date+ql.Period(2,ql.Days)
        ql.Settings.instance().evaluationDate = today;
        term_structure = ql.YieldTermStructureHandle(
            ql.FlatForward(settlement, 0.004875825, ql.Actual365Fixed())
        )

        index =ql.USDLibor(ql.Period(3, ql.Months),term_structure)

        model = ql.G2(term_structure);
        #engine = ql.TreeSwaptionEngine(model, 25)
        engine = ql.G2SwaptionEngine(model, 10, 400)

        swaptions = self.create_swaption_helpers(self.Sdata, index, term_structure, engine)

        optimization_method = ql.LevenbergMarquardt(1.0e-8, 1.0e-8, 1.0e-8)
        end_criteria = ql.EndCriteria(1000, 100, 1e-6, 1e-8, 1e-8)
        model.calibrate(swaptions, optimization_method, end_criteria)
        return model.params()

    def create_swaption_helpers(self, data, index, term_structure, engine):
        swaptions = []
        fixed_leg_tenor = ql.Period(1, ql.Years)
        fixed_leg_daycounter = ql.Actual360()
        floating_leg_daycounter = ql.Actual360()
        for d in data:
            vol_handle = ql.QuoteHandle(ql.SimpleQuote(d.volatility))
            helper = ql.SwaptionHelper(ql.Period(d.start, ql.Years),
                                       ql.Period(d.length, ql.Years),
                                       vol_handle,
                                       index,
                                       fixed_leg_tenor,
                                       fixed_leg_daycounter,
                                       floating_leg_daycounter,
                                       term_structure
                                       )
            helper.setPricingEngine(engine)
            swaptions.append(helper)
        return swaptions


class CDSCalibration(object):
    """Defines a class for the calibration of CDS objects.

    Input a CDS class, a DiscountCurve object, and a MarketData object,
    along with an initial guess for parameters, and the Calibration object
    can then be used to calibrate the CDS intensity process to market data.
    """

    def __init__(self, Today, DiscountCurve=FlatDiscountCurve(r=0.0), \
                 MarketData=None, Method="nm", CDS=None, Process=None, \
                 Guess=None):
        super(CDSCalibration, self).__init__()
        self.DiscountCurve = DiscountCurve
        self.R = 0.4
        self.CDS = CDS
        self.Process = Process
        self.Guess = Guess
        self.Method = Method
        self.Today=Today
        if MarketData is not None:
            self.SetMarketData(MarketData)

    def SetMarketData(self, MarketData):
        """docstring for SetMarketData"""
        self.MarketData = MarketData

    def ObjectiveFunction(self, gamma):
        """Calculates the error in estimation for use in our calibration
        routines.

        Currently we use the L^2 norm."""

        sum = 0
        for t, market_spread in self.MarketData.Data():

            CDS = self.CDS(self.Today, DiscountCurve=self.DiscountCurve(self.Today),
                           maturity=int(t))
            model_spread = CDS.ParSpread(gamma)
            sum += ((model_spread - market_spread)) ** 2
        # print(sum)
        return sum

    def Calibrate(self, method='nm'):
        """Performs the calibration and returns the optimal parameters.

        The built in Optimise method in SciPy uses Nelder-Mead optimisation."""
        # methods = {'nm': optimize.fmin,
        #            'powell': optimize.fmin_powell,
        #            'cg': optimize.fmin_cg,
        #            'bfgs': optimize.fmin_bfgs
        #            }
        #
        # if method == None:
        #     optimise = methods[self.Method]
        # else:
        #     optimise = methods[method]
        res =minimize(self.ObjectiveFunction,self.Guess, method='nelder-mead',tol=0.1)
        output = res.x
        # output=optimise(self.ObjectiveFunction,
        #                   self.Guess,
        #                   disp=0
        #                   )

        # if(2*output[0]*output[1]>=output[2]):
        self.calibrated_gamma = output
        return output
        # else:
        #     print('Try different guess')

    def RMSE(self):
        """Returns the RMSE for the calibrated parameters."""
        N = len(self.MarketData.Tenors())
        return sqrt(self.ObjectiveFunction(self.calibrated_gamma) / N)

class Default(object):
    def __init__(self,nsims,N):
        self.calendar = ql.TARGET()
        self.todays_date = self.calendar.adjust(ql.Date(12, 12, 2016))
        self.__nsims=nsims
        self.__N=N
        ql.Settings.instance().evaluation_date =self.todays_date
        self.__c = np.zeros([self.__nsims, self.__N + 1])

    def setENV(self):
# dummy curve
        self.ts_curve =ql.YieldTermStructureHandle(ql.FlatForward(self.todays_date, 0.01, ql.Actual365Fixed()))
        self.recovery_rate = 0.5
        self.quoted_spreads = [0.0046, 0.0053, 0.0071, 0.0094, 0.0116, 0.0137, 0.0158, 0.0173, 0.0176, 0.0179, 0.0178]
        self.tenors = [ql.Period(i, ql.Months) for i in [6, 12, 24, 36, 48, 60, 84, 120, 180, 240, 360]]
        self.maturities = [self.calendar.adjust(self.todays_date + self.tenors[i], ql.Following) for i in range(4)]

    def SurvivalFun(self):
        instruments = []
        for i in range(len(self.tenors)):
            helper = ql.SpreadCdsHelper(
            ql.QuoteHandle(ql.SimpleQuote(self.quoted_spreads[i])), self.tenors[i], 0, self.calendar, ql.Quarterly,
            ql.Following, ql.DateGeneration.TwentiethIMM, ql.Actual365Fixed(), self.recovery_rate, self.ts_curve)

            instruments.append(helper)

        hazard_rate_structure=ql.PiecewiseFlatHazardRate(self.todays_date,instruments,ql.Actual365Fixed())

        for i in range(self.__N):
            self.__c[:,i]=hazard_rate_structure.survivalProbability(self.todays_date + ql.Period(i,ql.Months))

        return self.__c