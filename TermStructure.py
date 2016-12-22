import QuantLib as ql
import matplotlib.pyplot as plt
import numpy as np
import datetime
'''
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
'''
'''
depo_maturities = [ql.Period(6,ql.Months), ql.Period(12, ql.Months)]
depo_rates = [5.25, 5.5]

# Bond rates
bond_maturities = [ql.Period(6*i, ql.Months) for i in range(3,21)]
bond_rates = [5.75, 6.0, 6.25, 6.5, 6.75, 6.80, 7.00, 7.1, 7.15,
              7.2, 7.3, 7.35, 7.4, 7.5, 7.6, 7.6, 7.7, 7.8]

calc_date = ql.Date(15, 1, 2015)
ql.Settings.instance().evaluationDate = calc_date

calendar = ql.UnitedStates()
bussiness_convention = ql.Unadjusted
day_count = ql.Thirty360()
end_of_month = True
settlement_days = 0
face_amount = 100
coupon_frequency = ql.Period(ql.Semiannual)
settlement_days = 0

depo_helpers = [ql.DepositRateHelper(ql.QuoteHandle(ql.SimpleQuote(r/100.0)),
                                     m,
                                     settlement_days,
                                     calendar,
                                     bussiness_convention,
                                     end_of_month,
                                     day_count )
                for r, m in zip(depo_rates, depo_maturities)]

bond_helpers = []
for r, m in zip(bond_rates, bond_maturities):
    termination_date = calc_date + m
    schedule = ql.Schedule(calc_date,
                   termination_date,
                   coupon_frequency,
                   calendar,
                   bussiness_convention,
                   bussiness_convention,
                   ql.DateGeneration.Backward,
                   end_of_month)

    helper = ql.FixedRateBondHelper(ql.QuoteHandle(ql.SimpleQuote(face_amount)),
                                        settlement_days,
                                        face_amount,
                                        schedule,
                                        [r/100.0],
                                        day_count,
                                        bussiness_convention,
                                        )
    bond_helpers.append(helper)
rate_helpers = depo_helpers + bond_helpers

#yieldcurve = ql.PiecewiseLogCubicDiscount(calc_date,rate_helpers,day_count)
ts_curve =ql.YieldTermStructureHandle(ql.PiecewiseCubicZero ( calc_date, rate_helpers,ql.Actual365Fixed()))

            self.Future=[(99.36,Date(31,1,2017)),(99.34
                        ,Date(28,2,2017)),(99.31,Date(31,3,2017)),(99.29,Date(28,4,2017)),(99.24,Date(31,5,2017))
                        ,(99.17,Date(30,6,2017)),(99.10,Date(31,7,2017)),(99.06,Date(31,8,2017)),(99.01
                        ,Date(29,9,2017)),(98.93,Date(31,10,2017)),(98.89,Date(30,11,2017)),(98.82,Date(29,12,2017)),(98.76
                        ,Date(31,1,2018)),(98.73,Date(28,2,2018)),(98.70,Date(29,3,2018)),(98.63,Date(30,4,2018)),(98.59
                        ,Date(31,5,2018)),(98.55,Date(29,6,2018)),(98.50,Date(31,7,2018)),(98.46,Date(31,8,2018)),(98.49
                        ,Date(28,9,2018)),(98.38,Date(31,10,2018)),(98.35,Date(30,11,2018)),(98.31,Date(31,12,2018))]
'''
import math




from QuantLib import *

class TermStructure:

    def __init__(self,date):
        self.today = date
        Settings.instance().evaluationDate = self.today


    def DataReader(self,file):
        # need implement of csv
        if(file==None):
            self.DepositRate=[(0.41, 0), (0.41, 1), (0.41, 2)]
            self.ShortTermOIS=[(0.6510, (1, Weeks)), (0.6201, (2, Weeks)),
                (0.626, (3, Weeks))]
            self.Future=[(99.31,Date(15,3,2017)),(99.17,Date(21,6,2017)),(99.01
                        ,Date(20,9,2017)),(98.82,Date(20,12,2017)),(98.70,Date(21,3,2018)),(98.55,Date(20,6,2018)),(98.49
                        ,Date(19,9,2018)),(98.31,Date(19,12,2018))]
            self.LongTermOIS=[(1.1536, (2,Years)),
                            (1.401, (3,Years)), (1.5885, (4,Years)),
                            (1.7130, (5,Years)), (1.84, (6,Years)),
                            (1.924, (7,Years)),  (2.105, (10,Years)),
                            (2.171, (12,Years)),(2.238, (15,Years)),
                            (2.287, (20,Years)),(2.3015, (25,Years)),
                            (2.295, (30,Years))]

    def OIS(self,file=None):
        self.DataReader(file)
        helpers = [ DepositRateHelper(QuoteHandle(SimpleQuote(rate/100)),Period(1,Days), fixingDays,TARGET(), Following, False, Actual360())
                    for rate, fixingDays in self.DepositRate ]

        eonia = Eonia()
        helpers += [ OISRateHelper(2, Period(*tenor),QuoteHandle(SimpleQuote(rate/100)), eonia)
                     for rate, tenor in self.ShortTermOIS ]

        '''helpers +=[ FuturesRateHelper(QuoteHandle(SimpleQuote(quote)),
                                     date, 1,
                                     TARGET(), ModifiedFollowing,
                                     True, Actual360(),
                                     QuoteHandle(SimpleQuote(0.0)))
                    for quote , date in self.Future]'''

        helpers += [ OISRateHelper(2, Period(*tenor),QuoteHandle(SimpleQuote(rate/100)), eonia)
                     for rate, tenor in self.LongTermOIS ]
        eonia_curve_c = PiecewiseLogCubicDiscount(0, TARGET(),
                                          helpers, Actual360())
        eonia_curve_c.enableExtrapolation()
        return eonia_curve_c


