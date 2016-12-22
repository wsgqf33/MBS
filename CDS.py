from numpy import arange, arctanh
from scipy.misc import derivative
from scipy.integrate import quad
from math import exp, sqrt, cosh, sinh, log, tanh
import QuantLib as ql
# from __future__ import division

#from CreditDerivativeCSVReader import *
from DiscountCurve import *


# ------------------------------------------------------------------------------

def coth(x):
    """docstring for coth"""
    return 1 / tanh(x)


# ------------------------------------------------------------------------------

class CreditDefaultSwap(object):
    """Abstract base class for CDS object.  Implements the PaymentDates(),
    ParSpread(), and ContinuousParSpead() methods.

    To use this, subclass and implement the SurivivalProbability() method."""

    def __init__(self, date,maturity=None, DiscountCurve=None, spread=None):
        super(CreditDefaultSwap, self).__init__()
        self.maturity = maturity
        self.spread = 0.1
        self.DiscountCurve = DiscountCurve
        self.R = 0.4
        self.Date=date

    def PaymentDates(self):
        """Returns array containing quarterly payment dates."""
        start_date = self.Date
        i=self.maturity
        end_date = self.Date+i*ql.Period(ql.Annual)-ql.Period(3,ql.Months)
        period = ql.Period(3, ql.Months)
        calendar = ql.UnitedStates()
        buss_convention = ql.ModifiedFollowing
        rule = ql.DateGeneration.Forward
        end_of_month = False

        schedule = ql.Schedule(start_date, end_date, period,
                               calendar, buss_convention, buss_convention,
                               rule, end_of_month)

        return schedule

    def ParSpread(self, parameters):
        """Returns the CDS par spreads"""
        dates = list(self.PaymentDates())
        #print(dates)
        prot_leg = 0
        loss_leg = 0
        # Calculate par spread using formula provided in thesis.
        for date in dates:
            t_start = date
            t_end = date + ql.Period(3, ql.Months)
            prot_leg += self.DiscountCurve.DF(t_end) * \
                       (self.SurvivalProbability(parameters, t_start) - \
                        self.SurvivalProbability(parameters, t_end))
            loss_leg += self.DiscountCurve.DF(t_end) * \
                        (self.SurvivalProbability(parameters,t_start)+self.SurvivalProbability(parameters, t_end)) * 0.25

        par_spread = (1 - self.R) * prot_leg / loss_leg/2
        # Par spread is expressed in basis points (bps)
        ps = par_spread * 10000
        #print(ps)
        return ps

    def SurvivalProbability(self):
        """Returns P(tau > t) - the probability that the entity survives past
        time t"""
        abstract

    '''def ContinuousParSpread(self, parameters):
        surv_prob = lambda t: self.SurvivalProbability(parameters, t)
        numerator = lambda t: self.DiscountCurve.DF(0, t) * \
                              derivative(surv_prob, t, dx=0.0001)
        denominator = lambda t: self.DiscountCurve.DF(0, t) * \
                                self.SurvivalProbability(parameters, t)

        cts_ps = (1 - self.R) * -quad(numerator, 0, self.maturity)[0] \
                 / quad(denominator, 0, self.maturity)[0] * 10000

        return cts_ps'''


# ------------------------------------------------------------------------------

class HPCreditDefaultSwap(CreditDefaultSwap):
    """Intensity is a constant gamma > 0."""

    def __init__(self, date,maturity=None, DiscountCurve=None, spread=0.100):
        super(HPCreditDefaultSwap, self).__init__(date,maturity,
                                                  DiscountCurve,
                                                  spread,
                                                  )

    def SurvivalProbability(self, gamma, t):
        return exp(- gamma * t)

    def ParSpread(self, gamma):
        """Returns the par spread for a CDS
        where the intensity of default is a Poisson process
        with parameter gamma"""

        par_spread = (1 - self.R) * gamma
        return par_spread * 10000





# ------------------------------------------------------------------------------

class CIRCreditDefaultSwap(CreditDefaultSwap):
    """docstring for OUCreditDefaultSwap"""

    def __init__(self, date,maturity=None, DiscountCurve=None, spread=0.100):
        super(CIRCreditDefaultSwap, self).__init__(date,maturity, DiscountCurve, \
                                                   spread=0.100)

    def SurvivalProbability(self, parameters, date):
        """Solves \phi_{CIR}(i, t; \kappa, \nu, \zeta, \gamma_0)
        Parameters: [ kappa, nu, zeta, gamma]"""
        assert (len(parameters) == 4)
        kappa, nu, sigma, lamb = parameters
        t=(date-self.Date)/365
        if t == 0.0:
            probability = 1
        else:
            gamma = sqrt(kappa ** 2 + 2 * sigma ** 2)
            A=((2*gamma*exp((kappa+gamma)*t/2))/(2*gamma+(kappa+gamma)*(exp(gamma*t)-1)))**(2*kappa*nu/sigma/sigma)
            B=(2*(exp(gamma*t)-1))/(2*gamma+(kappa+gamma)*(exp(gamma*t)-1))
            probability=A*exp(-B*lamb)
            # probability = (exp(kappa ** 2 * nu * t / vega ** 2) * \
            #                    exp(-2 * lamb / (kappa + gamma * coth(gamma * t / 2)))) \
            #                   / (coth(gamma * t / 2) \
            #                      + kappa * sinh(gamma * t / 2) / gamma) \
            #                     ** (2 * kappa * nu / vega ** 2)
        return probability


# ------------------------------------------------------------------------------

