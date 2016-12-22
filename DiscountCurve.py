import numpy as np
from TermStructure import *

class DiscountCurveBase(object):
    """Defines the interface for a DiscountCurve.  The key method is DF(),
    which returns the discount factor from t_1 to t_2.

    To use, subclass and implement the DF() method."""

    def __init__(self):
        super(DiscountCurveBase, self).__init__()

    def DF(self, t1, t2):
        """Subclass and implement.

        Returns the discount factor from t1 to t2."""
        raise NotImplementedError


# ------------------------------------------------------------------------------

class FlatDiscountCurve(DiscountCurveBase):
    """Implements a DiscountCurve where the yield curve is constant r."""

    def __init__(self, r):
        super(FlatDiscountCurve, self).__init__()
        self.r = r

    def DF(self, t1, t2):
        return (1 + self.r) ** -(t2 - t1)

class OISDisountCurve(DiscountCurveBase):

    def __init__(self,date):
        super(OISDisountCurve,self).__init__()
        self.ts=TermStructure(date).OIS()
    def DF(self,t):
        return self.ts.discount(t)




if __name__ == '__main__':
    x = FlatDiscountCurve(r=0.06)
    print
    x.DF(0, 10)