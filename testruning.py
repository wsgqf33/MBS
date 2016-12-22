import InterestrateModel as IRM
import PrepaymentModel as PPM
import numpy as np
import timeit
from numpy import*
import matplotlib.pyplot as plt
import QuantLib as ql
# from tabulate import tabulate
import random
from Calibration import *
from MarketData import *
from DiscountCurve import *
from CDS import *
from HazardRateSim import *
import math

# Q=array([[1,3,4,5,6,7,8,9],[1,2,4,3.3,2,3,4,9]])
# t=linspace(0,100,100)
# a=[]
# b=[]
# c=[]
# for i in range(101):
#     if i<=30:
#         a.append(0.002*i)
#         b.append(0.003 * i)
#         c.append(0.004 * i)
#     else:
#         a.append(0.06)
#         b.append(0.09)
#         c.append(0.12)
#
# plt.plot(a,label='100 PSA')
# plt.plot(b,label='150 PSA')
# plt.plot(c,label='200 PSA')
# # plt.suptitle('PSA Prepayment Model')
# plt.xlabel('Months')
# plt.ylabel('CPR')
# plt.savefig("PSA.eps")
# plt.legend()
# plt.show()
aa=np.array([1.1524968465880101, 1.1523798387811814, 1.1523944983050671, 1.1524606127562742, 1.1524116904741752,
             1.1524523929554222, 1.1524870534164149, 1.1525011726714554, 1.1525393068198284, 1.1524600107049259,
             1.1524872231948684, 1.1523576702481164, 1.1524632897106437, 1.1523868595399747, 1.1524093555625072,
             1.1523780732562217, 1.1524691000496126, 1.1524913479813768, 1.1524354223919142, 1.1524495527347196])
a=np.array([1.1524710133416749, 1.1524402294173564, 1.1524533296018602, 1.1524837587180927, 1.1524237802880466, 1.1524569035464334, 1.1524555567679391, 1.1524172251172431, 1.1524366445745755, 1.1524442820119003, 1.1524857857927608, 1.1524099947280164, 1.152464762344017, 1.1524899524257879, 1.1524197064180397, 1.1525471477119882, 1.15246243683476, 1.1524534955363632, 1.1524975270103353, 1.1524668953726422])
b=np.array([1.1524654820657447, 1.1524887459202877, 1.1524698997880527, 1.152482205501689, 1.152438002216738, 1.1524875905610144, 1.1524724160465707, 1.1524168992091424, 1.1524554264388955, 1.1524425943889187, 1.1524267702857702, 1.1524802982727282, 1.1524762363730672, 1.1524768466970756, 1.1524515275352929, 1.152470116415558, 1.1524540742224059, 1.1524671770718908, 1.1524338807661785, 1.152466705172724])
c=np.array([1.1524971304704361, 1.1524929973594391, 1.1524833082651849, 1.1524790902249307, 1.1524438391594993, 1.1524596570786152, 1.1524676075144735, 1.1524185635234507, 1.1524840758983783, 1.152468015799458, 1.1524729785642889, 1.1524515896020371, 1.1524590207217673, 1.1524426870872675, 1.1524549389638774, 1.152450742708617, 1.1524619681420141, 1.1524712098338208, 1.1524649909817757, 1.1524753432033252])
d=np.array([1.1524527831807698, 1.1524578406019996, 1.152468771982269, 1.1524604789217199, 1.1524619073748954, 1.1524556330050335, 1.1524638805352168, 1.1524463878311058, 1.1524714476631259, 1.152452632252474, 1.1524431696809256, 1.1524556329916462, 1.1524546205467705, 1.1524593105084682, 1.1524704601169908, 1.1524775866123083, 1.1524519649951299, 1.1524571479599353, 1.1524599010721472, 1.1524539202857353])


pp=[log10(aa.std()),log10(a.std()),log10(b.std()),log10(c.std()),log10(d.std())]
tt=[log10(100),log10(250),log10(500),log10(1000),log10(2000)]
plt.plot(tt,pp)
print(pp)
plt.suptitle('Monte Carlo Error Analysis')
plt.xlabel('Log(Npaths)')
plt.ylabel('Log(Error)')
plt.savefig('MCerror.eps')
plt.legend()
plt.show()