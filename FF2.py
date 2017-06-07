from __future__ import division
import numpy as np

def F2s(Q1,Q2,CP12):
    F2 = 5/7 + 1/2 * CP12 * (Q1/Q2 + Q2/Q1) + 2/7 * CP12**2
    return F2

def G2s(Q1,Q2,CP12):
    G2 = 3/7 + 1/2 * CP12 * (Q1/Q2 + Q2/Q1) + 4/7 * CP12**2
    return G2

def F3U(Q1,Q2,Q3,CP12,CP23,CP13):
    Q12 = np.sqrt(Q1**2 + Q2**2 + 2*Q1*Q2*CP12)
    Q23 = np.sqrt(Q2**2 + Q3**2 + 2*Q2*Q3*CP23)
    Q13 = np.sqrt(Q1**2 + Q2**2 + 2*Q1*Q2*CP13)

    alpha1 = (1/2) * (1 + Q2/Q1 * CP12 + Q3/Q1 * CP13 + (Q1*Q2*CP12 + Q1*Q3*CP13/(Q23**2)) + 1)
    alpha2 = (1/2) * (1 + (Q3*Q1*CP13 + Q3*Q2*CP23/(Q12**2)) + 1 + Q1/Q3 * CP13 + Q2/Q3 * CP23)

    beta1 = (Q1**2 + Q2**2 + Q3**2 + 2*Q1*Q2*CP12 + 2*Q1*Q3*CP13 + 2*Q2*Q3*CP23)/(2*Q23**2) * (Q2/Q1 * CP12 + Q3/Q1 * CP13)
    beta2 = (Q1**2 + Q2**2 + Q3**2 + 2*Q1*Q2*CP12 + 2*Q1*Q3*CP13 + 2*Q2*Q3*CP23)/(2*Q12**2) * (Q1/Q3 * CP13 + Q2/Q3 * CP23)

    F3 = 1/18 * (7*alpha1*F2s(Q2,Q3,CP23) + 2*beta1*G2s(Q2,Q3,CP23)) + G2s(Q1,Q2,CP12)/18 * (7*alpha2 + 2*beta2)
    return F3

def F3s(Q1, Q2, Q3, CP12, CP23, CP13):
    F3 = 1/6 * (F3U(Q1,Q2,Q3,CP12,CP23,CP13) + F3U(Q1,Q3,Q2,CP13,CP23,CP12) + F3U(Q2,Q1,Q3,CP12,CP13,CP23) + F3U(Q3,Q1,Q2,CP13,CP12,CP23) + F3U(Q2,Q3,Q1,CP23,CP13,CP12) + F3U(Q3,Q2,Q1,CP23,CP12,CP13))
    return F3
