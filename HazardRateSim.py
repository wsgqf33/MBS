import numpy as np
def HazardRateSim(para,N,T,delta):
    [kappa, mu, sigma, lamb0] = para
    div = int(T / delta)

    Q = np.zeros([N, (div + 1)])
    lamb = np.zeros([N, div + 1])
    c = sigma ** 2 * (1 - np.exp(-kappa * delta)) / 4 / kappa
    d = 4 * mu * kappa / sigma ** 2
    lamb[:, 0] = lamb0
    Q[:, 0] = 1
    num = 1
    while (num <= div):
        x = np.random.normal(0, 1, N)
        lamb[:, num] = lamb[:, num - 1] + kappa * (mu - lamb[:, num - 1]) * delta + sigma * np.sqrt(
            np.abs(lamb[:, num - 1] * delta)) * x
        Q[:, num] = np.exp(-lamb[:, (num - 1)] * delta) * Q[:, (num - 1)]
        num = num + 1
    return Q