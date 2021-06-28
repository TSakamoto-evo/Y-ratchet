import numerical
import matplotlib.pyplot as plt
import numpy as np

N = 10000
L1 = 0
L2 = 2000
u = 0.00001
s = 0.0025
c = 0.0

seps = 1001

fin_list = []

for i in range(seps):
    log_c = -7 + 4 * i / (seps - 1)
    c = 10 ** log_c
    ret_list = list(numerical.ret_rates_semi(N, s, u, c, L1, L2))

    fin_list.append([ret_list[0], ret_list[1], c])

np.savetxt('no_degeneration_additive.csv', fin_list, delimiter=',')
