import numerical
import matplotlib.pyplot as plt
import numpy as np

N = 10000
L1 = 1500
L2 = 1200
u = 0.00001
s = 0.007
c = 0.0

seps = 1001

fin_list = []

for i in range(seps):
    log_c = -7 + 4 * i / (seps - 1)
    c = 10 ** log_c
    ret_list = list(numerical.ret_rates(N, s, u, c, L1, L2))

    fin_list.append([ret_list[0], ret_list[1], c])

np.savetxt('no_degeneration_dominant.csv', fin_list, delimiter=',')
