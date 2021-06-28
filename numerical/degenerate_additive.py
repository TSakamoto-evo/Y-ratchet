import numerical
import matplotlib.pyplot as plt
import numpy as np


N = 10000
L1 = 1500
L2 = 1500
H = 0
u = 0.00001
s = 0.0025
c = 0.000001

seps = 100001
max_t = 100 * (seps - 1)

list = numerical.ret_degenerate_semi_process(N, s, u, c, L1, L2, H, max_t, seps)

if len(list) >= 1001:
    chosen_index = [int((len(list) - 1) * i / 1000) for i in range(1001)]
    sub_list = list[chosen_index, :]
else:
    sub_list = list

np.savetxt('degenerate_additive.csv', sub_list, delimiter=',')
