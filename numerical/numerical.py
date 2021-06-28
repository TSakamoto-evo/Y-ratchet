import numpy as np
import scipy.optimize as opt
import scipy.special as sc
from scipy.integrate import odeint
import warnings
import scipy.integrate.odepack

def for_ne(x, s, u, c, L1, L2, N):
    U = 2*x*u
    C = 2*x*c

    H = (L2 * (2*u)/(2*u+c) *
            sc.hyp1f1(2*U+1, 2*U+C+1, -(U+C)) / sc.hyp1f1(2*U, 2*U+C, -(U+C)))
    Ne = N * np.exp(-((c+u)*H + u*L1) / s)

    return(Ne - x)

def for_v(v, N, s, u):
    sigma = s/u
    left = sigma * np.log(N * u * sigma**(3/2))
    right = ((1 - v/2 * ((1 - np.log(v)) ** 2 + 1)) -
        sigma * np.log(np.sqrt(v**3 / (1-v)) * (1 - np.log(v)) /
        (1 - v + v * np.log(v) + 5 * sigma / 6)))
    return(left - right)

def ret_rates(N, s, u, c, L1, L2):
    if L2 != 0:
        Ne = opt.newton(for_ne, N, args=(s, u, c, L1, L2, N))

        U = 2*Ne*u
        C = 2*Ne*c
        H = (L2 * (2*u)/(2*u+c) *
                sc.hyp1f1(2*U+1, 2*U+C+1, -(U+C)) / sc.hyp1f1(2*U, 2*U+C, -(U+C)))
    else:
        Ne = N * np.exp(-L1 * u / s)
        H = 0

    if Ne * s > 1:
        alpha = np.sqrt(s / ((c+u) * H + u * L1))
        inv_T = alpha * s / np.exp(alpha * s * Ne) * np.sqrt(alpha * s * Ne / np.pi)

    else:
        v = opt.newton(for_v, 0.001, args=(N, s, (c+u) * H + u * L1))
        inv_T = ((c+u) * H + u * L1) * v

    return(u / ((c+u) * H + u * L1) * inv_T, (c+u) * H / ((c+u) * H + u * L1) / L2 * inv_T)

def return_diff(var, t, N, s, u, c):
    L1 = var[0]
    L2 = var[1]

    rates = ret_rates(N, s, u, c, L1, L2)
    dL1dt = -L1 * rates[0]
    dL2dt = -L2 * rates[1]
    dtdt = 1
    return [dL1dt, dL2dt, dtdt]

def ret_degenerate_process(N, s, u, c, L1, L2, max_t, seps):
    tmp_seps = seps

    error_list = [0]
    while len(error_list) == 1:
        print(tmp_seps, error_list)
        t_list = np.linspace(0.0, max_t, tmp_seps)
        var_init = [L1, L2, 0]

        with warnings.catch_warnings(record = True) as error_list:
            var_list = odeint(return_diff, var_init, t_list, args=(N, s, u, c))

        if len(error_list) >= 1:
            tmp_seps = int(tmp_seps / 2)

    var_init = [L1, L2, 0]
    t_list = np.linspace(0.0, max_t, tmp_seps)
    var_list = odeint(return_diff, var_init, t_list, args=(N, s, u, c))
    print(tmp_seps)

    return(var_list)

def equ_H1(H, N, s, u, c, L1, L2):
    U = 2*(L2-H)*u + (u+c)*H + L1*u

    alpha = np.sqrt(s/U)
    Ne = N*np.exp(-U/s)
    inv_Tf = alpha * s / np.exp(alpha * s * Ne) * np.sqrt(alpha * s * Ne / np.pi)

    r1 = u/U * inv_Tf
    r2 = (u+c)*H/U/L2 * inv_Tf

    Hf = 2*u*(L2-H) / U * inv_Tf - L2 * r2

    inv_Tb = 2*s*c*H*Ne
    Hb = -inv_Tb

    return( Hf + Hb + H * r2 )

def ret_H1(N, s, u, c, L1, L2):
    if L2 != 0:
        ret = opt.newton(equ_H1, 0.001, args=(N, s, u, c, L1, L2), full_output=True, disp=False)
        if ret[1].converged and ret[0] >= 0 and ret[0] <= L2:
            H = ret[0]
        else:
            H = opt.newton(equ_H1, L2, args=(N, s, u, c, L1, L2))
    else:
        H = 0

    return(H)

def equ_H2(H, N, s, u, c, L1, L2):
    U = 2*(L2-H)*u + (u+c)*H + L1*u
    v = opt.newton(for_v, 0.001, args=(N, s, U))
    inv_Tf = U*v

    r1 = u/U * inv_Tf
    r2 = (u+c)*H/U/L2 * inv_Tf

    Hf = 2*u*(L2-H) / U * inv_Tf - L2 * r2

    inv_Tb = (np.exp(1) / v - 1) * c * H / (1 - np.log(v))
    Hb = -inv_Tb

    return( Hf + Hb + H * r2 )

def ret_H2(N, s, u, c, L1, L2):
    if L2 != 0:
        ret = opt.newton(equ_H2, 0.001, args=(N, s, u, c, L1, L2), full_output=True, disp=False)
        if ret[1].converged and ret[0] >= 0 and ret[0] <= L2:
            H = ret[0]
        else:
            H = opt.newton(equ_H2, L2, args=(N, s, u, c, L1, L2))
    else:
        H = 0
    return(H)

def ret_rates_semi(N, s, u, c, L1, L2):
    H = ret_H1(N, s, u, c, L1, L2)
    U = 2*(L2-H)*u + (u+c)*H + L1*u

    if N*s*np.exp(-U/s) > 1:
        H = ret_H1(N, s, u, c, L1, L2)

        alpha = np.sqrt(s/U)
        Ne = N*np.exp(-U/s)
        inv_Tf = alpha * s / np.exp(alpha * s * Ne) * np.sqrt(alpha * s * Ne / np.pi)

        r1 = u/U * inv_Tf
        r2 = (u+c)*H/U/L2 * inv_Tf

        Hf = 2*u*(L2-H) / U * inv_Tf - L2 * r2

        inv_Tb = 2*s*c*H*Ne
        Hb = -inv_Tb

        return(r1, r2, Hf+Hb)

    else:
        H = ret_H2(N, s, u, c, L1, L2)
        U = 2*(L2-H)*u + (u+c)*H + L1*u

        v = opt.newton(for_v, 0.001, args=(N, s, U))
        inv_Tf = U*v

        r1 = u/U * inv_Tf
        r2 = (u+c)*H/U/L2 * inv_Tf

        Hf = 2*u*(L2-H) / U * inv_Tf - L2 * r2

        inv_Tb = (np.exp(1) / v - 1) * c * H / (1 - np.log(v))
        Hb = -inv_Tb

        return(r1, r2, Hf+Hb)

def ret_rates_semi_for_diff(N, s, u, c, L1, L2, H):
    U = 2*(L2-H)*u + (u+c)*H + L1*u

    if N*s*np.exp(-U/s) > 1:
        alpha = np.sqrt(s/U)
        Ne = N*np.exp(-U/s)
        inv_Tf = alpha * s / np.exp(alpha * s * Ne) * np.sqrt(alpha * s * Ne / np.pi)

        r1 = u/U * inv_Tf
        r2 = (u+c)*H/U/L2 * inv_Tf

        Hf = 2*u*(L2-H) / U * inv_Tf - L2 * r2

        inv_Tb = 2*s*c*H*Ne
        Hb = -inv_Tb

        return(r1, r2, Hf+Hb)

    else:
        U = 2*(L2-H)*u + (u+c)*H + L1*u

        v = opt.newton(for_v, 0.001, args=(N, s, U))
        inv_Tf = U*v

        r1 = u/U * inv_Tf
        r2 = (u+c)*H/U/L2 * inv_Tf

        Hf = 2*u*(L2-H) / U * inv_Tf - L2 * r2

        inv_Tb = (np.exp(1) / v - 1) * c * H / (1 - np.log(v))
        Hb = -inv_Tb

        return(r1, r2, Hf+Hb)

def return_semi_diff(var, t, N, s, u, c):
    L1 = var[0]
    L2 = var[1]
    H = var[2]

    rates = ret_rates_semi_for_diff(N, s, u, c, L1, L2, H)
    dL1dt = -L1*rates[0]
    dL2dt = -L2*rates[1]
    dHdt = rates[2]
    dtdt = 1
    return [dL1dt, dL2dt, dHdt, dtdt]

def ret_degenerate_semi_process(N, s, u, c, L1, L2, H, max_t, seps):
    tmp_seps = seps

    error_list = [0]
    while len(error_list) == 1:
        print(tmp_seps, error_list)
        for i in error_list:
            print(i)

        t_list = np.linspace(0.0, max_t, tmp_seps)
        var_init = [L1, L2, H, 0]

        with warnings.catch_warnings(record = True) as error_list:
            var_list = odeint(return_semi_diff, var_init, t_list, args=(N, s, u, c))

        if len(error_list) >= 1:
            tmp_seps = int(tmp_seps / 2)

    var_init = [L1, L2, H, 0]
    t_list = np.linspace(0.0, max_t, tmp_seps)
    var_list = odeint(return_semi_diff, var_init, t_list, args=(N, s, u, c))
    print(tmp_seps)

    return(var_list)
