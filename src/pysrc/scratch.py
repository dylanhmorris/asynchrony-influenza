import numpy as np
import matplotlib.pyplot as plt
import wh_popgen as whp
import legacy_wh_popgen as lwhp
from importlib import reload

reload(whp)
reload(lwhp)
f_target = 0.9
t = 3
R0 = 10
d_v = 4
bottleneck = 1
t_M = 2.75
k = 10
mu = .33e-4
c_w = 1
c_m = 0

delta = k * (c_w - c_m)

c_w = 0.8
c_m = 0.72
b = whp.p_repl(
    f_target,
    t,
    t_M,
    R0,
    d_v,
    bottleneck,
    k,
    c_w,
    c_m,
    mu)

c = lwhp.p_repl_old(
    f_target,
    t,
    t_M,
    R0,
    d_v,
    bottleneck,
    k,
    c_w,
    c_m,
    mu)

