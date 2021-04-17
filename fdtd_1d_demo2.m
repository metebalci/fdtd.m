% anti-reflection example #1 here: https://empossible.net/wp-content/uploads/2020/01/Lecture-Examples-of-1D-FDTD.pdf

clear;

c0 = 299793458; % m/s

e_r = 12
e_air = 1
e_ar = sqrt(e_r * e_air)
n_ar = sqrt(e_ar)
f0 = 2.4e9;
lambda_0 = c0 / f0
d_ar = floor(lambda_0 / (4 * n_ar) / 1e-3)

fdtd_1d_run([
    d_ar e_ar 1
    304.8 e_r 1
    d_ar e_ar 1
    ], 1e-3, 1e-4, 4e9, 10000)