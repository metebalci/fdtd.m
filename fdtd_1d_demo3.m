% bragg gratings example #2 here: https://empossible.net/wp-content/uploads/2020/01/Lecture-Examples-of-1D-FDTD.pdf

clear;

c0 = 299793458; % m/s

e_sin = 1.5*1.5;
e_sio2 = 2*2;
d1 = 163;
d2 = 122;

fdtd_1d_run([
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    ], 1e-9, 1e-9, 2*c0/980e-9, 10000)
