% example here: https://empossible.net/wp-content/uploads/2020/01/Lecture-Review-Walkthrough-of-1D-FDTD.pdf

clear;

fdtd_1d_run([
    30.48 2 6
    ], 1e-2, 1e-4, 1e9, 10000);
