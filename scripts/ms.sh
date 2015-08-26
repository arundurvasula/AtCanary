#!/bin/bash
# ms simulation to get the expected null distribution of haplotypes based on
# the demography of the Iberia-Morocco-Canary Island populations

ms 51 10000 -p 5 -t 280 -r 36 100000 -I 3 14 20 17 -n 1 8.78 -n 2 4.65 -n 3 0.08 -ma x 0.42 7.31 0.46 x 1.46 0 0 x -ej 2.42 3 2 -ej 3.38 2 1 -eg 3.38 1 0.5 -eg 2.42 2 0.3 -eg 2.42 3 -1

#51 1 -I 3 14 20 17 -n 1 8.78 -n 2 4.65 -n 3 0.08 -eg 3.38 -ej 2.42 3 2 -ej 3.38 2 1
# explanation:
# ms 51 100: 51 samples, 100 reps
# -p 5: print 5 decimal places for the positions
# -t 28: theta is 28 (see math below)
# -r 3.6 10000: rho is 3.6, number of base pairs is 10000 (see math below)
# -I 3 14 20 17: simulate an island demography with 3 subpopulations of size 14, 20, and 17
# -n 1 8.78: make the Iberian population size 8.78*N_0
# -n 2 4.65: make the Moroccan population size 4.65*N_0
# -n 3 0.08: make the Canary Island population size 0.08*N_0

# defined parameters:
# L = 100000
# N_0 = 10^5
# mu = 7e-9
# r = 1.8e-8/20

# math
# theta = 4*N_0*L*mu = 4*10^5*100000*7e-9 = 280
# rho = 4*N_0*L*r = 4*10^5*100000*1.8e-8/20 = 36
# IB growth rate: 1.64=8.78*e^(-x*3.38) = 0.5
# MO growth rate: 2.26=4.65*e^(-x*2.42) = 0.3
# CA growth rate: 0.91=0.08*e^(-x*2.42) = -1.0

# theta = 1
# Nref = 3.57e7
# theta = 4*N_0*L*mu = 4*3.57e7*10000*7e-9 = 10000
# rho = 4*N_0*L*r = 4*10^5*100000*1.8e-8/20 = 36


# dadi parameters:
# nIb (start):		1.64
# nIb (end):		8.78
# nMo (start):	2.26
# nMo (end):	4.65
# nCa (start):	0.91
# nCa (end):		0.08
# mIbMo:		0.20
# mMoIb:		0.23
# mIbCa:			3.66
# mMoCa:		0.73
# T1:				0.48
# T2:				1.21
# dadi -> ms:
# migration: multiply by 2
# 0.20945176,  0.22806955,  3.65669199,  0.7312348
# 0.4189035, 0.4561391 7.3133840 1.4624696
# time: divide by 2
# 0.47997895,  1.20825076
# 0.9599579 2.4165015

# switch canary demography
ms 51 10000 -p 5 -t 28 -r 3.6 10000 -I 3 14 20 17 -n 1 8.78 -n 2 4.65 -n 3 0.91 -ma x 0.42 7.31 0.46 x 1.46 0 0 x -ej 2.42 3 2 -ej 3.38 2 1 -eg 3.38 1 0.5 -eg 2.42 2 0.3 -eg 2.42 3 1

# previous simulations aren't scaled correctly
# switch canary demography, theta = 16.71
# need to change times
# nref = 46428
# l = 10000
# ib time: -log(1.64/8.78)/2.42 = 0.693
# mo time: -log(2.26/4.65)/2.42 = 0.298
# ca time: -log(0.08/0.91)/2.42 = 1.00
ms 51 1 -T -p 5 -t 16.71 -r 1.67 10000 -I 3 14 20 17 -n 1 8.78 -n 2 4.65 -n 3 0.91 -g 1 0.693 -g 2 0.298 -g 3 1 -ma x 0.46 0 0.42 x 0 7.31 1.46 x -ej 2.42 3 2 -ej 3.38 2 1 -eg 2.42 1 0.0 -eg 2.42 2 0.0 -eg 2.42 3 0.0
