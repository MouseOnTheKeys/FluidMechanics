# -*- coding: utf-8 -*-
"""
Spyder Editor

Graficko prikazivanja uniformnog strujanja, ponora i izvora, dvopola, Rankinovog Polutela...
"""

import math
import numpy
from matplotlib import pyplot

def uniformnaStruja(Uinf, beta, X, Y):
    cosBeta = math.cos(beta*math.pi/180)
    sinBeta = math.sin(beta*math.pi/180)
    u = Uinf*cosBeta
    v = Uinf*sinBeta
    phi = Uinf*(X*cosBeta + Y*sinBeta)
    psi = Uinf*(Y*cosBeta - X*sinBeta)
    return phi, psi, u, v

def izvor(izdasnost, x0, y0, X, Y):
    u = izdasnost/(2*numpy.pi)*(X-x0)/((X-x0)**2+(Y-y0)**2)
    v = izdasnost/(2*numpy.pi)*(Y-y0)/((X-x0)**2+(Y-y0)**2)
    phi = izdasnost/(2*numpy.pi)*numpy.log(numpy.sqrt((X-x0)**2 + (Y-y0)**2))
    psi = izdasnost/(2*numpy.pi)*numpy.arctan2((Y-y0), (X-x0))
    return phi, psi, u, v

def ponor(izdasnost, x0, y0, X, Y):
    u = -izdasnost/(2*numpy.pi)*(X-x0)/((X-x0)**2+(Y-y0)**2)
    v = -izdasnost/(2*numpy.pi)*(Y-y0)/((X-x0)**2+(Y-y0)**2)    
    phi = -izdasnost/(2*numpy.pi)*numpy.log(numpy.sqrt((X-x0)**2 + (Y-y0)**2))
    psi = -izdasnost/(2*numpy.pi)*numpy.arctan2((Y-y0), (X-x0))
    return phi, psi, u, v

def vrtlog(Gamma, x0, y0, X, Y):
    u = -Gamma/(2*math.pi) * (Y-y0)/((X-x0)**2 + (Y-y0)**2)
    v = Gamma/(2*math.pi) * (X-x0)/((X-x0)**2 + (Y-y0)**2)
    phi = Gamma/(2*numpy.pi)*numpy.arctan2((Y-y0), (X-x0))
    psi = -Gamma/(2*numpy.pi)*numpy.log(numpy.sqrt((X-x0)**2 + (Y-y0)**2))
    return phi, psi, u, v

def dvopol(M, x0, y0, X, Y):
    u = - M/(2*math.pi)*((X-x0)**2-(Y-y0)**2)/((X-x0)**2+(Y-y0)**2)**2
    v = - M/(2*math.pi)*2*(X-x0)*(Y-y0)/((X-x0)**2+(Y-y0)**2)**2
    psi =  -M/(2*math.pi)*(Y-y0)/((X-x0)**2+(Y-y0)**2)
    phi = M/(2*math.pi)*(X-x0)/((X-x0)**2+(Y-y0)**2)
    return phi, psi, u, v

N = 800                                  # broj tacaka u x i y pravcu
x_start, x_end = -2.0, 2.0               # xmin i xmax
y_start, y_end = -2.0, 2.0               # ymin i yman
x = numpy.linspace(x_start, x_end, N)    # 1D niz x koordinata
y = numpy.linspace(y_start, y_end, N)    # 1D niy y koordinata


X, Y = numpy.meshgrid(x, y)              # MREZA PRORACUNSKIH TACAKA

# phi1, psi1, u1, v1 = uniformnaStruja(5, 0, X, Y)
# phi2, psi2, u2, v2 = izvor(10, -1, 0, X, Y)
# phi3, psi3, u3, v3 = ponor(10, 1, 0, X, Y)
# u = u1+u2+u3
# v = v1 +v2+v3
# psi = psi1 + psi2+psi3

# # --------------   CiklicnoOPSTRUJAVANJE CILINDRA ---------------------------#
# phi1, psi1, u1, v1 = uniformnaStruja(5, 0, X, Y)
# phi2, psi2, u2, v2 = dvopol(5, 0, 0, X, Y)
# u = u1 + u2
# v = v1 + v2
# psi = psi1 + psi2




# ----------------------  ACIKLICNO OPSTRUJAVANJE CILINDRA ----------------------#
# R = 0.5
# Uinf = 1.0
# Gamma = 4.05*math.pi*R*Uinf
# M = R**2*math.pi*Uinf ## Poluprecnik kruga 1
# phi1, psi1, u1, v1 = uniformnaStruja(Uinf, 0, X, Y)
# phi2, psi2, u2, v2 = dvopol(M, 0, 0, X, Y)
# phi3, psi3, u3, v3 = vrtlog(-Gamma, 0, 0, X, Y)
# u = u1 + u2 + u3
# v = v1 + v2 + v3
# psi = psi1 + psi2 + psi3
# K = Gamma/(2*math.pi)*math.log(R)




# ----------         KELVINOVA OVALA  ---------------------------------------------
# Uinf = 10
# a = 1
# K = 1
# gamma = K*Uinf*a*2*math.pi
# phi_inf, psi_inf, u_inf, v_inf = uniformnaStruja(10, 0, X, Y) 
# phi_vrtlog, psi_vrtlog, u_vrtlog, v_vrtlog = vrtlog(15, 0, 0, X, Y)
# phi_vrtlog1, psi_vrtlog1, u_vrtlog1, v_vrtlog1 = vrtlog(-gamma, 0, 1, X, Y)
# phi_vrtlog2, psi_vrtlog2, u_vrtlog2, v_vrtlog2 = vrtlog(gamma, 0, -1, X, Y)
# u = u_inf + u_vrtlog1 + u_vrtlog2
# v = v_inf + v_vrtlog1 + v_vrtlog2
# psi = psi_inf + psi_vrtlog1 + psi_vrtlog2
#----------------------------------------------------------------------------------


# ---------------------- Cetiri IZVORA -------------------------------------
# phi_iz1, psi_iz1, u_iz1, v_iz1 = izvor(5, 1, 1, X, Y)
# phi_iz2, psi_iz2, u_iz2, v_iz2 = izvor(5, -1, 1, X, Y)
# phi_iz3, psi_iz3, u_iz3, v_iz3 = izvor(5, 1, -1, X, Y)
# phi_iz4, psi_iz4, u_iz4, v_iz4 = izvor(5, -1, -1, X, Y)
# u = u_iz1 + u_iz2 + u_iz3 + u_iz4
# v = v_iz1 + v_iz2 + v_iz3 + v_iz4
# psi = psi_iz1 + psi_iz2 + psi_iz3 + psi_iz4
# #-------------------------------------------------------------

#  --------------------- Rankinova OVALA ----------------------#
# Uinf = 10
# epsilon = 0.5*Uinf*math.pi
# phi_inf, psi_inf, u_inf, v_inf = uniformnaStruja(-Uinf, 0, X, Y) 
# phi_iz, psi_iz, u_iz, v_iz = izvor(epsilon, -1, 0, X, Y)
# phi_pon, psi_pon, u_pon, v_pon = ponor(epsilon, 1, 0, X, Y)
# u = u_inf + u_iz + u_pon
# v = v_inf + v_iz + v_pon
# psi = psi_inf + psi_iz + psi_pon
# -------------------------------------------------------------#

#------- Dva izvora i ponor (teorema o kruznici)---------     #
epsilon = 1
phi_iz1, psi_iz1, u_iz1, v_iz1 = izvor(epsilon, 1.5, 0, X, Y)
phi_iz2, psi_iz2, u_iz2, v_iz2 = izvor(epsilon, 0.5, 0, X, Y)
phi_pon, psi_pon, u_pon, v_pon = ponor(epsilon, 0, 0, X, Y)
u = u_iz1 + u_iz2 + u_pon
v = v_iz1 + v_iz2 + v_pon
psi = psi_iz1 + psi_iz2 + psi_pon
#-------------------------------------------------------------#


# -------------  DVA DVOPOLA --------------------------------#
#phi0, psi0, u0, v0 = uniformnaStruja(1., 0, X, Y)
# phi1, psi1, u1, v1 = dvopol( 10, 1, 0, X, Y)
# phi2, psi2, u2, v2 = dvopol(10, -1, 0, X, Y)
# u =  u1 + u2
# v =  v1 + v2
# psi =  psi1 + psi2


# ------------ Uinf + IZVOR + ZID ---------------------------#
# Uinf = 2*math.sqrt(2.0)
# Eps = 8*math.pi
# a = 2.
# phi1, psi1, u1, v1 = uniformnaStruja(Uinf, 0, X, Y)
# phi2, psi2, u2, v2 = izvor(Eps, 0, a, X, Y)
# phi3, psi3, u3, v3 = izvor(Eps, 0, -a, X, Y)
# u = u1 + u2 + u3
# v = v1 + v2 + v3
# psi = psi1 + psi2 + psi3


# Uinf = 5
# Gamma = 10*math.pi
# phi1, psi1, u1, v1 = uniformnaStruja(Uinf, 0, X, Y)
# phi2, psi2, u2, v2 = vrtlog(Gamma, 0, 0, X,Y)
# u = u1 + u2
# v = v1 + v2
# psi = psi1 + psi2


size = 8
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.axes().set_aspect('equal', 'datalim')
pyplot.streamplot(X, Y, u, v, density=2.2, linewidth=1, arrowsize=2, arrowstyle='->')


pyplot.contour(X, Y, psi, levels= [-0.5, 0, 0.5], colors='#CD2305', linewidths=2, linestyles='solid')

pyplot.show()

