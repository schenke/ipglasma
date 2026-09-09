#!/usr/bin/env python3
import numpy as np
#import collections
#from enum import Enum
import math
import pylab as pl
#from collections import Counter
#import pandas as pd
from numpy import *
from sys import argv, exit 
istep = 0#argv[1]
istep = int(istep)
ievent = 0#argv[1]
ievent = int(ievent)
aa = np.loadtxt("epsilon-u-Hydro-TauHydro-0.dat".format(ievent))
length =512

x = []
dx = 0.0585938
T00 = []
temp = aa[:,3]
largea = temp[:][ temp > 100
] 
print("number of large energy density point is ", np.max(temp), len(largea), len(largea)/length)
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
    x.append(-1*length/2.*dx + ii*dx)
T00 = np.array(T00)
x=np.array(x)
y=x
X, Y = np.meshgrid(x, y)
pl.figure(10020)
levels = np.linspace(0,200,50)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar()
    pl.title('Pb-Pb, ed, second shift, b={}fm'.format(istep), fontsize=15)# give plot a title
    pl.savefig("Png/ed_{}.png".format(istep))
    
T00 = []
temp = aa[:,5]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(100201)
levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, ux at second shift, b={}fm'.format(istep), fontsize=15)# give plot a title
    pl.savefig("Png/ux_{}.png".format(istep))
    
temp = aa[:,6]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(100202)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, uy at second shift, b={}fm'.format(istep), fontsize=15)# give plot a title
    pl.savefig("Png/uy_{}.png".format(istep))

temp = aa[:,8]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
print("pitautau ", np.min(T00),np.max(T00))
levels = np.linspace(np.min(T00),np.max(T00),21)
levels = np.linspace(-40,40,21)
pl.figure(1002022)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    name = '$\pi^{\\tau\\tau}$'
    pl.title('Pb-Pb, {} at second shift, b={}fm'.format(name,istep), fontsize=15)# give plot a title
    pl.savefig("Png/pitautau_{}.png".format(istep))
temp = aa[:,9]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
print(np.min(T00),np.max(T00))
#levels = np.linspace(-40,40,21)
pl.figure(10020212)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    name = '$\pi^{\\tau x}$'
    pl.title('Pb-Pb, {} at second shift, b={}fm'.format(name,istep), fontsize=15)# give plot a title
    pl.savefig("Png/pitaux_{}.png".format(istep))

temp = aa[:,10]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
print(np.min(T00),np.max(T00))
#levels = np.linspace(np.min(T00),np.max(T00),21)
pl.figure(1010202212)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    name = '$\pi^{\\tau y}$'
    pl.title('Pb-Pb, {} at second shift, b={}fm'.format(name,istep), fontsize=15)# give plot a title
    pl.savefig("Png/pitauy_{}.png".format(istep))
    
temp = aa[:,12]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
print("pixx = ",np.min(T00),np.max(T00))
#levels = np.linspace(-1,2000,21)
#levels = np.linspace(np.min(T00),np.max(T00),21)
pl.figure(101202022)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    name = '$\pi^{xx}$'
    pl.title('Pb-Pb, {} at second shift, b={}fm'.format(name,istep), fontsize=15)# give plot a title
    pl.savefig("Png/pixx_{}.png".format(istep))

temp = aa[:,15]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
print("piyy = ",np.min(T00),np.max(T00))
#levels = np.linspace(np.min(T00),np.max(T00),21)
pl.figure(10021022)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    name = '$\pi^{yy}$'
    pl.title('Pb-Pb, {} at second shift, b={}fm'.format(name,istep), fontsize=15)# give plot a title
    pl.savefig("Png/piyy_{}.png".format(istep))
    
temp = aa[:,13]
T00 = []
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
print("pixy = ", np.min(T00),np.max(T00))
#levels = np.linspace(np.min(T00),np.max(T00),21)
pl.figure(1002103422)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    #pl.xlim(-2.5,2.5)
    #pl.ylim(-2.5,2.5)
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    name = '$\pi^{xy}$'
    pl.title('Pb-Pb, {} at second shift, b={}fm'.format(name,istep), fontsize=15)# give plot a title
    pl.savefig("Png/pixy_{}.png".format(istep))

           
