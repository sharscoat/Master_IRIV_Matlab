#!/usr/bin/env python
#-*- coding:utf-8 -*



import numpy as np
import scipy as sp
from math import *
import scipy.optimize as optimization
from scipy.optimize import curve_fit
import pandas as pd
from matplotlib.mlab import find
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
# from matplotlib.mlab import *      # MAD ca c'est une mauvaise idee... et inutile ici
import sys
from tkFileDialog import askopenfilename




#############################################################################


# 3ère fonction : fonction "dérivée".

def derivee (t) : 
	dt = t[1:]-t[:-1]
	return dt




# calcule de la diagonale entre le point "(0,1)" et un point de la courbe
def diag(t0,t,b,d):
    dx = t - t0
    dy = d - b
    return np.sqrt(dy**2 + dx**2)




## formule de l'équation
def func(t, a, d, c, t0, t1, T1, T2): 
    """ la fonction 'funct' prend 8 parmetres :
    - t : 
    - a : 
    - d :
    - c :
    - t0:
    - t1:
    - T0:
    - T1:
    """
    base = c
    monte = ((d-c)-a*(t1-t))*(t>t0)*(t<=t1)*(1-np.exp(-(t-t0)/T1))
    descend = (d-c)*np.exp(-(t-t1)/T2)*(t>t1)
    return base + monte + descend


## formule de l'équation
def func2(t, a, d, c, t0, t1, T1, T2, t2): 
    """ la fonction 'funct' prend 8 parmetres :
    - t : 
    - a : 
    - d :
    - c :
    - t0:
    - t1:
    - T0:
    - T1:
    - t2:
    """
    base = c
    monte = ((d-c)-a*(t1-t))*(t>t0)*(t<=t1)*(1-np.exp(-(t-t0)/T1))  
    descend = (d-c)/(1+np.exp(T2*((t-t1)-(t2-t1))))*(t>t1)       	# 1/(1+np.exp(T2*((t-t1)-(t2-t1))))
    return base + monte + descend


 
# cadre
def cadre(t, a, d, c, t0, t1):
    """ dessine le cadre. Prend 6 parmetres :
    - t : 
    - a : 
    - d :
    - c :
    - t0:
    - t1:
    """
    x = [0, t0, t0, t1, t1, max(t)]
    y = [c, c,  d,  d,  c,  c]
    plt.plot(x, y)


# fonction : "summary"
def summary(valeur):
	minimum = round(np.min(valeur),2)
	q1 = round(np.percentile(valeur,25),2)
	mediane = round(np.median(valeur),2)
	q2 = round(np.percentile(valeur,75),2)
	maxiumum = round(np.max(valeur),2)
	moyenne = round(np.mean(valeur),2)
	ecart_type = round(np.std(valeur),2)
	variance = round(np.var(valeur),2)
	#resume = pd.DataFrame(np.array([minimum,q1,mediane,q2,moyenne,ecart_type,variance]), columns = ['min.','q1.','med.','q2','moy.','ec_ty.','var.'])
	return np.array([minimum,q1,mediane,q2,maxiumum,moyenne,ecart_type,variance])



# fonction de lissage de courbe par moyenne glissante : "lissage"
def lissage(Ly,p):
	'''Fonction qui débruite une courbe par une moyenne glissante sur 2P+1 points'''
	Lyout = []
	for i in range(p,len(Ly)-p):
		val=0
		for k in range(2*p):
			val+=Ly[i-p+k]
		Lyout.append(val/2/p)

	return Lyout




#############################################################################






#######################################################################################################################


# 1ère fonction : nettoyage de la courbe avec selection de la  portion utile
#				  nécessite l'utilisation de fonction 'find' du module 'matplotlib.mlab' (from matplotlib.mlab import find)

from matplotlib.mlab import find
 
def nettoyage_cap (yvalues,temps) :
	Y=max(yvalues)
	N=len(yvalues);
	i1=find(yvalues>Y/20)[0]-100
	i2=find(yvalues>Y/20)[-1]+10
	if i1<0:
	    i1=0
	if i2>N:
	    i2=N

	cap = {
			'yvalues' : yvalues[i1:i2],
			'temps' : temps[i1:i2]
			}
	return cap




# 2ère fonction : recherche les cycles de capnographie et les découpes en capnogramme.
#				  nécessite l'utilisation de fonction 'find' du module 'matplotlib.mlab' (from matplotlib.mlab import find)

def cycle_cap (yvalues,temps) :
	ym = (yvalues-min(yvalues))/( max(yvalues)-min(yvalues))  # Je normalise
	i=0
	b = []
	while True :
	    try:
	        imax =  i + find( (ym[i:]) > 0.4)[0]
	        #le premier index tel que y à 40% de developpement ou plus
	        imin1 = imax   + find( (ym[imax:]) < 0.2)[0]
	        #l'index suivant tel que y à 20% de developpement ou moins
	        imin2 = imin1  + find( (ym[imin1:]) > 0.25)[0]
	        #l'index encore suivant tel que y à 25% de developpement ou plus
	    except IndexError:
	        #print imax, imin1, imin2, "break"
	        break
	    i = np.floor((imin1+imin2)/2)    # i entre imin1 et imin2
	    b.append(i)
	    i = imin2
	B = np.array(b, dtype=int)
	#B1 = np.array(b, dtype=int)  # je force en entier pour pouvoir l'utiliser comme index : B = np.array(b, dtype=int)
	#B = np.concatenate(([0],B1,[len(ym)]),axis=0)    # permet de récuperer le premier cycle et le dernier cycle.
	NP = len(B)-1
	# Creation d'une liste de Capnogramme
	yp = []
	tpr = []
	for i in range(len(B)-1):
	    yp.append(yvalues[B[i]:B[i+1]])
	    tpr.append(temps[B[i]:B[i+1]])
	
	cap = {
			'yp' : yp,
			'tpr' : tpr
			}
	return cap





# 3ère fonction : permet une analyse des capnogrammes.

def AnalyseCapno(yvalues,temps):
	yp  = yvalues 		# pour avoir le capnogramme au plus juste : yp  = yvalues[10:-10]
	tpr = temps			# pour avoir le capnogramme au plus juste : tpr = temps[10:-10]
	MM  = max(yp)
	mm  = min(yp)
	y1  = yp[1:]-yp[:-1]
	y2  = y1[1:]-y1[:-1]               
	t0  = find(y1>max(y1)*0.1)[0]
	t1  = find(y1<min(y1)*0.1)[0]
	#t1  = find(yp==max(yp))       # un des enjeu majeur est de trouver la meilleur valeur de t1 ()
	t2  = find(y1==min(y1))
	c   = np.mean(yp[:t0])
	d   = MM
	a   = 0
	T1  = 5    #(t1-t0)/10
	T2  = T1/2
	max_dy = max(y1)   
	t = np.arange(len(yp))
	p01 = [a, d, c, t0, t1, T1, T2]
	popt1,pcov1 = sp.optimize.curve_fit(func, t, yp, p01)
	la1, ld1, lc1, lt01, lt11, lT11, lT21 = popt1
	ycalc1 = func(t,*popt1)
	ycalc_d1 = [ycalc1[1:]-ycalc1[:-1]]
	d1_max = find(ycalc_d1 == max(ycalc_d1))
	d1_mix = find(ycalc_d1 == min(ycalc_d1))
	p02 = [a, d, c, t0, t1, T1, T2, t2]
	popt2,pcov2 = sp.optimize.curve_fit(func2, t, yp, p02)
	la2, ld2, lc2, lt02, lt12, lT12, lT22, lt22 = popt2
	ycalc2 = func2(t,*popt2)
	ycalcm = (ycalc1+ycalc2)/2
	tn = (t-np.mean(t))/np.std(t)                      
	trn = (yp-np.mean(yp))/np.std(yp)                   
	t0n = (lt01-np.mean(t))/np.std(t)
	dn = (d-np.mean(yp))/np.std(yp)  
	diagonale = diag(t0n, tn, trn, dn)
	t2n = (t2-np.mean(t))/np.std(t)+0.5
	diagonale2 = diag(t2n, tn, trn, dn+0.5)
	tdiagmin = find(diagonale == min(diagonale))
	ydiagmin = yp[tdiagmin]
	tdiagmin2 = find(diagonale2 == min(diagonale2))
	ydiagmin2 = yp[tdiagmin2]
	xcad = [0, lt01, lt01, tdiagmin2, tdiagmin2, max(t)]
	ycad = [lc1, lc1, d, d, lc1, lc1]
	chi2 = np.sqrt( np.sum((yp-ycalc1**2)))
	aire = (tdiagmin2-lt01)*(d-lc1)            # surface du carre
 	onde = yp[int(lt01):int(tdiagmin2)]-lc1     # aire sous la courbe de capnographie   
	ratio = sum(onde) / aire         	# ratio entre sous la courbe et la surface du carre 
	ratio2 = ((aire-sum(onde)) / (aire/2))*100       # ratio entre sous la courbe et la surface du carre divisée par deux 
	# ratio20
	t20 = lt01 + 20
	if t20 <= lt11 :
		d20 = yp[t20]
		aire20 = (t20-lt01)*(d20-lc1)            # surface du carre
	 	onde20 = yp[int(lt01):int(t20)]-lc1     # aire sous la courbe de capnographie   
		ratio20 = ((aire20-sum(onde20)) / (aire20/2))*100
	else :
		ratio20 = 0
		d20 = 0
	xcad20 = [0, lt01, lt01, t20, t20, max(t)]
	ycad20 = [lc1, lc1, d20, d20, lc1, lc1]
	# ratio40
	t40 = lt01 + 40
	if t40 <= lt11 :
		d40 = yp[t40]
		aire40 = (t40-lt01)*(d40-lc1)            # surface du carre
	 	onde40 = yp[int(lt01):int(t40)]-lc1     # aire sous la courbe de capnographie   
		ratio40 = ((aire40-sum(onde40)) / (aire40/2))*100
	else :
		ratio40 = 0
		d40 = 0
	xcad40 = [0, lt01, lt01, t40, t40, max(t)]
	ycad40 = [lc1, lc1, d40, d40, lc1, lc1]
	# ratio60
	t60 = lt01 + 60
	if t60 <= lt11 :
		d60 = yp[t60]
		aire60 = (t60-lt01)*(d60-lc1)            # surface du carre
	 	onde60 = yp[int(lt01):int(t60)]-lc1     # aire sous la courbe de capnographie   
		ratio60 = ((aire60-sum(onde60)) / (aire60/2))*100
	else :
		ratio60 = 0
		d60 = 0
	xcad60 = [0, lt01, lt01, t60, t60, max(t)]
	ycad60 = [lc1, lc1, d60, d60, lc1, lc1]
	# ratio80
	t80 = lt01 + 80
	if t80 <= lt11 :
		d80 = yp[t80]
		aire80 = (t80-lt01)*(d80-lc1)            # surface du carre
	 	onde80 = yp[int(lt01):int(t80)]-lc1     # aire sous la courbe de capnographie   
		ratio80 = ((aire80-sum(onde80)) / (aire80/2))*100
	else :
		ratio80 = 0
		d80 = 0
	xcad80 = [0, lt01, lt01, t80, t80, max(t)]
	ycad80 = [lc1, lc1, d80, d80, lc1, lc1]
	# ratio100
	t100 = lt01 + 100
	if t100 <= lt11 :
		d100 = yp[t100]
		aire100 = (t100-lt01)*(d100-lc1)            # surface du carre
	 	onde100 = yp[int(lt01):int(t100)]-lc1     # aire sous la courbe de capnographie   
		ratio100 = ((aire100-sum(onde100)) / (aire100/2))*100
	else :
		ratio100 = 0
		d100 = 0
	xcad100 = [0, lt01, lt01, t100, t100, max(t)]
	ycad100 = [lc1, lc1, d100, d100, lc1, lc1]
	'''
	# ratio120
	t120 = lt01 + 120
	if t120 <= lt11 :
		d120 = yp[t120]
		aire120 = (t120-lt01)*(d120-lc1)            # surface du carre
	 	onde120 = yp[int(lt01):int(t120)]-lc1     # aire sous la courbe de capnographie   
		ratio120 = (((aire120-sum(onde120)) / (aire120/2))*100
	else :
		ratio120 = 0
	'''
	#######
	indice = (aire-sum(onde))/len(onde)
	itc = aire-sum(onde)
	residu1 = yp - ycalc1
	residu2 = yp - ycalc2
	residum = yp - ycalcm
	ybar = np.sum(yp)/len(yp)          # or sum(y)/len(y)
	ssreg1 = np.sum((ycalc1-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	sstot = np.sum((yp - ybar)**2)
	determination1 = ssreg1 / sstot
	ssreg2 = np.sum((ycalc2-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	sstot = np.sum((yp - ybar)**2)
	determination2 = ssreg2 / sstot
	ssregm = np.sum((ycalcm-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	sstot = np.sum((yp - ybar)**2)
	determinationm = ssregm / sstot
	ypu = yp[lt01:lt11]
	ycalcu1 = ycalc1[lt01:lt11]
	ybaru = np.sum(ypu)/len(ypu)          # or sum(y)/len(y)
	ssregu1 = np.sum((ycalcu1-ybaru)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	sstotu = np.sum((ypu - ybaru)**2)
	determinationu1 = ssregu1 / sstotu
	ycalcu2 = ycalc2[lt01:lt11]
	ssregu2 = np.sum((ycalcu2-ybaru)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	determinationu2 = ssregu2 / sstotu
	ycalcum = ycalcm[lt01:lt11]
	ssregum = np.sum((ycalcum-ybaru)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
	determinationum = ssregum / sstotu
	tpointeret = tdiagmin
	ypointeret = ydiagmin
	xdiag = [lt01,tdiagmin]
	ydiag = [d,ydiagmin]
	xdiag2 = [lt01,tdiagmin2]
	ydiag2 = [d,ydiagmin2]
	xtri = [lt01,tdiagmin2]
	ytri = [lc1,d]
	xpent = [tdiagmin,tdiagmin2]
	ypent = [ydiagmin,ydiagmin2]


	capno = {
				'yp'  : yp,
				'tpr' : tpr,
				'MM'  : MM,
				'mm'  : mm,
				'y1'  : y1,
				'y2'  : y2,          
				't0'  : t0,
				't1'  : t1,
				't2'  : t2,
				'c'   : c,
				'd'   : d,
				'a'   : a,
				'T1'  : T1,
				'T2'  : T2,
				'la1' : la1*1000,
				'lT11' : lT11,
				'la2' : la2*1000,
				'lT12' : lT12,
				'max_dy' : max_dy,
				't'   : t,
	    		'p01'  : p01,
	    		'popt1': popt1,
	    		'pcov1': pcov1,
	    		'ycalc1': ycalc1,
	    		'ycalc_d1': ycalc_d1,
	    		'd1_max' : d1_max,
			    'd1_mix' : d1_mix,
			    'p02'  : p02,
	    		'popt2': popt2,
	    		'pcov2': pcov2,
	    		'ycalc2': ycalc2,
	    		'ycalcm' : ycalcm,
			    'xcad' : xcad,
			    'ycad' : ycad,
			    'xcad20' : xcad20,
			    'ycad20' : ycad20,
			    'xcad40' : xcad40,
			    'ycad40' : ycad40,
			    'xcad60' : xcad60,
			    'ycad60' : ycad60,
			    'xcad80' : xcad80,
			    'ycad80' : ycad80,
			    'xcad100' : xcad100,
			    'ycad100' : ycad100,
			    'tn' : tn,                      
			    'trn' : trn,                 
			    't0n' : t0n, 
			    'dn' : dn,  
			    'diagonale': diagonale,
			    'diagonale2': diagonale2,
			    'tdiagmin': tdiagmin,
			    'ydiagmin': ydiagmin,
			    'chi2' : chi2, 
			    'aire' : aire,           
			    'onde' : onde,    
			    'ratio' : ratio,   
			    'ratio2' : ratio2,
			    'ratio20' : ratio20,
			    'ratio40' : ratio40, 
			    'ratio60' : ratio60, 
			    'ratio80' : ratio80, 
			    'ratio100' : ratio100, 
			    #'ratio120' : ratio120,   
			    'indice' : indice,
			    'itc' : itc,
			    'residu1' : residu1,
			    'residu2' : residu2,
			    'residum' : residum,
			    'r2' : determination1,
			    'r2u' : determinationu1,
			    'r22' : determination2,
			    'r22u' : determinationu2,
			    'r2m' : determinationm,
			    'r2mu' : determinationum,
			    'tpointeret':tpointeret,
			    'ypointeret':ypointeret,
			    'xdiag' : xdiag,
			    'ydiag' : ydiag,
			    'xdiag2' : xdiag2,
			    'ydiag2' : ydiag2,
			    'xtri' : xtri,
			    'ytri' : ytri,
			    'xpent' : xpent,
			    'ypent' : ypent,
			}
	return capno

