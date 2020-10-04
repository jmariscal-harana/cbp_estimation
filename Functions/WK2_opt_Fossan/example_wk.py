import numpy as np
import matplotlib.pylab as plt
import csv
from ODEschemes import euler
import scipy.optimize

def calcAverageError(refData, numData, data_type="P"):
    
    if data_type == "P":
        error = np.sum(np.abs((numData - refData)/refData))/len(refData)
    elif data_type == "Q":
        error = np.sum(np.abs((numData - refData)/np.amax(refData)))/len(refData)
        
    return error

def loadTPQfromCSV(type = 'young'):
    if type == 'old':
        reader = csv.DictReader(open('55arteryNetwork_age75.csv','rb'),delimiter=';')
    else:
        reader = csv.DictReader(open('55arteryNetwork_age19.csv','rb'),delimiter=';')
    rows = list(reader)
    numberOfPoints = len(rows)
    t = np.empty(numberOfPoints)
    P = np.empty(numberOfPoints)
    Q = np.empty(numberOfPoints)
    for i,row in enumerate(rows):
        t[i] = row['t [s]']
        P[i] = row['P [Pa]']
        Q[i] = row['Q [m^3/s]']
        
    return np.array(t),np.array(P),np.array(Q)


def WK2(p,t,q,prms):
    R=prms[0]
    C=prms[1]
    
    dpdt=q/C-p/(R*C) 
    return dpdt


def discrepancy(C):
    prms=np.zeros(2)
    prms[0]=R
    prms[1]=C
    p_wk=euler(WK2,p[0],t,Q,prms)
    
    method='ppm'
    
    if method=='ppm':
        error =  ((np.max(p_wk)-np.min(p_wk))-(np.max(p)-np.min(p)))**2   # Error in Pulse pressure
    else:
       error = calcAverageError(p, p_wk) # Average relativ error
       
    return error
    

# Main program
t, p, Q = loadTPQfromCSV()
print p.size
p = p/133.322365 # Convert to mmHg
Q = Q/1e-6 # Convert to ml/s

# Set parameters for the model
R = np.mean(p)/np.mean(Q) 
    
C0 = 1.0

optResult = scipy.optimize.minimize(discrepancy, C0, method='Nelder-Mead')

print "success: ", optResult['success']
C = optResult['x']
print 'C optimal:', C

prms = [R, C]
p_wk = euler(WK2, p[0], t, Q, prms)

fig = plt.figure()

plt.plot(t, p)             
plt.plot(t, p_wk)
plt.legend(['P','Wk2'])

plt.show()

