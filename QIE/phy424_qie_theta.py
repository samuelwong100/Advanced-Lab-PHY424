import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

filename = 'phy424_qie_theta.txt'
theta_data = np.loadtxt(filename, usecols=(0,1)).T
theta = theta_data[0]
corr = theta_data[1]

plt.figure(1, figsize=(8,6))
plt.scatter(theta,corr)
plt.xlabel(r'$\theta$')
plt.ylabel('Correlation Counts')
plt.savefig('phy424_qie_theta.png')
plt.clf()

mean_corr = []
a = theta[0]
sum_list = []

for i in range(len(theta)):
    if a == theta[i]:
        sum_list.append(corr[i])
    elif a != theta[i]:
        mean_corr.append(np.mean(sum_list))
        sum_list = []
        a = theta[i]
        sum_list.append(corr[i])
mean_corr.append(np.mean(sum_list))

theta_unique = np.unique(theta)

def model(x, a, b, c, d):
    #return a*np.sin(b*x + c) + d
    return a*np.sin(b*x + c)**2 + d

#popt, pcov = curve_fit(model, theta_unique, mean_corr, p0=(67,0.07,1,64))
popt, pcov = curve_fit(model, theta_unique, mean_corr, p0=(67,0.02,1,64))

x = np.arange(0, 90, 0.5)

plt.figure(1, figsize=(8,6))
plt.scatter(theta_unique, mean_corr)
#plt.plot(x, model(x,67,0.07,1,64))
plt.plot(x, model(x, *popt))
plt.xlabel(r'$\theta$')
plt.ylabel('Correlation Counts')
plt.savefig('phy424_qie_theta_mean.png')
