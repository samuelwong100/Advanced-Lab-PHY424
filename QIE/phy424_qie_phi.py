import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

filename = 'phy424_qie_phi_data.txt'
phi_data = np.loadtxt(filename, usecols=(0,1)).T
phi = phi_data[0]
corr = phi_data[1]

plt.figure(1, figsize=(8,6))
plt.scatter(phi,corr)
plt.xlabel(r'$\phi$')
plt.ylabel('Correlation Counts')
plt.savefig('phy424_qie_phi.png')
#plt.clf()

mean_corr = []
phi_unique = []
a = phi[0]
sum_list = []

for i in range(len(phi)):
    if a == phi[i]:
        sum_list.append(corr[i])
    elif a != phi[i]:
        mean_corr.append(np.mean(sum_list))
#        print(mean_corr)
        phi_unique.append(a)
        sum_list = []
        a = phi[i]
        sum_list.append(corr[i])
mean_corr.append(np.mean(sum_list))
phi_unique.append(a)


def model(x, a, b, c, d):
    #return a*np.sin(b*x + c) + d
    return a*np.sin(b*x + c)**2 + d

#plt.plot(x, model(x,67,0.07,1,64))
#plt.plot(x, model(x, *popt))
plt.xlabel(r'$\phi$')
#popt, pcov = curve_fit(model, theta_unique, mean_corr, p0=(67,0.07,1,64))
popt, pcov = curve_fit(model, phi_unique, mean_corr, p0=(30,0.07,1.2,15))

x = np.arange(0, 55, 0.5)

plt.figure(2, figsize=(8,6))
plt.scatter(phi_unique, mean_corr)
#plt.plot(x,model(x,30,0.07,1.2,15))
plt.plot(x, model(x, *popt))
plt.ylabel('Correlation Counts')#plt.savefig('phy424_qie_phi_mean.png')

