import numpy as np
import matplotlib.pyplot as plt

filename = 'phy424_qie_phi_data_raw.txt'
phi_data = np.loadtxt(filename, usecols=(0,1,2,3)).T
phi = phi_data[0]
#no need to atcually read dphi, since it is constant
#dphi = phi_data[1]
corr = phi_data[2]
dcorr = phi_data[3]

plt.figure(1, figsize=(8,6))
plt.scatter(phi,corr)
plt.xlabel(r'$\phi$')
plt.ylabel('Correlation Counts')
plt.savefig('phy424_qie_phi.png')
#plt.clf()

mean_corr = []
dmean_corr = []
phi_unique = []
a = phi[0]
sum_list = []
dcorr_sum_list = []

for i in range(len(phi)):
    if a == phi[i]:
        sum_list.append(corr[i])
        dcorr_sum_list.append(dcorr[i])
    elif a != phi[i]:
        mean_corr.append(np.mean(sum_list))
#        print(mean_corr)
        #calculate uncertainty of average
        dcorr_sum_array = np.array(dcorr_sum_list)
        quadrature = np.sqrt(np.sum(dcorr_sum_array**2)) #add in quadrature
        dmean_corr.append(quadrature/dcorr_sum_array.size)
        phi_unique.append(a)
        #re-initialize
        sum_list = []
        dcorr_sum_list = []
        a = phi[i]
        sum_list.append(corr[i])
        dcorr_sum_list.append(dcorr[i])
#do it one more time at the end
mean_corr.append(np.mean(sum_list))
phi_unique.append(a)
dcorr_sum_array = np.array(dcorr_sum_list)
quadrature = np.sqrt(np.sum(dcorr_sum_array**2)) #add in quadrature
dmean_corr.append(quadrature/dcorr_sum_array.size)

#since the error in phi is a constant of 0.5 degrees
dphi = 0.5*np.ones(shape=(len(phi_unique)))

#output
for i in range(len(phi_unique)):
    print("{} \t {} \t {} \t {}".format(round(phi_unique[i],2),round(dphi[i],1),
          round(mean_corr[i],2),round(dmean_corr[i],1)))

#def model(x, a, b, c, d):
#    #return a*np.sin(b*x + c) + d
#    return a*np.sin(b*x + c)**2 + d
#
##plt.plot(x, model(x,67,0.07,1,64))
##plt.plot(x, model(x, *popt))
#plt.xlabel(r'$\phi$')
##popt, pcov = curve_fit(model, theta_unique, mean_corr, p0=(67,0.07,1,64))
#popt, pcov = curve_fit(model, phi_unique, mean_corr, p0=(30,0.07,1.2,15))
#
#x = np.arange(0, 55, 0.5)
#
#plt.figure(2, figsize=(8,6))
#plt.scatter(phi_unique, mean_corr)
##plt.plot(x,model(x,30,0.07,1.2,15))
#plt.plot(x, model(x, *popt))
#plt.ylabel('Correlation Counts')#plt.savefig('phy424_qie_phi_mean.png')

