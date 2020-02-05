import numpy as np
import matplotlib.pyplot as plt

filename = 'phy424_qie_theta_data_raw.txt'
theta_data = np.loadtxt(filename, usecols=(0,1,2,3)).T
theta = theta_data[0]
#no need to atcually read dphi, since it is constant
#dphi = phi_data[1]
corr = theta_data[2]
dcorr = theta_data[3]

plt.figure(1, figsize=(8,6))
plt.scatter(theta,corr)
plt.xlabel(r'$\theta$')
plt.ylabel('Correlation Counts')
plt.savefig('phy424_qie_theta.png')
plt.clf()

mean_corr = []
dmean_corr = []
theta_unique = []
a = theta[0]
sum_list = []
dcorr_sum_list = []

for i in range(len(theta)):
    if a == theta[i]:
        sum_list.append(corr[i])
        dcorr_sum_list.append(dcorr[i])
    elif a != theta[i]:
        mean_corr.append(np.mean(sum_list))
#        print(mean_corr)
        #calculate uncertainty of average
        dcorr_sum_array = np.array(dcorr_sum_list)
        quadrature = np.sqrt(np.sum(dcorr_sum_array**2)) #add in quadrature
        dmean_corr.append(quadrature/dcorr_sum_array.size)
        theta_unique.append(a)
        #re-initialize
        sum_list = []
        dcorr_sum_list = []
        a = theta[i]
        sum_list.append(corr[i])
        dcorr_sum_list.append(dcorr[i])
#do it one more time at the end
mean_corr.append(np.mean(sum_list))
theta_unique.append(a)
dcorr_sum_array = np.array(dcorr_sum_list)
quadrature = np.sqrt(np.sum(dcorr_sum_array**2)) #add in quadrature
dmean_corr.append(quadrature/dcorr_sum_array.size)

#since the error in phi is a constant of 0.5 degrees
dtheta = 0.5*np.ones(shape=(len(theta_unique)))

#output
for i in range(len(theta_unique)):
    print("{} \t {} \t {} \t {}".format(round(theta_unique[i],2),
          round(dtheta[i],1),
          round(mean_corr[i],2),round(dmean_corr[i],1)))
