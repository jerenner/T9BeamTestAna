#This is a simple plotting code to fit the sigmaE/E against E to verify the proportionality constant. 

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

def linear(x, A, B):
    return x * A + B

def fitFunction(x, y, function):
    # Fit the Gaussian function to the data points
    x = np.array(x)
    y = np.array(y)
    params, covariance = curve_fit(linear, x, y, p0=[1, 0])
    return params, covariance

#using dict is so much better 
##############################
run432_el = {"momentum": 460, 
            "E": 6.249,
            "sigmaE": 0.654,
            "n" : 94629}

run432_mu = {"momentum": 460, 
            "E": 3.305,
            "sigmaE": 0.374}
##############################
run412_el = {"momentum": 560, 
            "E": 7.907,
            "sigmaE": 0.828,
            "n" : 95452}

run412_mu = {"momentum": 560, 
            "E": 3.474,
            "sigmaE": 0.383}
##############################
run438_el = {"momentum": 700, 
            "E": 10.558,
            "sigmaE": 1.107,
            "n": 74117}

run438_mu = {"momentum": 700, 
            "E": 3.806,
            "sigmaE": 0.427}
##############################
run398_el = {"momentum": 900, 
            "E": 13.959,
            "sigmaE": 1.518, 
            "n":55931}

run398_mu = {"momentum": 900, 
            "E": 3.901,
            "sigmaE": 0.554}
##############################
run393_el = {"momentum": 1000, 
            "E": 15.890,
            "sigmaE": 1.763,
            "n": 75116}

#careful, muon fit at 1000 is very poor, sigmaE is overestimated
run393_mu = {"momentum": 1000, 
            "E": 3.994,
            "sigmaE": 0.670}

run_numbers = [432, 412, 438, 398, 393]
all_electron_dict = [run432_el, run412_el,  run438_el,  run398_el, run393_el]
all_muon_dict = [run432_mu, run412_mu,  run438_mu,  run398_mu, run393_mu]


plt.figure(figsize = [10, 6])
for i, dict in enumerate(all_electron_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    plt.errorbar(E, (E/sigmaE)**2, yerr = sigmaE, fmt ='x', markersize = 10, label = "Run %i - p = %i MeV/c"%(run_numbers[i], dict["momentum"]))

plt.grid()
plt.xlabel("LG fitted mean E (a.u)", fontsize = 18)
plt.ylabel(r"LG fitted $(E/\sigma_E)^2$ (a.u)", fontsize = 18)
plt.title("LG electron WindowIntPE energy resolution", fontsize = 18)
plt.legend(fontsize = 15)
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()

plt.savefig("../%s/ElectronLGenergyRelationCheck.pdf"%("pdf_2plus2Coincidence_nosumTScut"))
plt.savefig("../%s/ElectronLGenergyRelationCheck.png"%("png_2plus2Coincidence_nosumTScut"))


###########################################################
plt.figure(figsize = [10, 6])
for i, dict in enumerate(all_muon_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    plt.errorbar(E, (E/sigmaE)**2, yerr = sigmaE, fmt ='x', markersize = 10, label = "Run %i - p = %i MeV/c"%(run_numbers[i], dict["momentum"]))

plt.grid()
plt.xlabel("LG fitted mean E (a.u)", fontsize = 18)
plt.ylabel(r"LG fitted $(E/\sigma_E)^2$ (a.u)", fontsize = 18)
plt.title("LG muon WindowIntPE energy resolution", fontsize = 18)
plt.legend(fontsize = 15)
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()

plt.savefig("../%s/MuonLGenergyRelationCheck.pdf"%("pdf_2plus2Coincidence_nosumTScut"))
plt.savefig("../%s/MuonLGenergyRelationCheck.png"%("png_2plus2Coincidence_nosumTScut"))

########################################

########################################

array_x = []
array_y = []
plt.figure(figsize = [10, 6])
for i, dict in enumerate(all_electron_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    mom = dict["momentum"]

    array_x.append(mom)
    array_y.append(sigmaE)
    
    
    plt.errorbar(mom, sigmaE, yerr = 0, fmt ='x', markersize = 10, label = "Run %i - p = %i MeV/c"%(run_numbers[i], dict["momentum"]))

print(array_x, array_y)
array_x = np.array(array_x)
array_y = np.array(array_y)
params, covariance = fitFunction(array_x, array_y, linear)
x = np.linspace(min(array_x)-10, max(array_x) + 10, 100)
plt.plot(x, linear(x, params[0], params[1]), 'k--', label = 'Linear Fit: %.2e x + %.2f'%(params[0], params[1]))


plt.grid()
plt.ylabel("LG fitted std E (a.u)", fontsize = 18)
plt.xlabel("Nominal beam momentum (MeV/c)", fontsize = 18)
plt.title("LG electron WindowIntPE std of energy against beam momentum", fontsize = 18)
plt.legend(fontsize = 15)
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()

plt.savefig("../%s/ElectronLGSigmaMomentumRelationCheck.pdf"%("pdf_2plus2Coincidence_nosumTScut"))
plt.savefig("../%s/ElectronLGSIgmaMomentumRelationCheck.png"%("png_2plus2Coincidence_nosumTScut"))
#######################################################
plt.figure(figsize = [10, 6])
for i, dict in enumerate(all_electron_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    mom = dict["momentum"]
    plt.errorbar(mom, (E/sigmaE)**2, yerr = 0, fmt ='x', markersize = 10, label = "Run %i - p = %i MeV/c"%(run_numbers[i], dict["momentum"]))

plt.grid()
plt.ylabel(r"LG fitted $(E/\sigma_E)^2$ (a.u)", fontsize = 18)
plt.xlabel("Nominal beam momentum (MeV/c)", fontsize = 18)
plt.title("LG electron WindowIntPE square of inverse energy resolution \n against beam momentum", fontsize = 18)
plt.legend(fontsize = 15)
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()

plt.savefig("../%s/ElectronInverseResolutionLGMomentumRelationCheck.pdf"%("pdf_2plus2Coincidence_nosumTScut"))
plt.savefig("../%s/ElectronInverseResolutionLGMomentumRelationCheck.png"%("png_2plus2Coincidence_nosumTScut"))

########################################

array_x = []
array_y = []
fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize = [15, 10])
for i, dict in enumerate(all_electron_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    mom = dict["momentum"]

    array_x.append(mom)
    array_y.append(E)
    ax0.errorbar(mom, E, yerr = sigmaE, fmt ='x', markersize = 10, label = "Run %i - p = %i MeV/c"%(run_numbers[i], dict["momentum"]))

print(array_x, array_y)
array_x = np.array(array_x)
array_y = np.array(array_y)
params, covariance = fitFunction(array_x, array_y, linear)
x = np.linspace(min(array_x)-10, max(array_x) + 10, 100)
ax0.plot(x, linear(x, params[0], params[1]), 'k--', label = 'Linear Fit: %.2e x + %.2f'%(params[0], params[1]))
for i, dict in enumerate(all_electron_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    mom = dict["momentum"]
    ax1.errorbar(mom, (linear(mom, params[0], params[1])-E)/E, yerr = sigmaE/E, fmt ='x', markersize = 10)

ax1.plot([min(array_x)-10, max(array_x)+10], [0, 0], "k--", label = r"Total abs(fit-E)/E: %.2e"%(sum((abs(linear(array_x, params[0], params[1])-array_y)/array_y))))



ax0.grid()
ax0.set_ylabel("LG fitted mean E (a.u)", fontsize = 18)
ax0.set_xlabel("Nominal beam momentum (MeV/c)", fontsize = 18)
fig.suptitle("LG electron WindowIntPE mean +/- std energy against beam momentum", fontsize = 20, weight = "bold")

#residuals
# ax1.grid()
ax1.set_ylabel("Fit fractional residuals", fontsize = 18)
ax1.set_xlabel("Nominal beam momentum (MeV/c)", fontsize = 18)
plt.grid()
ax1.legend()

ax0.legend(fontsize = 15)
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()

fig.savefig("../%s/ElectronLGMomentumRelationCheck.pdf"%("pdf_2plus2Coincidence_nosumTScut"))
fig.savefig("../%s/ElectronLGMomentumRelationCheck.png"%("png_2plus2Coincidence_nosumTScut"))
#################################################
plt.figure(figsize = [10, 6])
for i, dict in enumerate(all_electron_dict):
    E = dict["E"]
    sigmaE = dict["sigmaE"]
    mom = dict["momentum"]
    plt.errorbar(mom, sigmaE/E, yerr = 0, fmt ='x', markersize = 10, label = "Run %i - p = %i MeV/c"%(run_numbers[i], dict["momentum"]))

plt.grid()
plt.ylabel(r"LG fitted $(\sigma_E/E)$ (a.u)", fontsize = 18)
plt.xlabel("Nominal beam momentum (MeV/c)", fontsize = 18)
plt.title("LG electron WindowIntPE fractional energy \n resolution against beam momentum", fontsize = 18)
plt.legend(fontsize = 15)
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()

plt.savefig("../%s/ElectronesolutionLGMomentumRelationCheck.pdf"%("pdf_2plus2Coincidence_nosumTScut"))
plt.savefig("../%s/ElectronResolutionLGMomentumRelationCheck.png"%("png_2plus2Coincidence_nosumTScut"))

