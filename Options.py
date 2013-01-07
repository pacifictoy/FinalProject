import numpy as np;
import scipy as sp;
import math;
from ctypes import *;

class FDResult(Structure):
	_fields_ = [ ("S", c_float), ("Payoff", c_float), ("V", c_float) ]

def printFDResult(inputArray):
	for item in inputArray:
		print "%f, %f, %f" % (item.S, item.Payoff, item.V)

#NAS: Number of Asset Steps
#NTS: Number of Time Steps
def OptionsFDPricer( volHigh, volLow, r, optionType, strike, expiration, quantity, NAS ):
	ds = 2 * strike / NAS;
	dt = 0.9 / volHigh / volHigh / NAS / NAS;
	NTS = int(expiration/NAS) + 1;
	dt = expiration/NTS;

	S = np.array([]);
	Payoff = np.array([]);
	VOld = np.array([]);
	VNew = np.empty(NAS);  VNew.fill(0);
	RV = [];

	for i in range(0,NAS):
		S = np.append(S, [i * ds]);
		if optionType == "call":
			Payoff = np.append(Payoff, max(S[i]-strike,0) * quantity)
		elif optionType == "put":
			Payoff = np.append(Payoff, max(strike-S[i],0) * quantity)

		VOld = np.append(VOld,Payoff[i]);
		tempRV = FDResult(S[i],Payoff[i],0);
		RV.append(tempRV);

	for k in range(1, NTS):
		for i in range(1, NAS):

			delta = (VOld(i+1) - VOld(i-1))/2/ds; #central difference
			gamma = (VOld(i+1) - (2 * VOld(i)) + VOld(i-1)) / ds / ds;

			if gamma > 0:
				vol = volLow
			else:
				vol = volHigh

			theta = r * VOld(i) - math((0.5*vol*S[i]), 2) * gamma - r * S[i] * delta

			VNew[i] = VOld[i] - theta * dt;

		VNew[0] = VOld[0] * (1-r*dt);
		VNew[NAS] = 2*VNew[NAS-1] - VNew[NAS-2];

		for i in range(0,NAS):
			VOld[i] = VNew[i];

	for i in range(0, NAS):
		RV[i].V = VOld[i];

	return RV;



volHigh = 0.4
volLow = 0.1
r = 0.05
strike = 110
assetPrice = 100
NAS = 200
quantity =1
optionType = "call";
expiration = 1;

Portfolio = [

]

Result = OptionsFDPricer( volHigh, volLow, r, optionType, strike, expiration, quantity, NAS);
