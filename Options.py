import numpy as np;
import scipy as sp;
import math;
from ctypes import *;
from collections import *;

# class FDResult(Structure):
# 	_fields_ = [ ("S", c_float), ("Payoff", c_float), ("V", c_float) ]

FDResult = namedtuple('FDResult', ['S', 'Payoff', 'V']);
OptionContract = namedtuple('OptionContract', ['Type', 'Strike', 'Expiration' ])

def printFDResult(inputArray):
	for item in inputArray:
		print "%f, %f, %f" % (item.S, item.Payoff, item.V)

#NAS: Number of Asset Steps
#NTS: Number of Time Steps
def OptionsFDPricer( volHigh, volLow, r, optionType, strike, expiration, quantity, NAS ):
	ds = 2 * strike / NAS ;
	dt = 0.9 / volHigh / volHigh / NAS / NAS;
	NTS = int(expiration/dt) + 1;
	dt = float(expiration)/NTS;

	S = [];
	Payoff = [];
	VOld = [];
	VNew = [0]*NAS; 
	RV = [];

	for i in range(0,NAS):
		S.append(i * ds);
		if optionType == "call":
			Payoff.append(max(S[i]-strike,0) * quantity)
		elif optionType == "put":
			Payoff.append(max(strike-S[i],0) * quantity)

		VOld.append(Payoff[i]);
		tempRV = FDResult(S[i],Payoff[i],0);
		RV.append(tempRV);

	for k in range(0, NTS):
		print "K================== %d" % k
		for i in range(1, NAS-1):
			#print i
			delta = (VOld[i+1] - VOld[i-1])/2/ds; #central difference
			gamma = (VOld[i+1] - (2 * VOld[i]) + VOld[i-1]) / ds / ds;

			if gamma > 0:
				vol = volLow
			else:
				vol = volHigh

			temp = 0.5*vol*S[i];
			theta = ( r * VOld[i] ) - ( math.pow(temp, 2) * gamma ) - ( r * S[i] * delta )

			VNew[i] = VOld[i] - theta * dt;

		VNew[0] = VOld[0] * (1-r*dt);
		VNew[NAS-1] = 2*VNew[NAS-2] - VNew[NAS-3];
		
		VOld = VNew[:];

	for i in range(0, NAS):
		RV[i] = RV[i]._replace(V = VOld[i]);

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

LongCall = { "Contract": OptionContract( "call", 90, 1), "Quantity":1 };
LongPut = { "Contract": OptionContract("put", 100, 1), "Quantity": 1 };
Portfolio = [ 
	LongCall,
	LongPut,
]

Result = OptionsFDPricer( volHigh, volLow, r, optionType, strike, expiration, quantity, NAS);
