import numpy as np;
import scipy as sp;
import math;
from ctypes import *;
from collections import *;
from ImpliedVol import *;

# class FDResult(Structure):
# 	_fields_ = [ ("S", c_float), ("Payoff", c_float), ("V", c_float) ]

FDResult = namedtuple('FDResult', ['S', 'Payoff', 'V']);
OptionContract = namedtuple('OptionContract', ['Type', 'Strike', 'Expiration' ])

def printFDResult(inputArray):
	for item in inputArray:
		print "%f, %f, %f" % (item.S, item.Payoff, item.V)

def Heaviside( s, k):
	if s>k:
		return 1.0
	return 0.0

#NAS: Number of Asset Steps
#NTS: Number of Time Steps
def OptionsFDPricer( volHigh, volLow, r, optionType, strike, expiration, quantity, NAS ):
	ds = 2 * float( strike ) / NAS ;
	dt = 0.9 / volHigh / volHigh / NAS / NAS;
	NTS = int(expiration/dt) + 1;
	dt = float(expiration)/NTS;

	S = [];
	Payoff = [];
	VOld = [];
	VNew = [0]*NAS; #initialize with zero
	RV = [];

	for i in range(0,NAS):
		S.append(i * ds);
		if optionType == "call":
			Payoff.append(max(S[i]-strike,0) * quantity)
		elif optionType == "put":
			Payoff.append(max(strike-S[i],0) * quantity)
		elif optionType == "digital call":
			Payoff.append(Heaviside(S[i], strike))
		elif optionType == "digital put":
			Payoff.append(Heaviside(-S[i], -strike))


		VOld.append(float(Payoff[i]));
		tempRV = FDResult(S[i],Payoff[i],0);
		RV.append(tempRV);

	print "generating FD grid"
	for k in range(0, NTS):
		print "K=================%d" % k
		for i in range(1, NAS-1):

			delta = (VOld[i+1] - VOld[i-1])/2/ds; #central difference
			gamma = (VOld[i+1] - (2 * VOld[i]) + VOld[i-1]) / ds / ds;

			if gamma > 0:
				vol = volLow
			else:
				vol = volHigh

			temp = vol*S[i];
			theta = ( r * VOld[i] ) - (0.5* math.pow(temp, 2) * gamma ) - ( r * S[i] * delta )

			VNew[i] = VOld[i] - theta * dt;

		VNew[0] = VOld[0] * (1-r*dt);
		VNew[NAS-1] = 2*VNew[NAS-2] - VNew[NAS-3];
		
		VOld = VNew[:];

	for i in range(0, NAS):
		RV[i] = RV[i]._replace(V = VOld[i]);

	return RV;

def findClosest( inputFDResult, targetPrice):
	ds = abs(inputFDResult[0].S - inputFDResult[1].S);
	for i in range (0,len(inputFDResult)-1):
		if inputFDResult[i].S == targetPrice:
			return inputFDResult[i]
		elif (abs(inputFDResult[i].S - targetPrice) <= ds) and (inputFDResult[i+1].S-inputFDResult[i].S) <= ds:
			return FDResult((inputFDResult[i].S + inputFDResult[i+1].S)/2, (inputFDResult[i].Payoff + inputFDResult[i+1].Payoff)/2, (inputFDResult[i].V + inputFDResult[i+1].V)/2)

	return None;

assetPrice = 101
volHigh = 0.3
volLow = 0.3
r = 0.05
strike = 100 #at the money
NAS = 100
quantity =1
optionType = "call";
expiration = 1;


LongCall = { "Contract": OptionContract( "call", 100, expiration), "Quantity":1, "Name": "LongCall" };
LongPut = { "Contract": OptionContract("put", 100, expiration), "Quantity": 1, "Name": "LongPut" };
LongDC = { "Contract": OptionContract("digital call", 100, expiration), "Quantity": 1, "Name": "LongDC" };
Portfolio = [ 
	#LongCall,
	#LongPut,
	LongDC,
]



for item in Portfolio:
	asset = item["Contract"];
	q = item["Quantity"];
	priceFD  = OptionsFDPricer( volHigh, volLow, r, asset.Type, asset.Strike, asset.Expiration, q, NAS);
	price = findClosest( priceFD, assetPrice );
	if price != None:
		print price.V;


	# if asset.Type == "put" or asset.TYpe == "call":
	# 	EC = EuropeanPut( assetPrice, volHigh, 0, r, asset.Strike, asset.Expiration);

