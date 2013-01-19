import numpy as np;
import scipy as sp;
import math;
from ctypes import *;
from collections import *;
from ImpliedVol import *;
import matplotlib.pyplot as plt;

# class FDResult(Structure):
# 	_fields_ = [ ("S", c_float), ("Payoff", c_float), ("V", c_float) ]

FDResult = namedtuple('FDResult', ['S', 'Payoff', 'V']);
OptionContract = namedtuple('OptionContract', ['Type', 'Strike', 'Expiration' ])
OptionPrice = namedtuple('OptionPrice', ["ContractName", "FDGrid", "FDPrice", "BSPrice"]);

def printFDResult(inputArray):
	for item in inputArray:
		print "%f, %f, %f" % (item.S, item.Payoff, item.V)

def Heaviside( s, k):
	if s>k:
		return 1.0
	return 0.0

def findUniformNAS( strike, ds=1 ):
	NAS = (2 * float(strike)) / ds;
	return int(NAS);

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
			Payoff.append(Heaviside(S[i], strike) * quantity)
		elif optionType == "digital put":
			Payoff.append(Heaviside(-S[i], -strike) * quantity)


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

def drawDiagram( Result ):
	#prepare to draw
	maxXAxis = min([len(x.FDGrid) for x in Result]) #so we get uniform length
	xAxis = [ x.S for x in Result[0].FDGrid[:maxXAxis] ] 
	for i in range(0,len(Result)):
		if i == 0:
			yAxisV = np.array([ x.V for x in Result[i].FDGrid[:maxXAxis] ])
			yAxisPayoff = np.array([ x.Payoff for x in Result[i].FDGrid[:maxXAxis] ])
		else:
			yAxisV += np.array([ x.V for x in Result[i].FDGrid[:maxXAxis] ])
			yAxisPayoff +=np.array([ x.Payoff for x in Result[i].FDGrid[:maxXAxis] ])

	TotalPremium = sum([x.FDPrice for x in Result])
	yAxisV += TotalPremium
	yAxisPayoff += TotalPremium
	#Draw 
	plt.plot(xAxis,yAxisV, label='V');
	plt.plot(xAxis,yAxisPayoff, label='Payoff')
	plt.legend(loc='best');
	plt.title('Options Price and Payoff diagram');
	plt.show();


assetPrice = 100
volLow = 0.2
volHigh = 0.2
r = 0.02
strike = 100 #at the money
NAS = 200
expiration = 1;

LongDC = { "Contract": OptionContract("digital call", 100, expiration), "Quantity": 10, "Name": "LongDC" };
LongCall = { "Contract": OptionContract( "call", 90, expiration), "Quantity":0.008126, "Name": "LongCall" };
ShortCall = { "Contract": OptionContract("call", 110, expiration), "Quantity": -1, "Name": "ShortCall" };
Portfolio = [ 
	LongDC,
	LongCall,
	ShortCall,
]

Result = []

for item in Portfolio:
	asset = item["Contract"]
	q = item["Quantity"]
	NAS = findUniformNAS( asset.Strike )
	priceFD  = OptionsFDPricer( volHigh, volLow, r, asset.Type, asset.Strike, asset.Expiration, q, NAS)
	ResultPrice = OptionPrice( item["Name"],priceFD, 0,0)
	price = findClosest( priceFD, assetPrice )
	if price != None:
	 	print price.V
	 	ResultPrice = ResultPrice._replace(FDPrice=price.V)
	else:
	 	print "can't find price"
	 	
	Result.append(ResultPrice)

	# if asset.Type == "call":
	# 	BS = EuropeanCall( assetPrice, volHigh, 0, r, asset.Strike, asset.Expiration)
	# elif asset.Type == "digital call":
	# 	BS = DigitalCall( assetPrice, volHigh, 0, r, asset.Strike, asset.Expiration)

drawDiagram(Result)


