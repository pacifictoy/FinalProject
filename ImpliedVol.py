from scipy import stats
import math

def Heaviside( s, k):
	if s>k:
		return 1.0
	return 0.0

def dOne(asset, vol, divyld, intrate, strike, expiry):
	d1 = (math.log(float(asset)/float(strike)) + ((intrate - divyld + (0.5 * vol * vol)) * expiry)) / (float(vol) * math.sqrt(expiry))
	return d1	

def dTwo(asset, vol, divyld, intrate, strike, expiry):
	d2 = dOne(asset, vol, divyld, intrate, strike, expiry) - (vol * math.sqrt(expiry))
	return d2

def N1p(asset, vol, divyld, intrate, strike, expiry):
	N1p = stats.norm.cdf(dOne(asset, vol, divyld, intrate, strike, expiry))
	return N1p

def N2p(asset, vol, divyld, intrate, strike, expiry):
	N2p = stats.norm.cdf(dTwo(asset, vol, divyld, intrate, strike, expiry))
	return N2p

def EuropeanCall(asset, vol, divyld, intrate, strike, expiry):
	Call = (asset * math.exp(-divyld*expiry) * N1p(asset, vol, divyld, intrate, strike, expiry)) - (strike * math.exp(-intrate * expiry) * N2p(asset, vol, divyld, intrate, strike, expiry))
	return Call

def DigitalCall(asset, vol, divyld, intrate, strike, expiry):
	Payout = Heaviside( asset, strike )
	RV = Payout * N2p(asset, vol, divyld, intrate, strike, expiry) * math.exp(-intrate*expiry)
	return RV

def EuropeanPut(asset, vol, divyld, intrate, strike, expiry):
	Put = (strike * math.exp(-intrate*expiry) * N2p(asset, vol, divyld, intrate, strike, expiry)) - (asset * math.exp(-divyield * expiry) * N1p(asset, vol, divyld, intrate, strike, expiry))
	return Put

#use bruteforce method, accurate until 2nd decimal point (very poor accuracy)
def ImpliedVolCallBruteForce( optionprice, asset, divyld, intrate, strike, expiry, tolerance = 0.001):
	vol = 1;
	while vol>0: 
		call = EuropeanCall(asset, vol, divyld, intrate, strike, expiry)
		diff = optionprice - call
		if abs(diff) <= tolerance:
			return vol
		vol=vol-(tolerance/10)
	return -1 #error

def ImpliedVolCallBisection( optionprice, asset, divyld, intrate, strike, expiry, start, end, tolerance=0.001):
	vol = float(start+end)/2.0
	while count<10000:
		call = EuropeanCall(asset, vol, divyld, intrate, strike, expiry)
		diff = optionprice - call
		if abs(diff) <= tolerance:
				return vol
		

