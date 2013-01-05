import csv;
import scipy as sp;
import numpy as np;
import numpy.linalg as la;
import math;
import random;
import matplotlib.pyplot as plt;
import decimal;

#Trapezoid Integration rule for definite integral
#RV = dTau/2 * (f(x0)+2f(x1)+2f(x2)+f(xn))
#dTau = (b-a)/N
#xi = a+i(dTau)
def trapezoidIntegration(f, param, a, b, N):
	if b==0:
		return 0;
	dTau = (b-a)/N;
	M = 0.5 * f(param,a);
	for i in xrange(1, N):
		M += f(param, a + i*dTau);
	M = M + (0.5 * f(param, b));
	M = M * dTau;
	return M;

#param is a polynomial regression result (array)
#tau is an array of time column
#result is an array
def volMusiela(param, tau):
	if len(param) == 4:
		rv = param[3] + (tau*param[2]) + (tau*tau*param[1]) + (tau*tau*tau*param[0])
	if len(param)==1:
			#if( type(tau) == "int"):
			rv = 0.6431e-2; #param[0];
			#else:
			#	rv = np.array([0.6431e-2]);
			#	rv = np.resize(rv,(len(tau),))
	return rv;

def multiFactorMuMusiela(param1, param2, param3, tau):
	M1 = np.array([]);
	M2 = np.array([]);
	M3 = np.array([]);

	for i in xrange(0,len(tau)):
		tempM1 = trapezoidIntegration(volMusiela, param1, 0, tau[i],100);
		M1 = np.append( M1, tempM1 * volMusiela(param1, tau[i]) );
		
		tempM2 = trapezoidIntegration(volMusiela, param2, 0, tau[i],100);
		M2 = np.append( M2, tempM2 * volMusiela(param2, tau[i]) );
		
		tempM3 = trapezoidIntegration(volMusiela, param3, 0, tau[i],100);
		M3 = np.append( M3, tempM3 * volMusiela(param3, tau[i]) );

	M = M1+M2+M3;
	return M;

#initialRow is a matrix
def generateForwardCurve(initialRow, dt, tenor, header, mu, vol1, vol2, vol3):
	#assert that header and tenor is the same length
	i = 0 + dt; #because index 0 is already the initialRow
	prevRow = initialRow;
	maxI = tenor; #rows
	maxJ = initialRow.shape[0]; #columns
	F = initialRow;
	while i <= maxI:
		j=0;
		row = np.array([]);
		dx1 = random.normalvariate(0,1);
		dx2 = random.normalvariate(0,1);
		dx3 = random.normalvariate(0,1);

		while j < maxJ:
			if( j == maxJ - 1):
				dF = prevRow[j] - prevRow[j-1];
				dTau = header[j] - header[j-1];
			else:
				dF = prevRow[j+1] - prevRow[j];
				dTau = header[j+1] - header[j];

			newVal = prevRow[j] + (mu[j]*dt) + ( ( (vol1*dx1) + (vol2[j]*dx2) + (vol3[j]*dx3) ) * math.sqrt(dt) ) + ((dF/dTau)*dt);
			row = np.append(row, newVal);
			j += 1;
		
		prevRow = row;
		F = np.vstack((F, row));	
		i += dt;
	return F;

def generateRowLabel( dt, tenor):
	rv = np.array([]);
	i = 0;
	while i <= tenor:
		rv = np.append( rv, i );
		i += dt;
	return rv;

def linearSearch( inputArray, val ):
	for i in xrange(0,len(inputArray)):
		if abs(inputArray[i]-val) < 1e-7:
			return i;
	return(-1);

def ZCBRate( F, dt, startTenor, endTenor, rowLabels ):
	startIndex = linearSearch( rowLabels, startTenor );
	endIndex = linearSearch( rowLabels, endTenor );
	if endIndex ==-1 or startIndex == -1:
		print "can't find libor rates. Invalid Tenor!";
		return -1;
	total = sum(F[:,0][startIndex:endIndex]) * dt;
	return total;

def ZCBPricer( F, dt, startTenor, endTenor, rowLabels ):
	rate = ZCBRate( F, dt, startTenor, endTenor, rowLabels );
	if rate != -1:
		return math.exp(-rate);
	else:
		print "can't price ZCB"
		return -1;

def LiborRate( F, dt, startTenor, endTenor, rowLabels):
	startIndex = linearSearch( rowLabels, startTenor );
	endIndex = linearSearch( rowLabels, endTenor );
	rate = ZCBRate( F, dt, startTenor, endTenor, rowLabels)/dt;
	avg = rate/(endIndex-startIndex);
	return avg;

def capletsPricer( F, dt, tenor, rowLabels, strike ):
	i = 0+dt;
	totalCap = 0;
	while i<=tenor:
		ZCBPrice = ZCBPricer( F, dt, 0, i, rowLabels );
		LR = LiborRate( F, dt, i-dt, i, rowLabels );
		caplet = ZCBPrice * max(LR-strike,0) * dt;
		totalCap += caplet;
		i = i+dt;
	return totalCap;

def floorletsPricer( F, dt, tenor, rowLabels, strike ):
	i = 0+dt;
	totalCap = 0;
	while i<=tenor:
		ZCBPrice = ZCBPricer( F, dt, 0, i, rowLabels );
		LR = LiborRate( F, dt, i-dt, i, rowLabels );
		caplet = ZCBPrice * max(strike-LR,0) * dt;
		totalCap += caplet;
		i = i+dt;
	return totalCap;


myArray = [];  
CQFExampleHeaders = [  
0.0,     0.5,     1.0,     1.5,     2.0,
2.5,     3.0,     3.5,     4.0,     4.5,     5.0,     5.5,     6.0,     6.5,
7.0,     7.5,     8.0,     8.5,     9.0,     9.5,     10.0,    10.5,    11.0,
11.5,    12.0,    12.5,    13.0,    13.5,    14.0,    14.5,    15.0,    15.5,
16.0,    16.5,    17.0,    17.5,    18.0,    18.5,    19.0,    19.5,    20.0,
20.5,    21.0,    21.5,    22.0,    22.5,    23.0,    23.5,    24.0,    24.5,
25.0, ];

Headers = [ 
0.08,    0.17,    0.25,    0.33,    0.42,    0.50,    0.58,
0.67,    0.75,    0.83,    0.92,    1.00,    1.08,    1.17,    1.25,    1.33,
1.42,    1.50,    1.58,    1.67,    1.75,    1.83,    1.92,    2.00,    2.08,
2.17,    2.25,    2.33,    2.42,    2.50,    2.58,    2.67,    2.75,    2.83,
2.92,    3.00,    3.08,    3.17,    3.25,    3.33,    3.42,    3.50,    3.58,
3.67,    3.75,    3.83,    3.92,    4.00,    4.08,    4.17,    4.25,    4.33,
4.42,    4.50,    4.58,    4.67,    4.75,    4.83,    4.92,    5.00,    5.5,
6.0,     6.5,     7.0,     7.5,     8.0,     8.5,     9.0,     9.5,     10.0,
10.5,    11.0,    11.5,    12.0,    12.5,    13.0,    13.5,    14.0,    14.5,
15.0,    15.5,    16.0,    16.5,    17.0,    17.5,    18.0,    18.5,    19.0,
19.5,    20.0,    20.5,    21.0,    21.5,    22.0,    22.5,    23.0,    23.5,
24.0,    24.5,    25.0, ];

filename = "CQFExample.csv";
#filename = 'UKYieldCurveSmall.csv';
#filename = 'UKYieldCurveAll.csv';

print "Loading csv file";
with open(filename, 'rb') as f:
	reader = csv.reader(f)
 	for row in reader:
 		if len(row[1]) > 0: #skip bad row (empty row)
 			myRow = [];
 			for item in row:
 				myRow.append(float(item))
 			myArray.append(myRow);

myMatrix = sp.matrix(myArray);

print "Loading csv file completed"

#create difference matrix
print "creating difference matrix";
diffArray = [];
rowSize = myMatrix.shape[0];
columnSize = myMatrix.shape[1];
for i in range(0,rowSize-1):
	diffRow = [];
	for j in range(0, columnSize):
		diffRow.append(myArray[i+1][j]-myArray[i][j]);
	diffArray.append(diffRow);

myDiffMatrix = sp.matrix(diffArray);


#create covariance matrix
print "creating covariance matrix";
myCovMatrix = np.cov(myDiffMatrix,rowvar=0,bias=1);
myCovMatrix = myCovMatrix*252/10000; #normalize for 252 days per year

#Get eigenvector and eigenvalue
print "computing eigenvalues and eigenvectors";
myEigenValues, myEigenVectors = la.eigh(myCovMatrix); #input is symmetric matrix

l = len(myEigenValues);
#the sign is backward
PCAOne = myEigenVectors[:,l-1] * -1;
PCATwo = myEigenVectors[:,l-2] * -1;
PCAThree = myEigenVectors[:,l-3] ;

print "computing lambda1, lambda2, lambda3";
LambdaOne = myEigenValues[l-1];
LambdaTwo = myEigenValues[l-2];
LambdaThree = myEigenValues[l-3];

LambdaOneSq = math.sqrt(LambdaOne);
LambdaTwoSq = math.sqrt(LambdaTwo);
LambdaThreeSq = math.sqrt(LambdaThree);

print "computing vol1, vol2, vol3"
volOne = PCAOne * LambdaOneSq;
volTwo = PCATwo * LambdaTwoSq;
volThree = PCAThree * LambdaThreeSq;

print "computing linear regression factor";
if "CQF" in filename:
	H = CQFExampleHeaders
else:
	H = Headers;

la1 = np.polyfit(H, volOne, 0, full=True);
la2 = np.polyfit(H, volTwo, 3, full=True);
la3 = np.polyfit(H, volThree, 3, full=True)

print "computing 3 factors vol (musiela)";
vol1Musiela = volMusiela(la1[0],np.array(H));
vol2Musiela = volMusiela(la2[0],np.array(H));
vol3Musiela = volMusiela(la3[0],np.array(H));

print "computing mu (musiela)";
mu = multiFactorMuMusiela(la1[0],la2[0],la3[0],np.array(H));

#test only
bla1 = trapezoidIntegration(volMusiela,la1[0],0,0,100);
bla1 = bla1 * volMusiela(la1[0],0.5);
bla2 = trapezoidIntegration(volMusiela,la2[0],0,0,100);
bla2 = bla2 * volMusiela(la2[0],0.5);
bla3 = trapezoidIntegration(volMusiela,la3[0],0,0,100);
bla3 = bla3 * volMusiela(la3[0],0.5);

initialRow = np.array(myArray[len(myArray)-1]) * 1e-2; #get the last row as seed data (data is in percentage)

print "Please enter tenor to price ZCB and CapFloors (in year, for example, 1.5 for 1year 6months): ";
inputTenor = float(raw_input().split()[0]);

print "Please enter strike for CapFloor (in percent, for example: 5 for 5%): ";
inputStrike = float(raw_input().split()[0])/100;

dt = 0.01; #time interval
tenor = 10; #generate every dt until 10yr tenor
rowLabels = generateRowLabel( dt, tenor);
MCPaths = 1;
ZCBPrice = 0;
LiborR = 0;
Caplets = 0;
Floorlets = 0;

#MonteCarlo Simulations
for i in xrange(0, MCPaths):
	print " MC: %d" % i;
	F = generateForwardCurve( initialRow, dt, tenor, H, mu, vol1Musiela, vol2Musiela, vol3Musiela);
	

	# #series by column
	# plt.plot(rowLabels,F[:,0], label='spot');
	# #plt.show();

	# #series by row
	# plt.plot(H, F[0], label='0'); #0 year
	# plt.plot(H,F[200], label='2'); #2 year
	# plt.plot(H,F[500], label='5'); #5 year
	# plt.plot(H,F[700], label='7'); #7 year
	# plt.plot(H,F[1000], label='10'); #10 year
	# plt.legend(loc='best');
	# #plt.show();

	temp1 = ZCBPricer( F, dt, 0, inputTenor, rowLabels );
	ZCBPrice += temp1;
	print "ZCBPrice = %f" % temp1;

	temp2 = LiborRate( F, dt, 0, inputTenor, rowLabels );
	LiborR += temp2;
	print "LiborRate = %f" % temp2;

	Caplets += capletsPricer( F, dt, inputTenor, rowLabels, inputStrike );
	Floorlets += floorletsPricer( F, dt, inputTenor, rowLabels, inputStrike );

ZCBPrice = ZCBPrice / MCPaths;
LiborR = LiborR / MCPaths;
Cap = Caplets / MCPaths;
Floor = Floorlets / MCPaths;

print "=========================="
print "======FINAL RESULT========"
print "MCPaths= %d" % MCPaths;
print "inputTenor= %f" % inputTenor;
print "inputStrike= %f" % inputStrike;
print "ZCBPrice= %f" % ZCBPrice;
print "LiborR= %f" % LiborR;
print "Cap= %f" % Cap;
print "Floor= %f" % Floor;
print "=========================="


