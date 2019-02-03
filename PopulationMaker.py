
# This program will take a template file and flow cytometry files of the Unstimulated and 
# Stimulated FCS files listed therein. It will then calculate the three populations (True 
# Positives, Constitutively Low, and Constitutively High) that comprise this library. It will
# write an output text file that contains the proportion of each population, as well as histogram
# descriptions of all three populations. Finally, it will plot all histograms on a single set
# of axes. Plotting parameters may be adjusted as desired. 
#
# Script should be run in console directly as 'python SubtractHistograms.py'
#
# NOTE: The file SubtractHistogramTEMPLATE.txt must be in the same folder.
# NOTE: The user should adjust the Threshold parameter (line 69) to match their own data
#
# Template file 'SubtractHistogramTEMPLATE.txt' should be a plain text file containing:
# 1.  absolute path to directory containing the Unstimulated .fcs file
# 2.  the name of the curve (to be used in the legend)
# 3.  absolute path to directory containing the Stimulated .fcs file
# 4.  the name of the curve (to be used in the legend)
# 5.  the channel name (choose from Cy3-A, DAPI-A, PE-Texas Red-A, CFP-A, Alexa Fluor
#       488-A, etc. The channel names are listed in the FCS file header when it is
#       opened as a text file. No quotes. 
# 6.  the title of the plot, which also becomes the name of the output file

from __future__ import division
import numpy as np
import matplotlib.pylab as plt
from FlowCytometryTools import FCMeasurement, FCPlate, PolyGate

__author__ = 'Sarah Stainbrook'


# Parse template file.
with open('PopulationMakerTEMPLATE.txt','r') as f:
	NumFiles = 2
	FileLocs = [0 for x in range(NumFiles)]
	FileNames = [0 for x in range(NumFiles)]
	for file in range(NumFiles):
		FileLoc = f.readline()
		FileLocs[file] = FileLoc.strip()
		FileName = f.readline()
		FileNames[file] = FileName.strip()
	channel = f.readline()
	channel = str(channel.strip())
	title = f.readline()
	title = title.strip()

outfile = ''.join([title,'.txt'])

#Retrieve data from file and sort into bins
Frequencies = []
SetXvalues = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5,11.6,11.7,11.8,11.9,12.0,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13.0,13.1,13.2,13.3,13.4,13.5]
counter=0
for file in range(NumFiles):
	data = FCMeasurement(ID = 'Data', datafile = FileLocs[file])
	FP_data = data[channel]
	log_data = [np.log(x) for x in FP_data if x>1]
	scaledlog_data = [10*x for x in log_data]
	freqs = (np.bincount(scaledlog_data))/10               # one bin per every 0.1 log fluorescence unit
	scaledfreqs = [a/(len(scaledlog_data)/10) for a in freqs]  # scales the data by the total number of points such that files of different lengths can be compared.
	Frequencies.append(scaledfreqs)
	if len(scaledfreqs) < len(SetXvalues):
		for counter in range(0,(len(SetXvalues)-len(scaledfreqs))):
			scaledfreqs.append(0.0)
			counter += 1

IPTGarray=np.array(Frequencies[1])
MaxIPTGLocation = IPTGarray.argmax()

Threshold=6.6  #  adjust this cutoff based on the visual separation between populations

SubtractedHistograms = []
for item in range(0,len(SetXvalues)):
	Subtracted = Frequencies[1][item] - Frequencies[0][item]
	SubtractedHistograms.append(Subtracted)

TruePositives = []
for item in range(0,len(SetXvalues)):
	if SubtractedHistograms[item] >= 0:
		TruePositives.append(SubtractedHistograms[item])
	else:
		TruePositives.append(0)
ProportionTruePositives = sum(TruePositives) 

MovedTruePositives = []
for item in range(0,len(SetXvalues)):
	if SubtractedHistograms[item]<0:
		MovedTruePositives.append(SubtractedHistograms[item])
	else:
		MovedTruePositives.append(0)

ConstitutiveOff = []
for item in range(0,len(SetXvalues)):
	if SetXvalues[item]<=Threshold: 
		if MovedTruePositives[item] < 0:
			ConstitutiveOff.append(Frequencies[0][item]+MovedTruePositives[item])
		else:
			ConstitutiveOff.append(Frequencies[0][item])
	else:
		ConstitutiveOff.append(0)
ProportionConstitutiveOff = sum(ConstitutiveOff)

ConstitutiveOn = []
for item in range(0,len(SetXvalues)):
	if SetXvalues[item]>Threshold:
		ConstitutiveOn.append(Frequencies[0][item]+MovedTruePositives[item])
	else:
		ConstitutiveOn.append(0)
ProportionConstitutiveOn = sum(ConstitutiveOn)

OffTruePositives = [-1*x for x in MovedTruePositives]

## Write file containing population information
f = open(outfile, 'w')
f.write(str(len(SetXvalues)))
f.write('\n%.3f'%(ProportionTruePositives))
f.write('\n%.3f'%(ProportionConstitutiveOff))
f.write('\n%.3f'%(ProportionConstitutiveOn))
f.write('\n')
f.write(str(SetXvalues))
f.write('\n')
f.write(str(TruePositives))
f.write('\n')
f.write(str(ConstitutiveOff))
f.write('\n')
f.write(str(ConstitutiveOn))
f.write('\n')
f.write(str(OffTruePositives))
f.close()
fig = plt.figure(figsize=(5,4))

# Plot data points and fitted curve.
NoIPTG, = plt.plot(SetXvalues,Frequencies[0],lw=1,color='dimgray')#
IPTG, = plt.plot(SetXvalues,Frequencies[1],lw=1,color = 'darkorange')#marker='.',markersize=4

TruePosOff, = plt.plot(SetXvalues,MovedTruePositives,color = 'c',lw=1)
TruePos, = plt.plot(SetXvalues,TruePositives,color='royalblue',lw=1)
ConstOn, = plt.plot(SetXvalues,ConstitutiveOn,color='olivedrab',lw=1)
ConstOff, = plt.plot(SetXvalues,ConstitutiveOff,color='darkorchid',lw=1)

TruePosOffLabel = 'Unstimulated True Positives'
TruePosLabel = 'True Positives'
ConstOffLabel = 'Constitutively Low'
ConstOnLabel = 'Constitutively High'
#TruePosLabel = ''.join(['True Positives     Ratio: ','%.3f'%(ProportionTruePositives)])
#ConstOffLabel = ''.join(['\nConstitutively High   Ratio: ','%.3f'%(ProportionConstitutiveOff)])
#ConstOnLabel = ''.join(['\nConstitutively Low     Ratio: ','%.3f'%(ProportionConstitutiveOn)])

plt.plot((Threshold,Threshold),(-0.1,0.1),'k',linestyle=':',lw=1)
#plt.ylim(-0.001,0.1)

plt.xlim(2,10)
plt.xlabel('log GFP expression', size = 12)
plt.ylabel('Count', size = 12)
#plt.legend([NoIPTG, IPTG],['No IPTG','IPTG'],loc='upper left')
#plt.legend([NoIPTG, IPTG, TruePos, ConstOff, ConstOn],['No IPTG','IPTG',TruePosLabel,ConstOffLabel,ConstOnLabel])
#plt.legend([TruePos, ConstOff, ConstOn],[TruePosLabel,ConstOffLabel,ConstOnLabel])
plt.legend([NoIPTG, IPTG, TruePos, TruePosOff, ConstOff, ConstOn],['No IPTG','IPTG',TruePosLabel,TruePosOffLabel, ConstOffLabel,ConstOnLabel])
plt.title(title, size = 20)
plt.tight_layout()
plt.show()