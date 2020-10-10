import sys, csv, subprocess, os, argparse
from datetime import date
from collections import Counter

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python 8-7-flair-to-fusions-pipe.py -i flair.aligned.bam -o outputPrefix -b buffer -a path to annotations')
parser.add_argument('-o', '--output', action='store', dest='o', default=date.today().strftime("%d-%m-%Y"), help='output file name base (default: date)')
parser.add_argument('-b', '--fusionReads', action='store', dest='b', default="", help='.bed file')
parser.add_argument('-r', '--readMap', action='store', dest='r', default="", help='.txt file')
parser.add_argument('-i', '--isoforms', action='store', dest='i', default="", help='.bed file')
parser.add_argument('-f', '--fusions', action='store', dest='f', default="", help='.tsv file')
args = parser.parse_args()
prefix = args.o + args.b.split('Reads.')[0]
outfilename = args.o + prefix
readNames = []
readToInfo = {}
out = open(prefix + "IsoformReads.bed", 'w')
fusionsOut = open(prefix + 'IsoformFusions.tsv', 'w')
isoformSupport = {}
isoforms = []
readShortNames = []
allIsoformSupport = []
readsToIsoforms = {}
fusionsFound = {}
allReads = []
with open(args.r, 'r') as readMap:
	for line in readMap:
		line = line.rstrip().split('\t')
		isoformSupport[line[0]] = line[1].split(',')
		#allReads += line[1].split(',')
		for i in line[1].split(','):
			if i not in readsToIsoforms:
				readsToIsoforms[i] = []
			readsToIsoforms[i].append(line[0])
		isoforms.append(line[0])
		allIsoformSupport += line[1].split(',')
# countedReads = dict(Counter(allReads))
# for i in countedReads:
# 	if countedReads[i] > 1:
# 		print(i, countedReads[i])
# doubleMappedIsoforms = {}
# for i in readsToIsoforms:
# 	if len(readsToIsoforms[i]) > 1:
# 		for j in readsToIsoforms[i]:
# 			if j not in doubleMappedIsoforms: doubleMappedIsoforms[j] = {'reads':0, 'isos':[]}
# 			doubleMappedIsoforms[j]['reads'] += 1
# 			for k in readsToIsoforms[i]:
# 				if j != k and k not in doubleMappedIsoforms[j]['isos']: doubleMappedIsoforms[j]['isos'].append(k)
# dmIsos = set(doubleMappedIsoforms.keys())
with open(args.b, 'r') as reads:
	for line in reads:
		line = line.rstrip().split('\t')
		info = line[3].split('-.-')
		name = info[1] + '->' + str(round(int(line[1]), -5))
		readNames.append(name)
		readShortNames.append(info[1])
		readToInfo[name] = info
#print(dmIsos)
#print(list(set(allIsoformSupport) & set(readShortNames)))
with open(args.i, 'r') as theseIsoforms:
	for line in theseIsoforms:
		line = line.rstrip().split('\t')
		if line[3] in isoforms:
			count = 0
			currName = line[3]
			for read in isoformSupport[line[3]]:
				if read + '->' + str(round(int(line[1]), -5)) in readNames:
					if count == 0:
						currName = readToInfo[read + '->' + str(round(int(line[1]), -5))]
						if '-.-' in line[3]: line[3] = line[3].split('-.-')[1]
						currName[1] = line[3]
					if currName[0] not in fusionsFound: fusionsFound[currName[0]] = {'rs':0, 'isos':0, 'locs':[]}
					fusionsFound[currName[0]]['locs'].append(currName[2].split('/')[0])
					count += 1
			if count > 0:
				fusionsFound[currName[0]]['rs'] += count
				# if line[3] in dmIsos:
				# 	print(doubleMappedIsoforms[line[3]])
				# 	for i in doubleMappedIsoforms[line[3]]['isos']:
				# 		temp = currName.copy()
				# 		temp[1] = str(doubleMappedIsoforms[line[3]]['reads']) + '|' + currName[1][:int(len(currName)/2)] + '|' + i[:int(len(i)/2)]
				# 		line[3] = '-.-'.join(temp)
				# 		out.write('\t'.join(line) + '\n')
				# 		fusionsFound[currName[0]]['isos'] += 1
				# else:
				currName[1] = str(count) + '|' + currName[1]
				line[3] = "-.-".join(currName)
				out.write('\t'.join(line) + '\n')
				fusionsFound[currName[0]]['isos'] += 1
with open(args.f, 'r') as fusions:
	for line in fusions:
		if line[0] == '#': fusionsOut.write(line)
		line = line.rstrip().split('\t')
		if line[0] in fusionsFound:
			if len(list(dict.fromkeys(fusionsFound[line[0]]['locs']))) > 1:
				temp = list(dict.fromkeys(fusionsFound[line[0]]['locs']))
				for i in range(len(temp)):
					if temp[i][:3] == 'chr':
						temp[i] = '-'.join([temp[i].split('-')[0], str(round(float(temp[i].split('-')[-1]), -5))])
				if len(list(dict.fromkeys(temp))) > 1:
					line[1] = str(fusionsFound[line[0]]['rs']) + '/' + str(fusionsFound[line[0]]['isos'])
					print(line[0], line[1])
					fusionsOut.write('\t'.join(line) + '\n')
out.close()
fusionsOut.close()