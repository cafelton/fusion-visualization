import sys, csv, subprocess, os, argparse
from datetime import date
from statistics import median
from collections import Counter
import copy
import numpy as np

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python 8-7-flair-to-fusions-pipe.py -i flair.aligned.bam -o outputPrefix -b buffer -a path to annotations')
parser.add_argument('-o', '--output', action='store', dest='o', default=date.today().strftime("%d-%m-%Y"), help='output file name base (default: date)')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='.bed file')
parser.add_argument('-b', '--buffer', action='store', dest='b', default=50000, help='length of buffer for combining nearby regions')
parser.add_argument('-a', '--anno', action='store', dest='a', default='/private/groups/brookslab/cafelton/id-fusions-new/gencodeGenesShort.gtf', help='path to anno.gtf')
parser.add_argument('-p', '--bedProcess', action='store_true', dest='p', help='whether to take .bam and convert to .bed and process (True = assume existing processed .bam)')
parser.add_argument('-s', '--samConvert', action='store_true', dest='s', help='whether to convert .bam to .sam or (True = assume existing .sam)')
parser.add_argument('-m', '--includeMito', action='store_true', dest='m', help='whether to include fusions that are in the mitochondria (True=include)')

args = parser.parse_args()
prefix = args.r.rstrip('.bam')
if not args.p:
	myCommands = ['bedtools bamtobed -bed12 -i ' + args.r + ' > ' + prefix + '-bedtools.bed',
				  'sort-bed ' + prefix + '-bedtools.bed' + ' > ' + prefix + '-bedtools-sorted.bed',
				  'python ' + os.path.dirname(os.path.realpath(__file__)) + '/standardizeBed.py' + ' -i ' + prefix + '-bedtools-sorted.bed',
				  'bedtools intersect -wao -a ' + prefix + '-bedtools-sorted-std.bed' + ' -b ' + args.a + ' > ' + prefix + '-bedtools-genes.txt',
				  'python ' + os.path.dirname(os.path.realpath(__file__)) + '/bedtoolsGeneHelper.py' + ' -i ' + prefix + '-bedtools-genes.txt',
				  'bedtools intersect -wao -a ' + prefix + '-bedtools-genes-short.bed' + ' -b ' + os.path.dirname(os.path.realpath(__file__)) + '/repeatMaskerRegions.bed > ' + prefix + '-bedtools-genes-repeats.txt',
				  'python ' + os.path.dirname(os.path.realpath(__file__)) + '/repeatsHelper.py' + ' -i ' + prefix + '-bedtools-genes-repeats.txt',
				  'rm ' + prefix + '-bedtools-genes.txt ' + prefix + '-bedtools-sorted.bed ' + prefix + '-bedtools-genes-short.bed ' + prefix + '-bedtools-genes-repeats.txt ' + prefix + '-bedtools.bed']
	process = subprocess.Popen('; '.join(myCommands),stdout=subprocess.PIPE, shell=True)
	proc_stdout = process.communicate()[0].strip()
	print(proc_stdout)
	print('done with preprocessing')

outfilename = args.o + prefix.split('/')[-1]
bed = open(prefix + '-bedtools-genes-repeats-short.bed', 'r')
reads = open(outfilename + "Reads.bed", "w")
fusions = open(outfilename + "Fusions.tsv", "w")
potential_chimeric = {}  # {read name: [entries]}
print("finding potential fusions \n")
bedLines = []
buffer = args.b
maxMapQ, bedLineCount, avgMapQ = 0, 0, 0
geneInfo = {}
c = 0
#GET POTEINTIAL CHIMERIC (MULTIPLE MAPPING) READS FROM BED FILE
for line in bed:
	if len(line) > 20:
		bedLineCount += 1
		line = line.rstrip().split('\t')
		bedLines.append(line)
		readname, gene  = line[3].split('--')
		if '/' not in gene:
			gene = gene.replace('.0', '')
			geneInfo[gene] = None
		else:
			geneLoc = gene.split('/')
			if geneLoc[0] not in geneInfo:
				geneInfo[geneLoc[0]] = geneLoc[1:]
			gene = geneLoc[0]
		avgMapQ += int(line[4])
		if int(line[4]) > maxMapQ: maxMapQ = int(line[4])
		if readname in potential_chimeric:
			if gene in potential_chimeric[readname]:
				potential_chimeric[readname][gene].append(line)
			else:
				potential_chimeric[readname][gene] = [line]
		else:
			potential_chimeric[readname] = {gene: [line]}

fusions_found = {}  # {fused genes: count}
fusionReads = []
avgMapQ = avgMapQ/float(bedLineCount)
print("filtering potential fusions \n" + str(len(potential_chimeric.keys())))
c = 0
for read in potential_chimeric:
	c += 1
	#if c % 10000 == 0: print(c)
	if len(potential_chimeric[read]) == 1:
		continue
	elif len(potential_chimeric[read]) > 1:
		locs = list(potential_chimeric[read].keys())
		locs.sort()
		fusion_name = '--'.join(locs)
		if fusion_name not in fusions_found:
			fusions_found[fusion_name] = {'mapScore':0, 'readNames':[], 'repeatScore':0}
			for loc in locs:
				fusions_found[fusion_name][loc] = {'reads':[], 'left':[], 'right':[], 'strand':[], 'chr':potential_chimeric[read][loc][0][0]}
		for loc in locs:
			fusions_found[fusion_name][loc]['reads'] += potential_chimeric[read][loc]
			for i in potential_chimeric[read][loc]:
				fusions_found[fusion_name]['mapScore'] += int(i[4])
				fusions_found[fusion_name]['repeatScore'] += float(i[-1])
				fusions_found[fusion_name][loc]['left'].append(int(i[1]))
				fusions_found[fusion_name][loc]['right'].append(int(i[2]))
				fusions_found[fusion_name][loc]['strand'].append(i[5])
		fusions_found[fusion_name]['readNames'].append(read)

#AGGREGATE AND SORT NON-GENIC REGIONS IN FUSIONS
print('condensing fusions in non-genic regions')
leftLocs, rightLocs = {}, {}
c = 0
for i in fusions_found:
	locs = i.split('--')
	if len(locs[0].split('-')) > 1 and locs[0][:3] == 'chr':
		chr, loc = locs[0].split('-')
		if chr not in leftLocs:
			leftLocs[chr] = []
		leftLocs[chr].append(loc)
	if len(locs[-1].split('-')) > 1 and locs[-1][:3] == 'chr':
		chr, loc = locs[-1].split('-')
		if chr not in rightLocs:
			rightLocs[chr] = []
		rightLocs[chr].append(loc)

#JOIN CLOSE REGIONS TOGETHER AND MARK THEM FOR UPDATING
updateValues = [[leftLocs, {}], [rightLocs, {}]]
for k in range(2):
	for i in updateValues[k][0]:
		lastLoc = '0'
		lastKey = None
		updateValues[k][0][i].sort()
		for j in updateValues[k][0][i]:
			if int(j)-int(lastLoc) < buffer and lastKey == None:
				lastKey = '-'.join([i, lastLoc])
				updateValues[k][1]['-'.join([i, lastLoc])] = lastKey
				updateValues[k][1]['-'.join([i, j])] = lastKey
			elif float(j)-float(lastLoc) < buffer:
				updateValues[k][1]['-'.join([i, j])] = lastKey
			else:
				updateValues[k][1]['-'.join([i, j])] = '-'.join([i, str(int(j))])
				lastKey = None
			lastLoc = j

#CONDENSE NON-GENIC FUSION REGIONS INTO FEWER FUSIONS
new_fusions_found = {}
for i in fusions_found:
	locs = i.split('--')
	if len(locs[0].split('-')) > 1 and locs[0][:3] == 'chr':
		locs[0] = updateValues[0][1][locs[0]]
	if len(locs[-1].split('-')) > 1 and locs[-1][:3] == 'chr':
		locs[-1] = updateValues[1][1][locs[-1]]
	if '--'.join(locs) not in new_fusions_found.keys():
		new_fusions_found['--'.join(locs)] = {}
		new_fusions_found['--'.join(locs)]['mapScore'] = fusions_found[i]['mapScore']
		new_fusions_found['--'.join(locs)]['repeatScore'] = fusions_found[i]['repeatScore']
		new_fusions_found['--'.join(locs)]['readNames'] = fusions_found[i]['readNames']
		for j in range(len(locs)):
			new_fusions_found['--'.join(locs)][locs[j]] = fusions_found[i][i.split('--')[j]]
	else:
		new_fusions_found['--'.join(locs)]['mapScore'] += fusions_found[i]['mapScore']
		new_fusions_found['--'.join(locs)]['repeatScore'] += fusions_found[i]['repeatScore']
		new_fusions_found['--'.join(locs)]['readNames'] += fusions_found[i]['readNames']
		for j in range(len(locs)):
			if locs[j] not in new_fusions_found['--'.join(locs)]:
				new_fusions_found['--'.join(locs)][locs[j]] = fusions_found[i][i.split('--')[j]]
			elif locs[j] not in ['mapScore', 'readNames', 'repeatScore']:
				for key in ['reads', 'left', 'right', 'strand']:
					new_fusions_found['--'.join(locs)][locs[j]][key] += fusions_found[i][i.split('--')[j]][key]

orgFusions = []
allMatches = []
readNames = {}
avgQualScore = 0
print('filtering fusions and detecting breakpoints')
c = 0
for i in new_fusions_found:
	supportCount = len(new_fusions_found[i]['readNames'])
	mapScore = round((new_fusions_found[i]['mapScore']/float(supportCount * len(i.split('--'))))/maxMapQ, 3) #* len(i.split('-')))
	repeatScore = round((new_fusions_found[i]['repeatScore']/float(supportCount * len(i.split('--')))), 3) #* len(i.split('-')))
	avgBreakpointAgg = 0
	#print(i, supportCount, mapScore)
	#if mapScore > 0.5: print(i, supportCount, mapScore)
	if supportCount >= 3 and (mapScore >= avgMapQ/float(maxMapQ) or mapScore > 0.9) and mapScore > .5 \
			and (args.m or 'chrM' not in i):
		currFusion = [i, str(supportCount), str(mapScore), str(repeatScore)]
		distTo5 = []
		locInfo = {}
		for loc in new_fusions_found[i]:
			#print(loc)
			if loc not in ['mapScore', 'readNames', 'repeatScore']:
				#if i == 'CCR7--MALAT1' and loc == 'CCR7'
				leftMed1, rightMed1 = int(median(new_fusions_found[i][loc]['left'])), int(median(new_fusions_found[i][loc]['right']))
				#leftMed = new_fusions_found[i][loc]['left'][min(range(len(new_fusions_found[i][loc]['left'])), key = lambda i: abs(new_fusions_found[i][loc]['left'][i]-leftMed1))]
				#rightMed = new_fusions_found[i][loc]['right'][min(range(len(new_fusions_found[i][loc]['right'])), key = lambda i: abs(new_fusions_found[i][loc]['right'][i]-rightMed1))]
				leftSide, rightSide = np.asarray(new_fusions_found[i][loc]['left']), np.asarray(new_fusions_found[i][loc]['right'])
				leftMed = leftSide[(np.abs(leftSide - leftMed1)).argmin()]
				rightMed = rightSide[(np.abs(rightSide - rightMed1)).argmin()]
				locInfo[loc] = {'chr':new_fusions_found[i][loc]['chr'].strip('chr'), 'loc':leftMed, 'seq':'NNNNNNNNNNNNNNNNNNNN', 'cseq':'--------------------'}
				leftFracClose, rightFracClose = 0, 0
				avgDistToLeft, avgDistToRight = 0, 0
				for j in new_fusions_found[i][loc]['left']:
					if abs(j-leftMed) <= 10: leftFracClose += 1
					avgDistToLeft += abs(j-leftMed)
				for j in new_fusions_found[i][loc]['right']:
					if abs(j-rightMed) <= 10: rightFracClose += 1
					avgDistToRight += abs(j-rightMed)
				leftFracClose, rightFracClose = leftFracClose/float(supportCount), rightFracClose/float(supportCount)
				avgDistToLeft, avgDistToRight = avgDistToLeft/float(supportCount), avgDistToRight/float(supportCount)
				#print(i, loc, leftFracClose, rightFracClose)
				if geneInfo[loc] != None:
					geneStart = int(geneInfo[loc][1]) if geneInfo[loc][0] == '+' else int(geneInfo[loc][2])
				else: geneStart=None
				#s = 'm' if Counter(new_fusions_found[i][loc]['strand']).most_common(1)[0][0] == '-' else 'p'
				if leftFracClose > rightFracClose or (leftFracClose == rightFracClose and avgDistToLeft < avgDistToRight):
					currFusion.append(loc + '-' + new_fusions_found[i][loc]['chr'] + '-' + str(leftMed))
					avgBreakpointAgg += leftFracClose
					locInfo[loc]['side'] = 'l'
					if geneStart: distTo5.append(abs(geneStart-rightMed))
					else: distTo5.append(None)
				else:
					currFusion.append(loc + '-' + new_fusions_found[i][loc]['chr'] + '-' + str(rightMed))
					avgBreakpointAgg += rightFracClose
					locInfo[loc]['side'] = 'r'
					if geneStart: distTo5.append(abs(geneStart-leftMed))
					else: distTo5.append(None)
		for j in new_fusions_found[i]['readNames']:
			readNames[j] = {'fusion':i, **copy.deepcopy(locInfo)}
		avgBreakpointAgg = avgBreakpointAgg/float(len(new_fusions_found[i].keys()))
		currFusion.append(mapScore*supportCount*avgBreakpointAgg*(1-repeatScore))
		avgQualScore += mapScore*supportCount*avgBreakpointAgg*(1-repeatScore)
		if None not in distTo5:
			#print(distTo5)
			if distTo5[0] < distTo5[1]:
				temp = currFusion[4]
				currFusion[4] = "3'-" + currFusion[5]
				currFusion[5] = "5'-" + temp
			else:
				currFusion[4] = "3'-" + currFusion[4]
				currFusion[5] = "5'-" + currFusion[5]
				currFusion[-1] *= -1
		if len(currFusion)>6:
			if i == "IGKV4-1--IGKC" or i == "IGKC--IGKV4-1": print(currFusion, (currFusion[4].split('-')[-2] == currFusion[5].split('-')[-2] and
														abs(int(currFusion[4].split('-')[-1])-int(currFusion[5].split('-')[-1])) < args.b),
																   abs(int(currFusion[4].split('-')[-1])-int(currFusion[5].split('-')[-1])),
																   args.b, currFusion[4].split('-')[-2])
			if not(currFusion[4].split('-')[-2] == currFusion[5].split('-')[-2] and
					abs(int(currFusion[4].split('-')[-1])-int(currFusion[5].split('-')[-1])) < args.b):
				orgFusions.append(currFusion)
				allMatches += new_fusions_found[i]['readNames']

if not args.s:
	print('identifying breakpoint sequence - making sam file')
	process = subprocess.Popen('samtools view -h -o ' + prefix + '.sam ' + args.r,stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())
sam = open(prefix + '.sam', 'r')
myReads = set(readNames.keys())
print('identifying breakpoint sequence - searching sam file')
for line in sam:
	line = line.rstrip().split('\t')
	if line[0] in myReads:
		finalKey = None
		for key in readNames[line[0]].keys():
			if key != 'fusion':
				#print(line)
				if readNames[line[0]][key]['chr'] == line[2].strip('chr'):
					if finalKey == None: finalKey = key
					elif abs(readNames[line[0]][key]['loc'] - int(line[3])) < \
							abs(readNames[line[0]][finalKey]['loc'] - int(line[3])): finalKey = key
		#if len(line[9]) == 1:print('short', line[9])
		cigar, thisSeq = line[5], line[9]
		currString = ""
		if len(thisSeq) >= 20 and len(cigar) >= 2:
			#print(line[0], readNames[line[0]], finalKey)
			if readNames[line[0]][finalKey]['side'] == 'l':
				while len(currString) < 20 and len(cigar) >= 2:
					i = 0
					while cigar[i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']:
						i += 1
					if cigar[i] in ['M', 'I', 'X', 'S']:
						if cigar[i] == 'M' or cigar[i] == 'X':
							currString += thisSeq[:int(cigar[:i])]
						thisSeq = thisSeq[int(cigar[:i]):]
					if cigar[i] == 'D' or cigar[i] == 'N':
						currString += '-' * int(cigar[:i])
					cigar = cigar[i+1:]
				currString = currString[:20]
				readNames[line[0]][finalKey]['seq'] = line[9][:20]
			elif readNames[line[0]][finalKey]['side'] == 'r':
				while len(currString) < 20 and len(cigar) >= 2:
					i = -1
					temp = cigar.rstrip('MDISHXPN')
					while abs(i) < len(temp) and temp[i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']:
						i -= 1
					if cigar[-1] in ['M', 'I', 'X', 'S']:
						if cigar[-1] == 'M' or cigar[-1] == 'X':
							currString = thisSeq[(-1*int(temp[i+1:])):] + currString
						thisSeq = thisSeq[:(-1*int(temp[i+1:]))]
					if cigar[-1] == 'D' or cigar[-1] == 'N':
						currString = '-' * int(temp[i+1:]) + currString
					cigar = cigar[:i]
				currString = currString[-20:]
				readNames[line[0]][finalKey]['seq'] = line[9][-20:]
			readNames[line[0]][finalKey]['cseq'] = currString
	#if line[0] in names: print(readNames[line[0]])
sam.close()
#for i in readNames: print(i, readNames[i])
fusionNameFlip = []
for i in range(len(orgFusions)):
	locScores = []
	for loc in orgFusions[i][0].split('--'):
		endSeq = [{'A':0,'C':0,'G':0,'T':0,'N':0, '*':0, '-':0} for x in range(20)]
		for read in new_fusions_found[orgFusions[i][0]]['readNames']:
			thisSeq = readNames[read][loc]['cseq']
			if len(thisSeq)==20:
				for j in range(20):
					endSeq[j][thisSeq[j]] += 1
		if len(thisSeq)==20:
			seqCalc = []
			for j in range(20):
				thisTot, thisMax = 0, 0
				for k in ['A', 'C', 'G', 'T', 'N', '*', '-']:
					thisTot += endSeq[j][k]
					if endSeq[j][k] > thisMax and k != 'N' and k != '*': thisMax = endSeq[j][k]
				seqCalc.append(float(thisMax)/thisTot)
			locScores.append(sum(seqCalc)/len(seqCalc))
		else:
			locScores.append(0)
	orgFusions[i].insert(3, str(round(sum(locScores)/len(locScores), 3)))
	if orgFusions[i][-1] < 0:
		fusionNameFlip.append(orgFusions[i][0])
		orgFusions[i][0] = '--'.join(orgFusions[i][0].split('--')[::-1])

print('fusions filtered')
if len(orgFusions) > 0:
	avgQualScore = avgQualScore/len(orgFusions)
orgFusions.sort(key=lambda x:abs(x[-1]), reverse=True)
fusions.write("#name\tspanning reads\tmapping score(1 is good)\tseq agreement near breakpoint (1 is good)\tfraction repetitive (0 is good)\t3' breakpoint\t5' breakpoint\n")
for i in orgFusions:
	#print(i)
	if abs(i[-1]) > avgQualScore*.01:
		fusions.write('\t'.join(i[:-1]) + '\n')
#print(avgQualScore, avgMapQ/float(maxMapQ))
print('fusions written')
printNames = open(outfilename + "readNames.txt", "w")
for readName in allMatches:
	printNames.write(readName + '\n')
printNames.close()
for line in bedLines:
	thisName = line[3].split('--')[0]
	if thisName in allMatches:
		if readNames[thisName]['fusion'] in fusionNameFlip:fName = '--'.join(readNames[thisName]['fusion'].split('--')[::-1])
		else: fName = readNames[thisName]['fusion']
		currName = line[3].split('--')[1]
		if currName[:3] == 'chr':
			#currName = currName.split('.')[0]
			for loc in readNames[thisName]:
				if loc[:3] == 'chr':
					if abs(float(currName.split('-')[1]) - float(loc.split('-')[1])) <= args.b * 3:
						currName = loc
		#if currName.split('/')[0] == "IGKV4-1": print(readNames[thisName])
		if currName.split('/')[0] in readNames[thisName].keys():
			thisSeq = readNames[thisName][currName.split('/')[0]]['cseq']
			line[3] = '-.-'.join([fName, line[3].split('--')[0], currName, thisSeq])
			reads.write("\t".join(line[:-1]) + "\n")
reads.close()
fusions.close()
print('reads written')

# subprocess.run(['perl', os.path.dirname(os.path.realpath(__file__)) + '/bed12ToGTF.pl'],
# 			   stdin=open(args.o + prefix + 'Reads.bed', 'r'), stdout=open(args.o + prefix + 'Reads.gtf', 'w'))
