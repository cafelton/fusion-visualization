import sys, csv, subprocess, os, argparse
from datetime import date

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python 8-7-flair-to-fusions-pipe.py -i flair.aligned.bam -o outputPrefix -b buffer -a path to annotations')
parser.add_argument('-o', '--output', action='store', dest='o', default=date.today().strftime("%d-%m-%Y"), help='output file name base (default: date)')
parser.add_argument('-f', '--fastq', action='store', dest='f', default="", help='.fastq file')
parser.add_argument('-s', '--chrSize', action='store', dest='s', default=100, type=int, help='size of area to remap')
parser.add_argument('-l', '--fusions', action='store', dest='l', default="", help='.tsv file')
parser.add_argument('-r', '--readNames', action='store', dest='r', default="", help='.txt file')
parser.add_argument('-a', '--anno', action='store', dest='a', default='/private/groups/brookslab/reference_sequence/hg38.nobreaks.simple.fa', help='path to genome.fa')
#parser.add_argument('-p', '--bedProcess', action='store_true', dest='p', help='whether to filter fasta reads')
args = parser.parse_args()
prefix = args.r.split('readNames.')[0]
outfilename = args.o + prefix
readNames = []
faPrefix = args.o + '-' + ".".join(args.f.split('/')[-1].split('.')[:-1])
with open(args.r) as names:
	for line in names:
		readNames.append(line.rstrip())
#Filter reads to only the double mapped reads
with open(args.f, 'r') as reads, open(faPrefix + "Filtered.fa", "w") as faOut:
	writeRead = False
	for line in reads:
		if line[0] == '@':
			if line.rstrip().lstrip('@').split()[0] in readNames:
				writeRead = True
				faOut.write(">" + line.lstrip('@'))
			else:
				writeRead = False
		elif writeRead:
			faOut.write(line)
			writeRead = False
print('reads filtered')
fusions = {}
firstLine = []
with open(args.l, 'r') as thesefusions, open(faPrefix + 'Locs.bed', 'w') as bedFile:
	for line in thesefusions:
		line = line.rstrip().split('\t')
		if line[0][0] != '#' and len(line) > 6:
			chr1, center1 = line[5].split('-')[-2:]
			chr2, center2 = line[6].split('-')[-2:]
			if int(center1) > args.s and int(center2) > args.s:
				name1 = line[0] + '->' + '-'.join(line[5].lstrip("3'-").split('-')[:-2])
				name2 = line[0] + '->' + '-'.join(line[6].lstrip("5'-").split('-')[:-2])
				fusions[line[0]] = {'line':line}
				fusions[line[0]]['-'.join(line[5].lstrip("3'-").split('-')[:-2])] = {'side':"3'", 'chr':chr1, 'bp':int(center1), 'left':[], 'right':[], 'reads':[], 'mapQ':[]}
				fusions[line[0]]['-'.join(line[6].lstrip("5'-").split('-')[:-2])] = {'side':"5'", 'chr':chr2, 'bp':int(center2), 'left':[], 'right':[], 'reads':[], 'mapQ':[]}
				bedFile.write('\t'.join([chr1, str(int(center1)-args.s), str(int(center1) + args.s), name1]) + '\n')
				bedFile.write('\t'.join([chr2, str(int(center2)-args.s), str(int(center2) + args.s), name2]) + '\n')
		elif line[0][0] == '#':
			firstLine = line
process = subprocess.Popen('bedtools getfasta -fi ' + args.a + ' -bed ' + faPrefix + 'Locs.bed' + ' -fo ' + faPrefix + 'Genome.fa' + ' -name; ' +
						   'minimap2 -a ' + faPrefix + 'Genome.fa ' + faPrefix + "Filtered.fa" + ' > ' + faPrefix + 'Remapped.sam; ' +
						   "sam2bed < " + faPrefix + 'Remapped.sam' + ' > ' + faPrefix + 'Remapped.bed',stdout=subprocess.PIPE, shell=True)
print(process.communicate()[0].strip())
maxMapQ = 0
with open(faPrefix + 'Remapped.bed', 'r') as remapped:
	for line in remapped:
		line = line.rstrip().split('\t')
		if int(line[4]) > maxMapQ: maxMapQ = int(line[4])
		fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['mapQ'].append(int(line[4]))
		fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['reads'].append(line[3])
		fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['left'].append(int(line[1])-args.s)
		fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['right'].append(int(line[2])-args.s)
newFusions = open(faPrefix + 'FusionsRemapped.tsv', 'w')
firstLine.insert(2, 'confirmed reads')
newFusions.write('\t'.join(firstLine) + '\n')
doubleMappedReads = []
for fusion in fusions:
	good = 0
	fusionReads = []
	mapScores = []
	for loc in fusions[fusion]:
		if loc not in ['line', 'conf reads']:
			if len(fusions[fusion][loc]['left']) > 0:
				leftAvg, rightAvg = sum(fusions[fusion][loc]['left'])/len(fusions[fusion][loc]['left']), \
									sum(fusions[fusion][loc]['right'])/len(fusions[fusion][loc]['right'])
				if abs(leftAvg + rightAvg) > 50:
					good += 1
					# if abs(leftAvg) > abs(rightAvg):
					# 	fusions[fusion][loc]['bp'] += max(fusions[fusion][loc]['right'])
					# else:
					# 	fusions[fusion][loc]['bp'] += min(fusions[fusion][loc]['left'])
					fusionReads.append(fusions[fusion][loc]['reads'])
					mapScores += fusions[fusion][loc]['mapQ']
					if fusions[fusion][loc]['side'] == "3'":
						fusions[fusion]['line'][5] = '-'.join(fusions[fusion]['line'][5].split('-')[:-1] + [str(fusions[fusion][loc]['bp'])])
					else:
						fusions[fusion]['line'][6] = '-'.join(fusions[fusion]['line'][6].split('-')[:-1] + [str(fusions[fusion][loc]['bp'])])
	if good >= 2:
		shared = list(set(fusionReads[0]) & set(fusionReads[1]))
		if len(shared) > 0:#good >= 2:
			for i in shared: doubleMappedReads.append(fusion + i)
			fusions[fusion]['line'].insert(2, str(len(shared)))
			fusions[fusion]['line'][3] = str((sum(mapScores)/len(mapScores))/float(maxMapQ))
			newFusions.write('\t'.join(fusions[fusion]['line']) + '\n')
newFusions.close()
doubleMappedReads = set(doubleMappedReads)
with open(faPrefix + 'Remapped.bed', 'r') as remapped, open(faPrefix + 'Remapped-1.bed', 'w') as remapFilt:
	for line in remapped:
		line = line.rstrip().split('\t')
		if line[0].split('->')[0] + line[3] in doubleMappedReads:
			remapFilt.write('\t'.join(line) + '\n')

process = subprocess.Popen('python ' + os.path.dirname(os.path.realpath(__file__)) + '/makeAlnSeq.py -f ' + faPrefix +
						   'Genome.fa -r ' + faPrefix + 'Remapped-1.bed',stdout=subprocess.PIPE, shell=True)
print(process.communicate()[0].strip())