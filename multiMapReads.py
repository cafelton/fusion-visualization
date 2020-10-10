import sys, csv, subprocess, os, argparse
from datetime import date

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python 8-7-flair-to-fusions-pipe.py -i flair.aligned.bam -o outputPrefix -b buffer -a path to annotations')
parser.add_argument('-o', '--output', action='store', dest='o', default=date.today().strftime("%d-%m-%Y"), help='output file name base (default: date)')
parser.add_argument('-f', '--fastq', action='store', dest='f', default="", help='.fastq file')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='Remapp-seq.bed file')
args = parser.parse_args()
readNames = []
faPrefix = args.o + '-' + ".".join(args.f.split('/')[-1].split('.')[:-1])
path = os.getcwd() + "/" + faPrefix + "-temp"
if not os.path.exists(path):
	os.mkdir(path)
readDict = {}
with open(args.r) as names:
	for line in names:
		line = line.rstrip().split('\t')
		name = line[3]
		if line[3] in readNames:
			name = line[3] + '-' + str(len(readDict[line[3]]))
		else:
			readDict[line[3]] = []
			readNames.append(line[3])
		readDict[line[3]].append('>' + line[0] + '->' + name + '\n' + line[11] + '\n')
for uniqueRead in readDict:
	with open(path + "/" + uniqueRead + 'Reads.fa', 'w') as out:
		for i in readDict[uniqueRead]:
			out.write(i)
print('bed processed')
print(len(readNames))
readNames = set(readNames)
#Filter reads to only the double mapped reads
readsAddedToRef=  []
with open(args.f, 'r') as reads:
	faOut = None
	writeRead = False
	for line in reads:
		if line[0] == '@':
			lineName = line.rstrip().lstrip('@').split()[0]
			if lineName in readNames: #and lineName not in readsAddedToRef:
				writeRead = True
				faOut = open(path + "/" + lineName + "REF.fa", 'w')
				readsAddedToRef.append(lineName)
				faOut.write(">REF" + lineName + '\n')
		elif writeRead:
			faOut.write(line)
			faOut.close()
			writeRead = False
print('reads filtered')
c = 0
print(len(readsAddedToRef))
for i in readsAddedToRef:
	c += 1
	if i in readNames:
		process = subprocess.Popen('minimap2 -a ' + path + "/" + i + "REF.fa " + path + "/" + i + 'Reads.fa > ' + path + '/' + i + 'Mapped.sam; ' +
							   "sam2bed < " + path + '/' + i + 'Mapped.sam > ' + path + '/' + i + 'Mapped.bed; ' +
							   'python ' + os.path.dirname(os.path.realpath(__file__)) + '/makeAlnSeq.py -f ' + path + '/' +
								i + 'REF.fa -r ' + path + '/' + i + 'Mapped.bed',stdout=subprocess.PIPE, shell=True)
		print(c)
print('done')
process = subprocess.Popen('cat ' + path + '/*Mapp-seq.bed > ' + faPrefix + 'MappedToFasta.bed; rm -r ' + path,stdout=subprocess.PIPE, shell=True)
print(process.communicate()[0].strip())