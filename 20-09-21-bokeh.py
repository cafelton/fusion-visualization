from bokeh.models import ColumnDataSource, TableColumn, StringFormatter, NumberFormatter, DataTable, CustomJS, FixedTicker, \
	PrintfTickFormatter, HBar, LinearAxis, Range1d, FuncTickFormatter, PanTool,ResetTool,BoxZoomTool, CustomJSFilter, \
	CDSView, Span, BasicTickFormatter, SingleIntervalTicker, Dropdown, Text, Div, NumeralTickFormatter
from bokeh.layouts import column, row, gridplot
from bokeh.models.glyphs import Rect, Text
from bokeh.themes import Theme
from bokeh.colors import RGB
from bokeh.plotting import curdoc, figure
from collections import OrderedDict
from math import pi
import pandas as pd
import requests
from bokeh.models.widgets import FileInput
from pybase64 import b64decode
import io
import numpy as np

theme = Theme(json={
	'attrs': {
		'Plot': {
			'background_fill_color': 'white',
			'outline_line_color': None },
		'Axis': {
			'axis_line_color': None,
			'major_tick_line_color': 'grey',
			'minor_tick_line_color': None},
		'Grid': { 'grid_line_color': "#dddddd"}}})
curdoc().theme = theme
#Code framework for clicking on tables from https://stackoverflow.com/questions/55964945/chart-on-click-selection-from-data-table-in-bokeh/55970697#55970697
fusions = pd.read_csv('02-09-2020DRR-flair.alignedFusions.tsv', sep='\t', index_col=False)
reads = pd.read_csv('02-09-2020DRR-flair.alignedReads.bed', sep='\t', index_col=False,
					names=['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb',
						   'blockCount','blockSizes','blockStarts'])

currFusion = 0

def getReadOrder(myDF, fusionName, currChromPoints):
	seqLen = len(myDF.at[0, 'seq'])
	myDF = myDF.set_index('id')
	readsLeft = myDF.loc[(myDF['fusion']==fusionName) & (myDF['gene'] == currChromPoints[0][0]), :].copy()
	readsRight = myDF.loc[(myDF['fusion']==fusionName) & (myDF['gene'] == currChromPoints[1][0]), :].copy()
	readsRight = readsRight[~readsRight.index.duplicated(keep='first')]
	readsLeftOnly = readsLeft.loc[~readsLeft.index.isin(list(readsRight.index))].copy()
	readsRightOnly = readsRight.loc[~readsRight.index.isin(list(readsLeft.index))].copy()
	readsBoth = readsLeft.loc[readsLeft.index.isin(list(readsRight.index))].copy()
	readsBoth = readsBoth.rename(columns={'seq':'lseq'})
	temp = len(readsBoth)
	readsBoth['rseq'] = readsRight['seq']
	readsBoth=readsBoth.head(temp)
	readsRightOnly = readsRightOnly.head(25)
	readsLeftOnly = readsLeftOnly.head(25)
	readsBoth = readsBoth.head(50)
	readIDs = list(readsRightOnly.index) + list(readsLeftOnly.index) + list(readsBoth.index)
	leftSeq = [["-"]*seqLen for i in range(len(readsRightOnly))] + list(readsLeftOnly['seq']) + list(readsBoth['lseq'])
	rightSeq = list(readsRightOnly['seq']) + [["-"]*seqLen for i in range(len(readsLeftOnly))] + list(readsBoth['rseq'])
	return readIDs, leftSeq, rightSeq

def makeFilteredData(fusion_name, currChromPoints, reads_file):
	if isinstance(reads_file, pd.DataFrame):
		myReads = reads_file
	else: myReads = reads
	#temp = [x.split('-.-') for x in myReads['name']]#myReads['name'].str.split('-.-', expand=True)
	myReads[['fusionID','readID','geneName', 'seq']] = pd.DataFrame([x.split('-.-')[:4] for x in myReads['name']], index= myReads.index)#temp[temp.columns[0:4]]
	myReads['geneName'] = myReads['geneName'].str.split('/', expand = True)[0]
	#if loc != None: readsFiltered = myReads.loc[reads['fusionID']==fusions.at[loc, '#name'], :].copy()
	readsFiltered = myReads.loc[myReads['fusionID']==fusion_name, :].copy()
	#currChromPoints = [fusions.loc[fusions['#name']==fusion_name, "5' breakpoint"].item().split("-")[1:], fusions.loc[fusions['#name']==fusion_name, "3' breakpoint"].item().split("-")[1:]]
	# if currChromPoints[0][0][:3] == 'chr':
	# 	readsFiltered['chartLoc'] = readsFiltered['geneName'].split('-') == currChromPoints[0][0]
	# else:
	readsFiltered['chartLoc'] = readsFiltered['geneName'] == currChromPoints[0][0]
	readsFiltered['chartLoc'] = np.where(readsFiltered['chartLoc'], 'left', 'right')
	leftRF, rightRF = readsFiltered[readsFiltered['chartLoc'] == 'left'].fillna(method='bfill', axis=0), \
					  readsFiltered[readsFiltered['chartLoc'] == 'right'].fillna(method='bfill', axis=0)
	if len(leftRF) > 0 and len(rightRF) > 0:
		numForMean = int(len(leftRF)*0.05) if int(len(leftRF)*0.05) > 3 else 3
		if len(leftRF) <= 3:
			leftBounds, rightBounds = [min(list(leftRF['chromStart'])), max(list(leftRF['chromEnd']))], \
									  [min(list(rightRF['chromStart'])), max(list(rightRF['chromEnd']))]
		else:
			leftBounds, rightBounds = [int(np.partition(leftRF['chromStart'].to_numpy(), numForMean)[:numForMean].mean()),
									   int(np.partition(leftRF['chromEnd'].to_numpy(), (-1* numForMean))[(-1*numForMean):].mean())], \
									  [int(np.partition(rightRF['chromStart'].to_numpy(), numForMean)[:numForMean].mean()),
									   int(np.partition(rightRF['chromEnd'].to_numpy(), (-1* numForMean))[(-1*numForMean):].mean())]
		leftStart, rightStart = leftBounds[1]-30, rightBounds[0] - 10
		leftFlip, rightFlip = False, False
		currChrom = [currChromPoints[0][1], currChromPoints[1][1]]
		currPoints = [int(currChromPoints[0][2]), int(currChromPoints[1][2])]
		a = [abs(i-int(currChromPoints[0][2])) for i in leftBounds]
		if a[0] < a[1]:
			leftBounds = leftBounds[::-1]
			leftFlip = True
			leftStart = leftBounds[1] - 10
		a = [abs(i-int(currChromPoints[1][2])) for i in rightBounds]
		if a[1] < a[0]:
			rightBounds = rightBounds[::-1]
			rightFlip = True
			rightStart = rightBounds[0] - 30
		leftSeq, rightSeq = requests.get('http://togows.org/api/ucsc/hg38/chr' + str(currChromPoints[0][1]) + ':' + str(leftStart+1) + '-' + str(leftStart + 41) + '.fasta').text, \
							requests.get('http://togows.org/api/ucsc/hg38/chr' + str(currChromPoints[1][1]) + ':' + str(rightStart+1) + '-' + str(rightStart + 41) + '.fasta').text
		leftShortSeq, rightShortSeq = leftSeq.strip().split('\n')[1], rightSeq.strip().split('\n')[1]
		if leftFlip: leftShortSeq = leftShortSeq[::-1]
		if rightFlip: rightShortSeq = rightShortSeq[::-1]
		genomeRow = pd.DataFrame({'readID':['hg38 sequence', 'hg38 sequence'], 'chromStart': [leftStart, rightStart], 'chromEnd':[leftStart + 40, rightStart + 40], 'seq':[leftShortSeq, rightShortSeq]})
		tickSpaceL, tickSpaceR = int(abs(leftBounds[0] - leftBounds[1])/6.), int(abs(rightBounds[0] - rightBounds[1])/6.)
		if leftFlip: leftBounds[0] += tickSpaceL
		else: leftBounds[0] -= tickSpaceL
		if rightFlip: rightBounds[1] -= tickSpaceR
		else: rightBounds[1] += tickSpaceR
		readsFiltered = readsFiltered[readsFiltered['readID'].isin(list(set(readsFiltered['readID']))[:200])]
		sizes = readsFiltered['blockSizes'].str.split(',', expand=True).stack().str.strip().reset_index(level=1, drop=True)
		starts = readsFiltered['blockStarts'].str.split(',', expand=True).stack().str.strip().reset_index(level=1, drop=True)
		temp = pd.concat([starts,sizes], axis=1, keys=['starts','sizes'])
		readsExpanded = readsFiltered.join(temp).reset_index(drop=True)
		readsExpanded[['starts', 'sizes']] = readsExpanded[['starts', 'sizes']].apply(pd.to_numeric)
		readsExpanded['tStart'] = readsExpanded['starts'] + readsExpanded['chromStart']
		readsExpanded['tEnd'] = readsExpanded['sizes'] + readsExpanded['tStart']
		#print(readsExpanded.head())
		return readsExpanded, currChrom, currPoints, leftBounds, rightBounds, leftFlip, rightFlip, tickSpaceL, tickSpaceR, genomeRow
	else:
		print(len(readsFiltered), fusion_name)#, set(myReads['fusionID']))
		return readsFiltered, [], [], [], [], False, False, 0, 0, 0

def view_alignment(ids, seqs, chromPoints, side, fontsize="9pt", plot_width=800):
	"""Bokeh sequence alignment view"""
	text = [i for s in list(seqs) for i in s]
	clrs = {'T':RGB(153, 204, 153), 'A':RGB(255, 153, 153), 'G':RGB(255, 219, 153), 'C':RGB(153, 153, 255), '-':'white'}
	colors = [clrs[i] for i in text]
	N = len(seqs[0])/2
	center = int(chromPoints[2])
	x = np.arange(center-N,center+N)
	y = np.arange(0,len(seqs),1)
	#creates a 2D grid of coords from the 1D arrays
	xx, yy = np.meshgrid(x, y)
	#flattens the arrays
	gx = xx.ravel()
	gy = yy.flatten()
	#now we can create the ColumnDataSource with all the arrays
	#source = ColumnDataSource(dict(x=gx+0.5, y=gy, recty=gy+0.5, text=text, colors=colors))
	source = ColumnDataSource(dict(x=gx+0.5, y=gy, recty=gy+0.5, rectx1=gx, rectx2=gx+1, text=text, colors=colors))
	plot_height = len(seqs)*20+90
	view_range = (center-20, center+20)
	p1 = figure(title=' '.join(chromPoints), plot_width=plot_width, plot_height=plot_height,
				x_range=view_range, y_range=ids, tools="xpan,reset",
				min_border=0, toolbar_location='below')#, lod_factor=1)
	glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
				 text_font="monospace",text_font_size=fontsize)
	# p1.rect(x="x", y="recty",  width=0.9, height=1, fill_color="colors",
	# 			 line_color=None, fill_alpha=1)
	p1.segment(x0='rectx1', y0='recty', x1='rectx2',
			   y1='recty', color="colors", line_width=14, source=source)
	#p1.add_glyph(source, rects)
	p1.add_glyph(source, glyph)
	#print('break line at', chromPoints[2])
	breakLine = Span(location=int(chromPoints[2]), dimension='height', line_color='red', line_width=2)
	p1.renderers.extend([breakLine])
	p1.grid.visible = False
	p1.xaxis.major_label_text_font_style = "bold"
	p1.yaxis.minor_tick_line_width = 0
	p1.yaxis.major_tick_line_width = 0
	p1.below[0].formatter.use_scientific = False
	p1.xaxis.major_label_orientation = pi/4
	if side == 'left': p1.name = 'pl'
	else:
		p1.name = 'pr'
		p1.min_border_left = 10
		p1.yaxis.visible=False
	return p1

def makeFullPlot(fusion_name, currChromPoints, reads_file):
	readsFiltered, currChrom, currPoints, leftBounds, rightBounds, leftFlip, rightFlip, tickSpaceL, tickSpaceR, genomeRow = makeFilteredData(fusion_name, currChromPoints, reads_file)
	if len(readsFiltered) > 0:
		if '|' in readsFiltered.at[0, 'readID']:
			readsFiltered['isoSupport'] = readsFiltered['readID'].str.split('|', expand=True)[0].astype(int)
			readsFiltered = readsFiltered.sort_values(by='isoSupport')
		source = ColumnDataSource(readsFiltered[['readID','fusionID','geneName','chrom', 'tStart', 'tEnd']])
		readsFiltered = readsFiltered.reset_index(drop=True)
		tools=[BoxZoomTool(dimensions='width'), PanTool(), ResetTool()]
		pl = figure(name='pl', title=currChromPoints[0][0] + " " + str(currChrom[0]), y_range=list(OrderedDict.fromkeys(readsFiltered['readID'])),
					plot_width=(len(readsFiltered.loc[0, 'readID'])*6)+300, plot_height=len(list(set(readsFiltered['readID'])))*15+150,
					x_range=(leftBounds[0], leftBounds[1]), tools=tools, toolbar_location=None, min_border_bottom=100)
		pr = figure(name='pr', title=currChromPoints[1][0] + ' ' + str(currChrom[1]), y_range=pl.y_range, plot_width=300,
					plot_height=len(list(set(readsFiltered['readID'])))*15+150,
					x_range=(rightBounds[0], rightBounds[1]), tools=tools, toolbar_location='right', min_border_bottom=100)
		pl.segment(x0='tStart', y0='readID', x1='tEnd',
				   y1='readID', color="green", line_width=5, source=source)
		pr.segment(x0='tStart', y0='readID', x1='tEnd',
				   y1='readID', color="blue", line_width=5, source=source)

		breakLineL = Span(location=currPoints[0], dimension='height', line_color='red', line_width=3)
		breakLineR = Span(location=currPoints[1], dimension='height', line_color='red', line_width=3)
		pl.renderers.extend([breakLineL])
		pr.renderers.extend([breakLineR])

		pl.below[0].formatter.use_scientific = False
		#pl.add_layout(LinearAxis(ticker=SingleIntervalTicker(interval=tickSpaceL), formatter=BasicTickFormatter(use_scientific=False)), 'above')
		pl.xaxis.ticker = SingleIntervalTicker(interval=tickSpaceL)
		pl.xaxis.major_label_orientation = pi/4
		pl.ygrid.grid_line_color = None
		pl.xgrid.ticker =  SingleIntervalTicker(interval=tickSpaceL)
		pl.y_range.range_padding = 0

		pr.below[0].formatter.use_scientific = False
		#pr.add_layout(LinearAxis(ticker=SingleIntervalTicker(interval=tickSpaceR), formatter=BasicTickFormatter(use_scientific=False)), 'above')
		pr.xaxis.ticker = SingleIntervalTicker(interval=tickSpaceR)
		pr.xaxis.major_label_orientation = pi/4
		pr.ygrid.grid_line_color = None
		pr.xgrid.ticker = SingleIntervalTicker(interval=tickSpaceR)
		pr.y_range.range_padding = 0
		pr.yaxis.visible = False
		return pl, pr
	else:
		return figure(), figure()

hiddenSource = ColumnDataSource(reads)
hiddenShortSource = ColumnDataSource()
pl, pr = makeFullPlot("PIWIL4--GUCY1A2",[fusions.loc[currFusion, "5' breakpoint"].split("-")[1:],
										 fusions.loc[currFusion, "3' breakpoint"].split("-")[1:]], pd.DataFrame.from_dict(hiddenSource.data))
columns = [TableColumn(field = "#name", title = "Name", formatter=StringFormatter(), width=390),
		   TableColumn(field = "mapping score(1 is good)", title = "Map Score", formatter=NumberFormatter(format='0.000')),
		   TableColumn(field = "spanning reads", title = "Spanning Reads", formatter=NumberFormatter(format='0')),
		   TableColumn(field='seq agreement near breakpoint (1 is good)', title='Seq. Agreement', formatter=NumberFormatter(format='0.000')),
		   TableColumn(field='fraction repetitive (0 is good)', title='Frac. Repetitive', formatter=NumberFormatter(format='0.000'))]
tableSource = ColumnDataSource(fusions)
data_table = DataTable(source = tableSource, columns = columns, width = 420, height_policy ='auto', editable = False, reorderable = False, name='tableSource')

def py_callback(attr, old, new):
	a = 0
	for i in tableSource.selected.indices:
		a = i
	currLayout = curdoc().get_model_by_name('mainLayout').children
	currLayout.remove(curdoc().get_model_by_name('pl'))
	currLayout.remove(curdoc().get_model_by_name('pr'))
	#print('curr Breakpoints', tableSource.data["5' breakpoint"][a])
	currChromPoints = [['-'.join(tableSource.data["5' breakpoint"][a].strip("5'-").split("-")[:-2])] + tableSource.data["5' breakpoint"][a].split("-")[-2:],
					   ['-'.join(tableSource.data["3' breakpoint"][a].strip("3'-").split("-")[:-2])] + tableSource.data["3' breakpoint"][a].split("-")[-2:]]
	fusionName = tableSource.data["#name"][a]
	print(currChromPoints, fusionName)
	if "confirmed reads" in tableSource.data:
		readIDs, leftSeq, rightSeq = getReadOrder(pd.DataFrame.from_dict(hiddenShortSource.data), fusionName, currChromPoints)
		#print(readIDs)
		pl = view_alignment(readIDs, leftSeq, currChromPoints[0], 'left', plot_width=int(len(readIDs[0])*5.58) + 300)
		pr = view_alignment(readIDs, rightSeq, currChromPoints[1], 'right', plot_width=300)
	else:
		pl, pr = makeFullPlot(fusionName, currChromPoints, pd.DataFrame.from_dict(hiddenSource.data))
	currLayout.append(pl)
	currLayout.append(pr)

def upload_fit_data(attr, old, new):
	#currSubLayout = curdoc().get_model_by_name('subColumn').children
	decoded = b64decode(new)
	f = io.BytesIO(decoded)
	#firstLine = f.readline().strip().decode('utf-8')
	#f.seek(0)
	new_df = pd.read_csv(f, sep='\t', index_col=False)
	tableSource.data.update(ColumnDataSource(new_df).data)

def upload_reads_data(attr, old, new):
	#print("reads data upload succeeded")
	currLayout = curdoc().get_model_by_name('mainLayout').children
	decoded = b64decode(new)
	f = io.BytesIO(decoded)
	currChromPoints = [['-'.join(tableSource.data["5' breakpoint"][0].strip("5'-").split("-")[:-2])] + tableSource.data["5' breakpoint"][0].split("-")[-2:],
					   ['-'.join(tableSource.data["3' breakpoint"][0].strip("3'-").split("-")[:-2])] + tableSource.data["3' breakpoint"][0].split("-")[-2:]]
	print(currChromPoints)
	fusionName = tableSource.data["#name"][0]
	firstLine = f.readline().strip().decode("utf-8")
	f.seek(0)
	if (firstLine[:3] == 'chr'):
		new_df = pd.read_csv(f, sep='\t', index_col=False, names=['chrom','chromStart','chromEnd','name','score','strand',
															  'thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts'])
		pl, pr = makeFullPlot(fusionName, currChromPoints, new_df)
		hiddenSource.data.update(ColumnDataSource(new_df).data)
		if "confirmed reads" in tableSource.data:
			tableSource.remove("confirmed reads")
	else:
		readsList = []
		for line in f:
			line = line.decode('utf-8').strip().split('\t')
			readsList.append(line[0].split('->') + [line[3], line[11]])
		myDF = pd.DataFrame(readsList, columns=['fusion', 'gene', 'id', 'seq'])
		readIDs, leftSeq, rightSeq = getReadOrder(myDF, fusionName, currChromPoints)
		pl = view_alignment(readIDs, leftSeq, currChromPoints[0], 'left', plot_width=int(len(readIDs[0])*5.58) + 300)
		pr = view_alignment(readIDs, rightSeq, currChromPoints[1], 'right', plot_width=300)
		hiddenShortSource.data.update(ColumnDataSource(myDF).data)
		#print(hiddenSource.data['left'])
	currLayout.remove(curdoc().get_model_by_name('pl'))
	currLayout.remove(curdoc().get_model_by_name('pr'))
	currLayout.append(pl)
	currLayout.append(pr)
	#VERY IMPORTANT
	#hiddenSource.data.update(ColumnDataSource(new_df).data)

pageTitle = Div(text='Long read gene fusion visualization', style={'font-size': '24px', 'text-decoration':'underline'})
fileTitle = Div(text='Fusion file (.tsv)')
file2Title = Div(text='Reads file (.bed)')
file_input = FileInput(accept=".tsv", name='fusion')#title='Fusion file (.tsv)')
file_input.on_change('value', upload_fit_data)
file_input2 = FileInput(accept=".bed", name='reads')#title='Reads file (.bed)')
file_input2.on_change('value', upload_reads_data)


tableSource.selected.on_change('indices', py_callback)
#callback = CustomJS(args = dict(source = tableSource, filteredSource = tableSource), code = source_code)
#tableSource.selected.js_on_change('indices', callback)
#data_table.on_change('indices', data_callback)
subColumn = column([pageTitle, fileTitle, file_input, file2Title, file_input2, data_table], name='subColumn')
mainLayout = row([subColumn, pl, pr], name='mainLayout')
curdoc().add_root(mainLayout)

