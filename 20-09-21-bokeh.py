from bokeh.models import ColumnDataSource, TableColumn, StringFormatter, NumberFormatter, DataTable, PanTool,ResetTool, \
	BoxZoomTool, Span, SingleIntervalTicker, Div, TickFormatter, BasicTickFormatter, LinearAxis, WheelZoomTool, Range1d, FactorRange
from bokeh.models.glyphs import Text, Rect
from bokeh.models.widgets import FileInput
from bokeh.themes import Theme
from bokeh.colors import RGB
from bokeh.layouts import column, row, gridplot
from bokeh.plotting import curdoc, figure
from bokeh.events import Tap
from bokeh.io import show
from bokeh.util.compiler import TypeScript
from collections import OrderedDict
from math import pi
import pandas as pd
import requests
from pybase64 import b64decode
import io, subprocess
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
#Sequence aligner code from https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner
fusions = pd.read_csv('02-09-2020DRR-flair.alignedFusions.tsv', sep='\t', index_col=False)
reads = pd.read_csv('02-09-2020DRR-flair.alignedReads.bed', sep='\t', index_col=False,
					names=['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb',
						   'blockCount','blockSizes','blockStarts'])
reads[['fusionID','readID','geneName', 'seq']] = pd.DataFrame([x.split('-.-')[:4] for x in reads['name']], index= reads.index)#temp[temp.columns[0:4]]
reads['geneName'] = reads['geneName'].str.split('/', expand = True)[0]
currFusion = 0

def view_fasta(myDF, readName, fusionName, plot_width=750, fontsize='9pt'):
	myDF = myDF.set_index('id')
	inFusionReads = myDF.loc[(myDF['read']=='REF'+readName) & (myDF['fusion']==fusionName), :].copy()
	outFusionReads = myDF.loc[(myDF['read']=='REF'+readName) & (myDF['fusion']!=fusionName) & (myDF['fusion']!='reference'), :].copy()
	ref = myDF.loc[(myDF['read']=='REF'+readName) & (myDF['fusion']=='reference'), :].copy()
	print(len(ref))
	ids = list(outFusionReads.index) + list(inFusionReads.index) + ['RAW READ SEQUENCE']
	seqs = list(outFusionReads['seq']) + list(inFusionReads['seq']) + list(ref['seq'])
	t = [list(i) for i in seqs]
	text = [item for sublist in t for item in sublist]
	clrs = {'T':RGB(101,194,101), 'A':RGB(255, 84,84), 'G':RGB(255, 195,84), 'C':RGB(87,87, 255), '-':'white', 'N':'white'}
	colors = [clrs[i] for i in text]
	x = [list(np.arange(0, len(seqs[0]))) for i in range(len(seqs))]
	gx = [item for sublist in x for item in sublist]
	gx1 = [i + 0.5 for i in gx]
	gx2 = [i + 1 for i in gx]
	y = [[i]*len(seqs[0]) for i in ids]
	gy = [item for sublist in y for item in sublist]
	source = ColumnDataSource(dict(x=gx1, recty=gy, rectx1=gx, rectx2=gx2, text=text, colors=colors))# y=gy, recty=gy+0.5,
	plot_height = len(seqs)*20+100
	#entire sequence view (no text, with zoom)
	x_range = Range1d(0, len(seqs[0]), bounds='auto')
	p = figure(name='fullReadView', title=readName, plot_width= plot_width, plot_height=100,
			   x_range=x_range, y_range=list(OrderedDict.fromkeys(ids)), tools="xwheel_zoom, xpan",
			   min_border=0, toolbar_location='right')
	p.segment(x0='rectx1', y0='recty', x1='rectx2',
			   y1='recty', color="colors", line_width=3, source=source, alpha=0.5)
	p.yaxis.visible = False
	p.grid.visible = False
	p1 = figure(name='50bpView', plot_width=plot_width, plot_height=plot_height,
				x_range=(0, 50), y_range=list(OrderedDict.fromkeys(ids)), tools="xpan,reset",
				min_border_bottom=50, min_border_top=50, toolbar_location='right')#, lod_factor=1)
	p1.segment(x0='rectx1', y0='recty', x1='rectx2',
			   y1='recty', color="colors", line_width=14, source=source, alpha=0.5)
	glyph = Text(x="x", y="recty", text="text", text_align='center',text_color="black",
				 text_font="monospace",text_font_size=fontsize, text_baseline='middle')
	p1.add_glyph(source, glyph)
	p1.grid.visible = False
	p1.add_layout(LinearAxis(formatter=BasicTickFormatter(use_scientific=False), major_label_orientation = pi/4), 'above')
	p1.xaxis.major_label_text_font_style = "bold"
	p1.yaxis.minor_tick_line_width = 0
	p1.yaxis.major_tick_line_width = 0
	p1.below[0].formatter.use_scientific = False
	p1.xaxis.major_label_orientation = pi/4
	p.on_event(Tap, full_view_click_callback)
	print('done making plot')
	#p1.y_range = FactorRange(factors=ids)
	#p1.name='pl'
	return column([p, p1], name='pl')

def getReadOrder(myDF, fusionName, currChromPoints):
	seqLen = len(myDF.at[0, 'seq'])
	myDF = myDF.set_index('id')
	leftRefSeq = myDF.loc[(myDF.index =='reference') & (myDF['fusion']==fusionName) & (myDF['gene'] == currChromPoints[0][0]), :].copy()
	rightRefSeq = myDF.loc[(myDF.index =='reference') & (myDF['fusion']==fusionName) & (myDF['gene'] == currChromPoints[1][0]), :].copy()
	if len(leftRefSeq) < 1: leftRefSeq['ref'] = '-' * seqLen
	else: myDF.drop(myDF.loc[(myDF.index =='reference') & (myDF['fusion']==fusionName) & (myDF['gene'] == currChromPoints[0][0])].index, inplace=True)
	if len(rightRefSeq) < 1: rightRefSeq['ref'] = '-' * seqLen
	else: myDF.drop(myDF.loc[(myDF.index =='reference') & (myDF['fusion']==fusionName) & (myDF['gene'] == currChromPoints[1][0])].index, inplace=True)
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
	print("shared reads ", len(readsBoth))
	readsBoth = readsBoth.head(50)
	readIDs = list(readsRightOnly.index) + list(readsLeftOnly.index) + list(readsBoth.index) + ['REFERENCE SEQUENCE']
	leftSeq = [["-"]*seqLen for i in range(len(readsRightOnly))] + list(readsLeftOnly['seq']) + list(readsBoth['lseq']) + list(leftRefSeq['seq'])
	rightSeq = list(readsRightOnly['seq']) + [["-"]*seqLen for i in range(len(readsLeftOnly))] + list(readsBoth['rseq']) + list(rightRefSeq['seq'])
	return readIDs, leftSeq, rightSeq

def view_alignment(ids, seqs, chromPoints, side, fontsize="9pt", plot_width=800):
	"""Bokeh sequence alignment view"""
	text = [i for s in list(seqs) for i in s]
	clrs = {'T':RGB(153, 204, 153), 'A':RGB(255, 153, 153), 'G':RGB(255, 219, 153), 'C':RGB(153, 153, 255), '-':'white', 'N':'white'}
	colors = [clrs[i] for i in text]
	N = len(seqs[0])/2
	center = int(chromPoints[2])
	x = np.arange(center-N,center+N)
	y = np.arange(0,len(seqs),1)
	xx, yy = np.meshgrid(x, y)
	gx = xx.ravel()
	gy = yy.flatten()
	source = ColumnDataSource(dict(x=gx+0.5, y=gy, recty=gy+0.5, rectx1=gx, rectx2=gx+1, text=text, colors=colors)) #id=ids - too short, doesn't match
	plot_height = len(seqs)*20+200
	view_range = (center-20, center+20)
	p1 = figure(title=' '.join(chromPoints), plot_width=plot_width, plot_height=plot_height,
				x_range=view_range, y_range=ids, tools="xpan,reset",
				min_border_bottom=100, min_border_top=100, toolbar_location='below')#, lod_factor=1)
	glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
				 text_font="monospace",text_font_size=fontsize)
	p1.segment(x0='rectx1', y0='recty', x1='rectx2',
			   y1='recty', color="colors", line_width=14, source=source)
	p1.add_glyph(source, glyph)
	breakLine = Span(location=int(chromPoints[2]), dimension='height', line_color='red', line_width=2)
	p1.renderers.extend([breakLine])
	p1.grid.visible = False
	p1.add_layout(LinearAxis(formatter=BasicTickFormatter(use_scientific=False), major_label_orientation = pi/4), 'above')
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
	return p1, source

def makeFilteredData(fusion_name, currChromPoints, reads_file):
	myReads = reads_file
	# myReads[['fusionID','readID','geneName', 'seq']] = pd.DataFrame([x.split('-.-')[:4] for x in myReads['name']], index= myReads.index)#temp[temp.columns[0:4]]
	# myReads['geneName'] = myReads['geneName'].str.split('/', expand = True)[0]
	readsFiltered = myReads.loc[myReads['fusionID']==fusion_name, :].copy()
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
		leftFlip, rightFlip = False, False
		currChrom = [currChromPoints[0][1], currChromPoints[1][1]]
		currPoints = [int(currChromPoints[0][2]), int(currChromPoints[1][2])]
		a = [abs(i-int(currChromPoints[0][2])) for i in leftBounds]
		if a[0] < a[1]:
			leftBounds = leftBounds[::-1]
			leftFlip = True
		a = [abs(i-int(currChromPoints[1][2])) for i in rightBounds]
		if a[1] < a[0]:
			rightBounds = rightBounds[::-1]
			rightFlip = True
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
		return readsExpanded, currChrom, currPoints, leftBounds, rightBounds, leftFlip, rightFlip
	else:
		print(len(readsFiltered), fusion_name)#, set(myReads['fusionID']))
		return readsFiltered, [], [], [], [], False, False, 0, 0, 0

def makeFullPlot(fusion_name, currChromPoints, reads_file):
	readsFiltered, currChrom, currPoints, leftBounds, rightBounds, leftFlip, rightFlip = makeFilteredData(fusion_name, currChromPoints, reads_file)
	if len(readsFiltered) > 0:
		if '|' in readsFiltered.at[0, 'readID']:
			readsFiltered['isoSupport'] = readsFiltered['readID'].str.split('|', expand=True)[0].astype(int)
			readsFiltered = readsFiltered.sort_values(by='isoSupport')
		source = ColumnDataSource(readsFiltered[['readID','fusionID','geneName','chrom', 'tStart', 'tEnd']])
		readsFiltered = readsFiltered.reset_index(drop=True)
		tools=[BoxZoomTool(dimensions='width'), WheelZoomTool(dimensions='width'), PanTool(), ResetTool()]
		pl = figure(name='pl', title=currChromPoints[0][0] + " " + str(currChrom[0]), y_range=list(OrderedDict.fromkeys(readsFiltered['readID'])),
					plot_width=(len(readsFiltered.loc[0, 'readID'])*6)+300, plot_height=len(list(set(readsFiltered['readID'])))*15+220,
					x_range=(leftBounds[0], leftBounds[1]), tools=tools, toolbar_location=None, min_border_bottom=100, min_border_top=100)
		pr = figure(name='pr', title=currChromPoints[1][0] + ' ' + str(currChrom[1]), y_range=pl.y_range, plot_width=300,
					plot_height=len(list(set(readsFiltered['readID'])))*15+220,
					x_range=(rightBounds[0], rightBounds[1]), tools=tools, toolbar_location='right', min_border_bottom=100, min_border_top=100)
		pl.segment(x0='tStart', y0='readID', x1='tEnd', y1='readID', color="green", line_width=5, source=source)
		pr.segment(x0='tStart', y0='readID', x1='tEnd', y1='readID', color="blue", line_width=5, source=source)
		pl.renderers.extend([Span(location=currPoints[0], dimension='height', line_color='red', line_width=3)])
		pr.renderers.extend([Span(location=currPoints[1], dimension='height', line_color='red', line_width=3)])

		pl.below[0].formatter.use_scientific = False
		pl.add_layout(LinearAxis(formatter=BasicTickFormatter(use_scientific=False), major_label_orientation = pi/4), 'above')
		pl.xaxis.major_label_orientation = pi/4
		pl.ygrid.grid_line_color = None
		pl.y_range.range_padding = 0
		pr.below[0].formatter.use_scientific = False
		pr.add_layout(LinearAxis(formatter=BasicTickFormatter(use_scientific=False), major_label_orientation = pi/4), 'above')
		pr.xaxis.major_label_orientation = pi/4
		pr.ygrid.grid_line_color = None
		pr.y_range.range_padding = 0
		pr.yaxis.visible = False
		return pl, pr
	else:
		return figure(), figure()

hiddenSource = ColumnDataSource(reads)
hiddenShortSource = ColumnDataSource()
hiddenFilteredSource = ColumnDataSource()
hiddenFastaSource = ColumnDataSource()
pl, pr = makeFullPlot("PIWIL4--GUCY1A2",[fusions.loc[currFusion, "5' breakpoint"].split("-")[1:],
										 fusions.loc[currFusion, "3' breakpoint"].split("-")[1:]], pd.DataFrame.from_dict(hiddenSource.data))
columns = [TableColumn(field = "#name", title = "Name", formatter=StringFormatter(), width=390),
		   TableColumn(field = "mapping score(1 is good)", title = "Map Score", formatter=NumberFormatter(format='0.000')),
		   TableColumn(field = "spanning reads", title = "Spanning Reads"),#, formatter=NumberFormatter(format='0')),
		   TableColumn(field='seq agreement near breakpoint (1 is good)', title='Seq. Agreement', formatter=NumberFormatter(format='0.000')),
		   TableColumn(field='fraction repetitive (0 is good)', title='Frac. Repetitive', formatter=NumberFormatter(format='0.000'))]
tableSource = ColumnDataSource(fusions)
data_table = DataTable(source = tableSource, columns = columns, width = 420, height_policy ='auto', editable = False,
					   reorderable = False, name='tableSource', sortable=True)

def table_click_callback(attr, old, new):
	if len(tableSource.selected.indices) > 0:
		curdoc().get_model_by_name('informUser').text = 'fusion selection registered'
		a = tableSource.selected.indices[0]
		currLayout = curdoc().get_model_by_name('mainLayout').children
		currChromPoints = [['-'.join(tableSource.data["5' breakpoint"][a].strip("5'-").split("-")[:-2])] + tableSource.data["5' breakpoint"][a].split("-")[-2:],
						   ['-'.join(tableSource.data["3' breakpoint"][a].strip("3'-").split("-")[:-2])] + tableSource.data["3' breakpoint"][a].split("-")[-2:]]
		fusionName = tableSource.data["#name"][a]
		print(currChromPoints, fusionName)
		if "confirmed reads" in tableSource.data:
			readIDs, leftSeq, rightSeq = getReadOrder(pd.DataFrame.from_dict(hiddenShortSource.data), fusionName, currChromPoints)
			pl, myData = view_alignment(readIDs, leftSeq, currChromPoints[0], 'left', plot_width=int(len(readIDs[0])*5.58) + 300)
			pr, temp = view_alignment(readIDs, rightSeq, currChromPoints[1], 'right', plot_width=300)
			hiddenFilteredSource.data.update(myData.data)
			pl.on_event(Tap, chart_click_callback)
		else:
			pl, pr = makeFullPlot(fusionName, currChromPoints, pd.DataFrame.from_dict(hiddenSource.data))
		currLayout.remove(curdoc().get_model_by_name('pl'))
		currLayout.remove(curdoc().get_model_by_name('pr'))
		currLayout.append(pl)
		currLayout.append(pr)
		curdoc().get_model_by_name('informUser').text = 'fusion reads chart updated'

def chart_click_callback(event):
	curdoc().get_model_by_name('informUser').text = 'click registered'
	plTemp = curdoc().get_model_by_name('pl')
	print(plTemp.y_range.factors[round(event.y-0.5)])#(hiddenFilteredSource.data['id'][round(event.y-0.5)])
	if bool(hiddenFastaSource.data):
		currLayout = curdoc().get_model_by_name('mainLayout').children
		fusionName = tableSource.data["#name"][tableSource.selected.indices[0]]
		tableSource.selected.indices = []
		readName = plTemp.y_range.factors[round(event.y-0.5)]#hiddenFilteredSource.data['id'][round(event.y-0.5)]
		pl = view_fasta(pd.DataFrame.from_dict(hiddenFastaSource.data), readName, fusionName)
		pr = figure(name='pr')
		currLayout.remove(curdoc().get_model_by_name('pl'))
		currLayout.remove(curdoc().get_model_by_name('pr'))
		currLayout.append(pl)
		currLayout.append(pr)
		print(currLayout)
	curdoc().get_model_by_name('informUser').text = 'single read chart loaded'

def full_view_click_callback(event):
	curdoc().get_model_by_name('informUser').text = 'click registered'
	#fullView = curdoc().get_model_by_name('fullReadView')#50bpView
	center = int(event.x)#(hiddenFilteredSource.data['id'][round(event.y-0.5)])
	curdoc().get_model_by_name('informUser').text = 'click registered - ' + str(int(event.x))
	bpView = curdoc().get_model_by_name('50bpView')#'fullReadView')#50bpView
	bpView.x_range.start = center-25
	bpView.x_range.end = center+25

def upload_fasta_data(attr, old, new):
	curdoc().get_model_by_name('informUser').text = 'uploading fasta'
	decoded = b64decode(new)
	f = io.BytesIO(decoded)
	readsList = []
	for line in f:
		line = line.decode('utf-8').strip().split('\t')
		readsList.append([line[0], line[3].split('->')[0], ('->').join(line[3].split('->')[:-1]), line[11]])
	myDF = pd.DataFrame(readsList, columns=['read', 'fusion', 'id', 'seq'])
	hiddenFastaSource.data.update(ColumnDataSource(myDF).data)
	curdoc().get_model_by_name('informUser').text = 'fasta uploaded'

def upload_reads_data(attr, old, new):
	curdoc().get_model_by_name('informUser').text = 'uploading reads'
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
		new_df[['fusionID','readID','geneName', 'seq']] = pd.DataFrame([x.split('-.-')[:4] for x in new_df['name']], index= new_df.index)#temp[temp.columns[0:4]]
		new_df['geneName'] = new_df['geneName'].str.split('/', expand = True)[0]
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
		pl, myData = view_alignment(readIDs, leftSeq, currChromPoints[0], 'left', plot_width=int(len(readIDs[0])*5.58) + 300)
		pr, temp = view_alignment(readIDs, rightSeq, currChromPoints[1], 'right', plot_width=300)
		pl.on_event(Tap, chart_click_callback)
		hiddenFilteredSource.data.update(myData.data)
		hiddenShortSource.data.update(ColumnDataSource(myDF).data)
	currLayout.remove(curdoc().get_model_by_name('pl'))
	currLayout.remove(curdoc().get_model_by_name('pr'))
	currLayout.append(pl)
	currLayout.append(pr)
	curdoc().get_model_by_name('informUser').text = 'reads uploaded'

def upload_table_data(attr, old, new):
	curdoc().get_model_by_name('informUser').text = 'uploading fusions'
	decoded = b64decode(new)
	f = io.BytesIO(decoded)
	new_df = pd.read_csv(f, sep='\t', index_col=False, error_bad_lines=False) #DROPS 3-gene fusions completely
	if "confirmed reads" in new_df.columns:
		new_df['spanning reads'] = new_df['confirmed reads']
	elif "confirmed reads" in tableSource.data:
		tableSource.remove("confirmed reads")
	tableSource.data.update(ColumnDataSource(new_df).data)
	curdoc().get_model_by_name('informUser').text = 'fusions uploaded'

pageTitle = Div(text='Long read gene fusion visualization', style={'font-size': '24px', 'text-decoration':'underline'})
fusionUploadTitle = Div(text='Fusion file (.tsv)', name='fusionUploadText')
readsUploadTitle = Div(text='Reads file (.bed)')
fusionUpload = FileInput(accept=".tsv", name='fusion')
fusionUpload.on_change('value', upload_table_data)
readsUpload = FileInput(accept=".bed", name='reads')
readsUpload.on_change('value', upload_reads_data)
fastaUploadTitle = Div(text='Fusion fragments aligned to original reads file (.bed, optional)')
fastaUpload = FileInput(accept=".bed", name='fareads')
fastaUpload.on_change('value', upload_fasta_data)
informUser = Div(text='updates will appear here', name='informUser')
tableSource.selected.on_change('indices', table_click_callback)

subColumn = column([pageTitle, fusionUploadTitle, fusionUpload, readsUploadTitle, readsUpload, fastaUploadTitle, fastaUpload, informUser, data_table], name='subColumn')
mainLayout = row([subColumn, pl, pr], name='mainLayout')
curdoc().add_root(mainLayout)
