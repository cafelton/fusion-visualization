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
from bokeh.io import show
from bokeh.models import Button, CustomJS
from random import randint
from bokeh.io import output_file, show
from bokeh.layouts import widgetbox
from bokeh.models.widgets import RadioButtonGroup
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, HTMLTemplateFormatter
from bokeh.io import output_notebook, show
from bokeh.plotting import curdoc, figure
from random import randint
from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource
from bokeh.io import export_png
from PIL import Image
import base64

#comment just to edit and deploy lol
#next deploy comment lol
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
annotations = pd.read_csv('gencode.v37.annotation-edited.tsv', sep='\t', index_col=False, header=0) #names=['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb',
						   #'blockCount','blockSizes','blockStarts', 'starts', 'sizes', 'tStart', 'tEnd'])

# fusions = pd.read_csv('02-09-2020DRR-flair.alignedFusion-final.tsv', sep='\t', index_col=False)
# reads = pd.read_csv('02-09-2020DRR-flair.alignedReads-final.bed', sep='\t', index_col=False,
# 					names=['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb',
# 						   'blockCount','blockSizes','blockStarts', 'cigar', 'sequence'])
# reads[['fusionID','readID','geneName', 'seq']] = pd.DataFrame([x.split('-.-')[:4] for x in reads['name']], index= reads.index)#temp[temp.columns[0:4]]
# reads['geneName'] = reads['geneName'].str.split('/', expand = True)[0]
currFusion = 0

def getReadOrder(myDF, fusionName, currChromPoints):
	seqLen = len(myDF.at[0, 'seq'])
	myDF = myDF.set_index('id')
	leftRefSeq = myDF.loc[(myDF.index =='reference') & (myDF['fusion'].str.contains(fusionName)) & (myDF['gene'] == currChromPoints[0][0]), :].copy()
	rightRefSeq = myDF.loc[(myDF.index =='reference') & (myDF['fusion'].str.contains(fusionName)) & (myDF['gene'] == currChromPoints[1][0]), :].copy()
	if len(leftRefSeq) < 1: leftRefSeq['ref'] = '-' * seqLen
	else: myDF.drop(myDF.loc[(myDF.index =='reference') & (myDF['fusion'].str.contains(fusionName)) & (myDF['gene'] == currChromPoints[0][0])].index, inplace=True)
	if len(rightRefSeq) < 1: rightRefSeq['ref'] = '-' * seqLen
	else: myDF.drop(myDF.loc[(myDF.index =='reference') & (myDF['fusion'].str.contains(fusionName)) & (myDF['gene'] == currChromPoints[1][0])].index, inplace=True)
	readsLeft = myDF.loc[(myDF['fusion'].str.contains(fusionName)) & (myDF['gene'] == currChromPoints[0][0]), :].copy()
	readsRight = myDF.loc[(myDF['fusion'].str.contains(fusionName)) & (myDF['gene'] == currChromPoints[1][0]), :].copy()
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
	readIDs = list(readsRightOnly.index) + list(readsLeftOnly.index) + list(readsBoth.index) + ['REFERENCE SEQUENCE']
	leftSeq = [["-"]*seqLen for i in range(len(readsRightOnly))] + list(readsLeftOnly['seq']) + list(readsBoth['lseq']) + list(leftRefSeq['seq'])
	rightSeq = list(readsRightOnly['seq']) + [["-"]*seqLen for i in range(len(readsLeftOnly))] + list(readsBoth['rseq']) + list(rightRefSeq['seq'])
	return readIDs, leftSeq, rightSeq

def view_alignment(ids, seqs, chromPoints, side, fontsize="9pt", plot_width=800):
	"""Bokeh sequence alignment view"""
	text = [i for s in list(seqs) for i in s]
	clrs = {'T':RGB(153, 204, 153), 'A':RGB(255, 153, 153), 'G':RGB(255, 219, 153), 'C':RGB(153, 153, 255), '-':'white', 'N':'white', '*':'white'}
	colors = [clrs[i] for i in text]
	N = len(seqs[0])/2
	center = int(chromPoints[2])
	x = np.arange(center-N,center+N)
	y = np.arange(0,len(seqs),1)
	xx, yy = np.meshgrid(x, y)
	gx = xx.ravel()
	gy = yy.flatten()
	source = ColumnDataSource(dict(x=gx+0.5, y=gy, recty=gy+0.5, rectx1=gx, rectx2=gx+1, text=text, colors=colors))
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

def data_table_formatter(selected_opt, table_source):
	# https://stackoverflow.com/questions/50996875/how-to-color-rows-and-or-cells-in-a-bokeh-datatable
	if selected_opt == 0:
		template = """
					<p style="color:<%=
						(function colorfromint(){
							if( subgroup == " ")
								{return('green')}
							}()) %>;">
						<%= value %>
					</p>
					"""
	elif selected_opt == 1:
		template = """
					<p style="color:<%=
						(function colorfromint(){
							count = (name.split("--").length - 1)
							if ((remap_spanning_reads > 0 || subgroup.includes('--')) && count < 2)
								{return('red')}
							}()) %>;">
						<%= value %>
					</p>
					"""
	elif selected_opt == 2:
		template = """
					<p style="color:<%=
						(function colorfromint(){
							if ( num_of_isoform > 0)
								{return('blue')}
							}()) %>;">
						<%= value %>
					</p>
					"""
	formatter = HTMLTemplateFormatter(template=template)
	columns = [TableColumn(field="name", title="Name", formatter=formatter, width=390),
			   TableColumn(field="subgroup", title="Subgroup Name", formatter=formatter, width=450),  ##ADDED
			   TableColumn(field="mapping_score", title="Map Score", formatter=formatter),
			   TableColumn(field="spanning_reads", title="Spanning Reads", formatter=formatter),
			   TableColumn(field='seq_agreement_near_breakpoint', title='Seq. Agreement', formatter=formatter),
			   TableColumn(field='num_of_isoform', title='Isoforms', formatter=formatter),
			   TableColumn(field='fraction_repetitive', title='Frac. Repetitive',
						   formatter=formatter)]
	data_table = DataTable(source=table_source, columns=columns, fit_columns=True, width=420, height_policy='auto',
						   editable=False,
						   reorderable=False, name='tableSource', sortable=True, row_height=30)
	return data_table

hiddenSource = ColumnDataSource()#reads)
hiddenShortSource = ColumnDataSource()
hiddenFilteredSource = ColumnDataSource()
hiddenFastaSource = ColumnDataSource()
hiddenButtonSource = ColumnDataSource(data={"selected":[0]})
pl, pr = figure(name = "pl"), figure(name = "pr")
tableSource = ColumnDataSource()#fusions)
anoSource = ColumnDataSource(annotations)
data_table = data_table_formatter(0, tableSource)

def table_click_callback(attr, old, new):
	if len(tableSource.selected.indices) > 0:
		curdoc().get_model_by_name('informUser').text = 'fusion selection registered'
		a = tableSource.selected.indices[0] #<-- a = row that is selected
		currLayout = curdoc().get_model_by_name('mainLayout').children
		subLayout = curdoc().get_model_by_name('subColumn').children
		currChromPoints = [['-'.join(tableSource.data["5' breakpoint"][a].strip("5'-").split("-")[:-2])] + tableSource.data["5' breakpoint"][a].split("-")[-2:],
						   ['-'.join(tableSource.data["3' breakpoint"][a].strip("3'-").split("-")[:-2])] + tableSource.data["3' breakpoint"][a].split("-")[-2:]]
		if tableSource.data["name"][a] == " ":
			fusionName = tableSource.data["subgroup"][a]
		else:
			fusionName = tableSource.data["name"][a]
		if len(fusionName.split("--")) > 2:
			for i in range(len(fusionName.split("--")) - 2):
				currChromPoints += [['-'.join(tableSource.data[str(i + 2) + "-5'"][a].strip("5'-").split("-")[:-2])] +
									tableSource.data[str(i + 2) + "-5'"][a].split("-")[-2:],
									['-'.join(tableSource.data[str(i + 2) + "-3'"][a].strip("3'-").split("-")[:-2])] +
									tableSource.data[str(i + 2) + "-3'"][a].split("-")[-2:]]
		if(hiddenButtonSource.data["selected"][0] == 1):
			data_table = data_table_formatter(1, tableSource) ##
			readIDs, leftSeq, rightSeq = getReadOrder(pd.DataFrame.from_dict(hiddenShortSource.data), fusionName, currChromPoints)
			pl, myData = view_alignment(readIDs, leftSeq, currChromPoints[0], 'left', plot_width=int(len(readIDs[0])*5.58) + 400)
			pr, temp = view_alignment(readIDs, rightSeq, currChromPoints[1], 'right', plot_width=320)
			figList = [pl, pr]
			hiddenFilteredSource.data.update(myData.data)
		elif (hiddenButtonSource.data["selected"][0] == 0):
			data_table = data_table_formatter(0, tableSource) ##
			temp = pd.DataFrame.from_dict(hiddenSource.data)
			notisoformOnly = temp[temp["isoform support"].astype(int) == 0]
			if tableSource.data["name"][a] != " ":
				figList = makeFullPlot2(fusionName, currChromPoints, notisoformOnly, anoSource)
			else:
				figList = [figure(name="p0"), figure(name="p1")]
		elif (hiddenButtonSource.data["selected"][0] == 2):
			data_table = data_table_formatter(2, tableSource) ##
			isoformReads = float(tableSource.data["num_of_isoform"][a]) #"isoform spanning reads" ###changed from int to float
			if isoformReads > 0:
				temp = pd.DataFrame.from_dict(hiddenSource.data)
				isoformOnly = temp[temp["isoform support"].astype(int) > 0]
				if len(isoformOnly) > 0:
					figList = makeFullPlot2(fusionName, currChromPoints, isoformOnly, anoSource)
					# figList = makeFullPlot2(fusionName, currChromPoints, isoformOnly)
				else:
			  		figList = [figure(name="p0"), figure(name="p1")]
			else:
			  	figList = [figure(name="p0"), figure(name="p1")]
		else:
			figList = [figure(name = "p0"), figure(name = "p1")]
		if (buttonSource.data["press"][-1] == 1): ###
			currLayout.remove(curdoc().get_model_by_name('pl'))
			currLayout.remove(curdoc().get_model_by_name('pr'))
		else:
			for i in range(len(currLayout) -1):
				currLayout.remove(curdoc().get_model_by_name('p'+str(i)))
		subLayout.remove(curdoc().get_model_by_name("tableSource")) ##
		subLayout.append(data_table)
		for i in figList:
			currLayout.append(i)
		curdoc().get_model_by_name('informUser').text = 'fusion reads chart updated'
		buttonSource.stream(dict(press=[hiddenButtonSource.data["selected"][0]], figs=[len(figList)]))

def full_view_click_callback(event):
	curdoc().get_model_by_name('informUser').text = 'click registered'
	center = int(event.x)#(hiddenFilteredSource.data['id'][round(event.y-0.5)])
	curdoc().get_model_by_name('informUser').text = 'click registered - ' + str(int(event.x))
	bpView = curdoc().get_model_by_name('50bpView')#'fullReadView')#50bpView
	bpView.x_range.start = center-25
	bpView.x_range.end = center+25

def upload_reads_data(attr, old, new):
	curdoc().get_model_by_name('informUser').text = 'uploading reads'
	currLayout = curdoc().get_model_by_name('mainLayout').children
	decoded = b64decode(new)
	f = io.BytesIO(decoded)
	currChromPoints = [['-'.join(tableSource.data["5' breakpoint"][0].strip("5'-").split("-")[:-2])] + tableSource.data["5' breakpoint"][0].split("-")[-2:],
					   ['-'.join(tableSource.data["3' breakpoint"][0].strip("3'-").split("-")[:-2])] + tableSource.data["3' breakpoint"][0].split("-")[-2:]]
	fusionName = tableSource.data["name"][0]
	if len(fusionName.split("--")) > 2:
		for i in range(len(fusionName.split("--")) - 2):
			currChromPoints += [['-'.join(tableSource.data[str(i + 2) + "-5'"][0].strip("5'-").split("-")[:-2])] +
								tableSource.data[str(i + 2) + "-5'"][0].split("-")[-2:],
								['-'.join(tableSource.data[str(i + 2) + "-3'"][0].strip("3'-").split("-")[:-2])] +
								tableSource.data[str(i + 2) + "-3'"][0].split("-")[-2:]]

	new_df = pd.read_csv(f, sep='\t', index_col=False,
			names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts', 'cigar', 'sequence'])
	new_df[['fusionID', 'isoform support', 'remapped', 'readID', 'geneName', 'seq']] = pd.DataFrame([x.split('-.-')[:6] for x in new_df['name']],
																	index=new_df.index)  # temp[temp.columns[0:4]]
	new_df['geneName'] = new_df['geneName'].str.split('/', expand=True)[0]

	remappedOnly = new_df[new_df["remapped"] == "R"] #<-PUT BACK LATER after you fix bedfilecorrecter
	remappedOnly = remappedOnly[["fusionID", "geneName", "readID", "seq"]]
	remappedOnly.columns = ['fusion', 'gene', 'id', 'seq']
	figList = makeFullPlot2(fusionName, currChromPoints, new_df, anoSource)
	# figList = makeFullPlot2(fusionName, currChromPoints, new_df)
	hiddenShortSource.data.update(ColumnDataSource(remappedOnly).data)
	hiddenSource.data.update(ColumnDataSource(new_df).data)
	if buttonSource.data["press"][-1] == 1 or len(buttonSource.data["press"]) ==1:
		currLayout.remove(curdoc().get_model_by_name('pl'))
		currLayout.remove(curdoc().get_model_by_name('pr'))
	else:
		for i in range(buttonSource.data["figs"][-1]):
			currLayout.remove(curdoc().get_model_by_name('p' + str(i)))
	for i in figList:
		currLayout.append(i)
	curdoc().get_model_by_name('informUser').text = 'reads uploaded'

def upload_table_data(attr, old, new):
	curdoc().get_model_by_name('informUser').text = 'uploading fusions'
	decoded = b64decode(new)
	f = io.BytesIO(decoded)
	totalData = []
	count = 0
	c = 0
	for line in f:
		line = line.decode('utf-8').strip().split('\t')
		if c == 0:
			line.insert(1, "subgroup")
		else:
			line.insert(1, " ")
		c += 1
		genes = line[0].split('--')
		totalData.append(line)
		if len(genes) > 2:
			count += 1
			for i in range(len(genes)-1):
				temp = [" "]*len(line)
				temp[10] = line[-(((len(genes)-1)*2)-(2*i))]
				temp[11] = line[-(((len(genes)-1)*2)-(2*i)-1)]
				temp[1] = "--".join(genes[i:i+2])
				totalData.append(temp)
	if count == 0:
		curdoc().get_model_by_name('tableSource').columns = [curdoc().get_model_by_name('tableSource').columns[0]] + curdoc().get_model_by_name('tableSource').columns[2:]
	else:
		if len(curdoc().get_model_by_name('tableSource').columns) < 7:
			curdoc().get_model_by_name('tableSource').columns.insert(1, TableColumn(field = " ", title = "Subgroup Name", formatter=StringFormatter(), width=450))
	headers = totalData.pop(0)
	new_df = pd.DataFrame(totalData, columns=headers)
	tableSource.data.update(ColumnDataSource(new_df).data)
	data_table = data_table_formatter(0, tableSource)
	currLayout = curdoc().get_model_by_name('subColumn').children
	currLayout.remove(curdoc().get_model_by_name("tableSource"))
	currLayout.append(data_table)
	curdoc().get_model_by_name('informUser').text = 'fusions uploaded'

def makeFilteredData2(fusion_name, currChromPoints, reads_file):
	bounds = []
	flip = []
	myReads = reads_file
	readsFiltered = myReads.loc[(myReads['fusionID'] == fusion_name), :].copy() #filter
	readsFiltered["geneName"] = readsFiltered["geneName"].str.split(".").str[0]
	currName = fusion_name.split("--")
	currChrom = [i[1] for i in currChromPoints]
	currPoints = [int(i[2]) for i in currChromPoints]
	plotCount = 2 if len(currChromPoints) == 2 else int(2 + (len(currChromPoints) - 2) / 2)
	for i in range(plotCount):
		tempReads = readsFiltered[readsFiltered["geneName"] == currName[i]] #filter
		if (len(tempReads) > 0):
			if (i == 0  or i == plotCount-1):
				if i == plotCount-1: i = len(currChromPoints)-1
				tempMin = tempReads["chromStart"].min()
				tempMax = tempReads["chromEnd"].max()
				tickSpace = (tempMax-tempMin)/ 10
				if abs(tempMin-currPoints[i]) > abs(tempMax-currPoints[i]):
					bounds.append([int(tempMin-tickSpace), currPoints[i]])
					flip.append(False) if i == 0 else flip.append(True)
				else:
					bounds.append([currPoints[i], int(tempMax+tickSpace)])
					flip.append(True) if i == 0 else flip.append(False)
			else:
				bounds.append([currPoints[2*i-1], currPoints[2*i]])
				flip.append(False)
		else:
			bounds.append([0, 0])
			flip.append(False)
	if len(readsFiltered) > 0:
		readsFiltered = readsFiltered[readsFiltered['readID'].isin(list(set(readsFiltered['readID']))[:200])]
		sizes = readsFiltered['blockSizes'].str.split(',', expand=True).stack().str.strip().reset_index(level=1, drop=True)
		starts = readsFiltered['blockStarts'].str.split(',', expand=True).stack().str.strip().reset_index(level=1, drop=True)
		temp = pd.concat([starts,sizes], axis=1, keys=['starts','sizes'])
		readsExpanded = readsFiltered.join(temp).reset_index(drop=True)
		readsExpanded[['starts', 'sizes']] = readsExpanded[['starts', 'sizes']].apply(pd.to_numeric)
		readsExpanded['tStart'] = readsExpanded['starts'] + readsExpanded['chromStart']
		readsExpanded['tEnd'] = readsExpanded['sizes'] + readsExpanded['tStart']
		readsFiltered['tStart'] = readsFiltered['chromStart']
		readsFiltered['tEnd'] = readsFiltered['chromEnd']
		return readsExpanded, readsFiltered, currChrom, currPoints, bounds, flip
	else:
		return readsFiltered, readsFiltered, [], [], [], []

def makeFullPlot2(fusion_name, currChromPoints, reads_file, anoSource):
	readsFiltered, readsFilteredT, currChrom, currPoints, bounds,flip = makeFilteredData2(fusion_name, currChromPoints, reads_file)
	anots = pd.DataFrame.from_dict(anoSource.data)
	anotsNew = (anots[(anots['chrom'] == currChrom[0]) & (
			((anots['tStart'] > bounds[0][0]) & (anots['tStart'] < bounds[0][1])) | (
			(anots['tEnd'] > bounds[0][0]) & (anots['tEnd'] < bounds[0][1])))])
	currChromSet = list(set(currChrom))
	for i in range(1, len(currChromSet)):
		anotsNew = anotsNew.append(anots[(anots['chrom'] == currChromSet[i]) & (
				((anots['tStart'] > bounds[i][0]) & (anots['tStart'] < bounds[i][1])) | (
				(anots['tEnd'] > bounds[i][0]) & (anots['tEnd'] < bounds[i][1])))])
	anotsNew = anotsNew.rename(columns={'name': 'readID'})
	anotsT = anotsNew.loc[(anotsNew['type'] == 't'), :].copy()
	# anotsT = anotsT.append(readsFilteredT)
	anotsNew = anotsNew[anotsNew['type'] == 'e']

	cdsList = []
	cdsTList = []
	currChromPointsSet = list(dict.fromkeys([a[0] for a in currChromPoints]))
	for i in range(len(currChromPointsSet)):
		tempDF = readsFiltered[readsFiltered["geneName"] == currChromPointsSet[i]]
		print(tempDF.head(20))
		cdsList.append(ColumnDataSource(tempDF[['readID','fusionID','geneName','chrom', 'tStart', 'tEnd']]))
		tempDF = readsFilteredT[readsFilteredT["geneName"] == currChromPointsSet[i]]
		cdsTList.append(ColumnDataSource(tempDF[['readID', 'fusionID', 'geneName', 'chrom', 'tStart', 'tEnd']]))

	if '|' in readsFiltered.at[0, 'readID']:
		readsFiltered['isoSupport'] = readsFiltered['readID'].str.split('|', expand=True)[0].astype(int)
		readsFiltered = readsFiltered.sort_values(by='isoSupport')
	readsFiltered = readsFiltered.append(anotsNew) #ones with e
	if len(readsFiltered) > 0:
		source = ColumnDataSource(readsFiltered[['readID','fusionID','geneName','chrom', 'tStart', 'tEnd']])
		sourceT = ColumnDataSource(anotsT[['readID', 'chrom', 'tStart', 'tEnd']])
		readsFiltered = readsFiltered.reset_index(drop=True)
		tools=[BoxZoomTool(dimensions='width'), WheelZoomTool(dimensions='width'), PanTool(), ResetTool()]
		colorList = ["blue", "green", "magenta", "orange", "gold", "pink"]
		figList = []
		plotCount = 2 if len(currChromPoints) == 2 else int(2 + (len(currChromPoints) - 2) / 2)

		for i in range(plotCount):
			plotLoc = i
			if not (i == 0 or i == plotCount - 1):
				plotLoc = i * 2 - 1
			plotWidth = 300 if i > 0 else (len(readsFiltered.loc[0, 'readID']) * 6) + 300
			if i == plotCount - 1: plotLoc = len(currChromPoints) - 1
			if (bounds[i] != [0, 0]):
				if flip[i]:
					bounds[i] = bounds[i][::-1]
				fig = figure(name=('p'+str(i)), title=currChromPoints[plotLoc][0] + " " + str(currChrom[plotLoc]), y_range=list(OrderedDict.fromkeys(readsFiltered['readID'])),
							plot_width = plotWidth, plot_height=len(list(set(readsFiltered['readID'])))*15+220,
							x_range=(bounds[i][0], bounds[i][1]), tools=tools, toolbar_location=None, min_border_bottom=100, min_border_top=100)
				fig.segment(x0='tStart', y0='readID', x1='tEnd', y1='readID', color=colorList[i], line_width=6, source=cdsList[i])
				fig.segment(x0='tStart', y0='readID', x1='tEnd', y1='readID', color=colorList[i], line_width=1,
							source=sourceT)
				fig.segment(x0='tStart', y0='readID', x1='tEnd', y1='readID', color=colorList[i], line_width=1,
							source=cdsTList[i])
				fig.renderers.extend([Span(location=currPoints[plotLoc], dimension='height', line_color='red', line_width=3)])
				if not (i == 0 or i == plotCount - 1):
					fig.renderers.extend([Span(location=currPoints[plotLoc+1], dimension='height', line_color='red', line_width=3)])
				fig.below[0].formatter.use_scientific = False
				fig.add_layout(LinearAxis(formatter=BasicTickFormatter(use_scientific=False), major_label_orientation = pi/4), 'above')
				fig.xaxis.major_label_orientation = pi/4
				fig.ygrid.grid_line_color = None
				fig.y_range.range_padding = 0
				fig.yaxis.visible = False
				figList.append(fig)
			else:
				fig = figure(name=('p'+str(i)), title=currChromPoints[plotLoc][0] + " " + str(currChrom[plotLoc]), y_range=list(OrderedDict.fromkeys(readsFiltered['readID'])),
							plot_width = plotWidth, plot_height=len(list(set(readsFiltered['readID'])))*15+220,
							x_range=(0, 1), tools=tools, toolbar_location=None, min_border_bottom=100, min_border_top=100)
				fig.yaxis.visible = False
				figList.append(fig)
		figList[0].yaxis.visible = True
		figList[-1].toolbar_location = "right"
		return figList

	else:
		return [figure(), figure()]

pageTitle = Div(text='Long read gene fusion visualization', style={'font-size': '24px', 'text-decoration':'underline'})
fusionUploadTitle = Div(text='Fusion file (.tsv)', name='fusionUploadText')
readsUploadTitle = Div(text='Reads file (.bed)')
fusionUpload = FileInput(accept=".tsv", name='fusion')
fusionUpload.on_change('value', upload_table_data)
readsUpload = FileInput(accept=".bed", name='reads')
readsUpload.on_change('value', upload_reads_data)
fastaUploadTitle = Div(text='Fusion fragments aligned to original reads file (.bed, optional)')
fastaUpload = FileInput(accept=".bed", name='fareads')
informUser = Div(text='updates will appear here', name='informUser')
tableSource.selected.on_change('indices', table_click_callback)
buttonSource = ColumnDataSource(data=dict(press=[0], figs=[2]))

def buttonClick(attr):
	hiddenButtonSource.data.update(ColumnDataSource({"selected":[attr]}).data)
	table_click_callback(0, 0, 0)
def saveClick(attr):
	# currLayout = curdoc().get_model_by_name('mainLayout').children
	# curdoc().get_model_by_name()
	# for i in length of plots:
		# export_png(curdoc().get_model_by_name("p"+string(length)), filename="plot0.png")
	export_png(curdoc().get_model_by_name("p0"), filename="plot0.png")
	export_png(curdoc().get_model_by_name("p1"), filename="plot1.png")
	curdoc().get_model_by_name('informUser').text = 'files saved'

#https://stackoverflow.com/questions/34465697/python-bokeh-radio-button-group
buttons = RadioButtonGroup(labels=['Full View', 'Remapped Alignment View', 'Isoform Full Length View'], active=0)
buttons.on_click(buttonClick)


button = Button(label="Save")
button.on_click(saveClick)

subColumn = column([pageTitle, fusionUploadTitle, fusionUpload, readsUploadTitle, readsUpload, fastaUploadTitle, fastaUpload, buttons, button, informUser, data_table], name='subColumn')
mainLayout = row([subColumn, pl, pr], name='mainLayout')
curdoc().add_root(mainLayout)
