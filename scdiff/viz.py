#!/usr/bin/env python

import pdb,sys,os,json
import warnings
warnings.filterwarnings('ignore')
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
import pydiffmap.diffusion_map as pdm


#-----------------------------------------------------------------------
# export data ==>json
def GtoJson(G1,GL,dTD):
	#pdb.set_trace()
	xmatrix=[(item.__dict__)["E"] for item in G1.Cells]
	
	# pca
	pca=PCA(n_components=50)
	pca2=PCA(n_components=2)
	#pdb.set_trace()
	
	# tsne
	xpca2_matrix=pca2.fit_transform(xmatrix)
	xpca_matrix=pca.fit_transform(xmatrix)
	tnse=TSNE(n_components=2)
	xtsne_matrix=tnse.fit_transform(xpca_matrix)
	
	# isomap
	isomap=Isomap(n_components=2)
	isomap_matrix=isomap.fit_transform(xmatrix)
	
	# diffusion map 
	dmk=min(len(xmatrix),10)
	
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		dfmap=pdm.DiffusionMap(n_evecs = 2, epsilon ='bgh', alpha = 0.5,k=dmk)
		dfmap_matrix=dfmap.fit_transform(xmatrix)
		
	CL=[]
	for i in range(len(G1.Cells)):
		jci=(G1.Cells[i]).__dict__
		ixtsne=xtsne_matrix[i]
		ixpca=xpca2_matrix[i]
		ixisomap=isomap_matrix[i]
		idf=dfmap_matrix[i]
		
		jci={item:jci[item] for item in ["ID","T","typeLabel"]}
		jci["TE"]=list(ixtsne)
		jci["PE"]=list(ixpca)
		jci["IE"]=list(ixisomap)
		jci['ME']=list(idf)
		CL.append(jci)
			
	NL=[]
	
	MTL=[item.mT[0] for item in G1.Nodes]
	MIN_MTL=min(MTL)
    
	#DTL=[sum(item.DTA)/len(item.DTA) for item in G1.Nodes]
	#MIN_DTL=min(DTL)
    
	for j in G1.Nodes:
		ij=[G1.Cells.index(item) for item in j.cells]
		try:
			#pdb.set_trace()
			DT=max(0,(G1.Nodes[0].mT[0]-j.mT[0])/(G1.Nodes[0].mT[0]-MIN_MTL))
			#DT=max(0,(DTL[0]-sum(j.DTA)/len(j.DTA))/(DTL[0]-MIN_DTL))
		except:
			DT='NA'
		jcj={'ID':j.ID,'eTF':j.eTF,'T':j.T,'E':j.E,'D':DT,'CELL':ij,"parent": "null" if j.P==None else G1.Nodes.index(j.P),"children": "null" if j.C==[] else [G1.Nodes.index(item) for item in j.C]}
		NL.append(jcj)	
			
	EL=[]
	for j in G1.Edges:
		jfrom=G1.Nodes.index(j.fromNode)
		jto=G1.Nodes.index(j.toNode)
		je={'from':jfrom,'to':jto,'etf':j.etf,'de':j.diffG}
		EL.append(je)
	out=[[item.upper() for item in GL],CL,NL,EL,dTD]
	out="data="+str(out)+'\n'+"data=JSON.stringify(data);"
	return out
	
def getAvgEx(A):
		# get average experssion of Cluster A
		L=len(A.cells[0].E)         # dim of cell expression
		n=len(A.cells)              # number of cells for current cluster
		AE=[]
		for i in range(L):
			iAvg=sum([item.E[i] for item in A.cells])/n
			AE.append(iAvg)
		return AE

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# css template 
def viz(scg_name,G1,output):
	if os.path.exists(output)==False:
		os.mkdir(output)
		
	GL=G1.GL
	dTD=G1.dTD
	
	css_template="""
	
	
	body{
		width:1600px;
	}
	/*div svg */

	#div_svg{
			padding-left:10px;
			padding-top:10px;
			padding-right:10px;
			padding-bottom:10px;
			width:1200px;
			height:1200px;
			margin-left:320px;
			border:1px solid;
			overflow:auto;	
			/*background:#333333;*/
		}
		

	/*div config*/	


	#div_config{
		
		padding-left:10px;
		padding-top:10px;
		padding-right:10px;
		padding-bottom:10px;
		width:300px;
		height:1200px;
		border:1px solid;
		float:left;
		font-size:12px;
		background:silver;
	}

	#div_config button{
		height:15px;
		font-size:12px;
		margin-left:10px;
	}

	#div_config input{
		height:20px;
		font-size:10px;

	}

	#div_config input[type=submit]{
		border-radius:10px;
	}

	#div_config input[type=text]{
		height:15px;
		width:100px
	}
	/* vetrical toolbar starts*/

	nav {
	  font-family: Helvetica, Arial, "Lucida Grande";
	  font-size: 13px;
	  line-height: 1.5;
	  margin: 50px auto;
	  width: 250px;
	  -webkit-box-shadow: 2px 2px 5px rgba(0,0,0,0.2);
		 -moz-box-shadow: 2px 2px 5px rgba(0,0,0,0.2);
			  box-shadow: 2px 2px 5px rgba(0,0,0,0.2);
	}

	.menu-item {
	  background: silver;
	  width: 250px; 
	}

	/*Menu Header Styles*/
	.menu-item h4 {
	  border-bottom: 1px solid rgba(0,0,0,0.3);
	  border-top: 1px solid rgba(255,255,255,0.2);
	  color: #fff;
	  font-size: 12px;
	  font-weight: 500;
	  padding: 7px 12px;
	  background: steelblue;
	 
	}

	.menu-item h4:hover{  
	  background: black;ave javascript script
	}

	.menu-item h4 a {
	  color: white;
	  display: block;
	  text-decoration: none;
	  width: 200px;
	}

	/*ul Styles*/
	.menu-item ul {
	  background: #fff;
	  font-size: 12px;
	  line-height: 30px;
	  height: 0px;
	  list-style-type: none;
	  overflow: hidden;
	  padding: 0px;
	  
	  /*Animation*/
	  -webkit-transition: height 1s ease;
		 -moz-transition: height 1s ease;
		   -o-transition: height 1s ease;
		  -ms-transition: height 1s ease;
			  transition: height 1s ease;
	}


	.menu-item ul {
	  height: 93px;
	  
	}
	#globalconfig.menu-item ul{
		height:203px;
	}
	#downloadconfig.menu-item ul{
		height:240px;
	}
	
	#cellplot.menu-item  ul{
		height:100px;
	}
	#tfconfig.menu-item ul{
		height:140px;
	}
	#geneconfig.menu-item ul{
		height:160px;
	}

	/*li Styles*/

	.menu-item li {

	  border-bottom: 1px solid #eee;
	  margin-left: 10px;
	  text-decoration: none;
	  color: #black;
	  display: block;
	  width: 250px;

	}

	/*
	.menu-item li:hover {
	  background: #eee;
	}
	*/
	.menu-item li a:hover {
	  background: #eee;
	}

	/*First Item Styles*/
	.alpha p {
		padding: 8px 12px;
		color: #aaa;
	}

	input .checkbox{
	  width: 13px;
	  height: 13px;
	  padding: 0;
	  margin:0;
	  vertical-align:middle;
	  position: relative;
	  overflow: hidden;
	}

	/* vetrical toolbar ends*/

	/*bar chart */

	.axis path,
	.axis line {
	  fill: none;
	  stroke: #000;
	  shape-rendering: crispEdges;
	  color:red;
	}
	.axis{
		color: red;
	}
	
	
	"""

	#-----------------------------------------------------------------------
	# HTML template
	HTML_template="""

	<!DOCTYPE html>
	<html lang="en">
	  <head>
		<meta charset="utf-8">

		<title>Single Cell Data Visualization</title>
		<!--canvg library-->
		<link rel="stylesheet" type="text/css" href="style.css">
		<script src="https://code.jquery.com/jquery-1.11.3.min.js"></script>
		<script src="https://www.lactame.com/lib/ml/2.2.0/ml.min.js"></script>
		<script src="http://d3js.org/d3.v3.min.js"></script>
		<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.7.1/Chart.bundle.min.js"></script>
		<script src="parseJSON.js"></script>
		<script type="text/javascript" src="%s"></script>
		
		<script>
		var fold=function(divid){
			var x=document.getElementById(divid);
			if (x.style.display === 'none') {
				x.style.display = 'block';
			} else {
				x.style.display = 'none';
			}

		}

		</script>
		
	  </head>
	  <body onload="javascript:scviz.onload()">
		  <!-- This DIV contains static visualization configuration config-->
		  <div id="div_config">
			  <div class="tf">
				<nav>
					<div class="menu-item"><h4><a onclick="fold('globalconfig')">GLOBAL CONFIG</a></h4></div>  
					<div class="menu-item" id="globalconfig">
					  <ul>
						 <li>
						   <span>Reset: <input type="submit", onclick="scviz.config.resetconfig('zoomsliderbar','bgcolor')" value="RESET"></span><br>
							  <input type="range" id="zoomsliderbar" value="50" min="0" max="100" onchange="scviz.config.zoom('zoomsliderbar',this.value)"></input>
							  <span id="zoomslider">50</span>
							</li>
						 <li>
							 <label title="Enable/Disable the piechart visualization on mouseover">Enable/Disable tooltip:<input class="checkbox" id="tooltipcheck" type="checkbox"  checked> </input></label>
							 </li>
						 <li>
								<input type="submit" onclick="scviz.config.setbgcolor('bgcolor')" value="Set Background Color" title="change/set the background color the figure">
								<input type="color" id="bgcolor" value="#333333"><br>
								<input type="submit" onclick="scviz.config.settextcolor('textcolor')" value="Set Text Color" title="change/set the color of the text">
								<input type="color" id="textcolor" value="#ffffff"><br>
								<input type="submit" onclick="scviz.config.setpathcolor('pathcolor')" value="Set Edge Color" title="change/set the color of the edge">
								<input type="color" id="pathcolor" value="#ffffff">
							</li>
					  </ul>
					</div>
					
					<div class="menu-item"><h4><a onclick="fold('cellplot')">CELL PLOTS</a></h4></div>  
					<div class="menu-item" id="cellplot">
					  <ul>
						<li>
							<input type="submit" value="Plot Cells" onclick='scviz.vizCells("viztypep")'> <br>
							<input type="radio" name="viztype" value="tsne" id="viztypet" checked> T-SNE
							<input type="radio" name="viztype" value="pca" id="viztypep">  PCA <br>
							<input type="radio" name="viztype" value="dfmap" id="viztyped"> Diffusion Map
							<input type="radio" name="viztype" value="isomap" id="viztypem">  ISOMAP <br>
							</li>
						<!--<form action="/cgi-bin/colorCells.py" target="/cgi-bin/images/">
							<input type="submit" value="ColorCellsByGene" >
							<input type="text" name="cellGeneName" value="">
							</form>
							-->
					  </ul>
					</div>
					
					
					<div class="menu-item"><h4><a onclick="fold('tfconfig')">TF CONFIG</a></h4></div>   
					<div class="menu-item" id="tfconfig">
					  <ul>
						<li><label title="show top ranked 20 regulating TFs (if have) for each edge">Show/Hide TF for each path: <input class="checkbox" id="tfcheck" type="checkbox"  onchange="scviz.explore.showhideTF(checked)"> </input></label> </li>
						<li>
							<input type="submit" value="Explore TF" onclick="scviz.explore.exploretf('tfName')">
							<input type="text" id="tfName" value="" onkeydown="if (event.keyCode == 13) { scviz.explore.exploretf('tfName') }">
							</li>
						
						<li><label title="show top ranked 20 regulating eTFs (if have) for each edge">Show/Hide eTF for each path: <input class="checkbox" id="rtfcheck" type="checkbox"  onchange="scviz.explore.showhideRTF(checked)"> </input></label> </li>
						<li>
							<input type="submit" value="Explore eTF" onclick="scviz.explore.explorertf('rtfName')">
							<input type="text" id="rtfName" value="" onkeydown="if (event.keyCode == 13) {scviz.explore.explorertf('rtfName') }">
							</li>
							
					
					  </ul>
					</div>
					  
					<div class="menu-item"><h4><a onclick="fold('geneconfig')">GENE CONFIG</a></h4></div>  
					<div class="menu-item" id="geneconfig">
					  <ul>
						<li><label title="show top ranked 20 (both up and down regulated differentially expressed genes for each edge)">Show/Hide DE genes for each path:<input class="checkbox" id="genecheck" type="checkbox"  onchange="scviz.explore.showhideDE(checked)"> </input></label></li>
						<li>
							<input type="submit" value="Explore DE gene" onclick="scviz.explore.explorede('deName')">
							<input type="text" id="deName" value="" onkeydown="if (event.keyCode == 13) { scviz.explore.explorede('deName') }">
							</li>
						
						<li>
							<input type="submit" value="Find DE genes between:" onclick="scviz.explore.exploredebetween('node1dropdowndiv','node2dropdowndiv')"> <br>
							<div id="node1dropdowndiv">Node1: </div>
							<div id="node2dropdowndiv">Node2: </div>
							<!--<label >Node 1: </label><input type="text" id="denode1"><br>-->
							<!--<label >Node 2: </label><input type="text" id="denode2"> -->
							</li>
						
					  </ul>
					</div>
					  
					<div class="menu-item"><h4><a onclick="fold('downloadconfig')">DOWNLOAD</a></h4></div>  
						<div class="menu-item" id="downloadconfig">
						  <ul>
							<li><label><input type="submit" onclick="scviz.download.downloadfig('downloadlink')" value="Generate Figure:"></label></li>
							<li><a id="downloadlink"></a></li>
							<li><label><input type="submit" onclick="scviz.download.downloadjson('jsondownloadlink')" value="Generate Json download file:"></label></li>
							<li><a id="jsondownloadlink"></a></li>
							<li><label><input type="submit" onclick="scviz.download.downloadtf('tfdownloadlink')" value="Generate TF download file:"></label></li>
							<li><a id="tfdownloadlink"></a></li>
							<li><label><input type="submit" onclick="scviz.download.downloadde('dedownloadlink')" value="Generate DE download file:"></label></li>
							<li><a id="dedownloadlink"></a></li>
					  </ul>
					  </div>
					</div>
				</nav>
			  </div>
			  
			  <div class="gene">
				  </div>
		  </div>
		  <!-- This DIV contains dynamic data-driven svg elements -->
			<div id="div_svg"> </div>
		</body>
	</html>

	"""

	#-----------------------------------------------------------------------
	# javascript
	javascript="""


(function(scviz, $, undefined){
	//global parameters
	textfontsize=11;
	textfontcolor="#ffffff";
	svgWidth=1200;
	svgHeight=1200;
	defaultbgcolor="#333333";
	pathcolor="#ffffff";
	pathtextsize=8;
	
	margin = {top: 100, right: 20, bottom: 250, left: 20};
	colorList=getColors();
	colorList=Object.values(colorList);
	
	// public functions-------------------------------------------------
	
	scviz.onload=function(){
		RL = parseJSON(data);
		nodes = scviz.buildTree(RL[2]);
		root = nodes[0]; 
		GL = RL[0];    // list of genes
		cells = RL[1]; // list of cells
		edges = RL[3]; // list of edges
		dTD = RL[4];
		scviz.drawTree(root);
		updateNodes(nodes);
		//alert("done");
	};
	
	// visualize all the cells in 2D space projected by t-sne and pca
	scviz.vizCells=function(viztypep){
		//var plotType=document.getElementById(viztypep).checked;
		var plotType=$("input[name=viztype]:checked").val();
		var xData=[];
		var xLabels=[];
		var xTypeLabels=[];
		for (var inode of nodes){
			var inode_id=inode["ID"];
			var inode_cells=inode["CELL"];
			var inode_cells=inode_cells.map(function(d){return cells[d]});
			for (var cell of inode_cells){
				if (plotType=="pca"){
					xData.push(cell.PE);
				}else if (plotType=="isomap"){
					xData.push(cell.IE);
				}else if (plotType=="dfmap"){
					xData.push(cell.ME)
				}else {
					xData.push(cell.TE);
				}
				xLabels.push(inode_id);
				xTypeLabels.push(cell.typeLabel);
			}
		}
		scviz.plots.scatterplot(xData,xLabels,xTypeLabels);
	};
	
	// scviz plots : scatter, line, bar 
	scviz.plots={
		scatterplot: function(xData,xLabels,xTypeLabels){
			var dlabels={};
			for (var ix in xData){
				var di=xData[ix];
				var li=xLabels[ix];
				if (!(li in dlabels)){
					dlabels[li]=[{"x": di[0],"y": di[1]}];
				}else{
					dlabels[li].push({"x": di[0],"y": di[1]});
				}
			}
				
			var CombinedDataSet=[];
			var colorCounter=0;
			for (var di in dlabels){
				var xlabel=di;
				var xList=dlabels[di];
				CombinedDataSet.push({"label": xlabel, "data": xList, "backgroundColor": colorList[colorCounter]});
				colorCounter+=1;
			}
			
			var plotdata={"datasets":CombinedDataSet};
			
			var newW3 = open('','_blank','height=1200,width=1400')
			newW3.document.write('<head><title>scatter plot</title> </head><body></body>');
			d3.select(newW3.document.body).append("canvas")
						.attr("id","scatterplot")
						.attr("width","1200px")
						.attr("height","1000px");
						
			var ctx=newW3.document.getElementById("scatterplot").getContext('2d');
			
			
			var scatterChart = new Chart(ctx, {
				type: 'scatter',
				data: plotdata,
				options: {
				responsive : false,
				scales: {
					xAxes: [{
						type: 'linear',
						position: 'bottom'
					}]
				},
				tooltips: {
					callbacks:{
						label: function(tooltipItem,data){
							var pos=[tooltipItem.xLabel,tooltipItem.yLabel];
							for (var ix in  xData){
								if (xData[ix].toString()==pos.toString()){
									return xTypeLabels[ix];
								}							
							}
						}
					}
				}
			}
				
			});	
		},
		//line plot
		lineplot:function(X,Y,Z){
			//X: x-axis values
			//Y: y-axis values
			//Z: legend
			
			var newW3 = open('','_blank','height=600,width=800');
			newW3.document.write('<head><title>line plot</title> </head><body></body>');
			d3.select(newW3.document.body).append("canvas")
						.attr("id","lineplot")
						.attr("width","700px")
						.attr("height","500px");
						
			var ctx=newW3.document.getElementById("lineplot").getContext('2d');
			
			var lineChart = new Chart(ctx, {
				type: 'line',
				data:{
					labels: X,
					datasets:  [{
						data: Y,
						fill: false,
						label: ""
					}]
				},
				options: {
					responsive : false,
					title: {
						display: true,
						text: Z
					}
				
				}
			});
		},
		barplot: function(X,Y,Z,P){
			//X: x-axis values
			//Y: y-axis values
			//Z: legend
			
			var newW3 = open('','_blank','height=600,width=800,left='+P);
			newW3.document.write('<head><title></title> </head><body></body>');
			d3.select(newW3.document.body).append("canvas")
						.attr("id","lineplot")
						.attr("width","700px")
						.attr("height","500px");
						
			var ctx=newW3.document.getElementById("lineplot").getContext('2d');
			
			var lineChart = new Chart(ctx, {
				type: 'bar',
				data:{
					labels: X,
					datasets:  [{
						data: Y,
						backgroundColor: "blue",
						label: ""
					}]
				},
				options: {
					responsive : false,
					title: {
						display: true,
						text: Z
					}
				
				}
			});
		},
		
	}
	
	// prepare the data to build the tree 
	scviz.buildTree = function(nodes){
		var xc;
		for (x of nodes){
			if(x.parent!="null"){x.parent=nodes[x.parent];}
			xc=[];
			if (x.children!="null"){
				for (y of x.children){
					xc.push(nodes[y]);
				}
			}else{
				xc=null;
			}
			x.children=xc;
		}
		return nodes;	
	};
	
	// build svg tree diagram
	scviz.drawTree=function (root){
		width = svgWidth - margin.right - margin.left;
		height =svgHeight - margin.top - margin.bottom;
		var i = 0;
		//var tree = d3.layout.tree();
		tree=d3.layout.tree();
		tree.size([width,height]);
		
		var diagonal = d3.svg.diagonal()
			.projection(function(d) { return [d.x, d.y]; });
		
		//default bgcolor
		var default_bgcolor=defaultbgcolor;
		var bgcolor;
		var input_bgcolor=d3.select("#bgcolor").value;
		colorisOK  = /(^#[0-9A-F]{6}$)|(^#[0-9A-F]{3}$)/i.test(input_bgcolor)
		if (colorisOK){
			bgcolor=input_bgcolor;
		}else{
			bgcolor=default_bgcolor;
		}
		
		// create svg for visualization
		svg = d3.select("#div_svg").append("svg")
			.attr("width", svgWidth)
			.attr("height",svgHeight)
			.attr("id","svg")
			.style("background",bgcolor)
			.attr("viewBox","0 0 "+svgWidth+" "+svgHeight)
			.attr("preserveAspectRatio", "none")
			.append("g")
			.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
			
		var allnodes=tree.nodes(root);
		var links=tree.links(allnodes);
		 
		// Declare the nodes
		var node = svg.selectAll("g.node")
		  .data(nodes)
		  .data(nodes, function(d) { return d.id = ++i; });

		// Enter the nodes.
		var nodeEnter = node.enter().append("g")
		  .attr("class", "node")
		  .attr("transform", function(d) { 	 
			  return "translate(" + d.x + "," + d.y + ")"; })

		nodeEnter.append("circle")
		  .attr("r", 16)
		  .attr("fill","#fff")
		  .attr("stroke","steelblue")
		  .attr("stroke-width","3px")
		  .on("mouseover",inNode)  
		  .on("mouseout",outNode)
		  .on("click",onClickfunction);

		var textEnter=nodeEnter.append("text")
		  .attr("dy", "0em")
		  .attr("class","nodetext")
		  .attr("fill",textfontcolor)
		  .style("font-size",textfontsize)
		  .attr("text-anchor", function(d) { 
			  return  "start"; })
		  .style("fill-opacity", 1)
		  .call(addwraptext);

		//creat TF display svg text 
		createEdgeTF();
		
		//create RTF display svg text
		if (nodes[0].eTF!=undefined){
			createRTFs();
		}
		
		//create DE display SVG text
		createDE();
		//create colorbar 
		createColorBar();
		
		
		// Declare the links
		var link = svg.selectAll("path.link")
		  .data(links, function(d) { return d.target.id; })
		  .enter();

		// Enter the links.
		link.insert("path", "g")
		  .attr("class", "link")
		  .attr("stroke",pathcolor)
		  .attr("stroke-width","2px")
		  //.fill("none")
		  .attr("fill","none")
		  .attr("d", diagonal)
		  .attr("id",function(d,i){return "s"+i;})
		  .attr("color","black")
		  //.on("click",testClick);
		  //.on("mouseout", mouseoutfunction)

		// text along the path
		link.append("text")
		.append("textPath")
		.attr("class","linktext")
		.attr("xlink:href",function(d,i){return "#s"+i;})
		.attr("startOffset","20%")
		.attr("dy","-1em");

		// Define the div for the tooltip
		tooltipsvg = d3.select("body").append("div")	
			.attr("class", "tooltip")
			.append("svg");
			
		var NodeIDList=nodes.map(function(d){return "E"+d.T+"_"+d.ID;});
		// create dropdown
		createDropDown(NodeIDList[0],"#node1dropdowndiv","node1dropdown",NodeIDList,null);
		createDropDown(NodeIDList[1],"#node2dropdowndiv","node2dropdown",NodeIDList,null);
	};
	
	//config functions
	scviz.config={
		resetconfig: function(zoomsliderbarid,bgcolorid){
			scviz.config.zoom(zoomsliderbarid,50);
			document.getElementById(zoomsliderbarid).value=50;
			document.getElementById(bgcolorid).value="#333333";
			scviz.config.setbgcolor(bgcolorid);
			resetPath();
			location.reload();
		},
		setbgcolor: function(bgcolorid){
			var color=document.getElementById(bgcolorid).value;
			d3.select("svg").style("background", color);
		},
		settextcolor: function(textcolorid){
			var color=document.getElementById(textcolorid).value;
			textfontcolor=color;
			d3.selectAll(".nodetext")
			.attr("fill",textfontcolor);
			
		},
		setpathcolor: function(pathcolorid){
			var color=document.getElementById(pathcolorid).value;
			pathcolor=color;
			d3.selectAll(".link")
			.attr("stroke",pathcolor);
		},
		//zoom slider bar action
		zoom: function(zoomsliderid,newValue){
			document.getElementById(zoomsliderid).innerHTML=newValue;
			var wd=svgWidth;
			var ht=svgHeight;
			var sv=50;
			var zx=newValue/sv;
			var newwd=wd*zx;
			var newht=ht*zx;
			
			d3.select("#div_svg").select("svg")
			.attr("width",newwd)
			.attr("height",newht)
			.attr("preserveAspectRatio", "none");
		}

	};
	
	//scviz explore functions
	scviz.explore={
		//handle explore tf
		//handle explore tf
		exploretf: function(tfnameid){
			resetPath();
			var tfnameid="#"+tfnameid;
			var tfinput=d3.select(tfnameid).property("value").toUpperCase();
			var paths=d3.select("svg").selectAll("path")[0];
			var selectedpaths=[];
			var tnodeList=[];
			var etfList=[];
			for (var edge of edges){
				var etf=edge.etf.map(function(d){return d[1];});
				if(etf.indexOf(tfinput)!=-1){
					tnodeList.push(nodes[edge.to]);
					etfList.push(etf.indexOf(tfinput));
				}
			}
			
			
			
			var tpathList=[];
			for (var path of paths){
				var pt=path.__data__.target;
				if (tnodeList.indexOf(pt)!=-1){
					index_pt=tnodeList.indexOf(pt);
					tpathList.push([path,etfList[index_pt]]);
				}
			}
				
			
			for (var path of tpathList){
				var colorpath=hsl_col_perc(path[1]/60,0,120);
				d3.select(path[0]).attr("stroke",colorpath);
			}	
			
			
			//path text
			var textpaths=d3.select("svg").selectAll(".linktext")[0];
			var tpathIndex=tpathList.map(function(d){return paths.indexOf(d[0]);});
			var ttextpath=tpathIndex.map(function(d){return textpaths[d];});
			for (i in ttextpath){
				var tfi=etfList[i];
				var tpath=ttextpath[i];
				d3.select(tpath).text("p-value rank:"+tfi)
				.attr("fill","#fff");
			}
			plottf(tfinput);
		},
		
		//handle explore rtf
		explorertf: function(rtfnameid){
			var searchCut=40;
			resetPath();
			rtfnameid="#"+rtfnameid;
			var rtfinput=d3.select(rtfnameid).property("value").toUpperCase();
			var paths=d3.select("svg").selectAll("path")[0];
			
			
			
			var selectedpaths=[];
			var tnodeList=[];
			var etfList=[];
			for (var edge of edges){
				var rtf=nodes[edge.to].eTF
				rtf=rtf.slice(0,searchCut);
				var rtfList=rtf.map(function(d){return d[1].toUpperCase();});
				if (rtfList.indexOf(rtfinput)!=-1){
					tnodeList.push(nodes[edge.to]);
					etfList.push(rtfList.indexOf(rtfinput));
				}
			}
			var tpathList=[];
			for (var path of paths){
				var pt=path.__data__.target;
				if (tnodeList.indexOf(pt)!=-1){
					index_pt=tnodeList.indexOf(pt);
					tpathList.push([path,etfList[index_pt]]);
				}
			}
			
			for (var path of tpathList){
				var colorpath=hsl_col_perc(path[1]/60,0,120);
				d3.select(path[0]).attr("stroke",colorpath);
			}	
			
			//path text
			var textpaths=d3.select("svg").selectAll(".linktext")[0];
			var tpathIndex=tpathList.map(function(d){return paths.indexOf(d[0]);});
			var ttextpath=tpathIndex.map(function(d){return textpaths[d];});
			for (i in ttextpath){
				var tfi=etfList[i];
				var tpath=ttextpath[i];
				d3.select(tpath).text("RTF p-value rank:"+tfi)
				.attr("fill","#fff");
			}
			if (GL.indexOf(rtfinput)!=-1){		
				plottf(rtfinput);
			}
			
			//console.log(rtfinput);
		},
		//explore de
		explorede: function(denameid){
			resetPath();
			denameid="#"+denameid;
			var deinput=d3.select(denameid).property("value").toUpperCase();
			var paths=d3.select("svg").selectAll("path")[0];
			var selectedpaths=[];
			var tnodeList=[];
			var deList=[];
			for (var edge of edges){
				var de=edge.de.map(function(d){return d.toUpperCase();});
				if(de.indexOf(deinput)!=-1){
					tnodeList.push(nodes[edge.to]);
					deList.push(de.indexOf(deinput)/de.length);
				}
			}
			
			var tpathList=[];
			for (var path of paths){
				var pt=path.__data__.target;
				if (tnodeList.indexOf(pt)!=-1){
					index_pt=tnodeList.indexOf(pt);
					tpathList.push([path,deList[index_pt]]);
				}
			}
				
			
			
			for (var path of tpathList){
				var colorpath=hsl_col_perc(path[1],0,120);
				d3.select(path[0]).attr("stroke",colorpath);
			}	
			
			//path text
			var textpaths=d3.select("svg").selectAll(".linktext")[0];
			var tpathIndex=tpathList.map(function(d){return paths.indexOf(d[0]);});
			var ttextpath=tpathIndex.map(function(d){return textpaths[d];});
			for (i in ttextpath){
				//var dei=deList[i];
				var dei=tpathList[i][1];
				var tpath=ttextpath[i];
				d3.select(tpath).text("DE gene top "+((dei*100).toFixed(1))+"%")
				.attr("fill","#fff");
			}
			if (GL.indexOf(deinput)!=-1){		
				plottf(deinput);
			}
			//alert("explore de");
		},			
		//public hanlde exploredebetween
		exploredebetween: function(node1dropdownid,node2dropdownid){
			//var xdenode1=d3.select("#denode1").property("value").toUpperCase();
			//var xdenode2=d3.select("#denode2").property("value").toUpperCase();
			
			var xdenode1=document.getElementById('node1dropdown').value.toUpperCase();
			var xdenode2=document.getElementById('node2dropdown').value.toUpperCase();
			
			var denode1=xdenode1.split("_").slice(1);
			denode1=denode1.join("_");
			
			var denode2=xdenode2.split("_").slice(1);
			denode2=denode2.join("_");
			
			//var denode1=xdenode1;
			//var denode2=xdenode2;
			
			var node1;
			var node2;
			for (var node of nodes){
				if (node.ID==denode1){
					node1=node;
				}
				if (node.ID==denode2){
					node2=node;
				}
			}
			
			var DEGList=[]
			var decut=0.6;
			var topcut=800;
			for (var gi in node1.E){
				var fc=node2.E[gi]-node1.E[gi];
				fc=fc.toFixed(3);
				if (Math.abs(fc)>decut){
					DEGList.push([GL[gi],node1.E[gi].toFixed(3),node2.E[gi].toFixed(3),fc]);	
				}
			} 
			DEGList.sort(function(a,b){
				return b[3]-a[3];
			});
			var upList=DEGList.slice(0,topcut);
			var upList=upList.filter(function(d){return d[3]>0});
			upList.unshift(["Gene","E"+node1.T+"_"+node1.ID,"E"+node2.T+"_"+node2.ID,"log2 fold change"]);
			
			
			DEGList.sort(function(a,b){
				return a[3]-b[3];
			});
			var downList=DEGList.slice(0,topcut);
			var downList=downList.filter(function(d){return d[3]<0});
			downList.unshift(["Gene","E"+node1.T+"_"+node1.ID,"E"+node2.T+"_"+node2.ID,"log2 fold change"]);
			
			//create table for DEG between 2 given ndoes
			//console.log(DEGList);
			
			var newW = open('','_blank','height=600,width=800,left=200,top=200,scrollbars=yes')
			newW.document.write("<head><title>Differentially expressed genes between 2 given nodes</title></head><body></body>");
			newW.document.title="loading...";

			//adding TF details
			var tdiv=d3.select(newW.document.body)
			.style("background","white")
			.append("div")
			.style("padding-left","50px")
			.style("padding-top","50px")
			.attr("width",400)
			.attr("height",600)
			.attr("class","deg_between");
			
		
			tdiv
			.append("p")
			.text("Table: Differentially expressed genes (UP-regulated) between : "+xdenode1+" and "+xdenode2);
			
			var species="MOUSE"
			//add GO analysis
			tdiv
			.append("button")
			.text("functional analysis")
			.on("click",function(){
				if (d3.event.shiftKey){
					toppgenegoInput(upList.map(function(d){return d[0];}));
				}else{
					panthergoInput(upList.map(function(d){return d[0];}), species);
				}
			});
				
			//add download link
			tdiv.append("button")
			.text("click to create the excel table")
			.on("click",function(){
					table2XLS(newW,"deg_between","degupdlink");
				});
			
			tdiv.append("a")
			.attr("href","#")
			.attr("id","degupdlink");
			
			tdiv.append("p");
			//add table main body		
			createTable(tdiv,"deg_between",upList);
			
			tdiv
			.append("p")
			.text("Table: Differentially expressed genes (Down-regulated) between : "+xdenode1+" and "+xdenode2);
			
			//add GO analysis
			tdiv
			.append("button")
			.text("functional analysis")
			.on("click",function(){
				if (d3.event.shiftKey){
					toppgenegoInput(downList.map(function(d){return d[0];}));
				}else{
					panthergoInput(downList.map(function(d){return d[0];}), species);
				}
			});
			
			//add download link
			tdiv.append("button")
			.text("click to create the excel table")
			.on("click",function(){
					table2XLS(newW,"deg_between","degdowndlink");
				});
			
			tdiv.append("a")
			.attr("href","#")
			.attr("id","degdowndlink");
			
			tdiv.append("p");
			createTable(tdiv,"deg_between",downList);
			
				
			
			
		},
		
		//show/hide TF
		showhideTF: function(checked){
			
			if (checked){
				//pre-check
				var genecheck=d3.select("#div_config").select("#genecheck").property("checked");
				if(genecheck){
					hideDE();
					d3.select("#div_config").select("#genecheck").property("checked",false)
				}
				showTF();
			}else{
				hideTF();
			}
		},
		//show hide DE
		showhideDE: function(checked){
			if (checked){
				//pre-check
				var tfcheck=d3.select("#div_config").select("#tfcheck").property("checked");
				if (tfcheck){
					hideTF();
					d3.select("#div_config").select("#tfcheck").property("checked",false);
				}
				showDE();
			}else{
				hideDE();
			}
		},
		//show hide RTF
		showhideRTF:function(checked){
			
			if (checked){
				//pre-check
				var rtfcheck=d3.select("#div_config").select("#rtfcheck").property("checked");
				if(rtfcheck){
					hideRTF();
				}
				showRTF();
			}else{
				hideRTF();
			}
		}
	};

	//scviz download functions
	scviz.download={
		//download json
		downloadjson: function(jsondownloadlinkid){
			plswait(jsondownloadlinkid);
			window.setTimeout(function(){createjsondownload(jsondownloadlinkid);},10);
		},
		//download fig
		downloadfig: function(downloadlinkid){
			svgToCanvas(downloadlinkid);

		},
		//download de
		downloadde: function(dedownloadlinkid){
			plswait(dedownloadlinkid);
			window.setTimeout(function(){creatededownload(dedownloadlinkid);},10);
		},
		//download tf
		downloadtf: function(tfdownloadlinkid){
			plswait(tfdownloadlinkid);
			window.setTimeout(function(){createtfdownload(tfdownloadlinkid)},10);
		}

		
	};

	
	//private functions-------------------------------------------------
	function plottf(tfinput){
		if (GL.indexOf(tfinput)!=-1){		
			//plottf(tfinput);
			scviz.plots.barplot(nodes.map(function(d){return "E"+d.T+"_"+d.ID;}),nodes.map(function(d){return (d.E[GL.indexOf(tfinput)]).toFixed(3);}),"Expression of "+tfinput);
			var targets=dTD[tfinput];
			var targets=targets.filter(function(d){
				if (GL.indexOf(d)!=-1){
					return true;
				} return false;
			});
			
			var targetIndex=targets.map(function(d){
					return GL.indexOf(d);
				})
			
			var nodeIDList=[];
			var nodeValueList=[];
			var fcut=0.5;
			for (var node of nodes){
				if (node.endEdge){
					nodeIDList.push("E"+node.T+"_"+node.ID);
					var detargets=[];
					var tfde=(node.fc)[tfinput];
					for (var dt of targets){
						if (Math.abs(node.fc[dt])>fcut){
							detargets.push([dt,(node.fc)[dt]]);	
						}
					}
					if (detargets.length>0){
						var deE=detargets.map(function(d){return d[1];});
						var sumdeE=deE.reduce(function(a,b){return a+b;});
						var sumdeE=(sumdeE/detargets.length).toFixed(3);
					}else{
						var sumdeE=0;
					}
					nodeValueList.push(sumdeE);
					
				}
				
				
			}
			scviz.plots.barplot(nodeIDList,nodeValueList,"Average targets log2 fold change of "+tfinput+"(only considered differential targets, the log2 fold change is calculated along the edge ending at given node)",500);
			//console.log(targets);
		}
	}
	
	//function on testClick
	function onClickfunction(){
		var newW = open('','_blank','height=600,width=800,left=200,top=200,scrollbars=yes')
		newW.document.write("<head><title>Edge TF/Gene</title></head><body></body>");
		newW.document.title="loading...";

		//adding TF details
		var tdiv=d3.select(newW.document.body)
		.style("background","white")
		.append("div")
		.style("padding-left","50px")
		.style("padding-top","50px")
		.attr("width",400)
		.attr("height",600)
		.attr("class","div_table_edge");
		
		
		cnode=this.__data__;
		// cell details at selected node 
		var CellIDs=[["ID","Label"]];
		for (var icell of cnode.CELL){
			var cid=cells[icell].ID;
			var clabel=cells[icell].typeLabel;
			CellIDs.push([cid,clabel]);
		}
		
		tdiv
		.append("p")
		.text("Table 0. Cells at selected node: "+cnode.T+"_"+cnode.ID);
		
		
		createTable(tdiv,"celltable",CellIDs);
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"celltable","celldlink");
			});
			
		tdiv.append("a")
		.attr("href","#")
		.attr("id","celldlink");
		
		
		//TF details for edge ending at selected node
		pnode=cnode.parent;
		var chosenEdge;
		for (var edge of edges){
			if (nodes[edge.to]==cnode){
				chosenEdge=JSON.parse(JSON.stringify(edge));
			}
		}
		
		var resList=[]; //TF details
		resList.push(["TF","p-value","TF log2 fold change","Mean target log2 fold change","Mean DE target log2 fold change","edgeID"]);

		resList=exportTFEdge(resList,chosenEdge);
		
		tdiv
		.append("p")
		.text("Table 1. TFs(Predicted using the expression of the TF targets) for the edge ending at selected node:"+cnode.T+"_"+cnode.ID);
		
		
		createTable(tdiv,"tftable",resList);
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"tftable","tfdlink");
			});
			
		tdiv.append("a")
		.attr("href","#")
		.attr("id","tfdlink");
		
		
		////TF details for path ending at selected node
		cnode=this.__data__;
		
		/*
		pnode=cnode.parent;
		var chosenPath=[]
		
		while(pnode!="null"){
			var ec=cnode.endEdge;
			cnode=pnode;
			pnode=pnode.parent;
			chosenPath.push(ec);
		}
		*/
		
		
		var AllTFList=[];
		AllTFList.push(["TF","p-value","log2 fold change to parent","log2 fold change to siblings ","edgeID"]);
		
		for (var etf of cnode.eTF){
			AllTFList.push([etf[1],etf[0].toExponential(2),etf[2].toFixed(2),etf[3].toFixed(2),'E'+pnode.T+"_"+pnode.ID+"->E"+cnode.T+"_"+cnode.ID])
		}
	
		AllTFList.sort(function(a,b){return a[1]-b[1];});
		tdiv
		.append("p")
		.text("Table 2. eTFs(Predicted using the expression of TFs) for the edge ending at selected node:"+cnode.T+"_"+cnode.ID);
		
		createTable(tdiv,"alltftable",AllTFList);
		
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"alltftable","alltfdlink");
			});
			
		tdiv.append("a")
		.attr("href","#")
		.attr("id","alltfdlink")
		.text("");
		
		//-------------------------------------------------------------
		//adding DE details
		cnode=this.__data__;
		pnode=cnode.parent;
		var resDEList=[]; //DE details
		//resDEList.push(["DE gene","Expression_E"+pnode.T+"_"+pnode.ID,"Expression_"+cnode.T+"_"+cnode.ID,"Fold change","edgeID"]);
		
		var resDEList=exportDEEdge(resDEList,chosenEdge);
		
		
		var species="MOUSE"
		var functionalCut=800;
		
		// up-regulated 
		tdiv.append("p")
		.text("Table 3. DE gene (up-regulated) details for edge ending at selected node:"+cnode.T+"_"+cnode.ID)
	
		var upresDEList=resDEList.filter(function(d){
			if (d[3]>0){
				return true;
			} return false;
		});
		upresDEList.unshift(['DE gene','E'+pnode.T+"_"+pnode.ID+" expression",'E'+cnode.T+"_"+cnode.ID+" expression",'log2 fold change','edgeID'])
		
		createTable(tdiv,"detable",upresDEList);
		tdiv.append("button")
		.text("functional analysis")
		.on("click",function(){
			if (d3.event.shiftKey){
					toppgenegoInput(upresDEList.slice(0,functionalCut).map(function(d){return d[0];}));
				}else{
					panthergoInput(upresDEList.slice(0,functionalCut).map(function(d){return d[0];}), species);
				}
		});
		
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"detable","dedlink");
			});
		tdiv.append("a")
		.attr("href","#")
		.attr("id","dedlink");
		
		
		// down-regulated
		
		tdiv.append("p")
		.text("Table 4. DE gene (down-regulated) details for edge ending at selected node:"+cnode.T+"_"+cnode.ID)
		
		var downresDEList=resDEList.filter(function(d){
			if (d[3]<=0){
				return true;
			} return false;
		});
		
		downresDEList.unshift(['DE gene','E'+pnode.T+"_"+pnode.ID+" expression",'E'+cnode.T+"_"+cnode.ID+" expression",'log2 fold change','edgeID'])
		
		createTable(tdiv,"detable",downresDEList);
		tdiv.append("button")
		.text("functional analysis")
		.on("click",function(){
			if (d3.event.shiftKey){
					toppgenegoInput(downresDEList.slice(0,functionalCut).map(function(d){return d[0];}));
				}else{
					panthergoInput(downresDEList.slice(0,functionalCut).map(function(d){return d[0];}), species);
				}
		});
		
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"ddetable","ddedlink");
			});
		tdiv.append("a")
		.attr("href","#")
		.attr("id","ddedlink");
		
		//loading complete
		newW.document.title="loading complete";
		
	}
	
	// function on mouse out
	function outNode(){
		tooltipsvg
		.style("opacity", 0)
		.attr("width",0)
		.attr("height",0)
		.selectAll("*").remove();;
	}
	
	//private function on mouse over
	function inNode(){
		if (d3.select("#tooltipcheck").property("checked")==false){
			return false;
		}
		var CELL=this.__data__.CELL;
		var ctdict={}
		for(x of CELL){
			xlabel=cells[x].typeLabel;
			if (xlabel in ctdict==false){
				ctdict[xlabel]=1;
			}else{
				ctdict[xlabel]+=1;
			}
		}

		var ltdict=[]
		var asum=0;
		for (key in ctdict){
			asum+=ctdict[key];
			ltdict.push([key,ctdict[key]])
		}
		
		var cellCounter=0;
		
		var arc=d3.svg.arc()
				.innerRadius(100)
				.outerRadius(150)
				.startAngle(function(d,i){
					var sa=2*Math.PI*(cellCounter/asum);
					return sa;
					})
				.endAngle(function(d,i){
					cellCounter+=d[1];
					var ea=2*Math.PI*(cellCounter/asum);
					return ea;
					});
			
		tooltipsvg
			.style("left",(d3.event.pageX+30)+"px")
			.style("top",(d3.event.pageY-28)+"px")
			.style("position","absolute")
			.style("opacity",1)
			.style("background","gray")
			//.style("border","1px solid")
			.attr("class","tooltipsvg")
			.attr("width",300)
			.attr("height",300)
		
		var colorCodes=["red","blue","purple","green","lime","magenta","orange","olive","cyan","yellow"];
		var colorDex=0;
		tooltipsvg.selectAll("path")
			.data(ltdict)
			.enter()
			//.append("circle")
			//.attr("r",10);
			.append("path")
			.attr("d",arc)
			.attr("transform","translate (150,150)")
			.attr("fill",function(d,i) { 
				colorDex=colorDex % 10;
				var colorChoosen=colorCodes[colorDex];
				colorDex+=1;
				return colorChoosen;}
				)
				
		cellCounter=0;
		colorDex=0;
		tooltipsvg.selectAll("text")
			.data(ltdict)
			.enter()
			.append("text")
			.style("border","1px solid")
			.attr("transform",function(d,i){
				var tx=60;
				var ty=100+i*20;
				return "translate ("+tx+","+ty+")";
				})
			.attr("dy",20)
			.text(function(d,i){return d[0]+':'+d[1];})
			.style("font-size","10px")
			.style("fill",function(d,i) { 
				colorDex=colorDex % 10;
				var colorChoosen=colorCodes[colorDex];
				colorDex+=1;
				return colorChoosen;}
				);
	}
		
	// private wrap multi- line text 
	function addwraptext(text){
		text.each(function(d){
			CELL=d.CELL;
			var ctdict={};
			for(x of CELL){
				xlabel=cells[x].typeLabel;
				if (xlabel in ctdict==false){
					ctdict[xlabel]=1;
				}else{
					ctdict[xlabel]+=1;
				}
			}
			d3.select(this).append("tspan")
				.attr("x","1em")
				.attr("dy","1em")
				.text(function(d){
					return "ID: E"+d.T+"_"+d.ID;
					})
			// add diff stage
			d3.select(this).append("tspan")
				.attr("x","1em")
				.attr("dy","1em")
				.text(function(d){
					return "Differentiation: "+(d.D*100).toFixed(1)+"%";
				})
				
			for ( var key in ctdict){
				d3.select(this).append("tspan")
				.attr("x",function(d){ return 20;})
				.attr("dy","1em")
				.text(key+","+ctdict[key]);
			}
		});
	}
	
	//this function is used to pre-compute all needed data for nodes (e.g. target fold change).
	function updateNodes(nodes){
		for (var node of nodes){
			nodeUpdate(node);
		}
		
		
	}
	//calculate fold change for all gene for given node (parent->node)
	//calculate edge ending at chosen node 

	function nodeUpdate(node){
		var pnode=node.parent;
		var chosenEdge;
		var fc={};
		for (var edge of edges){
			if (nodes[edge.to]==node){
				chosenEdge=JSON.parse(JSON.stringify(edge));
			}
		}
		if (chosenEdge!=null){
			for (var i in GL){
				fcg=node.E[i]-pnode.E[i];
				gi=GL[i];
				fc[gi]=fcg;
			}
		}
		node["fc"]=fc;
		node["endEdge"]=chosenEdge;
	}

	//handle DIV_CONFIG section
	//reset
	
	//reset path
	function resetPath(){
		d3.select("svg").selectAll("path")
		.attr("stroke","#fff");
		
		d3.selectAll(".linktext")
		.text("");
	}

	//reset color bar
	function resetcolorbar(){
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
		
		d3.select(".colorbar").selectAll("text")
		.data([-4,4])
		.text(function(d){return d;});
	}

	//----------------------------------------------------------------------
	//hide TF
	function hideTF(){
		//hide TF text
		d3.select("svg").selectAll(".TFText")
		.style("opacity",0);
		
		d3.select("svg").selectAll(".TFText")
		.selectAll("tspan")
		//.style("font-size",0)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		//hide TF color bar
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
	}

	function hideDE(){
		//HIDE DE text
		d3.select("svg").selectAll(".DEText")
		.style("opacity",0);
		
		d3.select("svg").selectAll(".DEText")
		.selectAll("tspan")
		//.style("font-size",0)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		//HIDE color bar
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
	}

	//hide DE

	//---------------------------------------------------------------------
	//show controls
	//show RTF
	function showRTF(){
		//change the opacity
		d3.select("svg").selectAll(".RTFText")
		.style("opacity",1);
		
		//change the text size
		d3.select("svg").selectAll(".RTFText")
		.selectAll("tspan");
		
		resetcolorbar();
		d3.select("svg").select(".colorbar")
		.style("opacity",1);
		
	}
	
	function hideRTF(){
		//change the opacity
		d3.select("svg").selectAll(".RTFText")
		.style("opacity",0);
		
		//change the text size
		d3.select("svg").selectAll(".RTFText")
		.selectAll("tspan");
		
		resetcolorbar();
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
		
	}
	
	
	//show TF 
	function showTF(){
		
		//change the opacity
		d3.select("svg").selectAll(".TFText")
		.style("opacity",1);
		//change the text size
		d3.select("svg").selectAll(".TFText")
		.selectAll("tspan");
		
		
		resetcolorbar();
		d3.select("svg").select(".colorbar")
		.style("opacity",1);
		
	}

	//show DE
	function showDE(){
		d3.select("svg").selectAll(".DEText")
		.style("opacity",1);
		
		//change the text size
		d3.select("svg").selectAll(".DEText")
		.selectAll("tspan");
		
		
		//color bar
		resetcolorbar();
		d3.select("svg").select(".colorbar")
		.style("opacity",1);
		
	}

	//----------------------------------------------------------------------
	//create svg text for showing DE gens
	function createDE(){
		var gc_enter=svg.selectAll("g")
		.append("text")
		.attr("class","DEText")
		.attr("opacity",0)
		.attr("x","-8em")
		.attr("y","0em")
		gc_enter.call(addDEText,20);
	}

	//add DE texts 
	function addDEText(text,showcutoff){
		text.each(function(d){
			var chosenEdge=null;
			for (var edge of edges){
				if (nodes[edge.to]==d){
					chosenEdge=JSON.parse(JSON.stringify(edge));
				}
			}
			if (chosenEdge!=null){
				var de=chosenEdge.de;
				de=de.map(function(x){return x.toUpperCase();});
				var de=de.splice(0,20).reverse();
				for (var gd of de){
					var pnodeE=nodes[chosenEdge.from].E;
					var cnodeE=d.E;
					var tfindex=GL.indexOf(gd);
					var fc=cnodeE[tfindex]-pnodeE[tfindex];
					fc=fc.toFixed(2);
					nfc=Math.min(Math.max(fc,-4),4);
					var per=(4-nfc)/(4+4);
					var cl=hsl_col_perc(per,0,120);
					d3.select(this).append("tspan")
					.style("font-size",pathtextsize)
					.attr("x","-8em")
					.attr("dy","-1em")
					.style("fill",cl)
					.text(gd+' ('+fc+')');
				}
			}
		});
	}

	//create svg text for showing TFs
	function createEdgeTF(){
		var gc_enter=svg.selectAll("g")
		.append("text")
		.attr("class","TFText")
		.attr("opacity",0)
		.attr("x","-8em")
		.attr("y","0em")
		gc_enter.call(addTFText,20);
		
	}
	
	// creat svg text for showing repressing TFs
	function createRTFs(){
		var gc_enter=svg.selectAll("g")
		.append("text")
		.attr("class","RTFText")
		.attr("opacity",0)
		.attr("x","0em")
		.attr("y","0em")
		gc_enter.call(addRTFText,20);
	}
	
	//add RTF texts
	function addRTFText(text,showcutoff){
		text.each(function(d){
			var dRTF=d.eTF;
			if (dRTF.length>0){
				var etf=dRTF
				var etf=etf.slice(0,showcutoff).reverse();
				for (var tf of etf){
					var color_tf;
					var fc=tf[2];
					var per=(4-fc)/(4+4);
					var cl=hsl_col_perc(per,0,120);
					d3.select(this).append("tspan")
					.style("font-size",pathtextsize)
					.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
					.attr("x","2em")
					.attr("dy","-1em")
					.attr("text-anchor","start")
					.style("fill",cl)
					.text(tf[1]+' ('+tf[0].toExponential(2)+')');
				}
			}
		});
	}

	// add TF texts 
	function addTFText(text,showcutoff){
		text.each(function(d){
			var chosenEdge=null;
			for (var edge of edges){
				if (nodes[edge.to]==d){
					chosenEdge=JSON.parse(JSON.stringify(edge));
				}
			}
			if (chosenEdge!=null){
				var etf=chosenEdge.etf;
				var etf=etf.splice(0,20).reverse();
				for (var tf of etf){
					var color_tf;
					var pnodeE=nodes[chosenEdge.from].E;
					var cnodeE=d.E;
					var tfindex=GL.indexOf(tf[1]);
					var fc=cnodeE[tfindex]-pnodeE[tfindex];
					fc=Math.min(Math.max(fc,-4),4);
					var per=(4-fc)/(4+4);
					var cl=hsl_col_perc(per,0,120);
					d3.select(this).append("tspan")
					.style("font-size",pathtextsize)
					.attr("x","-8em")
					.attr("dy","-1em")
					.style("fill",cl)
					.text(tf[1]+' ('+tf[0].toExponential(2)+')');
				}
			}
		});
	}

	// create color bar 
	function createColorBar(){
		var cx=[];
		for (var ci=1;ci>0;ci-=0.005){
			cx.push(ci);
		}
		var xstart=10;
		var ystart=15;
		var xspan=1;
		var xspan_height=10;
		var xspan_width=4;
		d3.select("svg").append("g")
		.attr("opacity",0)
		.attr("class","colorbar")
		.selectAll("rect")
		.data(cx)
		.enter()
		.append("rect")
		.attr("x",function(d,i){
			var xx=xstart+i*xspan;
			return xx;
			})
		.attr("y",ystart)
		.attr("width",xspan_width)
		.attr("height",xspan_height)
		.attr("fill",function(d){
			var cl=hsl_col_perc(d,0,120);
			return cl;
			});
		xmm=[0,1]
		d3.select("svg").select(".colorbar").append("g")
		.selectAll("text")
		.data(xmm)
		.enter()
		.append("text")
		.attr("x",function(d,i){
				xx=xstart+d*xspan*cx.length;
				return xx;
			})
		.attr("y",ystart-2)
		.text(function(d){
			if(d==0){
				return -4;
			}else{
				return 4;
				}})
		.attr("text-anchor",function(d,i){
				if(d==0){
					return "start";
				}else{
					return "end";
				}
			})
		.attr("fill","white");	
	}

	//color convertion
	function hsl_col_perc(percent,start,end) {

		 var a=percent;
		 b = end*a;
		 c = b+start;

		//Return a CSS HSL string
		return 'hsl('+c+',100%,50%)';
	}

	// wrap multi- line text 
	function addwraptext(text){
		text.each(function(d){
			CELL=d.CELL;
			var ctdict={};
			for(x of CELL){
				xlabel=cells[x].typeLabel;
				if (xlabel in ctdict==false){
					ctdict[xlabel]=1;
				}else{
					ctdict[xlabel]+=1;
				}
			}
			d3.select(this).append("tspan")
				.attr("x",20)
				.attr("dy","1em")
				.text(function(d){
					return "ID: E"+d.T+"_"+d.ID;
					})
			// add diff stage
			d3.select(this).append("tspan")
				.attr("x",20)
				.attr("dy","1em")
				.text(function(d){
					return "Differentiation: "+(d.D*100).toFixed(1)+"%";
				})
				
			for ( var key in ctdict){
				d3.select(this).append("tspan")
				.attr("x",function(d){ return 20;})
				.attr("dy","1em")
				.text(key+","+ctdict[key]);
			}
		});
	}

	//function export tf details of a edge
	function exportTFEdge(resList,chosenEdge){
		var tf,tfID,pv,rank,fc,tfc;
		var etf=chosenEdge.etf;
		var de=chosenEdge.de;
		de=de.map(function(d){return d.toUpperCase();});
		cnode=nodes[chosenEdge.to];
		pnode=nodes[chosenEdge.from];
		var edgeid='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;
		for (i in etf){
			tf=etf[i];
			var tfID=tf[1].toUpperCase();
			var tfp=tf[0].toExponential(2);
			var tfrank=i;
			var fc=cnode.fc[tfID].toFixed(3);
			var target=dTD[tfID].filter(function(d){return d in cnode.fc;});
			var tfc=target.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
			tfc=(tfc.reduce(function(a,b){return a+b;})/tfc.length).toFixed(3);
			var dtfc=target.filter(function(d){
					if (de.indexOf(d)!=-1){
						return true;
					}else{
						return false;
					}
				});
			dtfc=dtfc.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
			dtfc=(dtfc.reduce(function(a,b){return a+b;})/dtfc.length).toFixed(3);
			var tfresList=[tfID,tfp,fc,tfc,dtfc,edgeid];
			resList.push(tfresList);
		}
		return resList;
	}

	function exportDEEdge(resDEList,chosenEdge){
		var cnode=nodes[chosenEdge.to];
		var pnode=nodes[chosenEdge.from];
		
		var edgeID='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;


		var etf=chosenEdge.etf;
		var de=chosenEdge.de.map(function(d){return d.toUpperCase();});
		
		for (i in de){
			g=de[i];
			ef=pnode.E[GL.indexOf(g)];
			et=cnode.E[GL.indexOf(g)];
			fc=et-ef;
			resDEList.push([g,ef.toFixed(3),et.toFixed(3),fc.toFixed(3),edgeID]);
		}
		return resDEList;
	}

	//append a table to the cant
	function createTable(cant,tableid,data){
		cant.append("table")
		.attr("id",tableid)
		.style("border","1px solid")
		.style("border-collapse","collapse")
		.selectAll("tr")
		.data(data)
		.enter()
		.append("tr")
		.style("border","1px solid")
		.selectAll("td")
		.data(function(d){return d;})
		.enter()
		.append("td")
		.style("border","1px solid")
		.text(function(d){return d;})
		.style("font-weight",function(d){
			if (d in dTD ){
				return "bold";
			}
		});
	}
	
	//create drop down
	function createDropDown(FirstRow,tfdropdowndiv,dropdownid,TFs,onChange){
		var Keys=[FirstRow];
		TFs.sort();
		for (var tf of TFs){
			Keys.push(tf);
		}
		
		d3.select(tfdropdowndiv)
		.select("select")
		.remove();
		
		d3.select(tfdropdowndiv)
		.append("select")
		.attr("id",dropdownid)
		.on("change",onChange)
		.selectAll("option")
		.data(Keys)
		.enter()
		.append("option")
		.text(function(d){
			return d;
		})
		.attr("value",function(d){
			return d;
			});
	
	}
	
	//action on mouse out
	function outNode(){
		tooltipsvg
		.style("opacity", 0)
		.attr("width",0)
		.attr("height",0)
		.selectAll("*").remove();;
	}
	
	// action on mouse over
	function inNode(){
		if (d3.select("#tooltipcheck").property("checked")==false){
			return false;
		}
		var CELL=this.__data__.CELL;
		var ctdict={}
		for(x of CELL){
			xlabel=cells[x].typeLabel;
			if (xlabel in ctdict==false){
				ctdict[xlabel]=1;
			}else{
				ctdict[xlabel]+=1;
			}
		}

		var ltdict=[]
		var asum=0;
		for (key in ctdict){
			asum+=ctdict[key];
			ltdict.push([key,ctdict[key]])
		}
		
		var cellCounter=0;
		
		var arc=d3.svg.arc()
				.innerRadius(100)
				.outerRadius(150)
				.startAngle(function(d,i){
					var sa=2*Math.PI*(cellCounter/asum);
					return sa;
					})
				.endAngle(function(d,i){
					cellCounter+=d[1];
					var ea=2*Math.PI*(cellCounter/asum);
					return ea;
					});
			
		tooltipsvg
			.style("left",(d3.event.pageX+30)+"px")
			.style("top",(d3.event.pageY-28)+"px")
			.style("position","absolute")
			.style("opacity",1)
			.style("background","gray")
			//.style("border","1px solid")
			.attr("class","tooltipsvg")
			.attr("width",300)
			.attr("height",300)
		
		var colorCodes=["red","blue","purple","green","lime","magenta","orange","olive","cyan","yellow"];
		var colorDex=0;
		tooltipsvg.selectAll("path")
			.data(ltdict)
			.enter()
			//.append("circle")
			//.attr("r",10);
			.append("path")
			.attr("d",arc)
			.attr("transform","translate (150,150)")
			.attr("fill",function(d,i) { 
				colorDex=colorDex % 10;
				var colorChoosen=colorCodes[colorDex];
				colorDex+=1;
				return colorChoosen;}
				)
				
		cellCounter=0;
		colorDex=0;
		tooltipsvg.selectAll("text")
			.data(ltdict)
			.enter()
			.append("text")
			.style("border","1px solid")
			.attr("transform",function(d,i){
				var tx=60;
				var ty=100+i*20;
				return "translate ("+tx+","+ty+")";
				})
			.attr("dy",20)
			.text(function(d,i){return d[0]+':'+d[1];})
			.style("font-size","10px")
			.style("fill",function(d,i) { 
				colorDex=colorDex % 10;
				var colorChoosen=colorCodes[colorDex];
				colorDex+=1;
				return colorChoosen;}
				);
	}
	
	
	// parse json file
	function parseJSON(data){
		var data=JSON.parse(data);
		var GL=data[0];
		var CellList=data[1];
		var NodeList=data[2];
		var EdgeList=data[3];
		var dTD=data[4];	
		return [GL,CellList,NodeList,EdgeList,dTD];
	};

	//-------------------------------------------------------------------------
	//download functions here
	
	function createjsondownload(jsondownloadlinkid){
		 var dataObj=JSON.parse(data);
		 var outObj={"GeneList":dataObj[0], "CellList":dataObj[1], "NodeList":dataObj[2],"EdgeList":dataObj[3]};
		 var outObjStr=JSON.stringify(outObj,null,"	");
		 var blob=new Blob([outObjStr],{type:'application/json;charset=utf-8'});
		 outurl=window.URL.createObjectURL(blob);
		 d3.select("#"+jsondownloadlinkid)
		.attr("href",outurl)
		.attr("download","download.json")
		.text("Ready,Click to download");
	}
	
	

	function creatededownload(dedownloadlinkid){
		var resList=[];
		resList.push(["DE gene","Expression_from","Expression_to","log2 Fold change","edgeID"]);

		for (var edge of edges){
			var etf=edge.etf;
			var de=edge.de;
			de=de.map(function(d){return d.toUpperCase();});
			var cnode=nodes[edge.to];
			var pnode=nodes[edge.from];
			var edgeid='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;
			for (var g of de){
				var gindex=GL.indexOf(g);
				var efrom=pnode.E[gindex];
				var eto=cnode.E[gindex];
				var fc=eto-efrom;
				var deresList=[g,efrom,eto,fc,edgeid];
				resList.push(deresList);
			}
		}
		//resList=resList.splice(0,10000);
		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#"+dedownloadlinkid)
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}


	
	function createtfdownload(tfdownloadlinkid){
		var resList=[];
		//target tfs
		resList.push(["TF","p-value","TF log2 fold change","Mean target log2 fold change","Mean DE target log2 fold change","edgeID"]);
		for (var edge of edges){
			var etf=edge.etf;
			var de=edge.de;
			de=de.map(function(d){return d.toUpperCase();});
			var cnode=nodes[edge.to];
			var pnode=nodes[edge.from];
			var edgeid='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;
			for (var i in etf){
				tf=etf[i];
				var tfID=tf[1].toUpperCase();
				var tfp=tf[0].toExponential(2);
				var tfrank=i;
				var fc=cnode.fc[tfID].toFixed(3);
				var target=dTD[tfID].filter(function(d){return d in cnode.fc;});
				var tfc=target.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
				tfc=(tfc.reduce(function(a,b){return a+b;})/tfc.length).toFixed(3);
				var dtfc=target.filter(function(d){
						if (de.indexOf(d)!=-1){
							return true;
						}else{
							return false;
						}
					});
				dtfc=dtfc.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
				dtfc=(dtfc.reduce(function(a,b){return a+b;})/dtfc.length).toFixed(3);
				var tfresList=[tfID,tfp,fc,tfc,dtfc,edgeid];
				resList.push(tfresList);
			}
		}
		
		resList.push([]) //two empty lines to separate eTFs from TFs
 		resList.push([])
		//expression TFs
		resList.push(['eTF','p-value','log2 fold change to parent','log2 fold change to siblings','edgeID'])
		for (var node of nodes){
			var etfs=node.eTF;
			var pnode=node.parent
			if (pnode!="null" && etf!=undefined){
				for (var etf of etfs){
					var tfresList=[etf[1],etf[0].toExponential(2),etf[2],etf[3],'E'+pnode.T+"_"+pnode.ID+'->E'+node.T+"_"+node.ID]
					resList.push(tfresList)
				}
			}
			
		}
		
		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#"+tfdownloadlinkid)
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}

	function plswait(id){
		document.getElementById(id).innerHTML="wait...";
		d3.select("#"+id)
		.text("Generating file,please wait...");	
	}

	//convert 2D-list to  TSV
	function List2TSV(LST){
		var out="";
		for (var i in LST){
			var LSTI=LST[i].join("	");
			out+=LSTI+'\\n';
		}
		return out;
	}
	//write download link

	//svgToCanvas
	function svgToCanvas(downloadlinkid){
		//get svg element.
		var svg = document.getElementById("svg");

		//get svg source.
		var serializer = new XMLSerializer();
		var source = serializer.serializeToString(svg);

		//add name spaces.
		if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)){
			source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
		}
		if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)){
			source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');
		}

		//add xml declaration
		source = '<?xml version="1.0" standalone="no"?>\\n' + source;

		//convert svg source to URI data scheme.
		var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(source);

		//set url value to a element's href attribute.
		document.getElementById(downloadlinkid).href = url;
		document.getElementById(downloadlinkid).innerHTML="Figure ready,Right click me to save!"
		//you can download svg file by right click menu	
	}

	//convert table to excel download link
	//dlink: link id
	function table2XLS(newW,tlink,dlink){
		var datatype="data:application/vnd.ms-excel";
		var table_div=newW.document.getElementById(tlink).outerHTML;
		var url='data:application/vnd.ms-excel,'+encodeURIComponent(table_div);
		newW.document.getElementById(dlink).href=url;
		newW.document.getElementById(dlink).download="download.xls";
		newW.document.getElementById(dlink).innerHTML="Table ready!click to save! (please add the right file extension .xls if not prompted)"
	}
	
	// panther go analysis
	
	function panthergoInput(keys,species){
		var species=`
		<body>
		<span> choose your species: </span> <br>
		<select
		   id="rte_species"
		   class="form-control"
		   name="species">
		  <option value="HUMAN">Homo sapiens</option>
		  <option value="MOUSE">Mus musculus</option>
		  <option value="RAT">Rattus norvegicus</option>
		  <option value="CHICK">Gallus gallus</option>
		  <option value="DANRE">Danio rerio</option>
		  <option value="DROME">Drosophila melanogaster</option>
		  <option value="CAEEL">Caenorhabditis elegans</option>
		  <option value="YEAST">Saccharomyces cerevisiae</option>
		  <option value="SCHPO">Schizosaccharomyces pombe</option>
		  <option value="DICDI">Dictyostelium discoideum</option>
		  <option value="ARATH">Arabidopsis thaliana</option>
		  <option value="ECOLI">Escherichia coli</option>
		  <option value="EMENI">Emericella nidulans</option>
		  <option value="ANOCA">Anolis carolinensis</option>
		  <option value="ANOGA">Anopheles gambiae</option>
		  <option value="AQUAE">Aquifex aeolicus</option>
		  <option value="ASHGO">Ashbya gossypii</option>
		  <option value="BACCR">Bacillus cereus</option>
		  <option value="BACSU">Bacillus subtilis</option>
		  <option value="BACTN">Bacteroides thetaiotaomicron</option>
		  <option value="BATDJ">Batrachochytrium dendrobatidis</option>
		  <option value="BOVIN">Bos taurus</option>
		  <option value="BRADI">Brachypodium distachyon</option>
		  <option value="BRAJA">Bradyrhizobium japonicum</option>
		  <option value="BRAFL">Branchiostoma floridae</option>
		  <option value="CAEBR">Caenorhabditis briggsae</option>
		  <option value="CANAL">Candida albicans</option>
		  <option value="CANFA">Canis familiaris</option>
		  <option value="CHLTR">Chlamydia trachomatis</option>
		  <option value="CHLRE">Chlamydomonas reinhardtii</option>
		  <option value="CHLAA">Chloroflexus aurantiacus</option>
		  <option value="CIOIN">Ciona intestinalis</option>
		  <option value="CLOBH">Clostridium botulinum</option>
		  <option value="COXBU">Coxiella burnetii</option>
		  <option value="CRYNJ">Cryptococcus neoformans</option>
		  <option value="DAPPU">Daphnia pulex</option>
		  <option value="DEIRA">Deinococcus radiodurans</option>
		  <option value="DICTD">Dictyoglomus turgidum</option>
		  <option value="DICPU">Dictyostelium purpureum</option>
		  <option value="ENTHI">Entamoeba histolytica</option>
		  <option value="HORSE">Equus caballus</option>
		  <option value="FELCA">Felis catus</option>
		  <option value="FUSNN">Fusobacterium nucleatum</option>
		  <option value="GEOSL">Geobacter sulfurreducens</option>
		  <option value="GIAIC">Giardia intestinalis</option>
		  <option value="GLOVI">Gloeobacter violaceus</option>
		  <option value="SOYBN">Glycine max</option>
		  <option value="HAEIN">Haemophilus influenzae</option>
		  <option value="HALSA">Halobacterium salinarum</option>
		  <option value="IXOSC">Ixodes scapularis</option>
		  <option value="KORCO">Korarchaeum cryptofilum</option>
		  <option value="LEIMA">Leishmania major</option>
		  <option value="LEPIN">Leptospira interrogans</option>
		  <option value="LISMO">Listeria monocytogenes</option>
		  <option value="MACMU">Macaca mulatta</option>
		  <option value="METJA">Methanocaldococcus jannaschii</option>
		  <option value="METAC">Methanosarcina acetivorans</option>
		  <option value="MONDO">Monodelphis domestica</option>
		  <option value="MONBE">Monosiga brevicollis</option>
		  <option value="MYCTU">Mycobacterium tuberculosis</option>
		  <option value="NEMVE">Nematostella vectensis</option>
		  <option value="NEUCR">Neurospora crassa</option>
		  <option value="ORNAN">Ornithorhynchus anatinus</option>
		  <option value="ORYSJ">Oryza sativa</option>
		  <option value="PANTR">Pan troglodytes</option>
		  <option value="PHANO">Phaeosphaeria nodorum</option>
		  <option value="PHYPA">Physcomitrella patens</option>
		  <option value="PHYIT">Phytophthora infestans</option>
		  <option value="PLAF7">Plasmodium falciparum</option>
		  <option value="POPTR">Populus trichocarpa</option>
		  <option value="PRIPA">Pristionchus pacificus</option>
		  <option value="PSEAE">Pseudomonas aeruginosa</option>
		  <option value="PUCGT">Puccinia graminis</option>
		  <option value="PYRAE">Pyrobaculum aerophilum</option>
		  <option value="PYRKO">Pyrococcus kodakaraensis</option>
		  <option value="RHOBA">Rhodopirellula baltica</option>
		  <option value="SALTY">Salmonella typhimurium</option>
		  <option value="SCHMA">Schistosoma mansoni</option>
		  <option value="SCLS1">Sclerotinia sclerotiorum</option>
		  <option value="SHEON">Shewanella oneidensis</option>
		  <option value="SOLLC">Solanum lycopersicum</option>
		  <option value="SORBI">Sorghum bicolor</option>
		  <option value="STAA8">Staphylococcus aureus</option>
		  <option value="STRR6">Streptococcus pneumoniae</option>
		  <option value="STRCO">Streptomyces coelicolor</option>
		  <option value="STRPU">Strongylocentrotus purpuratus</option>
		  <option value="SULSO">Sulfolobus solfataricus</option>
		  <option value="PIG">Sus scrofa</option>
		  <option value="SYNY3">Synechocystis</option>
		  <option value="TAKRU">Takifugu rubripes</option>
		  <option value="TETTS">Tetrahymena thermophila</option>
		  <option value="THAPS">Thalassiosira pseudonana</option>
		  <option value="THEYD">Thermodesulfovibrio yellowstonii</option>
		  <option value="THEMA">Thermotoga maritima</option>
		  <option value="TRIVA">Trichomonas vaginalis</option>
		  <option value="TRIAD">Trichoplax adhaerens</option>
		  <option value="brucei">TRYB2 Trypanosoma brucei</option>
		  <option value="USTMA">Ustilago maydis</option>
		  <option value="VIBCH">Vibrio cholerae</option>
		  <option value="VITVI">Vitis vinifera</option>
		  <option value="XANCP">Xanthomonas campestris</option>
		  <option value="XENTR">Xenopus tropicalis</option>
		  <option value="YARLI">Yarrowia lipolytica</option>
		  <option value="YERPE">Yersinia pestis</option>
		</select>
		<input type="submit" id="gospecies">
		</body>
		`
		var ww=open("",'_blank','height=300,width=400,scrollbars=yes');
		ww.document.write(species);
		d3.select(ww.document.body).select("#gospecies")
		.on("click",function(){
			var ss=ww.document.getElementById("rte_species").value;
			panthergoInputSubmit(keys,ss);
		});		
	}
	
	function panthergoInputSubmit(keys,species){
		var link_pre="http://pantherdb.org/webservices/go/overrep.jsp?input=";
		if (species==undefined){
			var link_suffix="&species=MOUSE";
		}else{
			var link_suffix="&species="+species;
		}
		
		keys=keys.join("\\n");
		key=encodeURIComponent(keys);
		var link=link_pre+keys+link_suffix;
		link=encodeURI(link);
		var ww=open(link,'_blank','height=600,width=800,left=1200,top=200,scrollbars=yes');
	}

	
	// toppgene go analysis 
	function toppgenegoInput(keys){
		var link_pre="https://toppgene.cchmc.org/CheckInput.action?query=TOPPFUN&type=HGNC_SYNONYMS&training_set="
		var link_suffix="";
		keys=keys.join("+");
		key=encodeURIComponent(keys);
		var link=link_pre+keys+link_suffix;
		link=encodeURI(link);
		open(link,'_blank','height=600,width=800,left=1200,top=200,scrollbars=yes');
	}
	
	// colors
	function getColors(){
			var colorObject={
		  "red": "#ff0000",
		  "orange": "#ffa500",
		  "yellow": "#ffff00",
		  "green": "#008000",
		  "blue": "#0000ff",
		  "cyan": "#00ffff",
		  "purple": "#800080",
		  "pink": "#ffc0cb",
		  "black": "#000000",
		  "brown": "#a52a2a", 
		  "burlywood": "#deb887",
		  "cadetblue": "#5f9ea0",
		  "chartreuse": "#7fff00",
		  "chocolate": "#d2691e",
		  "coral": "#ff7f50",
		  "cornflowerblue": "#6495ed",
		  "cornsilk": "#fff8dc",
		  "crimson": "#dc143c",
		  "darkblue": "#00008b",
		  "darkcyan": "#008b8b",
		  "darkgoldenrod": "#b8860b",
		  "darkgray": "#a9a9a9",
		  "darkgreen": "#006400",
		  "darkgrey": "#a9a9a9",
		  "darkkhaki": "#bdb76b",
		  "darkmagenta": "#8b008b",
		  "darkolivegreen": "#556b2f",
		  "darkorange": "#ff8c00",
		  "darkorchid": "#9932cc",
		  "darkred": "#8b0000",
		  "darksalmon": "#e9967a",
		  "darkseagreen": "#8fbc8f",
		  "darkslateblue": "#483d8b",
		  "darkslategray": "#2f4f4f",
		  "darkslategrey": "#2f4f4f",
		  "darkturquoise": "#00ced1",
		  "darkviolet": "#9400d3",
		  "deeppink": "#ff1493",
		  "deepskyblue": "#00bfff",
		  "dimgray": "#696969",
		  "dimgrey": "#696969",
		  "dodgerblue": "#1e90ff",
		  "firebrick": "#b22222",
		  "floralwhite": "#fffaf0",
		  "forestgreen": "#228b22",
		  "fuchsia": "#ff00ff",
		  "gainsboro": "#dcdcdc",
		  "ghostwhite": "#f8f8ff",
		  "gold": "#ffd700",
		  "goldenrod": "#daa520",
		  "gray": "#808080",
		  "greenyellow": "#adff2f",
		  "grey": "#808080",
		  "honeydew": "#f0fff0",
		  "hotpink": "#ff69b4",
		  "indianred": "#cd5c5c",
		  "indigo": "#4b0082",
		  "ivory": "#fffff0",
		  "khaki": "#f0e68c",
		  "lavender": "#e6e6fa",
		  "lavenderblush": "#fff0f5",
		  "lawngreen": "#7cfc00",
		  "lemonchiffon": "#fffacd",
		  "lightblue": "#add8e6",
		  "lightcoral": "#f08080",
		  "lightcyan": "#e0ffff",
		  "lightgoldenrodyellow": "#fafad2",
		  "lightgray": "#d3d3d3",
		  "lightgreen": "#90ee90",
		  "lightgrey": "#d3d3d3",
		  "lightpink": "#ffb6c1",
		  "lightsalmon": "#ffa07a",
		  "lightseagreen": "#20b2aa",
		  "lightskyblue": "#87cefa",
		  "lightslategray": "#778899",
		  "lightslategrey": "#778899",
		  "lightsteelblue": "#b0c4de",
		  "lightyellow": "#ffffe0",
		  "lime": "#00ff00",
		  "limegreen": "#32cd32",
		  "linen": "#faf0e6",
		  "magenta": "#ff00ff",
		  "maroon": "#800000",
		  "mediumaquamarine": "#66cdaa",
		  "mediumblue": "#0000cd",
		  "mediumorchid": "#ba55d3",
		  "mediumpurple": "#9370db",
		  "mediumseagreen": "#3cb371",
		  "mediumslateblue": "#7b68ee",
		  "mediumspringgreen": "#00fa9a",
		  "mediumturquoise": "#48d1cc",
		  "mediumvioletred": "#c71585",
		  "midnightblue": "#191970",
		  "mintcream": "#f5fffa",
		  "mistyrose": "#ffe4e1",
		  "moccasin": "#ffe4b5",
		  "navajowhite": "#ffdead",
		  "navy": "#000080",
		  "oldlace": "#fdf5e6",
		  "olive": "#808000",
		  "olivedrab": "#6b8e23",
		  "orangered": "#ff4500",
		  "orchid": "#da70d6",
		  "palegoldenrod": "#eee8aa",
		  "palegreen": "#98fb98",
		  "paleturquoise": "#afeeee",
		  "palevioletred": "#db7093",
		  "papayawhip": "#ffefd5",
		  "peachpuff": "#ffdab9",
		  "peru": "#cd853f",
		  "plum": "#dda0dd",
		  "powderblue": "#b0e0e6",
		  "rebeccapurple": "#663399",
		  "rosybrown": "#bc8f8f",
		  "royalblue": "#4169e1",
		  "saddlebrown": "#8b4513",
		  "salmon": "#fa8072",
		  "sandybrown": "#f4a460",
		  "seagreen": "#2e8b57",
		  "seashell": "#fff5ee",
		  "sienna": "#a0522d",
		  "silver": "#c0c0c0",
		  "skyblue": "#87ceeb",
		  "slateblue": "#6a5acd",
		  "slategray": "#708090",
		  "slategrey": "#708090",
		  "snow": "#fffafa",
		  "springgreen": "#00ff7f",
		  "steelblue": "#4682b4",
		  "tan": "#d2b48c",
		  "teal": "#008080",
		  "thistle": "#d8bfd8",
		  "tomato": "#ff6347",
		  "turquoise": "#40e0d0",
		  "violet": "#ee82ee",
		  "wheat": "#f5deb3",
		  "white": "#ffffff",
		  "whitesmoke": "#f5f5f5",
		  "aliceblue": "#f0f8ff",
		  "antiquewhite": "#faebd7",
		  "aqua": "#00ffff",
		  "aquamarine": "#7fffd4",
		  "azure": "#f0ffff",
		  "beige": "#f5f5dc",
		  "bisque": "#ffe4c4",
		  "blanchedalmond": "#ffebcd",
		  "blueviolet": "#8a2be2",
		  "yellowgreen": "#9acd32"
		};
		return colorObject;
	}

	
}(window.scviz =window.scviz ||{}, jQuery));
	
	"""
	#-----------------------------------------------------------------------
	
	for i in G1.Nodes:
		i.E=getAvgEx(i)
			
	GJ=GtoJson(G1,GL,dTD)	
	f=open(output+'/'+scg_name+'.json','w')
	f.write(GJ)
	f.close()
	
	HTML_template=HTML_template%(scg_name+'.json')
	f=open(output+'/'+scg_name+'.html','w')
	f.write(HTML_template)
	f.close()

	f=open(output+'/style.css','w')
	f.write(css_template)
	f.close()


	f=open(output+'/parseJSON.js','w')
	f.write(javascript)
	f.close()
	
