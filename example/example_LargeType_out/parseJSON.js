
	function onload(){
		RL=parseJSON(data);
		nodes=buildTree(RL[2]);
		root=nodes[0];
		cells=RL[1];
		edges=RL[3];
		GL=RL[0];
		dTD=RL[4];
		drawTree(root);
		updateNodes(nodes);
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
	/*
	function nodeUpdate(node){
		var pnode=node.parent;
		var chosenEdge;
		var tf,tfID,pv,rank,fc,tfc,tfcde;
		for (var edge of edges){
			if (nodes[edge.to]==node){
				chosenEdge=edge;
			}
		}
		if (chosenEdge!=null){
			var etf=chosenEdge.etf;
			var de=chosenEdge.de;
			var de=de.map(function(d){
					return d.toUpperCase();
				})
			for (var i in etf){
				tf=etf[i];
				tfID=tf[1];
				pv=tf[0];
				rank=i;
				tfindex=GL.indexOf(tfID.toUpperCase());
				fc=node.E[tfindex]-pnode.E[tfindex];
				var tftarget=dTD[tfID];
				var target=tftarget.filter(function(d){
					if(GL.indexOf(d)!=-1){
						return true;
					}else{
						return false;
					}});
				//fold change for all targets
				tfc=target.map(function(d){
						var di=GL.indexOf(d);
						fci=node.E[di]-pnode.E[di];
						return fci;
					});
				// fold change for only de targets
				var detarget=target.filter(function(d){
						if (de.indexOf(d)!=-1){
							return true;
						}else{
							return false;
						}
					});
				tfcde=detarget.map(function(d){
						var di=GL.indexOf(d);
						fci=node.E[di]-pnode.E[di];
						return fci;
				});
				
				tfc=tfc.reduce(function(a,b){
						return a+b;
					})/tfc.length;
				tfcde=tfcde.reduce(function(a,b){
						return a+b;
					})/tfcde.length;

				//
			}
		}
		
	}

	*/
	//----------------------------------------------------------------------
	//handle DIV_CONFIG section
	//reset

	function resetconfig(){
		zoom(50);
		document.getElementById('zoomsliderbar').value=50;
		document.getElementById("bgcolor").value="#333333";
		setbgcolor();
		resetPath();
	}
	//set bgcolor
	function setbgcolor(){
		var color=document.getElementById("bgcolor").value;
		d3.select("svg").style("background", color);
	}

	//handle explore tf

	function exploretf(){
		resetPath();
		var tfinput=d3.select("#tfName").property("value").toUpperCase();
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
		if (GL.indexOf(tfinput)!=-1){		
			plottf(tfinput);
		}
		//alert("implements explore tf:(1) highlight tf re-regualting path; (2): plot tf expression (3) plot tf target expression, (4) more");
	}


	//plot tf
	function plottf(tfinput){
		var newW2 = open('','_blank','height=600,width=600,left=1400,top=200,scrollbars=yes')
		newW2.document.write("<head><title>Plot TF</title> <link rel='stylesheet' type='text/css' href='style.css'></head><body></body>");
		
		var chart=d3.select(newW2.document.body).append("div")
					.attr("class","plottf")
					.append("svg")
					.attr("fill","black")
					.attr("width",600)
					.attr("height",700)
					//.style("border","1px solid")
					.style("margin-left","10px")
					.style("padding-elft","100px");
					
		var width = 400;
		var height = 200;
		var yoffset=300;
		
		var x = d3.scale.ordinal()
			.rangeRoundBands([0, width], .1);

		var y = d3.scale.linear()
			.range([height, 0]);
		
		var xAxis = d3.svg.axis()
		.scale(x)
		.orient("bottom");
		
		var yAxis = d3.svg.axis()
			.scale(y)
			.orient("left");
			
		var tfindex=GL.indexOf(tfinput);
		var data=[]
		for (node of nodes){
			var ne=node.E[tfindex].toFixed(1);
			var idata={name: 'E'+node.T+'_'+node.ID, value: ne};
			data.push(idata);
		}

		x.domain(data.map(function(d) { return d.name; }));
		y.domain([0, d3.max(data, function(d) { return d.value; })]);
		
		


		var bar = chart.selectAll("g")
		  .data(data)
		  .enter().append("g")
		  .attr("transform", function(d) { return "translate(" + (x(d.name)+100) + ","+(yoffset)+")"; });

		bar.append("rect")
		  .attr("y", function(d) { return y(d.value); })
		  .attr("height", function(d) { return height - y(d.value); })
		  .attr("width", x.rangeBand())
		  .attr("fill","steelblue");

		bar.append("text")
		  .attr("x", x.rangeBand() / 2-12)
		  //.attr("y",100)
		  .attr("y", function(d) { return y(d.value) -13; })
		  .attr("dy", ".75em")
		  .text(function(d) { return d.value; });
		
		chart.append("g")
		  .attr("class", "x axis")
		  .attr("transform", "translate(100," + (height+yoffset) + ")")
		  .call(xAxis)
		  .selectAll("text")
		  .attr("transform","rotate(90)")
		  .style("stroke","black")
		  .attr("y",0)
		  .attr("x",6)
		  .attr("fill","black")
		  .style("text-anchor", "start");
		  //attr("dy","3em")


		chart.append("g")
		  .attr("class", "y axis")
		  .attr("transform", "translate(100," + (0+yoffset) + ")")
		  .call(yAxis)
		  .append("text")
		  .style("stroke","black")
		  .attr("dx","0em")
		  .attr("y",50)
		  .attr("transform","rotate(90)")
		  .attr("dy","0.71em")
		  .attr("fill","black")
		  .style("text-anchor","start")
		  .text("Expression of "+tfinput);

	}
	//handle explorede

	function explorede(){
		resetPath();
		var deinput=d3.select("#deName").property("value").toUpperCase();
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
	}

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


	//zoom slider bar action
	function zoom(newValue){
		document.getElementById("zoomslider").innerHTML=newValue;
		var wd=1000;
		var ht=1200;
		var sv=50;
		var zx=newValue/sv;
		var newwd=wd*zx;
		var newht=ht*zx;
		
		d3.select("#div_svg").select("svg")
		.attr("width",newwd)
		.attr("height",newht)
		.attr("preserveAspectRatio", "none");
	}

	// plot the SVG figure
	function drawTree(root){
		var margin = {top: 50, right: 20, bottom: 50, left: 20},
		width = 1000 - margin.right - margin.left,
		height =1200 - margin.top - margin.bottom;
		var i = 0;
		var tree = d3.layout.tree();
		tree.size([width,height]);
		
		var diagonal = d3.svg.diagonal()
			.projection(function(d) { return [d.x, d.y]; });
		
		//default bgcolor
		var default_bgcolor="#333";
		var bgcolor;
		var input_bgcolor=d3.select("#bgcolor").value;
		colorisOK  = /(^#[0-9A-F]{6}$)|(^#[0-9A-F]{3}$)/i.test(input_bgcolor)
		if (colorisOK){
			bgcolor=input_bgcolor;
		}else{
			bgcolor=default_bgcolor;
		}
		svg = d3.select("#div_svg").append("svg")
			.attr("width", width + margin.right + margin.left)
			.attr("height", height + margin.top + margin.bottom)
			.attr("id","svg")
			.style("background",bgcolor)
			.attr("viewBox","0 0 1000 1200")
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
		  .style("font","8px Arial")
		  .attr("transform", function(d) { 	 
			  return "translate(" + d.x + "," + d.y + ")"; })

		nodeEnter.append("circle")
		  .attr("r", 16)
		  .attr("fill","#fff")
		  .attr("stroke","steelblue")
		  .attr("stroke-width","3px")
		  .on("mouseover",inNode)  
		  .on("mouseout",outNode)
		  .on("click",testClick);

		var textEnter=nodeEnter.append("text")
		  .attr("dy", ".35em")
		  .attr("class","nodetext")
		  .attr("fill","white")
		  //.style("font","8px Arial")
		  .attr("text-anchor", function(d) { 
			  return  "start"; })
		  .style("fill-opacity", 1)

		textEnter.call(addwraptext);
		  
		//creat TF display svg text 
		createEdgeTF();Object
		
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
		  .attr("stroke","#aaa")
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
	}
	//----------------------------------------------------------------------
	//show/hide controls

	//show/hide TF/GE

	function showhideTF(checked){
		
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
	}
	//show/hide GE
	function showhideDE(checked){Object
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
	}
	//----------------------------------------------------------------------
	//hide controls

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

	//show TF 
	function showTF(){
		
		//change the opacity
		d3.select("svg").selectAll(".TFText")
		.style("opacity",1);
		//change the text size
		d3.select("svg").selectAll(".TFText")
		.selectAll("tspan")
		.style("font-size",8)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		
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
		.selectAll("tspan")
		.style("font-size",8)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		
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
		.attr("x","-10em")
		.attr("y","-3em")
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
					.style("font-size",8)
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
		.attr("x","-10em")
		.attr("y","-3em")
		gc_enter.call(addTFText,20);
		
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
					.style("font-size",8)
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

	//get the expression for given node
	function getNodeEx(node){
		var ne=new Array(cells[node.CELL[0]].E.length);
		ne.fill(0);
		
		for(var cell of node.CELL){
			var ce=cells[cell].E;
			for (var x in ce){
				ne[x]+=ce[x];
			}
		}
		for (x in ne){
			ne[x]/=node.CELL.length;
		}
		return ne;
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


	// get Node text
	function getNodeText(d){
		return d.ID;
	}

	//node click function
	function testClick(){
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
		resList.push(["TF","p-value","TF fold change","Mean target fold change","Mean DE target fold change","edgeID"]);

		resList=exportTFEdge(resList,chosenEdge);
		
		tdiv
		.append("p")
		.text("Table 1. TF details for the edge ending at selected node:"+cnode.T+"_"+cnode.ID);
		
		
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
		pnode=cnode.parent;
		var chosenPath=[]
		
		while(pnode!="null"){
			var ec=cnode.endEdge;
			cnode=pnode;
			pnode=pnode.parent;
			chosenPath.push(ec);
		}
		

		var AllTFList=[];
		AllTFList.push(["TF","p-value","TF fold change","Mean target fold change","Mean DE target fold change","edgeID"]);

		for (var cedge of chosenPath){
			AllTFList=exportTFEdge(AllTFList,JSON.parse(JSON.stringify(cedge)));
			
		}

		AllTFList.sort(function(a,b){return a[1]-b[1];});
		tdiv
		.append("p")
		.text("Table 2. TFs along the path ending at selected node:"+cnode.T+"_"+cnode.ID);
		
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
		
		//adding DE details
		cnode=this.__data__;
		pnode=cnode.parent;
		var resDEList=[]; //DE details
		resDEList.push(["DE gene","Expression_E"+pnode.T+"_"+pnode.ID,"Expression_"+cnode.T+"_"+cnode.ID,"Fold change","edgeID"]);

		var resDEList=exportDEEdge(resDEList,chosenEdge);
		tdiv.append("p")
		.text("Table 3. DE gene details for edge ending at selected node:"+cnode.T+"_"+cnode.ID)
		
		createTable(tdiv,"detable",resDEList);
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"detable","dedlink");
			});
		tdiv.append("a")
		.attr("href","#")
		.attr("id","dedlink");
		
		//loading complete
		newW.document.title="loading complete";
		
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
		.text(function(d){return d;});
	}


	//action on mouse out nodeID: E1_16.0_0Proliferative AT2 Early Precursor,6Proliferative Bi-potential Precursor,1

	function outNode(){
		tooltipsvg
		.style("opacity", 0)
		.attr("width",0)
		.attr("height",0)
		.selectAll("*").remove();;
	}

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



	//action on mouse over node
	function oldnodemouseoverfunction(){
		
		var CELL=this.__data__.CELL;
		ctdict={}
		for(x of CELL){
			xlabel=cells[x].typeLabel;
			if (xlabel in ctdict==false){
				ctdict[xlabel]=1;
			}else{
				ctdict[xlabel]+=1;
			}
		}
		
		tooltipdiv
		.html(JSON.stringify(ctdict))
		.style("left",(d3.event.pageX)+"px")
		.style("top",(d3.event.pageY-28)+"px")
		.style("position","absolute")
		.style("opacity",0.9)
		.transition()
		.duration(200);
		
	}

	function nodeClick(){
		var newW=window.open("","","width=300,height=250");
		newW.document.write("<html><body></body></html>");
		var svg = d3.select(newW.document.body)
			.append("svg")
			.attr("width",400)
			.attr("height",200)
			.append("circle")
			.attr("cx",200)
			.attr("cy",100)
			.attr("r",20)
	}

	//create tree structure
	function buildTree(nodes){
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
	}

	// parse json file
	function parseJSON(data){
		data=JSON.parse(data);
		var GL=data[0];
		var CellList=data[1];
		var NodeList=data[2];
		var EdgeList=data[3];
		var dTD=data[4];	
		return [GL,CellList,NodeList,EdgeList,dTD];
	};


	//---------------------------------------------------------------------
	//download functions here

	//download tf
	function downloadde(){
		plswait("dedownloadlink");
		window.setTimeout(creatededownload,10);
	}

	function creatededownload(){
		var resList=[];
		resList.push(["DE gene","Expression_from","Expression_to","Fold change","edgeID"]);

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
				var eto=pnode.E[gindex];
				var fc=eto-efrom;
				var deresList=[g,efrom,eto,fc,edgeid];
				resList.push(deresList);
			}
		}
		//resList=resList.splice(0,10000);
		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#dedownloadlink")
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}


	//download tf
	function downloadtf(){
		plswait("tfdownloadlink");
		window.setTimeout(createtfdownload,10);
	}

	//download tf
	function createtfdownload(){
		var resList=[];
		resList.push(["TF","p-value","TF fold change","Mean target fold change","Mean DE target fold change","edgeID"]);
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
		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#tfdownloadlink")
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
			out+=LSTI+'\n';
		}
		return out;
	}
	//write download link


	//download figure

	function downloadfig(){
		svgToCanvas();

	}

	//svgToCanvas
	function svgToCanvas(){
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
		source = '<?xml version="1.0" standalone="no"?>\n' + source;

		//convert svg source to URI data scheme.
		var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(source);

		//set url value to a element's href attribute.
		document.getElementById("downloadlink").href = url;
		document.getElementById("downloadlink").innerHTML="Figure ready,Right click me to save!"
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

	