#!/usr/bin/env python

"""
author: jun ding
date: Jul.26,2016, 15:04
function: a GUI for Lung Single cell tools 
"""

#-----------------------------------------------------------------------

import pdb,sys,os
if sys.version_info[0]<3:
	import Tkinter as tk
	import tkFileDialog
	from Tkinter import *
	import ttk
	from ScrolledText import *
else:
	import tkinter as tk
	from tkinter import *
	import tkinter.filedialog as tkFileDialog
	import tkinter.ttk as ttk
	import tkinter.scrolledtext as ScrolledText
	from tkinter.scrolledtext import *

import subprocess
import multiprocessing
import threading
import pkg_resources 
#----------------------------------------------------------------------

class App:
	def __init__(self, parent):
		self.myParent = parent
		#---------------------------------------------------------------
		# frme 0: Software Logo
		self.frameLogo=Frame(parent,background='white')
		imgpath=pkg_resources.resource_filename(__name__,"img/logo.gif")
		#pdb.set_trace()
		img=PhotoImage(file=imgpath)
		imgLabel=Label(self.frameLogo,image=img)
		imgLabel.image=img
		imgLabel.pack()
		self.frameLogo.grid(row=0,column=1,sticky='ew')
		#---------------------------------------------------------------
		# frame 1: expression frame
		self.frameEx=Frame(parent)

		self.lex1=ttk.Label(self.frameEx,text="Please Read in the Single Cell Expression data: ")

		self.fileName=""
		self.bex1=ttk.Button(self.frameEx,text='Read in Expression',command=self.readEx)
		self.vex1=StringVar()
		self.lex2=ttk.Label(self.frameEx,textvariable=self.vex1)

		self.lex1.pack(side='left',fill=None,expand=False,padx=6,pady=20)
		self.bex1.pack(side='left',fill=None,expand=False,padx=20,pady=20)
		self.lex2.pack(side='left',fill=None,expand=False,padx=40,pady=20)

		# frame grid
		row_ct=1
		self.frameEx.grid(row=row_ct,column=1,sticky='ew')
		row_ct+=1
		#---------------------------------------------------------------
		# frame 2: TF-DNA interaction file
		self.frameTFDNA=Frame(parent)
		self.ltf1=ttk.Label(self.frameTFDNA,text="Please read in the TF-DNA interaction file: ")
		self.btf1=ttk.Button(self.frameTFDNA,text='Read in TF-DNA',command=self.readTF)
		self.vtf1=StringVar()
		self.ltf2=ttk.Label(self.frameTFDNA,textvariable=self.vtf1)

		self.ltf1.pack(side='left',fill=None,expand=False,padx=6)
		self.btf1.pack(side='left',fill=None,expand=False,padx=40)
		self.ltf2.pack(side='left',fill=None,expand=False,padx=20)

		# frame grid
		self.frameTFDNA.grid(row=row_ct,column=1,sticky='ew')
		row_ct+=1

		#---------------------------------------------------------------
		# frame 3: K frame (optimal number of clusters for each time point)
		self.frameK=Frame(parent)

		self.KName="auto"
		self.lk1=ttk.Label(self.frameK,text="Please Specificy the Optimal Number of Cluster K for Each Time Point: ")
		self.vk1=IntVar()
		self.rbk1=ttk.Radiobutton(self.frameK,text='user-defined',variable=self.vk1,value=1,command=self.readK)
		self.rbk2=ttk.Radiobutton(self.frameK,text='auto',variable=self.vk1,value=0)

		self.lk1.pack(side='left',padx=6,pady=20)
		self.rbk1.pack(side='left',padx=29)
		self.rbk2.pack(side='left',padx=20)

		self.frameK.grid(row=row_ct,column=1,sticky='ew')
		row_ct+=1
		#---------------------------------------------------------------
		# frame 4: Output
		self.frameoutfolder=Frame(parent)

		self.lo1=ttk.Label(self.frameoutfolder,text="Please Specify the Output Folder Name: ")

		self.ev1=StringVar()
		self.eo1=ttk.Entry(self.frameoutfolder,textvariable=self.ev1)

		self.lo1.pack(side='left',padx=6)
		self.eo1.pack(side='left',padx=20)

		self.frameoutfolder.grid(row=row_ct,column=1,sticky='ew')
		row_ct+=1
		#--------------------------------------------------------------

		# frame 5: run
		self.framerun=Frame(parent)
		self.br1=ttk.Button(self.framerun,text='Click to Run!', command=self.run)
		self.br1.pack(side='left',padx=6,pady=20)
		self.framerun.grid(row=row_ct,column=1,sticky='ew')
		row_ct+=1

		#---------------------------------------------------------------


		# seperator of input and output
		s1 = ttk.Separator(parent, orient=HORIZONTAL)
		s1.grid(row=row_ct,column=1,sticky='ew',pady=10)
		row_ct+=1

		#---------------------------------------------------------------
		# frame 6: progress bar

		self.framePb=Frame(parent)

		self.lpb1=ttk.Label(self.framePb,text='Running Progress: ')
		self.pb=ttk.Progressbar(self.framePb,orient='horizontal',mode="determinate",length=580)
		self.lpb1.pack(side='left',padx=6,pady=10)
		self.pb.pack(side='left',padx=20)


		self.framePb.grid(row=row_ct,column=1,sticky='ew')
		row_ct+=1

		#---------------------------------------------------------------


		#---------------------------------------------------------------
		# frame 7: output display area

		self.framedisplay=Frame(parent)
		self.ld1=ttk.Label(self.framedisplay,text='Running log: ')

		self.textarea1=ScrolledText(self.framedisplay)

		self.ld1.pack(side='left',padx=6,pady=20)
		self.textarea1.pack(side='left',padx=50)

		self.framedisplay.grid(row=row_ct,column=1,sticky='ew',pady=20)
		row_ct+=1


	def readEx(self):
		self.fileName=''
		self.fileName=tkFileDialog.askopenfilename()
		self.vex1.set(self.fileName.split('/')[-1])

	def readK(self):
		self.KName=''
		self.KName=tkFileDialog.askopenfilename()

	def readTF(self):
		self.TFName=''
		self.TFName=tkFileDialog.askopenfilename()
		self.vtf1.set(self.TFName.split('/')[-1])

	def run(self):
		#check whether input valid
		if self.isInputValid():
			# run
			self.p=threading.Thread(target=self.trun)
			self.p.start()
		
		
	def trun(self):
		#pdb.set_trace()
		self.o=self.ev1.get()
		#pdb.set_trace()
		python=sys.executable
		proc = subprocess.Popen([python,'scdiff.py','-i',self.fileName,'-t',self.TFName,'-k',self.KName,'-o',self.o],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		
		ct=0
		maxCT=10000
		self.textarea1.insert(tk.INSERT,"starting...")
		self.textarea1.see(tk.END)
		for line in iter(proc.stdout.readline,''):
			self.textarea1.insert(tk.INSERT,line)
			self.textarea1.see(tk.END)
			self.updateProgress(maxCT,ct)
			ct+=1
			self.textarea1.update_idletasks()
		self.updateProgress(maxCT,maxCT)
		
		proc.stdout.close()
		proc.wait()
		self.textarea1.insert(tk.INSERT,"end!")
		self.textarea1.see(tk.END)
		
	def updateProgress(self,MAX,VAL):
		self.pb['maximum']=int(MAX)
		self.pb['value']=int(VAL)

	def isInputValid(self):
		flag=0
		# check if expression file exists
		if os.path.isfile(self.fileName):
			flag+=1
		else:
			self.textarea1.insert(tk.INSERT,'\nError: Input single Cell Expression data not found!\n')
			self.textarea1.see(tk.END)
			self.textarea1.update_idletasks()

		# check if K valid
		if self.vk1.get()==0:
			flag+=1
		elif self.vk1.get()==1:
			if os.path.isfile(self.KName):
				flag+=1
			else:
				self.textarea1.insert(tk.INSERT,'\nError: User-defined Optimal Number of Clusters K file not found!\n')
				self.textarea1.see(tk.END)
				self.textarea1.update_idletasks()
		
		# check if output name valid
		try: 
			if os.path.exists(self.ev1.get())==False:
				os.mkdir(self.ev1.get())
			flag+=1
		except:
			self.textarea1.insert(tk.INSERT,'\nError: Output folder name invalid!\n')
			self.textarea1.see(tk.END)
			self.textarea1.update_idletasks()
				
		#pdb.set_trace()
		if flag>=3:
			return True
		return False 


def main():
	root=tk.Tk()
	#root.columnconfigure(0, weight=1)
	app=App(root)
	root.geometry("800x900")
	root.title("SCDIFF")
	root.mainloop()
		 	 
if __name__=='__main__':
	main()
	
	
