import abc

class File:
	def __init__(self,fname):
		self.name=fname
	def read(self):
		raise NotImplementedError('Please implement this method')
	def Linux2Win(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip()+'\r\n' for item in lf]
		lf=''.join(lf)
		f=open(self.name+'.win.txt','w')
		f.write(lf)
		f.close()


class TabFile(File):
	def read(self,sep):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip().split(sep) for item in lf]
		return lf

class FastaFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split('>')[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=[lf[i][0],''.join(lf[i][1:])]
		return lf

class TFBSFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split('>')[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=[item for item in lf[i] if item!='']
		return lf
	
class LineFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip() for item in lf]
		return lf
			

class MotifFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split('>')[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=[item for item in lf[i] if item!='']
			lf[i]=[item.split() for item in lf[i]]
			lf[i][0]=lf[i][0][0]
		return lf

class MicroArrayFile(File):
	def __init__(self,FileList):
		self.FileList=FileList
		self.E=self.read()
		
	def read(self):
		out=[]
		for i in self.FileList:
			fi=TabFile(i).read('\t')
			ID=[item[0] for item in fi]
			fi=[item[1] for item in fi]
			out.append(fi)
		L=len(ID)
		ME=[]
		for i in range(L):
			ei=[item[i] for item in out]
			ME.append([ID[i]]+ei)
		return ME
		
	def toGeneSymbol(self,geneMap):
		gm=TabFile(geneMap).read('\t')
		dgm={}
		for i in gm:
			if len(i)>9:
				mi=i[9].split('//')
				if len(mi)>1:
					mi=mi[1].strip()
					dgm[i[0]]=mi
		for i in range(len(self.E)):
			if self.E[i][0] in dgm:
				self.E[i][0]=dgm[self.E[i][0]]
			else:
				self.E[i]=[]
		self.E=[item for item in self.E if item!=[]]
		return self.E
		
		
		
		
