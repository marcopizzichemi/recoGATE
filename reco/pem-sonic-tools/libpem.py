import shutil
import os
from os.path import join, isdir, isfile, basename
from os import listdir
import re
from struct import Struct
from gzip import GzipFile
from math import pi

structELM2 = Struct("=dbffffffbffffbf")

def mkdir2(d):
	try:
		os.mkdir(d)
	except OSError, e:
		if e.errno != 17:
			raise e

def move(src, dst):
	shutil.move(src, dst)		

def scanDistance(fn):
	if fn[len(fn)-2:] == 'gz':
		f = GzipFile(fn)
	else:
		f = open(fn)

	data = f.read(structELM2.size)
	tt, random, distance, yozRot, x1, y1, z1, e1, n1, x2, y2, z2, e2, n2, dt = structELM2.unpack(data)
	f.close()
	return distance



class Exam:
	def __init__(self, exam):
		self.exam = exam
		self.path = join(exam.path, "pem")

	def storeELM2(self, fileName):
		move(fileName, join(self.path, "data.elm2"))
		return None

	def getELM2Name(self):
		if isfile(join(self.path, "data.elm2.gz")):
			return join(self.path, "data.elm2.gz")
		elif isfile(join(self.path, "data.elm2")):
			return join(self.path, "data.elm2")
		else:
			return None

	def getDistance(self):
		if isfile(join(self.path, "distance.cache")):
			f = open(join(self.path, "distance.cache"))
			s = f.read().strip('\n')
			f.close()
			return float(s)

		
		fn = self.getELM2Name();
		distance = scanDistance(fn)
		
		
		f = open(join(self.path, "distance.cache"), "w")
		f.write(str(distance))
		f.close()
		return distance

	def getAngles(self):
		if isfile(join(self.path, "angles.cache")):
			f = open(join(self.path, "angles.cache"))
			angles = [ float(x) for x in f ]
			f.close()

			return set(angles)

		fn = self.getELM2Name();
                if fn[len(fn)-2:] == 'gz':
                        f = GzipFile(fn)
                else:
                        f = open(fn)

		angles = []
		while True:
			data = f.read(structELM2.size)
			if len(data) != structELM2.size:
				break;
			tt, random, distance, yozRot, x1, y1, z1, e1, n1, x2, y2, z2, e2, n2, dt = structELM2.unpack(data)
			if yozRot not in angles:
				angles.append(yozRot)
		f.close()

		angles = [ x * 180 / pi for x in angles ]

		f = open(join(self.path, "angles.cache"), "w")
		for angle in angles:
			f.write("%f\n" % angle)
		f.close()

		return self.getAngles()
		

	def setInjectedDose(self, dose):
		f = open(join(self.path, "injectedDose"), 'w')
		f.write(str(dose))
		f.close()

	def getInjectedDose(self):
		try:
			f = open(join(self.path, "injectedDose"))
			s = f.read().strip('\n')
			f.close()
			return float(s)
		
		except IOError, e:
			return None


	def setInjectedDoseTime(self, time):
		f = open(join(self.path, "injectedDoseTime"), "w")
		f.write(str(time))
		f.close()
		return None

	def getInjectedDoseTime(self):
		f = open(join(self.path, "injectedDoseTime"), "r")
		s = f.read().strip('\n');
		f.close();
		return datetime(*strptime(s, "%Y-%m-%d %H:%M:%S")[0:6])

	def setIsotopeHalfLifeTime(self, time):
		f = open(join(self.path, "isotopeHalfLifeTime"), "w")
		f.write(str(time))
		f.close()
		return None

	def getIsotopeHalfLifeTime(self):
		f = open(join(self.path, "injectedDoseTime"), "r")
		s = f.read().strip('\n');
		f.close();
		return float(s)
		
## 	def storePreview(self, hFile, dFile):
## 		path = join(self.path, "preview")
## 		mkdir2(path)
## 		dstHeader = join(path, "data.hv")
## 		dstData   = join(path, "data.v")

## 		hIn = open(hFile, "r");
## 		hOut = open(dstHeader, "w");
## 		for line in hIn:
## 			if re.match('name of data file', line, re.I):
## 				hOut.write("name of data file := %s\n" % 'data.v')
## 			else:
## 				hOut.write(line)
## 		hIn.close()
## 		hOut.close()
## 		move(dFile, dstData)
## 		return None

## 	def getPreviewName(self):
## 		return (join(self.path, "preview", "data.hv"), join(self.path, "preview", "data.v"))

	def getPreview(self):
		path = join(self.path, "preview")
		mkdir2(path)
		return Processing(path)

	def getProcessings(self):
		rn = re.compile("[0-9]+")
		return sorted([ Processing(join(self.path, dirName)) 
			for dirName  in listdir(join(self.path)) \
			 if rn.match(dirName) and isdir(join(self.path, dirName)) ], key=Processing.getID)

	def getProcessingByID(self, ppid):
		processings = self.getProcessings();
		ppids = [ pp.getID() for pp in processings ]
		return processings[ppids.index(ppid)]

	def addProcessing(self):
		maxID = max([0]+[ max([0]+[ max([0]+[ pr.getID() for pr in e.pem().getProcessings() ]) for e in p.getExams() ]) for p in self.exam.patient.db.getPatients() ])
		newID = maxID + 1
		path = join(self.path, "%05d" % newID)
		mkdir2(path)
		return Processing(path)
		

class Processing:
	def __init__(self, path):		
		self.path = path
		return None

	def getID(self):
		return int(basename(self.path))

## 	def setFlags(self, flags):
## 		f = open(join(self.path, "flags"), "w")
## 		for f in flags:
## 			f.write("%s\n", f)
## 		f.close()
## 		return None

## 	def getFlags(self):
## 		f = open(join(self.path, "flags"), "r");
## 		flags = [ x[:-1] for x in f ]
## 		f.close()
## 		return flags

## 	def getSinogramName(self):
## 		return (join(self.path, "data.hs"), join(self.path, "data.s"))

## 	def storeSinogram(self, headerFile, dataFile):
## 		dstHeader = join(self.path, "data.hs")
## 		dstData   = join(self.path, "data.s")

## 		hIn = open(headerFile, "r");
## 		hOut = open(dstHeader, "w");
## 		for line in hIn:
## 			if re.match('name of data file', line, re.I):
## 				hOut.write("name of data file := %s\n" % "data.s")
## 			else:
## 				hOut.write(line)
## 		hIn.close()
## 		hOut.close()
## 		move(dataFile, dstData)
## 		return None

	def getImageName(self):
		return (join(self.path, "data.hv"), join(self.path, "data.v"))

	def storeImage(self, headerFile, dataFile):
		dstHeader = join(self.path, "data.hv")
		dstData   = join(self.path, "data.v")

		hIn = open(headerFile, "r");
		hOut = open(dstHeader, "w");
		for line in hIn:
			if re.match('name of data file', line, re.I):
				hOut.write("name of data file := %s\n" % 'data.v')
			else:
				hOut.write(line)
		hIn.close()
		hOut.close()
		move(dataFile, dstData)
		return None
	
	def setDescription(self, desc):
		f = open(join(self.path, "recon_par"), "w")
		f.write(desc);
		f.close()
		return None

	def getDescription(self):
		f = open(join(self.path, "recon_par"), "r")
		desc = f.read()
		f.close()
		return desc
