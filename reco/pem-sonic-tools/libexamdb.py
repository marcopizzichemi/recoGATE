

from os.path import join, isdir, basename
from os import listdir
import os
import shutil

from time import strptime, strftime
from datetime import date, datetime

import libpem

def mkdir2(d):
	try:
		os.mkdir(d)
	except OSError, e:
		if e.errno != 17:
			raise e

def move(src, dst):
	shutil.move(src, dst)		

class ExamDB:
	def __init__(self, path):
		self.path = path
		mkdir2(path)
		return None

	def getPatients(self):
		return sorted([ Patient(self, join(self.path, dirName))
			 for dirName  in listdir(self.path) \
			 if isdir(join(self.path, dirName)) ], key=Patient.getID)

	def getPatientByID(self, pid):
		patients = self.getPatients();
		patientIDs = [ p.getID() for p in patients ]
		return patients[patientIDs.index(pid)]

	def addPatient(self):
		maxID = max([0]+[ p.getID() for p in self.getPatients() ])
		newID = maxID + 1
		path = join(self.path, "%05d" % newID)
		mkdir2(path)
		patient = Patient(self, path)
		return patient
	
class Patient:
	def __init__(self, db, path):
		self.db = db
		self.path = path
		return None

	def getID(self):
		return int(basename(self.path))

	def getName(self):
		f = open(join(self.path, "patientName"), "r")
		name = f.read().strip('\n')
		f.close();
		return name

	def setName(self, patientName):
		f = open(join(self.path, "patientName"), "w")
		f.write(patientName)
		f.close()
		return None

	def setBirthDate(self, birthDate):
		f = open(join(self.path, "birthDate"), "w")
		f.write(str(birthDate))
		f.close()
		return None

	def getBirthDate(self):
		f = open(join(self.path, "birthDate"), "r")
		s = f.read().strip('\n')
		f.close();
		return date(*strptime(s, "%Y-%m-%d")[0:3])

	def getExams(self):
		return sorted([ Exam(self, join(self.path,dirName))
			 for dirName  in listdir(self.path) \
			 if isdir(join(self.path,dirName)) ], key=Exam.getID)

	def getExamByID(self, eid):
		exams = self.getExams()
		examIDs = [ e.getID() for e in exams ]
		return exams[examIDs.index(eid)]

	def addExam(self):
		maxID = max([0]+[ max([0]+[ e.getID() for e in p.getExams() ]) for p in self.db.getPatients() ])
		newID = maxID + 1
		path = join(self.path, "%05d" % newID)
		mkdir2(path)
		mkdir2(join(path, "pem"))
		return Exam(self, path)



class Exam:
	def __init__(self, patient, path):
		self.patient = patient
		self.path = path

		mkdir2(join(self.path, "pem"))
		mkdir2(join(self.path, "us"))
		return None

	def getID(self):
		return int(basename(self.path))

	def setDescription(self, description):
		f = open(join(self.path, "description"), "w")
		f.write(description)
		f.close()
		return None

	def getDescription(self):
		f = open(join(self.path, "description"), "r")
		d = f.read().strip('\n')
		f.close()
		return d

	def setTime(self, time):
		f = open(join(self.path, "time"), "w")
		f.write( strftime("%Y-%m-%d %H:%M:%S", time.timetuple()))
		f.close()
		return None

	def getTime(self):
		f = open(join(self.path, "time"), "r")
		s = f.read().strip('\n');
		f.close();
		return datetime(*strptime(s, "%Y-%m-%d %H:%M:%S")[0:6])

	def setPatientWeight(self, weight):
		f = open(join(self.path, "patientWeight"), "w")
		f.write(str(weight))
		f.close()
		return None

	def getPatientWeight(self):
		try:
			f = open(join(self.path, "patientWeight"), "r")
			s = f.read().strip('\n');
			f.close();
			return float(s)
		except IOError, e:
			return None

	def pem(self):
		return libpem.Exam(self)

