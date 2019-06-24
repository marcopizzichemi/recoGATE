import sys
from optparse import OptionParser
from os.path import isfile, join, expanduser, dirname
import os
from tempfile import mkstemp, mkdtemp, NamedTemporaryFile
import pexpect
import struct
from time import strptime
from datetime import date, datetime
import md5

numberOfSegments = 127

def writeReconParFile(fName, src, prefix, sensitivity, parameters):
	f = open(fName, "w")
	nIterations = parameters['iterations']
	nSubsets = parameters['subsets']
	nSubIterations = nIterations * nSubsets
	ffwhm = parameters['ffwhm']
	fmp = parameters['fmp']

	template = \
	"""OSMAPOSLParameters :=
	input file := %(src)s
	output filename prefix := %(prefix)s
	zoom := 1

	sensitivity image := %(sensitivity)s

	Output file format := Interfile
	Interfile Output File Format Parameters :=
	        byte order := littleendian     
	        number format := float
        	number of bytes per pixel := 4
  	End Interfile Output File Format Parameters :=

	number of subsets: = %(nSubsets)d
	start at subset := 0
	number of subiterations := %(nSubIterations)d
	save images at subiteration intervals := %(nSubIterations)d
	start at subiteration number := 1

	projector pair type := Matrix
	Projector Pair Using Matrix Parameters :=
		Matrix type := pem
		PEM matrix parameters :=
			disable caching := 1
			store only basic bins in cache := 1
			number of rays in tangential direction to trace for each bin := 10
			do symmetry 180degrees min phi := 0
			do symmetry swap segment := 1
			do symmetry swap s := 1
			do symmetry shift z := 1
		End PEM matrix parameters :=
	End Projector Pair Using Matrix Parameters :=

	maximum relative change := 100
	minimum relative change := 0
	zero end planes of segment 0 := 0

	inter-update filter type := Separable Cartesian Metz
	Separable Cartesian Metz Filter Parameters :=
		x-dir filter FWHM (in mm):= %(ffwhm)f 
		y-dir filter FWHM (in mm):= %(ffwhm)f
		z-dir filter FWHM (in mm):= %(ffwhm)f
		x-dir filter Metz power:= %(fmp)d
		y-dir filter Metz power:= %(fmp)d
		z-dir filter Metz power:= %(fmp)d
	END Separable Cartesian Metz Filter Parameters :=

	END :=
	"""

	f.write(template % locals())
	f.flush()
	f.close()
	return ("%(prefix)s_%(nSubIterations)d.hv" % locals(), "%(prefix)s_%(nSubIterations)d.v" % locals())


def makeSinogram(binPath, tmpDir, elm2Name, distance, angleList, parameters, normalization):
	if normalization is True:
		catFifo = join(tmpDir, "norm.cat.fifo")
	else:
		catFifo = join(tmpDir, "data.cat.fifo")
	os.mkfifo(catFifo)

	if elm2Name[len(elm2Name)-2:] == 'gz':
                catCmd = "pigz -d -c %(elm2Name)s > %(catFifo)s" % locals()
        else:
                catCmd = "cat %(elm2Name)s > %(catFifo)s" % locals()
        cat = pexpect.spawn("sh -c '%s'" % catCmd)


	if normalization is False:
		prefix = join(tmpDir, "data")
		rotateSource = None
		lmrotate = None
		mergerSource = catFifo
	else:
		prefix = join(tmpDir, "norm")
		rotateSource = catFifo
		mergerSource = prefix + '.merger.fifo'; os.mkfifo(mergerSource);	

	lmfName = prefix + '.lm.fifo'; os.mkfifo(lmfName);

	if normalization is True:
		anglesFileName = prefix + '.angles.txt'

		anglesFile = open(anglesFileName, "w")
		for angle in angleList:
			anglesFile.write("%f\n" % angle);
		anglesFile.close();
		lmrotate = pexpect.spawn("%(binPath)s/lmrotate2 %(anglesFileName)s %(rotateSource)s %(mergerSource)s" % locals(), timeout=1E6)

	else:
		lmrotate = None

	fmtParams = dict(locals())
	fmtParams.update(parameters)
	lmmerger = pexpect.spawn("%(binPath)s/lmmerger2 --energy-low %(eMin)f --energy-high %(eMax)f "\
		  "--time-window %(dMax)f "\
		  "%(lmfName)s %(mergerSource)s" % fmtParams, timeout=1E6)

	rebinner = pexpect.spawn("%(binPath)s/PlanarRebinner %(prefix)s %(lmfName)s" % locals(), timeout=1E6)
	rebinner.expect("How many angular positions were used");
	rebinner.sendline("%d" % 4)
	rebinner.expect("What is the distance between")
	rebinner.sendline("%f" % (distance - 63))
	rebinner.expect("What is the size of the pixel")
	rebinner.sendline("%f" % parameters['voxelSize'])
	rebinner.expect("Output coords");
	rebinner.sendline("n");


	rebinner.expect(pexpect.EOF)
	lmmerger.expect(pexpect.EOF)
	if lmrotate is not None: lmrotate.expect(pexpect.EOF)
	cat.expect(pexpect.EOF)
	
	hsName = prefix + '.hs'
	dsName = prefix + '.s'
	tsName = prefix + '.s.inv'

	if normalization:
		os.system("mv -f %(dsName)s %(tsName)s" % locals());		
		inverter = pexpect.spawn("%(binPath)s/invertsino %(tsName)s %(dsName)s" % locals(), timeout=1E6);
		inverter.expect(pexpect.EOF)

	return hsName, dsName
	

def getDigest(elm2Name, distance, angleList, parameters, binPath):
	digest = md5.new()
	digest.update("version serial 4");
	digest.update(elm2Name)
        nearDistance = round(distance)
        nv = round(distance / parameters['voxelSize'])*parameters['voxelSize']
        digest.update(str(nearDistance))
        digest.update(str(nv))
	digest.update(str(angleList))
	for key, value in parameters.items():
		if key in ['eMin', 'eMax', 'dMax', 'noNorm', 'voxelSize', 'subsets']:
			digest.update(str((key, value)))
			
	for app in ['lmrotate2', 'lmmerger2', 'PlanarRebinner', 'sensitivity']:
		f = open(join(binPath, app), "r")
		digest.update(f.read())
		f.close()

	return digest.hexdigest()


def getSensitivity(examTime, distance, angleList, parameters, basePath, binPath, tmpDir, quiet):
	
	if parameters['noNorm']:
		elm2Name = '/dev/null'
	else:
		nearDistance = int(round(distance))	
		dDir = join(basePath, "pemNormalization", "listMode", str(nearDistance))
#		print "debug: ", dDir
		try:
			byDateFolders = os.listdir(dDir)
		except OSError,e:
			sys.stderr.write("Error: Normalization data for %dmm not found!\n" % nearDistance)
			return None, None

		normalizationDates = [ date(*strptime(s, "%Y-%m-%d")[0:3]) for s in byDateFolders ]

		deltas = [ abs(examTime.date() - d) for d in normalizationDates ]
		minDelta = min(deltas)
		index = deltas.index(minDelta)

		elm2Name = join(dDir, byDateFolders[index], 'data.elm2')
		if isfile(elm2Name + '.gz'):
			elm2Name = elm2Name + '.gz'

	digest =  getDigest(elm2Name, distance, angleList, parameters, binPath)
	sensDir = join(basePath, "cache", "stir-sensitivities")
	prefix = join(sensDir, digest)
	hvSens = prefix + '.hv'
	dvSens = prefix + '.v'


	if os.path.isfile(hvSens) and os.path.isfile(dvSens):
		return hvSens, dvSens

	if not quiet:
		sys.stdout.write("No suitable sensitivity image found, creating a new one...\n")

	if not os.path.isdir(sensDir): os.mkdir(sensDir)

	if not quiet: 
		sys.stdout.write("Creating normalisation sinogram, please wait...\n")

	hsName, dsName = makeSinogram(binPath, tmpDir, elm2Name, distance, angleList, parameters, normalization=True)
	
	parFileName = join(tmpDir, "sens.par");
	hvName, dvName = writeReconParFile(parFileName, hsName, join(tmpDir, "none"), hvSens, parameters)


	sensitivity = pexpect.spawn("%(binPath)s/sensitivity %(parFileName)s" % locals(), timeout=1E6)
	sensitivity.expect("Get attenuation")
	sensitivity.sendline("0")
	sensitivity.expect("Get normalisation")
	sensitivity.sendline(hsName)

	n = numberOfSegments/2
	t0 = datetime.now()
	for i in range(n):
		sensitivity.expect("Starting to process segment");
		progress = 100.0 * (i+1) / n
		delta = datetime.now() - t0
		p = int(1000 * progress)
		eta = (100000 - p) * delta / p
		seconds = eta.days * 24 * 3600 + eta.seconds
		minutes = seconds / 60
		seconds = seconds % 60
		if not quiet: 
			sys.stdout.write("Sensitivity progress: %4.1f%% ETA: %02d:%02d\r" % (progress, minutes, seconds))
			sys.stdout.flush()		
	sensitivity.expect(pexpect.EOF)
	
	if not quiet:
		sys.stdout.write("\n")
	
	return hvSens, dvSens
	
