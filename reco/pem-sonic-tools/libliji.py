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

def getDigest(elm2Name, distance, angleList, parameters, binPath):
        digest = md5.new()
        digest.update("version serial 5");
        digest.update(elm2Name)
        nearDistance = round(distance)
        nv = round(distance / parameters['voxelSize'])*parameters['voxelSize']
        digest.update(str(nearDistance))
        digest.update(str(nv))
        digest.update(str(angleList))
        for key, value in parameters.items():
                if key in ['eMin', 'eMax', 'dMax', 'noNorm', 'voxelSize', 'rays']:
                        digest.update(str((key, value)))
        for app in ['elm2todkfz', 'Norm_Total_Gen']:
                f = open(join(binPath, app), "r")
                digest.update(f.read())
                f.close()

        return digest.hexdigest()


def getNormalization(examTime, distance, angleList, parameters, basePath, binPath, tmpDir, quiet):

        if parameters['noNorm'] is True:
                elm2Name = '/dev/null'
        else:
                nearDistance = int(round(distance))     
                dDir = join(basePath, "pemNormalization", "listMode", str(nearDistance))
#               print "debug: ", dDir
                try:
                        byDateFolders = os.listdir(dDir)
                except OSError,e:
                        sys.stderr.write("Error: Normalization data for %dmm not found!\n" % nearDistance)
                        return None, None

                normalizationDates = [ date(*strptime(s, "%Y-%m-%d")[0:3]) for s in byDateFolders ]

                deltas = [ abs(examTime.date() - d) for d in normalizationDates ]
                minDelta = min(deltas)
                index = deltas.index(minDelta)

                print "Using normalization from ", str(normalizationDates[index])

                elm2Name = join(dDir, byDateFolders[index], 'data.elm2')
                if isfile(elm2Name + '.gz'):
                        elm2Name = elm2Name + '.gz'

        digest =  getDigest(elm2Name, distance, angleList, parameters, binPath)
        sensDir = join(basePath, "cache", "dkfz-norm")
        normPrefix = join(sensDir, digest)

        crystalNorm = normPrefix + '_crystals'
        sensImage = normPrefix + '_total'

        if os.path.isfile(sensImage):
                return normPrefix

        if not quiet:
                sys.stdout.write("No suitable sensitivity image found, creating a new one...\n")

        if not os.path.isdir(sensDir): os.mkdir(sensDir)

        if not quiet: 
                sys.stdout.write("Processing normalization data, please wait...\n")


        catFifo = join(tmpDir, "norm.cat.fifo"); os.mkfifo(catFifo)
        if elm2Name[len(elm2Name)-2:] == 'gz':
                catCmd = "pigz -d -c %(elm2Name)s > %(catFifo)s" % locals()
        else:
                catCmd = "cat %(elm2Name)s > %(catFifo)s" % locals()
        cat = pexpect.spawn("sh -c '%s'" % catCmd)


        prefix = join(tmpDir, "norm")           
        lmfName = prefix + '.lm.fifo'; os.mkfifo(lmfName);
        anglesFileName = prefix + '.angles.txt'

        anglesFile = open(anglesFileName, "w")
        for angle in angleList:
                anglesFile.write("%f\n" % angle);
        anglesFile.close();
                                                                                     
        fmtParams = dict(locals())
        fmtParams.update(parameters)

        elm2todkfzCmd = "%(binPath)s/elm2todkfz --energy-low %(eMin)f --energy-high %(eMax)f "\
                        "--time-window %(dMax)f --use-direct-correction "\
                        "-o %(lmfName)s %(catFifo)s" % fmtParams

        normGenCmd = "nice -n 19 %(binPath)s/Norm_Total_Gen %(distance)f %(voxelSize)f %(anglesFileName)s "\
                        "--rays %(rays)d "\
                        "%(lmfName)s %(normPrefix)s" % fmtParams
        #print elm2todkfzCmd
        #print normGenCmd

        elm2todkfz = pexpect.spawn(elm2todkfzCmd, timeout=1E6)
        normGen = pexpect.spawn(normGenCmd, timeout=1E6)

        normGen.expect(pexpect.EOF)
        elm2todkfz.expect(pexpect.EOF);
        cat.expect(pexpect.EOF)

        if os.path.isfile(sensImage):
                return normPrefix
        else:
                return None
        
        
