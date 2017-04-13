#!/usr/bin/python


import sys
from optparse import OptionParser
from ConfigParser import SafeConfigParser, NoOptionError, NoSectionError
from os.path import join, expanduser, dirname
import os.path
import os
from tempfile import mkdtemp
import pexpect
from datetime import datetime

import libexamdb
from libstir import writeReconParFile, makeSinogram, getSensitivity, numberOfSegments
import libliji


def main(argv):
    config = SafeConfigParser()
    config.read(expanduser('~/.pem.rc'))
    try:
        basePath = config.get('global', 'base')
    except    NoOptionError, e:
        basePath = None
    except NoSectionError, e:
        basePath = None

    cmd = argv[0]
    parser = OptionParser()
    parser.add_option("--base", dest="base", default=basePath)
    parser.add_option("--batch", dest="batch", default=False, action="store_true")
    parser.add_option("--patient", dest="patient", default=1)
    parser.add_option("--exam", dest="exam")
    parser.add_option("--preview", dest="preview", default=False, action="store_true")
    parser.add_option("--energy-low", dest="eMin", default=400)
    parser.add_option("--energy-high", dest="eMax", default=650)
    parser.add_option("--time-window", dest="dMax", default=4)
    parser.add_option("--voxel-size", dest="voxelSize", default=2)
    parser.add_option("--iterations", dest="iterations", default=6)
    parser.add_option("--no-norm", dest="noNorm", default=False, action="store_true")
    parser.add_option("--algo", dest="algo", default="stir")
    parser.add_option("--filter-fwhm", dest="ffwhm", default=1.5)
    parser.add_option("--filter-power", dest="fmp", default=0)
    parser.add_option("--no-osem", dest="osem", default=True, action="store_false")
    parser.add_option("--dw-correction", dest="dwCorrection", default=False, action="store_true")
    parser.add_option("--rays", dest="rays", default=0)
    parser.add_option("--corrections", dest="corrections", default=False)
    parser.add_option("--segmentation", dest="segmentation", default="auto")
    parser.add_option("--segmentation-threshold", dest="segmentationThreshold", default=1)

    (options, argv) = parser.parse_args(argv)

    for o in ["patient", "exam"]:
        if not options.__dict__[o]:
            print 'Usage: %(cmd)s [ --base <base path> ] --patient <patient ID> --exam <exam id>' % locals()
            return None

    parameters = {'eMin': float(options.eMin), 'eMax': float(options.eMax), 'dMax': float(options.dMax),
                  'voxelSize': float(options.voxelSize), 'iterations': int(options.iterations),
                  'noNorm': options.noNorm, 'algo': options.algo.upper(), 'ffwhm': float(options.ffwhm),
                  'fmp': float(options.fmp), 'dwCorrection': options.dwCorrection, 'rays': int(options.rays)}

    if parameters['algo'] == 'STIR':
        parameters['subsets'] = 4
        parameters['noNorm'] = options.noNorm

        if options.preview is True:
            parameters['voxelSize'] = 4.0
            parameters['iterations'] = 1
            parameters['subsets'] = 4

    elif parameters['algo'] == 'DKFZ':
        parameters['osem'] = options.osem
        parameters['osem'] = False

        if options.preview == True:
            parameters['voxelSize'] = 4.0

    else:
        print 'Algorithm must be "stir" or "dkfz"'
        return 1

    pid = int(options.patient)
    eid = int(options.exam)

    basePath = options.base
    if basePath is None:
        print "Need to specify a base path in %s or with --base!" % expanduser('~/.pem.rc')

    basePath = expanduser(basePath)
    db = libexamdb.ExamDB(join(basePath, "exam"))

    tmpDir = join(basePath, "tmp")
    tmpDir = mkdtemp(dir=tmpDir)
    srcPath = dirname(__file__)
    if srcPath is "": srcPath = '.'
    binPath = srcPath

    patient = db.getPatientByID(pid)
    exam = patient.getExamByID(eid)
    elm2Name = exam.pem().getELM2Name()

    if not options.batch:
        sys.stdout.write("Patient: %s\n" % patient.getName())
        sys.stdout.write('Exam description: %s\n' % exam.getDescription())
        sys.stdout.write('Exam time: %s\n' % str(exam.getTime()))
        sys.stdout.write('Reconstruction parameters: %s\n' % str(parameters))
        sys.stdout.write("Scaning data...\r")
        sys.stdout.flush()

    distance = exam.pem().getDistance()
    angleList = exam.pem().getAngles()

    # Normalize anglelist...
    angleList = [round(x * 2) / 2 for x in angleList]
    angleFrequencies = [angleList.count(x) for x in frozenset(angleList)]
    if max(angleFrequencies) == min(angleFrequencies):
        angleList = [x for x in frozenset(angleList)]

    sys.stdout.write("Acquisition distance: %f mm\n" % distance)
    sys.stdout.write("Acquisition angles: %s\n" % (', ').join([str(x) for x in angleList]))
    sys.stdout.flush()

    if parameters['algo'] == 'STIR':
        hvSens, dv_sens = getSensitivity(exam.getTime(), distance, angleList, parameters, basePath, binPath, tmpDir,
                                        options.batch)

        if hvSens is None:
            os.system("rm -rf %s" % tmpDir)
            return None

        if not options.batch:
            sys.stdout.write("Creating data sinogram, please wait...\n")
            hsName, ds_name = makeSinogram(binPath, tmpDir, elm2Name, distance, angleList, parameters,
                                          normalization=False)

            parFileName = join(tmpDir, "data.par")
            hv_name, dvName = writeReconParFile(parFileName, hsName, join(tmpDir, "data"), hvSens, parameters)
            osmapol = pexpect.spawn("%(binPath)s/OSMAPOSL %(parFileName)s" % locals(), timeout=1E6)
            n = parameters['iterations'] * parameters['subsets'] * (numberOfSegments / 2)
            t0 = datetime.now()
            for i in range(n):
                osmapol.expect("Starting to process segment");
                progress = 100.0 * (i + 1) / n
                delta = datetime.now() - t0
                p = int(1000 * progress)
                eta = (100000 - p) * delta / p
                seconds = eta.days * 24 * 3600 + eta.seconds
                minutes = seconds / 60
                seconds %= 60
                if options.batch is False:
                    sys.stdout.write("Reconstruction progress: %4.1f%% ETA: %02d:%02d\r" % (progress, minutes, seconds))
                    sys.stdout.flush()
            osmapol.expect(pexpect.EOF)

        if options.batch is False:
            sys.stdout.write("\n")

        descString = "STIR; %(voxelSize)1.1f mm voxel; %(iterations)d iterations; %(eMin)d-%(eMax)d keV; %(dMax)1.1f ns; filter %(ffwhm)2.1f FWHM, MP %(fmp)d" % parameters
        if parameters['noNorm'] is True:
            descString += "; no normalization"
        if parameters['subsets'] != 4:
            descString += "; %(subsets)d subsets" % parameters

    elif parameters['algo'] == 'DKFZ':

        ### Files

        cat_fifo = join(tmpDir, "data.cat.fifo")
        os.mkfifo(cat_fifo)
        lmf_name = join(tmpDir, "data.lm")
        #os.mkfifo(lmf_name)
        scatter_total_events_file = join(tmpDir, "scatter_total.txt")
        rand_acq = join(tmpDir, "rand_acq.lm")
        rand_name = join(tmpDir, "rand_data.lm")
        rand_combined = join(tmpDir, "rand_comb.lm")
        lmf_name_corr = join(tmpDir, "data_corr.lm")
        rand_combined_corr = join(tmpDir, "rand_corr.lm")
        mesh_name = "tfinal_mesh.vtk" # mesh_name is hardcoded in MeshCreator

        output_prefix = join(tmpDir, "data_img")

        if elm2Name[len(elm2Name) - 2:] == 'gz':
            cat_cmd = "pigz -d -c %(elm2Name)s > %(cat_fifo)s" % locals()
        else:
            cat_cmd = "cat %(elm2Name)s > %(cat_fifo)s" % locals()

        norm_prefix = libliji.getNormalization(exam.getTime(), distance, angleList, parameters, basePath, binPath,
                                              tmpDir, options.batch)

        fmt_params = parameters
        if fmt_params['osem']:
            fmt_params['osem'] = 1
        else:
            fmt_params['osem'] = 0

        fmt_params.update(locals())

        elm2todkfz_cmd = "%(binPath)s/elm2todkfz -o %(lmf_name)s --randoms-file %(rand_acq)s " \
                         "--energy-low %(eMin)f --energy-high %(eMax)f " \
                         "--time-window %(dMax)f --scatter-total-events %(scatter_total_events_file)s " \
                         "%(cat_fifo)s" % fmt_params

        print elm2todkfz_cmd

        if fmt_params['dwCorrection']:
            elm2todkfz_cmd += " --dw-correction true"

        recon_cmd = "nice -n 19 %(binPath)s/ClearPEM_LMRec -i %(lmf_name)s -o %(output_prefix)s " \
                    " -d %(distance)f -n %(norm_prefix)s --pixel-length %(voxelSize)f " \
                    " --rays %(rays)d " \
                    " --iterations %(iterations)d --osem %(osem)d " \
                    " --filter-fwhm %(ffwhm)f --filter-power %(fmp)d" % fmt_params

        # print fmt_params ## DEBUG info
        hv_name = '%s_%03d.hv' % (output_prefix, parameters['iterations'] - 1)

        print ">>> Reading data..."
        # Extract data
        cat1 = pexpect.spawn("sh -c '%s'" % cat_cmd)
        elm2todkfz = pexpect.spawn(elm2todkfz_cmd, timeout=1E6)
        elm2todkfz.expect(pexpect.EOF)
        
        recon = pexpect.spawn(recon_cmd, timeout=1E6, logfile=sys.stdout)

        # elm2todkfz.expect(pexpect.EOF)

        print ">>> Basic reconstruction..."
        n = parameters['iterations']
        t0 = datetime.now()
        for i in range(n):
            recon.expect('Finished iteration')
            if options.batch is False:
                progress = 100.0 * (i + 1) / n
                delta = datetime.now() - t0
                p = int(1000 * progress)
                eta = (100000 - p) * delta / p
                seconds = eta.days * 24 * 3600 + eta.seconds
                minutes = seconds / 60
                seconds %= 60
                sys.stdout.write("Reconstruction progress: %4.1f%% ETA: %02d:%02d\r" % (progress, minutes, seconds))
                sys.stdout.flush()

        recon.expect(pexpect.EOF)

        # Corrections #########
        if options.corrections:
            mask_type = options.segmentation
            segmentation_threshold = options.segmentationThreshold
            # This number shall never be lower than 100e6
            scatter_min_photons = 100e6
            scatter_total_events = int(open(scatter_total_events_file).readlines()[0])
            
            os.chdir(binPath)
            
            # HACK allow for a custom segmentation.
            segmentedImage = "teste_segmented_holesFilled"
            bettaOverride = "betta_override"
            if os.path.isfile("%s.nii" % bettaOverride):
                print ">>>>>> Python vede e usa betta override! "
                segmentedImage = bettaOverride
            
            fmt_params.update(locals())

            #mask_cmd = "%(binPath)s/MaskCreator %(hv_name)s %(mask_type)s %(segmentation_threshold)s" % fmt_params
            #print mask_cmd
            
            mesh_cmd = "%(binPath)s/MeshCreator %(hv_name)s %(mask_type)s %(segmentation_threshold)s " % fmt_params
            
            true_corr_cmd = "%(binPath)s/CorrectAttenuation_SV %(mesh_name)s %(lmf_name)s "\
                            " %(lmf_name_corr)s LM" % fmt_params
            
            scatter_cmd = "%(binPath)s/CorrectScatter_multiCore_LMSW %(segmentedImage)s.nii "\
                          " Scatter_out %(scatter_min_photons)s %(distance)d %(eMin)d %(scatter_total_events)s" % fmt_params
            print scatter_cmd
            random_events_cat = "cat %(rand_acq)s Scatter_out_comptonEvSc.lm > %(rand_combined)s " % fmt_params
            rand_corr_cmd = "%(binPath)s/CorrectAttenuation_SV %(mesh_name)s %(rand_combined)s"\
                            " %(rand_combined_corr)s LM" % fmt_params

            random_cmd = "nice -n19 %(binPath)s/Random_Gen %(distance)d %(voxelSize)f %(norm_prefix)s "\
                         "%(rand_combined_corr)s %(rand_name)s" % fmt_params
            
            recon_cmd = "nice -n 19 %(binPath)s/ClearPEM_LMRec -i %(lmf_name_corr)s -o %(output_prefix)s " \
                        " -d %(distance)f -n %(norm_prefix)s --pixel-length %(voxelSize)f " \
                        " --rays %(rays)d " \
                        " --iterations %(iterations)d --osem %(osem)d " \
                        " --filter-fwhm %(ffwhm)f --filter-power %(fmp)d" \
                        " --rand-image %(rand_name)s" % fmt_params

            print ">>> Calculating mesh..."
            os.system(mesh_cmd)
            
            print ">>> Attenuation correction (true events)..."
            #cat2 = pexpect.spawn("sh -c '%s'" % cat_cmd)
            #elm2todkfz2 = pexpect.spawn(elm2todkfz_cmd, timeout=1E6)
            #true_corr2 = pexpect.spawn(true_corr_cmd, timeout=1E6)
            #elm2todkfz2.expect(pexpect.EOF)
            #true_corr2.expect(pexpect.EOF)
            
            os.system(true_corr_cmd)
            print ">>> Scatter correction..."
            os.system(scatter_cmd)
            print ">>> Attenuation correction (random events)..."
            os.system(random_events_cat)
            os.system(rand_corr_cmd)
            print ">>> Random correction..."
            os.system(random_cmd)
            print ">>> Final reconstruction..."
            #os.system(recon_cmd)
            #cat = pexpect.spawn("sh -c '%s'" % cat_cmd)
            #elm2todkfz = pexpect.spawn(elm2todkfz_cmd, timeout=1E6)
            recon = pexpect.spawn(recon_cmd, timeout=1E6)

            n = parameters['iterations']
            t0 = datetime.now()
            for i in range(n):
                recon.expect('Finished iteration')
                if options.batch is False:
                    progress = 100.0 * (i + 1) / n
                    delta = datetime.now() - t0
                    p = int(1000 * progress)
                    eta = (100000 - p) * delta / p
                    seconds = eta.days * 24 * 3600 + eta.seconds
                    minutes = seconds / 60
                    seconds = seconds % 60
                    sys.stdout.write("Reconstruction progress: %4.1f%% ETA: %02d:%02d\r" % (progress, minutes, seconds))
                    sys.stdout.flush()

            #elm2todkfz.expect(pexpect.EOF)
            recon.expect(pexpect.EOF)

        descString = "DKFZ; %(voxelSize)1.1f mm voxel; %(iterations)d iterations; %(eMin)d-%(eMax)d keV; %(dMax)1.1f ns; filter %(ffwhm)2.1f FFWM %(fmp)d MP; %(rays)d rays" % fmt_params

        hv_name = '%s_%03d.hv' % (output_prefix, parameters['iterations'] - 1)
        dvName = '%s_%03d.v' % (output_prefix, parameters['iterations'] - 1)


    # Segmentate image!
    # hv_name

    # print hv_name, dvName

    if options.preview is False:
        pp = exam.pem().addProcessing();
    else:
        pp = exam.pem().getPreview()

    pp.storeImage(hv_name, dvName)
    pp.setDescription(descString)
    hv_name, dvName = pp.getImageName()
    if options.batch is False:
        print hv_name
    else:
        if options.preview is False:
            print "%05d" % pp.getID()

    os.system("rm -rf %s" % tmpDir)
    return None


if __name__ == '__main__':
    main(sys.argv)
