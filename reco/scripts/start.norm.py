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

# import libexamdb
# from libstir import writeReconParFile, makeSinogram, getSensitivity, numberOfSegments
# import libliji

def main(argv):

    # basePath = "./"
    # basePath = expanduser(basePath)
    # tmpDir = join(basePath, "tmp")
    # tmpDir = mkdtemp()  #dir=tmpDir
    tmpDir = "./"
    cat_fifo = join(tmpDir, "data.cat.fifo")
    os.mkfifo(cat_fifo)
    lmf_name = join(tmpDir, "data.lmf.fifo")
    os.mkfifo(lmf_name)

    parameters = {'eMin': 450.0, \
                  'eMax': 600.0, \
                  'dMax' : 4.0, \
                  'distance' : 250, \
                  'voxelSize': 2.0, \
                  'norm_prefix' : "/home/marco/cernbox/Universita/Ideas/ComptonRecovery/testPipe/sensitivity.lm",\
                  'rays': 0, \
                  'binPath': "/home/marco/cernbox/Universita/Ideas/ComptonRecovery/recoGATE/reco/lmrec-trunk/bin", \
                  'lmf_name' : lmf_name, \
                  'scatter_total_events_file' :"scatter.txt", \
                  'cat_fifo' : cat_fifo, \
                  'output_prefix' : "output",\
                  'iterations': 8,\
                  'osem': 0,\
                  'ffwhm' : 1.5, \
                  'fmp' : 0, \
                  'inputelm2file' : "norm.gz"}



    cat_cmd = "pigz -d -c %(inputelm2file)s > %(cat_fifo)s" % parameters
    # binPath  = "~/cernbox/Universita/Ideas/ComptonRecovery/recoGATE/reco/pem-sonic-tools/bin"
    # eMin     = 450.0
    # eMax     = 600.0
    # dMax     = 4.0
    # scatter_total_events_file ="scatter.txt"



    elm2todkfz_cmd = "%(binPath)s/elm2todkfz -o %(lmf_name)s " \
                     "--energy-low %(eMin)f --energy-high %(eMax)f " \
                      "--time-window %(dMax)f --scatter-total-events %(scatter_total_events_file)s" \
                      " %(cat_fifo)s" % parameters

    recon_cmd = "nice -n 19 %(binPath)s/ClearPEM_LMRec -i %(lmf_name)s -o %(output_prefix)s " \
                " -d %(distance)f -n %(norm_prefix)s --pixel-length %(voxelSize)f " \
                " --rays %(rays)d " \
                " --iterations %(iterations)d --osem %(osem)d " \
                " --filter-fwhm %(ffwhm)f --filter-power %(fmp)d" % parameters


    print cat_cmd
    print elm2todkfz_cmd

    cat1 = pexpect.spawn("sh -c '%s'" % cat_cmd)
    elm2todkfz = pexpect.spawn(elm2todkfz_cmd, timeout=1E6)
    elm2todkfz.expect(pexpect.EOF)

    # recon = pexpect.spawn(recon_cmd, timeout=1E6, logfile=sys.stdout)


    return None

# def main(argv):
#
#     basePath = "./"
#     config = SafeConfigParser()
#
#     parser = OptionParser()
#     # parser.add_option("--base", dest="base", default=basePath)
#     parser.add_option("--input", dest="elm2Name")
#     parser.add_option("--binPath", dest="binPath")
#     # parser.add_option("--patient", dest="patient", default=1)
#     # parser.add_option("--exam", dest="exam")
#     # parser.add_option("--preview", dest="preview", default=False, action="store_true")
#     parser.add_option("--distance", dest="distance", default=250)
#     parser.add_option("--energy-low", dest="eMin", default=450)
#     parser.add_option("--energy-high", dest="eMax", default=600)
#     parser.add_option("--time-window", dest="dMax", default=4)
#     parser.add_option("--voxel-size", dest="voxelSize", default=2)
#     parser.add_option("--iterations", dest="iterations", default=8)
#     parser.add_option("--no-norm", dest="noNorm", default=False, action="store_true")
#     # parser.add_option("--algo", dest="algo", default="dkfz")
#     parser.add_option("--filter-fwhm", dest="ffwhm", default=1.5)
#     parser.add_option("--filter-power", dest="fmp", default=0)
#     # parser.add_option("--no-osem", dest="osem", default=True, action="store_false")
#     # parser.add_option("--dw-correction", dest="dwCorrection", default=False, action="store_true")
#     parser.add_option("--rays", dest="rays", default=0)
#     # parser.add_option("--corrections", dest="corrections", default=False)
#     # parser.add_option("--segmentation", dest="segmentation", default="auto")
#     # parser.add_option("--segmentation-threshold", dest="segmentationThreshold", default=1)
#
#     (options, argv) = parser.parse_args(argv)
#
#     parameters = {'eMin': float(options.eMin), 'eMax': float(options.eMax), 'dMax': float(options.dMax),
#                   'voxelSize': float(options.voxelSize), 'iterations': int(options.iterations),
#                   'noNorm': options.noNorm, 'ffwhm': float(options.ffwhm),
#                   'fmp': float(options.fmp), 'rays': int(options.rays) , 'distance': float(options.distance), 'binPath': options.binPath }
#
#
#     # sys.stdout.write("Patient: %s\n" % patient.getName())
#     # sys.stdout.write('Exam description: %s\n' % exam.getDescription())
#     # sys.stdout.write('Exam time: %s\n' % str(exam.getTime()))
#     sys.stdout.write('Reconstruction parameters: %s\n' % str(parameters))
#     sys.stdout.write("Scaning data...\n")
#     sys.stdout.flush()
#
#     basePath = expanduser(basePath)
#     tmpDir = join(basePath, "tmp")
#     sys.stdout.write(tmpDir)
#     sys.stdout.write("\n")
#     sys.stdout.flush()
#     tmpDir = mkdtemp()  #dir=tmpDir
#
#     cat_fifo = join(tmpDir, "data.cat.fifo")
#     os.mkfifo(cat_fifo)
#     lmf_name = join(tmpDir, "data.lm")
#     output_prefix = join(tmpDir, "data_img")
#     scatter_total_events_file = join(tmpDir, "scatter_total.txt")
#
#     elm2Name = options.elm2Name
#     fmt_params = parameters
#
#
#     if elm2Name[len(elm2Name) - 2:] == 'gz':
#         cat_cmd = "pigz -d -c %(elm2Name)s > %(cat_fifo)s" % locals()
#     else:
#         cat_cmd = "cat %(elm2Name)s > %(cat_fifo)s" % locals()
#
#     elm2todkfz_cmd = "%(binPath)s/elm2todkfz -o %(lmf_name)s " \
#                      "--energy-low %(eMin)f --energy-high %(eMax)f " \
#                      "--time-window %(dMax)f --scatter-total-events %(scatter_total_events_file)s" \
#                      "%(cat_fifo)s" % fmt_params
#     print elm2todkfz_cmd
#     recon_cmd = "nice -n 19 %(binPath)s/ClearPEM_LMRec -i %(lmf_name)s -o %(output_prefix)s " \
#                 " -d %(distance)f -n %(norm_prefix)s --pixel-length %(voxelSize)f " \
#                 " --rays %(rays)d " \
#                 " --iterations %(iterations)d " \
#                 " --filter-fwhm %(ffwhm)f --filter-power %(fmp)d" % fmt_params
#     print recon_cmd
#
#
#
#     # elif parameters['algo'] == 'DKFZ':
#     #
#     #     ### Files
#     #
#     #     cat_fifo = join(tmpDir, "data.cat.fifo")
#     #     os.mkfifo(cat_fifo)
#     #     lmf_name = join(tmpDir, "data.lm")
#     #     #os.mkfifo(lmf_name)
#     #     scatter_total_events_file = join(tmpDir, "scatter_total.txt")
#     #     rand_acq = join(tmpDir, "rand_acq.lm")
#     #     rand_name = join(tmpDir, "rand_data.lm")
#     #     rand_combined = join(tmpDir, "rand_comb.lm")
#     #     lmf_name_corr = join(tmpDir, "data_corr.lm")
#     #     rand_combined_corr = join(tmpDir, "rand_corr.lm")
#     #     mesh_name = "tfinal_mesh.vtk" # mesh_name is hardcoded in MeshCreator
#     #
#     #     output_prefix = join(tmpDir, "data_img")
#     #
#     #     if elm2Name[len(elm2Name) - 2:] == 'gz':
#     #         cat_cmd = "pigz -d -c %(elm2Name)s > %(cat_fifo)s" % locals()
#     #     else:
#     #         cat_cmd = "cat %(elm2Name)s > %(cat_fifo)s" % locals()
#     #
#     #     norm_prefix = libliji.getNormalization(exam.getTime(), distance, angleList, parameters, basePath, binPath,
#     #                                           tmpDir, options.batch)
#     #
#     #     fmt_params = parameters
#     #     if fmt_params['osem']:
#     #         fmt_params['osem'] = 1
#     #     else:
#     #         fmt_params['osem'] = 0
#     #
#     #     fmt_params.update(locals())
#     #
#     #     elm2todkfz_cmd = "%(binPath)s/elm2todkfz -o %(lmf_name)s --randoms-file %(rand_acq)s " \
#     #                      "--energy-low %(eMin)f --energy-high %(eMax)f " \
#     #                      "--time-window %(dMax)f --scatter-total-events %(scatter_total_events_file)s " \
#     #                      "%(cat_fifo)s" % fmt_params
#     #
#     #     print elm2todkfz_cmd
#     #
#     #     if fmt_params['dwCorrection']:
#     #         elm2todkfz_cmd += " --dw-correction true"
#     #
#     #     recon_cmd = "nice -n 19 %(binPath)s/ClearPEM_LMRec -i %(lmf_name)s -o %(output_prefix)s " \
#     #                 " -d %(distance)f -n %(norm_prefix)s --pixel-length %(voxelSize)f " \
#     #                 " --rays %(rays)d " \
#     #                 " --iterations %(iterations)d --osem %(osem)d " \
#     #                 " --filter-fwhm %(ffwhm)f --filter-power %(fmp)d" % fmt_params
#     #
#     #     # print fmt_params ## DEBUG info
#     #     hv_name = '%s_%03d.hv' % (output_prefix, parameters['iterations'] - 1)
#     #
#     #     print ">>> Reading data..."
#     #     # Extract data
#     #     cat1 = pexpect.spawn("sh -c '%s'" % cat_cmd)
#     #     elm2todkfz = pexpect.spawn(elm2todkfz_cmd, timeout=1E6)
#     #     elm2todkfz.expect(pexpect.EOF)
#     #
#     #     recon = pexpect.spawn(recon_cmd, timeout=1E6, logfile=sys.stdout)
#     #
#     #     # elm2todkfz.expect(pexpect.EOF)
#     #
#     #     print ">>> Basic reconstruction..."
#     #     n = parameters['iterations']
#     #     t0 = datetime.now()
#     #     for i in range(n):
#     #         recon.expect('Finished iteration')
#     #         if options.batch is False:
#     #             progress = 100.0 * (i + 1) / n
#     #             delta = datetime.now() - t0
#     #             p = int(1000 * progress)
#     #             eta = (100000 - p) * delta / p
#     #             seconds = eta.days * 24 * 3600 + eta.seconds
#     #             minutes = seconds / 60
#     #             seconds %= 60
#     #             sys.stdout.write("Reconstruction progress: %4.1f%% ETA: %02d:%02d\r" % (progress, minutes, seconds))
#     #             sys.stdout.flush()
#     #
#     #     recon.expect(pexpect.EOF)
#
#     os.system("rm -rf %s" % tmpDir)
#     return None


if __name__ == '__main__':
    main(sys.argv)
