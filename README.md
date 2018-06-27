# recoGATE

--------------------
|  SIM TRANSLATE   |
--------------------

execute ./deploy.sh from code translate folder

the star.translate script can be run from any folder (with full path to the build folder of course) with input from any folder

example:

/path/to/build/start.translate --input /path/to/data/data.root --output-dir ./

will produce the relevant output in the folder from where you run it

----------------------------
|  COMPILE RECONSTRUCTION  |
----------------------------


==> Step 0: compile reconstruction code

0.1) go to folder reco

cd /path/to/recoGATE/reco

0.2) compile

./deploy 

If the ITK and VTK paths are not set or not correct, the compiling of reconstruction correction will fail. Standard reconstruction will work anyway. Set the ITK and VTK paths in the deploy script.


--------------------
|  COPY WITH GZ    |
--------------------

Concatenate, compress and copy the ouput files .elm2 from somewhere (for example cernbox) in a single .gz file in the location where the reco is performed. Then start the first step of the reco (elmtodkfz) without uncompressing the file:

1) Concatenate, compress and copy:

cat `ls out*_2*` | gzip -c | ssh gatemaster@pcepcmasimu3 "cat > /home/gatemaster/Compton/simulations/cylindricalPET/data_full/full.gz"

this will take all the files matching the 'out*_2*' format (for example out0_2cry.elm2 , out1_2cry.elm2 , etc.), cat them into a pipe sent to gzip, which will compress it and send via ssh to the target PC, where the pipe will be redirected to a single file.

2) start step 1.1a) and 2.1) of RECONSTRUCTION with a pipe from .gz files:

zcat full.gz | ~/Compton/recoGATE/reco/lmrec-trunk/elm2todkfz -o full.lm --energy-low 400 --energy-high 650 --time-window 4 --scatter-total-events scatter.txt -

see point 1.1a) and 2.1) of next paraghaph for details on the command


-------------------
| RECONSTRUCTION  |
-------------------


==> Step 1: produce sensitivity

1.1a) elm2todkfz to transform normalization acq in format dkfz

example for elm2todkfz (data.elm2 is the data of normalization acquisition, scatter.txt a flie produced by the scatter correction, but since it's not used for now, you can just create an empty file):

/path/to/recoGATE/reco/pem-sonic-tools/elm2todkfz -o data.lm --energy-low 400 --energy-high 650 --time-window 4 --scatter-total-events scatter.txt  data.elm2

1.2a) Norm_Total_Gen to produce sensitivity image

/path/to/recoGATE/reco/pem-sonic-tools/Norm_Total_Gen 200 182 2 data.lm norm.lm

OR:

1.1b) or generate a unit image, if you don't want to use sensitivity (image of NxMxL voxels all filled with 1) using:

/path/to/recoGATE/reco/lmrec-trunk/UnitImage 28 402 402

where N = 28, M = 402, L = 402 in this case. The output sile will be unit.lm.


==> Step 2: recostruct image

2.1) elm2todkfz to transform patient acq in format dkfz 

/path/to/recoGATE/reco/pem-sonic-tools/elm2todkfz -o data.image.lm --energy-low 400 --energy-high 650 --time-window 4 --scatter-total-events scatter.txt  unif.elm2

in this example unif.elm2 is the patient acquisition data. scatter.txt is a text file where elm2todkfz write some data.

2.2) ClearPEM_LMRec to produce the image

/path/to/recoGATE/reco/pem-sonic-tools/ClearPEM_LMRec -i data.image.lm -d 200 --axial-size 182  --pixel-length 2 --threads 1 -n norm.lm -o image

where data.image.lm is produced in the step before, -d is the distance between detector/heads or the scanner diameter, pixel-length is the voxel size in mm, threads is the number of threads for calculating the reconstruction, norm.lm is the sensitivity image produced by Norm_Total_Gen.
