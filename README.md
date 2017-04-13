# recoGATE


execute ./deploy.sh from code translate folder

the star.translate script can be run from any folder (with full path to the build folder of course) with input from any folder

example:

/path/to/build/start.translate --input /path/to/data/data.root --output-dir ./

will produce the relevant output in the folder from where you run it



RECO

Step 1: produce sensitivity

- elm2todkfz to transform normalization acq in format dkfz

example for elm2todkfz:

/home/marco/Universita/Ideas/ComptonRecovery/recoGATE/reco/pem-sonic-tools/elm2todkfz -o data.lm --energy-low 400 --energy-high 650 --time-window 4 --scatter-total-events scatter.txt  data.elm2

- Norm_Total_Gen to produce sensitivity image

/home/marco/Universita/Ideas/ComptonRecovery/recoGATE/reco/pem-sonic-tools/Norm_Total_Gen 200 182 2 data.lm norm.lm


Step 2: recostruct image

- elm2todkfz to transform patient acq in format dkfz

home/marco/Universita/Ideas/ComptonRecovery/recoGATE/reco/pem-sonic-tools/elm2todkfz -o data.image.lm --energy-low 400 --energy-high 650 --time-window 4 --scatter-total-events scatter.txt  unif.elm2

- ClearPEM_LMRec to produce the image

/home/marco/Universita/Ideas/ComptonRecovery/recoGATE/reco/pem-sonic-tools/ClearPEM_LMRec -i data.image.lm -d 200 --axial-size 182  --pixel-length 2 --threads 1 -n norm.lm -o image
