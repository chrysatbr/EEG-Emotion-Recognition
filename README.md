# EEG-Emotion-Recognition
Emotion recognition through EEG by using HOS method

# Directory structure

```
lib
└── hosa
pdfs
├── DEAP_Histo_Coupling_Frequencies.pdf
├── DEAP_Participant_1.pdf
├── DEAP_Participant_3.pdf
├── DEAP_Participant_4.pdf
└── DEAP_Participant_5.pdf
utils
├── loadAndDecomposeDEAP.m
├── loadAndDecomposeSEED.m
└── maxMatrix.m
simulationDEAP.m
simulationSEED.m
```
- In the `lib` directory, there is the HOSA matlab toolbox that is slightly modified by the original to fit our plotting requirements.
- In the `utils` directory, there is a list of helping functions that we use in our simulation profiles.
- In the `pdfs` directory, there is a list with the output of the bulk visualization process of our analysis using bispectrum. The files `DEAP_Participant_<id>.pdf` contain all the bispectrum plots per participant for all the channels and videos. For the DEAP dataset we have 32 channels and 40 videos, so we created a canvas of 40x32 bispectrum plots. In the title of each plot, there is the information of the peak location. Observing the aforementioned canvas for some of the participants we conclude that there is a sharp peak (with some exceptions/outliers). So later on we tried to observe the frequency pairs that these peaks are located. This information is depicted with histograms in the file `DEAP_Histo_Coupling_Frequencies.pdf`.

# Set up

- In the simulations files, the local path of the datasets should be specified.
- The simulations are divided to code sections in order to give us flexibility on what we want to plot.