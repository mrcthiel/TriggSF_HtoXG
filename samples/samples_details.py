from glob import glob
from typing import TypedDict


class Sample(TypedDict):
    files: list[str]
    year: str
    data_or_mc: str


# samples to process
samples: dict[str, Sample] = {
    # MC
    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_2016": {
        "files": glob(
            "/afs/cern.ch/user/m/mthiel/public/TOP-RunIISummer20UL16NanoAODv9-00036_2927.root"
#            "/eos/cms/store/user/mthiel/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/DY16/221026_231315/0002/TOP-RunIISummer20UL16NanoAODv9-00036_2927.root"
        ),
        "year": "2016",
        "data_or_mc": "mc",
    }
}

# OLD
#samples: dict[str, Sample] = {
    # Data
#    "Run2018D_2018": {
#        "files": glob(
#            "/eos/cms/store/user/mthiel/SingleMuon/Run2018D/220808_100446/100*/*.root"
##            "/eos/cms/store/user/mthiel/SingleMuon/Run2018A/220804_234515/000*/ReReco-Run2018A-MuonEG-UL2018_MiniAODv2_NanoAODv9-00001_4*.root"
#        ),
#        "year": "2018",
#       "data_or_mc": "data",
#    },
#    #MC
#    "DY_2016": {
#        "files": glob(
#            "/afs/cern.ch/work/m/mthiel/private/HiggsToUpsilonGamma/analyzer/new_analyzer_hlt_test/HZUpsilonPhotonRun2NanoAOD/DY18.root"
#            "/afs/cern.ch/work/m/mthiel/private/HiggsToUpsilonGamma/CMSSW_10_6_26/src/*.root"
#        ),
#        "year": "2016",
#        "data_or_mc": "mc",
#    },
#}

# build samples files and descriptions
samples_files = {}
for sample in samples:
    samples_files[sample] = samples[sample]["files"]

data_samples_files = {}
for sample in samples:
    if samples[sample]["data_or_mc"] == "data":
        data_samples_files[sample] = samples[sample]["files"]

mc_samples_files = {}
for sample in samples:
    if samples[sample]["data_or_mc"] == "mc":
        mc_samples_files[sample] = samples[sample]["files"]
