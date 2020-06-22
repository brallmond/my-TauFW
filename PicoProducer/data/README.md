# Correction Tools
Several tools to get corrections, efficiencies, scale factors (SFs), event weights, etc.

[Pileup reweighting](#pileup-reweighting)  
[Lepton efficiencies](#lepton-efficiencies)
[B tagging tools](#b-tagging-tools)
[Test SFs](test-sfs)


## Pileup reweighting

[`PileupTool.py`](../python/corrections/PileupTool.py) provides the pileup event weight based on the data and MC profiles in [`pileup/`](pileup).

The data profile can be computed with the `pileupCalc.py` tool.
The MC profile can be taken from the distribution of the `Pileup_nTrueInt` variable in nanoAOD, for each MC event:
```
    self.out.pileup.Fill(event.Pileup_nTrueInt)
```
and then extracted with [`pileup/getPileupProfiles.py`](pileup/getPileupProfiles.py). Comparisons are shown [here for 2016](https://ineuteli.web.cern.ch/ineuteli/pileup/2016/), [here for 2017](https://ineuteli.web.cern.ch/ineuteli/pileup/2017/) and [here for 2018](https://ineuteli.web.cern.ch/ineuteli/pileup/2018/).

Please note that some 2017 samples had a buggy pileup module, and need to be treated separately, see [this](https://hypernews.cern.ch/HyperNews/CMS/get/generators/4060.html?inline=-1) or [this HyperNews post](https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3128.html).
`PileupTool.py` manually checks for some of these samples with the `hasBuggyPU` help function.
Also, [`pileup/getPileupProfiles.py`](pileup/getPileupProfiles.py) splits the 2017 MC pileup profiles into those of the buggy (`old_pmx`) and fixed (`new_pmx`) samples.
You can find out if your favorite samples has a buggy pileup profile by passing its DAS path to [`pileup/checkBuggyPileup2017.sh`](pileup/checkBuggyPileup2017.sh), which makes use of `dasgoclient`:
```
./pileup/checkBuggyPileup2017.sh /DY*JetsToLL_M-50_TuneCP5*mad*/RunIIFall17*/NANOAOD*
```


## Lepton efficiencies

Several classes are available to get corrections for electrons, muons and hadronically-decayed tau leptons:

* `ScaleFactorTool.py`
  * `ScaleFactor`: general class to get SFs from histograms
  * `ScaleFactorHTT`: class to get SFs from histograms, as measured by the [HTT group](https://github.com/CMS-HTT/LeptonEfficiencies)
* `MuonSFs.py`: class to get muon trigger / identification / isolation SFs
* `ElectronSFs.py` class to get electron trigger / identification / isolation SFs

`ROOT` files with efficiencies and SFs are saved in [`lepton`](lepton) and [`tau`](tau). 
Scale factors can be found here:
* muon efficiencies and SFs: [Muon POG Run-II Recommendations](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceSelectionAndCalibrationsRun2) ([2016 Legacy](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2016LegacyRereco), [2017](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017), [2018](https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2018))
* electron efficiencies and SFs: [Electron POG](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2) ([2017](https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations))
* tau triggers SFs (ditau, mutau, etau; [2016](https://github.com/rmanzoni/triggerSF/tree/moriond17), [2017](https://github.com/truggles/TauTriggerSFs/tree/final_2017_MCv2))

In case you use lepton scale factors and efficiencies as measured by the HTT group, you need to make sure you get them with
```
cd lepton
git clone https://github.com/CMS-HTT/LeptonEfficiencies HTT
```



## B tagging tools

[`BTagTool.py`](../python/corrections/BTagTool.py) provides two classes: `BTagWPs` for saving the working points (WPs) per year and type of tagger, and `BTagWeightTool` to provide b tagging weights. These can be called during the initialization of you analysis module, e.g.:
```
class ModuleMuTau(Module):
    
    def __init__(self, ... ):
        
        ...
        
        if not self.isData:
          self.btagTool = BTagWeightTool('DeepCSV','medium',channel=channel,year=year)
        self.deepcsv_wp = BTagWPs('DeepCSV',year=year)
        
    
    def analyze(self, event):
        
        nbtag  = 0
        jetIds = [ ]
        for ijet in range(event.nJet):
            ...
            jetIds.append(ijet)
            if event.Jet_btagDeepB[ijet] > self.deepcsv_wp.medium:
              nbtag += 1
        
        if not self.isData:
          self.out.btagweight[0] = self.btagTool.getWeight(event,jetIds)
```

`BTagWeightTool` calculates b tagging reweighting based on the [SFs provided from the BTagging group](https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation#Recommendation_for_13_TeV_Data) and analysis-dependent efficiencies measured in MC. These are saved in `ROOT` files in [`btag/`](btag).
The event weight is calculated according to [this method](https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1a_Event_reweighting_using_scale).

### Computing the b tag efficiencies

The b tag efficiencies are analysis-dependent. They can be computed from the analysis output run on MC samples. For each MC event, fill the numerator and denominator histograms with `BTagWeightTool.fillEfficiencies`, after removing overlap with other selected objects, e.g. the muon and tau object in [`ModuleMuTau.py`](../modules/ModuleMuTau.py):
<pre>
    def analyze(self event):
    
        # select isolated muon and tau
        ...
        
        for ijet in range(event.nJet):
            if event.Jet_pt[ijet] < 30: continue
            if abs(event.Jet_eta[ijet]) > 4.7: continue
            <b>if muon.DeltaR(jets[ijet].p4()) < 0.5: continue
            if tau.DeltaR(jets[ijet].p4()) < 0.5: continue</b>
            jetIds.append(ijet)
        
        if not self.isData:
          self.btagTool.fillEfficiencies(event,jetIds)
        
        ...
</pre>
Do this for as many MC samples as possible, to gain as much statistics as possible (also note that jets in Drell-Yan, W+jets and ttbar events typically have different jet flavor content). Then use [`btag/getBTagEfficiencies.py`](btag/getBTagEfficiencies.py) to extract all histograms from analysis output, add them together for maximum statistics, and compute the efficiencies. (You should edit this script to read in your analysis output.)
Examples of efficiency maps per jet flavor, and as a function of jet pT versus jet eta for the mutau analysis in 2017 are shown [here](https://ineuteli.web.cern.ch/ineuteli/btag/2017/?match=mutau).





## Test SFs

`testSFs.py` provides a simple and direct way of testing the correction tool classes, without running the whole framework.

