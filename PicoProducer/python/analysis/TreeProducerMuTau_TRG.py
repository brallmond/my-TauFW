# Author: Izaak Neutelings (June 2020)
# Extended: Braden Allmond (May 2024)
# Sources:
#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#Synchronisation
#   https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html
from TauFW.PicoProducer.analysis.TreeProducerTauPair import TreeProducerTauPair

# TagAndProbe should just be a subset of TRG, so just make this TRG, and let the TagAndProbe module use it
class TreeProducerMuTau_TRG(TreeProducerTauPair):
  """Class to create and prepare a custom output file & tree."""
  
  def __init__(self, filename, module, **kwargs):
    print("Loading TreeProducerMuTau_TRG for %r"%(filename))
    super(TreeProducerMuTau_TRG,self).__init__(filename,module,**kwargs)
   
    ############
    #  TRIGGER #
    ############
    # 1 is muon, 2 is tau (TauFW convention)
    self.addBranch('trig_match_1', '?')
    self.addBranch('trig_match_2', '?')
    # shouldn't these be higher level?
    #self.addBranch('tag_pass', '?')
    #self.addBranch('tag_fail', '?')
    #self.addBranch('probe_pass', '?')
    #self.addBranch('probe_fail', '?')
    # L1
    self.addBranch('trig_L1_lowest_singlemuon', '?')
    self.addBranch('trig_L1_lowest_singleele', '?')
    self.addBranch('trig_L1_lowest_singletau', '?')
    self.addBranch('trig_L1_lowest_ditau',     '?')
    self.addBranch('trig_L1_lowest_ditaujet',  '?')
    self.addBranch('trig_L1_lowest_mutau',     '?')
    self.addBranch('trig_L1_lowest_etau',      '?')
    self.addBranch('trig_L1_lowest_VBF',       '?')
    self.addBranch('trig_L1_lowest_VBFTau',    '?')
    # SingleLepton
    self.addBranch('trig_lowest_singlemuon', '?')
    self.addBranch('trig_lowest_singleele',  '?')
    # SingleTau
    self.addBranch('trig_singletau_deeptau',     '?')
    self.addBranch('trig_singletau_pnet_loose',  '?')
    self.addBranch('trig_singletau_pnet_medium', '?')
    self.addBranch('trig_singletau_pnet_tight',  '?')
    self.addBranch('trig_singletau_deeptau_mon',     '?')
    self.addBranch('trig_singletau_pnet_loose_mon',  '?')
    self.addBranch('trig_singletau_pnet_medium_mon', '?')
    self.addBranch('trig_singletau_pnet_tight_mon',  '?')
    # DiTau
    self.addBranch('trig_ditau_deeptau',     '?')
    self.addBranch('trig_ditau_pnet_medium', '?')
    self.addBranch('trig_ditau_pnet_tight',  '?')
    self.addBranch('trig_ditau_deeptau_mon',     '?')
    self.addBranch('trig_ditau_pnet_medium_mon', '?')
    self.addBranch('trig_ditau_pnet_tight_mon',  '?')
    # DiTau + Jet
    self.addBranch('trig_ditaujet_deeptau',        '?')
    self.addBranch('trig_ditaujet_deeptau_backup', '?')
    self.addBranch('trig_ditaujet_pnet',           '?')
    self.addBranch('trig_ditaujet_pnet_backup',    '?')
    self.addBranch('trig_ditaujet_deeptau_mon',        '?')
    self.addBranch('trig_ditaujet_pnet_mon',           '?')
    # MuTau # is its own monitoring path
    self.addBranch('trig_mutau_deeptau',     '?')
    self.addBranch('trig_mutau_pnet_loose',  '?')
    self.addBranch('trig_mutau_pnet_medium', '?')
    self.addBranch('trig_mutau_pnet_tight',  '?')
    # ETau
    self.addBranch('trig_etau_deeptau',     '?')
    self.addBranch('trig_etau_pnet_loose',  '?')
    self.addBranch('trig_etau_pnet_medium', '?')
    self.addBranch('trig_etau_pnet_tight',  '?')
    self.addBranch('trig_etau_deeptau_mon',     '?')
    self.addBranch('trig_etau_pnet_loose_mon',  '?')
    self.addBranch('trig_etau_pnet_medium_mon', '?')
    self.addBranch('trig_etau_pnet_tight_mon',  '?')
    # VBF SingleTau
    self.addBranch('trig_VBFsingletau_deeptau',        '?')
    self.addBranch('trig_VBFsingletau_deeptau_backup', '?')
    self.addBranch('trig_VBFsingletau_pnet',           '?')
    self.addBranch('trig_VBFsingletau_pnet_backup',    '?')
    self.addBranch('trig_VBFsingletau_deeptau_mon',    '?')
    self.addBranch('trig_VBFsingletau_pnet_mon',       '?')
    # VBF DiTau
    self.addBranch('trig_VBFditau_deeptau',        '?')
    self.addBranch('trig_VBFditau_pnet',           '?')
    self.addBranch('trig_VBFditau_deeptau_mon',    '?')
    self.addBranch('trig_VBFditau_pnet_mon',       '?')
 

    ############
    #   MUON   #
    ############
    
    self.addBranch('pt_1',       'f')
    self.addBranch('eta_1',      'f')
    self.addBranch('phi_1',      'f')
    self.addBranch('m_1',        'f')
    self.addBranch('y_1',        'f')
    self.addBranch('dxy_1',      'f')
    self.addBranch('dz_1',       'f')
    self.addBranch('q_1',        'i')
    self.addBranch('iso_1',      'f', title="relative isolation, pfRelIso04_all")
    self.addBranch('tkRelIso_1', 'f')
    self.addBranch('idMedium_1', '?')
    self.addBranch('idTight_1',  '?')
    self.addBranch('idHighPt_1', 'i')

    ###########
    #   TAU   #
    ###########
    
    self.addBranch('pt_2',                       'f')
    self.addBranch('eta_2',                      'f')
    self.addBranch('phi_2',                      'f')
    self.addBranch('m_2',                        'f')
    self.addBranch('y_2',                        'f')
    self.addBranch('dxy_2',                      'f')
    self.addBranch('dz_2',                       'f')
    self.addBranch('q_2',                        'i')
    self.addBranch('dm_2',                       'i')
    self.addBranch('iso_2',                      'f', title="rawIso")
    self.addBranch('idiso_2',                    'i', title="rawIso WPs")
    #self.addBranch('rawAntiEle_2',               'f') # not available anymore in nanoAODv9
    #self.addBranch('rawMVAoldDM2017v2_2',        'f')
    #self.addBranch('rawMVAnewDM2017v2_2',        'f')
    self.addBranch('rawDeepTau2017v2p1VSe_2',    'f')
    self.addBranch('rawDeepTau2017v2p1VSmu_2',   'f')
    self.addBranch('rawDeepTau2017v2p1VSjet_2',  'f')

    self.addBranch('rawDeepTau2018v2p5VSe_2',    'f')
    self.addBranch('rawDeepTau2018v2p5VSmu_2',   'f')
    self.addBranch('rawDeepTau2018v2p5VSjet_2',  'f')


    #self.addBranch('idAntiEle_2',                'i')
    #self.addBranch('idAntiMu_2',                 'i')
    self.addBranch('idDecayMode_2',              '?', title="oldDecayModeFinding")
    self.addBranch('idDecayModeNewDMs_2',        '?', title="newDecayModeFinding")
    #self.addBranch('idMVAoldDM2017v2_2',         'i')
    #self.addBranch('idMVAnewDM2017v2_2',         'i')
    self.addBranch('idDeepTau2017v2p1VSe_2',     'i')
    self.addBranch('idDeepTau2017v2p1VSmu_2',    'i')
    self.addBranch('idDeepTau2017v2p1VSjet_2',   'i')

    self.addBranch('idDeepTau2018v2p5VSe_2',     'i')
    self.addBranch('idDeepTau2018v2p5VSmu_2',    'i')
    self.addBranch('idDeepTau2018v2p5VSjet_2',   'i')


    self.addBranch('leadTkPtOverTauPt_2',        'f')
    self.addBranch('chargedIso_2',               'f')
    self.addBranch('neutralIso_2',               'f')
    self.addBranch('photonsOutsideSignalCone_2', 'f')
    self.addBranch('puCorr_2',                   'f')
    self.addBranch('jpt_match_2',                'f', -1, title="pt of jet matching tau")

    if self.module.ismc:
      self.addBranch('jpt_genmatch_2',      'f', -1, title="pt of gen jet matching tau")
      self.addBranch('genmatch_1',          'i', -1)
      self.addBranch('genmatch_2',          'i', -1)
      self.addBranch('genvistaupt_2',       'f', -1)
      self.addBranch('genvistaueta_2',      'f', -9)
      self.addBranch('genvistauphi_2',      'f', -9)
      self.addBranch('gendm_2',             'i', -1)
      self.addBranch('idisoweight_1',       'f', 1., title="muon ID/iso efficiency SF")
      self.addBranch('idweight_2',          'f', 1., title="tau ID efficiency SF, Tight")
      self.addBranch('idweight_dm_2',       'f', 1., title="tau ID efficiency SF, Tight, DM-dependent")
      self.addBranch('idweight_medium_2',   'f', 1., title="tau ID efficiency SF, Medium")
      self.addBranch('ltfweight_2',         'f', 1., title="lepton -> tau fake rate SF")
      if self.module.dosys: # systematic variation (only for nominal tree)
        self.addBranch('idweightUp_2',      'f', 1.)
        self.addBranch('idweightDown_2',    'f', 1.)
        self.addBranch('idweightUp_dm_2',   'f', 1.)
        self.addBranch('idweightDown_dm_2', 'f', 1.)
        self.addBranch('ltfweightUp_2',     'f', 1.)
        self.addBranch('ltfweightDown_2',   'f', 1.)
      if self.module.domutau:
        self.addBranch('mutaufilter',       '?', title="has tautau -> mutau, pT>18, |eta|<2.5")
    
