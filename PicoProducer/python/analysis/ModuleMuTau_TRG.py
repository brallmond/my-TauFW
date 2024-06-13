# Author: Izaak Neutelings (June 2020)
# Description: Simple module to pre-select mutau events
import sys
import numpy as np
from TauFW.PicoProducer.analysis.TreeProducerMuTau_TRG import *
from TauFW.PicoProducer.analysis.ModuleTauPair import *
from TauFW.PicoProducer.analysis.utils import LeptonTauPair, loosestIso, idIso, matchgenvistau, matchtaujet, filtermutau
from TauFW.PicoProducer.corrections.MuonSFs import *
#from TauFW.PicoProducer.corrections.TrigObjMatcher import loadTriggerDataFromJSON, TrigObjMatcher
#from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool, TauESTool, campaigns


class ModuleMuTau_TRG(ModuleTauPair):
  
  def __init__(self, fname, **kwargs):
    kwargs['channel'] = 'mutau'
    super(ModuleMuTau_TRG,self).__init__(fname,**kwargs)
    self.out = TreeProducerMuTau_TRG(fname,self)
    
    # TRIGGERS -- for Tau TRG capture the mutau cross-triggers also
    if self.year == 2024:
      self.trigger    = lambda e: e.HLT_IsoMu24 or e.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 or e.HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1
    else:
      self.trigger    = lambda e: e.HLT_IsoMu24 or e.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
    self.muonCutPt  = lambda e: 24
    self.muonCutEta = lambda e: 2.4
    self.tauCutPt     = 20
    self.tauCutEta    = 2.5
    
    # CORRECTIONS
    if self.ismc:
      self.muSFs      = MuonSFs(era=self.era,verb=self.verbosity) # muon id/iso/trigger SFs
      #self.tesTool    = TauESTool(tauSFVersion[self.year]) # real tau energy scale corrections
      #self.fesTool    = TauFESTool(tauSFVersion[self.year]) # e -> tau fake negligible
      #self.tauSFsT    = TauIDSFTool(tauSFVersion[self.year],'DeepTau2018v2p5VSjet','Tight')
      #self.tauSFsM    = TauIDSFTool(tauSFVersion[self.year],'DeepTau2018v2p5VSjet','Medium')
      #self.tauSFsT_dm = TauIDSFTool(tauSFVersion[self.year],'DeepTau2018v2p5VSjet','Tight', dm=True)
      #self.etfSFs     = TauIDSFTool(tauSFVersion[self.year],'DeepTau2018v2p5VSe',  'VLoose')
      #self.mtfSFs     = TauIDSFTool(tauSFVersion[self.year],'DeepTau2018v2p5VSmu', 'Tight')
    
    # CUTFLOW
    self.out.cutflow.addcut('none',         "no cut"                     )
    self.out.cutflow.addcut('trig',         "trigger"                    )
    self.out.cutflow.addcut('muon',         "muon"                       )
    self.out.cutflow.addcut('tau',          "tau"                        )
    self.out.cutflow.addcut('pair',         "pair"                       )
    self.out.cutflow.addcut('can_trig_match',         "can trig match"   )
    self.out.cutflow.addcut('weight',       "no cut, weighted", 15       )
    self.out.cutflow.addcut('weight_no0PU', "no cut, weighted, PU>0", 16 ) # use for normalization; bug in pre-UL 2018 caused small fraction of events with nPU<=0
    # Important cutflow entries to make stitching with exclusive mutauh sample (DYJetsToTauTauToMuTauh)
    self.out.cutflow.addcut('weight_mutaufilter', "no cut, mutaufilter", 17 )    
    self.out.cutflow.addcut('weight_mutaufilter_NUP0orp4', "no cut, weighted, mutau, 0 or >4 jets", 18 )
    self.out.cutflow.addcut('weight_mutaufilter_NUP1', "no cut, weighted, mutau, 1 jet", 19 )
    self.out.cutflow.addcut('weight_mutaufilter_NUP2', "no cut, weighted, mutau, 2 jets", 20 )
    self.out.cutflow.addcut('weight_mutaufilter_NUP3', "no cut, weighted, mutau, 3 jets", 21 )
    self.out.cutflow.addcut('weight_mutaufilter_NUP4', "no cut, weighted, mutau, 4 jets", 22 )
    
  
  def beginJob(self):
    """Before processing any events or files."""
    super(ModuleMuTau_TRG,self).beginJob()
    print(">>> %-12s = %s"%('tauwp',      self.tauwp))
    print(">>> %-12s = %s"%('muonCutPt',  self.muonCutPt))
    print(">>> %-12s = %s"%('muonCutEta', self.muonCutEta))
    print(">>> %-12s = %s"%('tauCutPt',   self.tauCutPt))
    print(">>> %-12s = %s"%('tauCutEta',  self.tauCutEta))
    print(">>> %-12s = %s"%('trigger',    self.trigger))
    pass
    
  
  def analyze(self, event):
    """Process and pre-select events; fill branches and return True if the events passes,
    return False otherwise."""
    sys.stdout.flush()
    
    
    ##### NO CUT #####################################
    if not self.fillhists(event):
      return False

    
    ##### TRIGGER ####################################
    if not self.trigger(event):
      return False
    self.out.cutflow.fill('trig')
    
    
    ##### MUON #######################################
    muons = [ ]
    for muon in Collection(event,'Muon'):
      if muon.pt<self.muonCutPt(event): continue
      if abs(muon.eta)>self.muonCutEta(event): continue
      if abs(muon.dz)>0.2: continue
      if abs(muon.dxy)>0.045: continue
      if not muon.mediumId: continue
      if muon.pfRelIso04_all>0.50: continue
      muons.append(muon)
    if len(muons)==0:
      return False
    self.out.cutflow.fill('muon')
    
    
    ##### TAU ########################################
    taus = [ ]
    for tau in Collection(event,'Tau'):
      if abs(tau.eta)>self.tauCutEta: continue
      if abs(tau.dz)>0.2: continue
      if tau.decayMode not in [0,1,10,11]: continue
      if abs(tau.charge)!=1: continue
      if tau.idDeepTau2018v2p5VSe<1: continue # VVVLoose
      if tau.idDeepTau2018v2p5VSmu<1: continue # VLoose
      if tau.idDeepTau2018v2p5VSjet<self.tauwp: continue
      if self.ismc:
        tau.es   = 1 # store energy scale for propagating to MET
        genmatch = tau.genPartFlav
        if genmatch==5: # real tau
          if self.tes!=None: # user-defined energy scale (for TES studies)
            tes = self.tes
          else: # recommended energy scale (apply by default)
            tes = self.tesTool.getTES(tau.pt,tau.decayMode,unc=self.tessys)
          if tes!=1:
            tau.pt   *= tes
            tau.mass *= tes
            tau.es    = tes # store for later reuse
        elif self.ltf and 0<genmatch<5: # lepton -> tau fake
          tau.pt   *= self.ltf
          tau.mass *= self.ltf
          tau.es    = self.ltf # store for later reuse
        #elif genmatch in [1,3]: # electron -> tau fake (apply by default, override with 'ltf=1.0')
        #  fes = self.fesTool.getFES(tau.eta,tau.decayMode,unc=self.fes)
        #  tau.pt   *= fes
        #  tau.mass *= fes
        #  tau.es    = fes
        elif self.jtf!=1.0 and genmatch==0: # jet -> tau fake
          tau.pt   *= self.jtf
          tau.mass *= self.jtf
          tau.es    = self.jtf
      if tau.pt<self.tauCutPt: continue
      taus.append(tau)
    if len(taus)==0:
      return False
    self.out.cutflow.fill('tau')
    
    
    ##### MUTAU PAIR #################################
    ltaus = [ ]
    for muon in muons:
      for tau in taus:
        if tau.DeltaR(muon)<0.5: continue
        ltau = LeptonTauPair(muon,muon.pfRelIso04_all,tau,tau.rawDeepTau2018v2p5VSjet)
        ltaus.append(ltau)
    if len(ltaus)==0:
      return False
    muon, tau = max(ltaus).pair
    muon.tlv  = muon.p4()
    tau.tlv   = tau.p4()
    genmatch  = -1 if self.isdata else tau.genPartFlav
    self.out.cutflow.fill('pair')

    # VETOES
    extramuon_veto, extraelec_veto, dilepton_veto = getlepvetoes(event,[ ],[muon],[tau],self.channel, self.era)
    self.out.extramuon_veto[0], self.out.extraelec_veto[0], self.out.dilepton_veto[0] = getlepvetoes(event,[ ],[muon],[ ],self.channel,self.era)
    self.out.lepton_vetoes[0]       = self.out.extramuon_veto[0] or self.out.extraelec_veto[0] or self.out.dilepton_veto[0]
    self.out.lepton_vetoes_notau[0] = extramuon_veto or extraelec_veto or dilepton_veto
    
    # TIGHTEN PRE-SELECTION
    if self.dotight: # do not save all events to reduce disk space
      fail = (self.out.lepton_vetoes[0] and self.out.lepton_vetoes_notau[0]) or\
             tau.idDeepTau2018v2p5VSjet<1 or tau.idDeepTau2018v2p5VSmu<2 or tau.idDeepTau2018v2p5VSe<1
      if (self.tes not in [1,None] or self.tessys!=None) and (fail or tau.genPartFlav!=5):
        return False
      if (self.ltf not in [1,None] or self.fes!=None) and (tau.genPartFlav<1 or tau.genPartFlav>4):
        return False
      ###if self.jtf!=1 and tau.genPartFlav!=0:
      ###  return False

    ##### TRIGGER MATCHED ###########################
    keep_ids = [13, 15]
    trigobjs = [ ]
    trig_match_muon = False
    trig_match_tau  = False
    for i,trigobj in enumerate(Collection(event,'TrigObj')):
      if abs(trigobj.id not in keep_ids): continue # keep muons and taus
      if trigobj.DeltaR(muon) < 0.5: 
        trig_match_muon = True
        trig_idx = i
      if trigobj.DeltaR(tau) < 0.5:  
        trig_match_tau = True
        trig_idx = i
      trigobjs.append(trigobj)
    self.out.trig_match_1[0] = trig_match_muon
    self.out.trig_match_2[0] = trig_match_tau
    if len(trigobjs)==0:
      return False

    # now for all trigger matched objects, save their filterbit, l1pt, l1iso, pt, eta, phi
    
    # EVENT
    self.fillEventBranches(event)
   
    # TRIGGER
    # L1
    #self.out.trig_L1_lowest_singlemuon[0] = e.
    #self.out.trig_L1_lowest_singleele[0]  = e.
    #self.out.trig_L1_lowest_singletau[0]  = e.
    #self.out.trig_L1_lowest_ditau[0]      = e.
    #self.out.trig_L1_lowest_ditaujet[0]   = e.
    #self.out.trig_L1_lowest_mutau[0]      = e.
    #self.out.trig_L1_lowest_etau[0]       = e.
    #self.out.trig_L1_lowest_VBF[0]        = e.
    #self.out.trig_L1_lowest_VBFTau[0]     = e.
    # SingleLepton
    self.out.trig_lowest_singlemuon[0]    = lambda e: e.HLT_IsoMu24
    self.out.trig_lowest_singleele[0]     = lambda e: e.HLT_Ele32_WPTight_Gsf # is it actually 30? 
    # SingleTau
    self.out.trig_singletau_deeptau[0]     = lambda e: e.HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1
    self.out.trig_singletau_deeptau_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1
    # DiTau
    self.out.trig_ditau_deeptau[0]     = lambda e: e.HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1
    self.out.trig_ditau_deeptau_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1
    # DiTau + Jet
    self.out.trig_ditaujet_deeptau[0]            = lambda e: e.HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60
    self.out.trig_ditaujet_deeptau_backup[0]     = lambda e: e.HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75
    self.out.trig_ditaujet_deeptau_mon[0]        = lambda e: e.HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1
    # MuTau
    self.out.trig_mutau_deeptau[0]    = lambda e: e.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
    # ETau
    self.out.trig_etau_deeptau[0]     = lambda e: e.HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1
    self.out.trig_etau_deeptau_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1
    # VBF SingleTau
    #self.out.trig_VBFsingletau_deeptau[0]        = e.
    #self.out.trig_VBFsingletau_deeptau_backup[0] = e.
    self.out.trig_VBFsingletau_deeptau_mon[0]    = lambda e: e.HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1
    # VBF DiTau
    self.out.trig_VBFditau_deeptau[0]     = lambda e: e.HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1
    self.out.trig_VBFditau_deeptau_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1
 
    if ("2024" in self.era): # store PNet paths
      # SingleTau
      self.out.trig_singletau_pnet_loose[0]      = lambda e: e.HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3
      self.out.trig_singletau_pnet_medium[0]     = lambda e: e.HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3
      self.out.trig_singletau_pnet_tight[0]      = lambda e: e.HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3
      self.out.trig_singletau_pnet_loose_mon[0]  = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1
      self.out.trig_singletau_pnet_medium_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1
      self.out.trig_singletau_pnet_tight_mon[0]  = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1
      # DiTau
      self.out.trig_ditau_pnet_medium[0]     = lambda e: e.HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3 
      self.out.trig_ditau_pnet_tight[0]      = lambda e: e.HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3 
      self.out.trig_ditau_pnet_medium_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1
      self.out.trig_ditau_pnet_tight_mon[0]  = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1
      # DiTau + Jet
      self.out.trig_ditaujet_pnet[0]            = lambda e: e.HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60
      self.out.trig_ditaujet_pnet_backup[0]     = lambda e: e.HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75
      self.out.trig_ditaujet_pnet_mon[0]        = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1
      # MuTau
      self.out.trig_mutau_pnet_loose[0]      = lambda e: e.HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1
      self.out.trig_mutau_pnet_medium[0]     = lambda e: e.HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1
      self.out.trig_mutau_pnet_tight[0]      = lambda e: e.HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1
      # ETau
      self.out.trig_etau_pnet_loose[0]      = lambda e: e.HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1
      self.out.trig_etau_pnet_medium[0]     = lambda e: e.HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1
      self.out.trig_etau_pnet_tight[0]      = lambda e: e.HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1
      self.out.trig_etau_pnet_loose_mon[0]  = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring
      self.out.trig_etau_pnet_medium_mon[0] = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring
      self.out.trig_etau_pnet_tight_mon[0]  = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring
      # VBF SingleTau
      self.out.trig_VBFsingletau_pnet[0]        = lambda e: e.HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3
      self.out.trig_VBFsingletau_pnet_backup[0] = lambda e: e.HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3
      self.out.trig_VBFsingletau_pnet_mon[0]    = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1
      # VBF DiTau
      self.out.trig_VBFditau_pnet[0]           = lambda e: e.HLT_VBF_DoublePNetTauhPFJet20_eta2p2
      self.out.trig_VBFditau_pnet_mon[0]       = lambda e: e.HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1
 
    # MUON
    self.out.pt_1[0]                       = muon.pt
    self.out.eta_1[0]                      = muon.eta
    self.out.phi_1[0]                      = muon.phi
    self.out.m_1[0]                        = muon.mass
    self.out.y_1[0]                        = muon.tlv.Rapidity()
    self.out.dxy_1[0]                      = muon.dxy
    self.out.dz_1[0]                       = muon.dz
    self.out.q_1[0]                        = muon.charge
    self.out.iso_1[0]                      = muon.pfRelIso04_all # relative isolation
    self.out.tkRelIso_1[0]                 = muon.tkRelIso
    self.out.idMedium_1[0]                 = muon.mediumId
    self.out.idTight_1[0]                  = muon.tightId
    self.out.idHighPt_1[0]                 = muon.highPtId
    
    
    # TAU
    self.out.pt_2[0]                       = tau.pt
    self.out.eta_2[0]                      = tau.eta
    self.out.phi_2[0]                      = tau.phi
    self.out.m_2[0]                        = tau.mass
    self.out.y_2[0]                        = tau.tlv.Rapidity()
    self.out.dxy_2[0]                      = tau.dxy
    self.out.dz_2[0]                       = tau.dz
    self.out.q_2[0]                        = tau.charge
    self.out.dm_2[0]                       = tau.decayMode
    self.out.iso_2[0]                      = tau.rawIso
    self.out.idiso_2[0]                    = idIso(tau) # cut-based tau isolation (rawIso)
    #self.out.rawAntiEle_2[0]               = tau.rawAntiEle # not available anymore in nanoAODv9
    #self.out.rawMVAoldDM2018v2_2[0]        = tau.rawMVAoldDM2018v2
    #self.out.rawMVAnewDM2018v2_2[0]        = tau.rawMVAnewDM2018v2
    self.out.rawDeepTau2018v2p5VSe_2[0]    = tau.rawDeepTau2018v2p5VSe
    self.out.rawDeepTau2018v2p5VSmu_2[0]   = tau.rawDeepTau2018v2p5VSmu
    self.out.rawDeepTau2018v2p5VSjet_2[0]  = tau.rawDeepTau2018v2p5VSjet
    #self.out.idAntiEle_2[0]                = tau.idAntiEle
    #self.out.idAntiMu_2[0]                 = tau.idAntiMu
    self.out.idDecayMode_2[0]              = tau.idDecayMode
    self.out.idDecayModeNewDMs_2[0]        = tau.idDecayModeNewDMs
    #self.out.idMVAoldDM2018v2_2[0]         = tau.idMVAoldDM2018v2
    #self.out.idMVAnewDM2018v2_2[0]         = tau.idMVAnewDM2018v2
    self.out.idDeepTau2018v2p5VSe_2[0]     = tau.idDeepTau2018v2p5VSe
    self.out.idDeepTau2018v2p5VSmu_2[0]    = tau.idDeepTau2018v2p5VSmu
    self.out.idDeepTau2018v2p5VSjet_2[0]   = tau.idDeepTau2018v2p5VSjet
    self.out.chargedIso_2[0]               = tau.chargedIso
    self.out.neutralIso_2[0]               = tau.neutralIso
    self.out.leadTkPtOverTauPt_2[0]        = tau.leadTkPtOverTauPt
    self.out.photonsOutsideSignalCone_2[0] = tau.photonsOutsideSignalCone
    self.out.puCorr_2[0]                   = tau.puCorr
    
    
    # GENERATOR
    if self.ismc:
      self.out.genmatch_1[0]     = muon.genPartFlav
      self.out.genmatch_2[0]     = tau.genPartFlav
      pt, phi, eta, status       = matchgenvistau(event,tau)
      self.out.genvistaupt_2[0]  = pt
      self.out.genvistaueta_2[0] = eta
      self.out.genvistauphi_2[0] = phi
      self.out.gendm_2[0]        = status
      if self.domutau:
        self.out.mutaufilter[0]  = self.ismutau #filtermutau(event) # for stitching DYJetsToTauTauToMuTauh
    
    
    # JETS
    jets, met, njets_vars, met_vars = self.fillJetBranches(event,muon,tau)
    if self.ismc:
      self.out.jpt_match_2[0], self.out.jpt_genmatch_2[0] = matchtaujet(event,tau,self.ismc)
    else:
      self.out.jpt_match_2[0] = matchtaujet(event,tau,self.ismc)[0]
    
    
    # WEIGHTS
    if self.ismc:
      self.fillCommonCorrBranches(event,jets,met,njets_vars,met_vars)
      if muon.pfRelIso04_all<0.50 and tau.idDeepTau2018v2p5VSjet>=2:
        self.btagTool.fillEffMaps(jets,usejec=self.dojec)
      
      # MUON WEIGHTS
      self.out.trigweight[0]          = self.muSFs.getTriggerSF(muon.pt,muon.eta) # assume leading muon was triggered on
      self.out.idisoweight_1[0]       = self.muSFs.getIdIsoSF(muon.pt,muon.eta)
      
      # DEFAULTS
      self.out.idweight_2[0]          = 1.
      self.out.idweight_dm_2[0]       = 1.
      self.out.idweight_medium_2[0]   = 1.
      
      self.out.ltfweight_2[0]         = 1.
      if self.dosys:
        self.out.idweightUp_2[0]      = 1.
        self.out.idweightDown_2[0]    = 1.
        self.out.idweightUp_dm_2[0]   = 1.
        self.out.idweightDown_dm_2[0] = 1.
        self.out.ltfweightUp_2[0]     = 1.
        self.out.ltfweightDown_2[0]   = 1.
      
      # TAU WEIGHTS
      if tau.genPartFlav==5: # real tau
        self.out.idweight_2[0]        = self.tauSFsT.getSFvsPT(tau.pt)
        self.out.idweight_medium_2[0] = self.tauSFsM.getSFvsPT(tau.pt)
        self.out.idweight_dm_2[0]     = self.tauSFsT_dm.getSFvsDM(tau.pt,tau.decayMode)
        if self.dosys:
          self.out.idweightUp_2[0]    = self.tauSFsT.getSFvsPT(tau.pt,unc='Up')
          self.out.idweightDown_2[0]  = self.tauSFsT.getSFvsPT(tau.pt,unc='Down')
          self.out.idweightUp_dm_2[0]   = self.tauSFsT_dm.getSFvsDM(tau.pt,tau.decayMode,unc='Up')
          self.out.idweightDown_dm_2[0] = self.tauSFsT_dm.getSFvsDM(tau.pt,tau.decayMode,unc='Down')
      elif tau.genPartFlav in [1,3]: # muon -> tau fake
        self.out.ltfweight_2[0]       = self.etfSFs.getSFvsEta(tau.eta,tau.genPartFlav)
        if self.dosys:
          self.out.ltfweightUp_2[0]   = self.etfSFs.getSFvsEta(tau.eta,tau.genPartFlav,unc='Up')
          self.out.ltfweightDown_2[0] = self.etfSFs.getSFvsEta(tau.eta,tau.genPartFlav,unc='Down')
      elif tau.genPartFlav in [2,4]: # electron -> tau fake
        self.out.ltfweight_2[0]       = self.mtfSFs.getSFvsEta(tau.eta,tau.genPartFlav)
        if self.dosys:
          self.out.ltfweightUp_2[0]   = self.mtfSFs.getSFvsEta(tau.eta,tau.genPartFlav,unc='Up')
          self.out.ltfweightDown_2[0] = self.mtfSFs.getSFvsEta(tau.eta,tau.genPartFlav,unc='Down')
      self.out.weight[0]              = self.out.genweight[0]*self.out.puweight[0]*self.out.trigweight[0]*self.out.idisoweight_1[0] #*self.out.idisoweight_2[0]
    elif self.isembed:
      ###self.applyCommonEmbdedCorrections(event,jets,jetIds50,met,njets_vars,met_vars)
      self.out.genweight[0]           = event.genWeight
      self.out.trackweight[0]         = 0.975 if tau.decayMode==0 else 1.0247 if tau.decayMode==1 else 0.927 if tau.decayMode==10 else 0.974 if tau.decayMode==11 else 1.0
    
    
    # MET & DILEPTON VARIABLES
    self.fillMETAndDiLeptonBranches(event,muon,tau,met,met_vars)
    
    
    self.out.fill()
    return True
    
