void momentum(){
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({"filelist.txt"}), std::vector<std::string>({"rTree"}));

  auto* eve_header = new AnalysisTree::EventHeader();
  auto* rec_header = new AnalysisTree::EventHeader();
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* trd_tracks = new AnalysisTree::TrackDetector();
  auto* tof_hits = new AnalysisTree::HitDetector();
  auto* vtx_tof_matching = new AnalysisTree::Matching();

  const int sp = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("p"); //simulated
  const int sx = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("x");
  const int sy = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("y");
  const int sz = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("z");
  const int smother_id = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("mother_id");
  const int rp = treeIn->GetConfiguration()->GetBranchConfig("VtxTracks").GetFieldId("p"); //reconstructed
  const int mass2 = treeIn->GetConfiguration()->GetBranchConfig("TofHits").GetFieldId("mass2"); //tofhits
  const int qp_tof = treeIn->GetConfiguration()->GetBranchConfig("TofHits").GetFieldId("qp_tof");
  const int dE_over_dx = treeIn->GetConfiguration()->GetBranchConfig("VtxTracks").GetFieldId("dE_over_dx");
  
  //TRD
  const int trd_p = treeIn->GetConfiguration()->GetBranchConfig("TrdTracks").GetFieldId("p");
  const int trd_eloss0 = treeIn->GetConfiguration()->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_0");
  const int trd_eloss1 = treeIn->GetConfiguration()->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_1");
  const int trd_eloss2 = treeIn->GetConfiguration()->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_2");
  const int trd_eloss3 = treeIn->GetConfiguration()->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_3");


  treeIn->SetBranchAddress("SimParticles.", &sim_tracks);
  treeIn->SetBranchAddress("VtxTracks.", &vtx_tracks);
  treeIn->SetBranchAddress("TrdTracks.", &trd_tracks);
  treeIn->SetBranchAddress("RecEventHeader.", &rec_header);
  treeIn->SetBranchAddress("TofHits.", &tof_hits);
  treeIn->SetBranchAddress("VtxTracks2TofHits.", &vtx_tof_matching);
  int N = 0;

  const int Nevents = treeIn->GetEntries();

  TFile* fileOut = TFile::Open("momentum.root", "recreate");
  TH2F hc_qp_mass2("hc_qp_mass2", "TOF qp mass2; sign(q)*p ;mass^2", 500, -16, 16, 7, -2, 5); //tof hits
  TH2F dEdx_p("dEdx_p", "STS dE/dx p;dE/dx ;p", 400, 0, 20, 400, 0, 100000); //tof hits
  TH2F eloss0_p("eloss0_p", " TRD deltaE_0 p ;eloss0 ;p", 400, -2, 20, 400, -4, 100); 
  TH2F eloss1_p("eloss1_p", " TRD deltaE_1 p ;eloss1 ;p", 400, -2, 20, 400, -4, 100); 
  TH2F eloss2_p("eloss2_p", " TRD deltaE_2 p ;eloss2 ;p", 400, -2, 20, 400, -4, 100); 
  TH2F eloss3_p("eloss3_p", " TRD deltaE_3 p ;eloss3 ;p", 400, -2, 20, 400, -4, 100); 
  TH2F p_mass2("p_mass2", "Matching STS p TOF mass2 ;p ;mass^2", 400, -2, 20, 7, -2, 5); 

  for(int i=0; i<Nevents; i++){
    treeIn -> GetEntry(i);

    //tof hits
    for (const auto& tof_hit : *tof_hits){
      float tof_mass2 = tof_hit.GetField<float>(mass2);
      const float tof_qp_tof = tof_hit.GetField<float>(qp_tof);
      if(tof_mass2 != -1)
      cout << "mass: " << tof_mass2 << '\t';
      
      // cout << "qp: " << tof_qp_tof << '\n';
      hc_qp_mass2.Fill(tof_qp_tof, tof_mass2);
      
      const int tof_id = vtx_tof_matching->GetMatch(tof_hit.GetId()); 
      // cout << "tofid: " << tof_id << '\t';
      // cout << "track id: " << tof_hit.GetId() << '\n';
      if(tof_id<0) continue;
      const float vtx_p = vtx_tracks->GetChannel(tof_id).GetField<float>(rp);
      // cout << "p: " << vtx_p << '\t';
      // cout << "m2: " << tof_mass2 << '\n';

      p_mass2.Fill(vtx_p, tof_mass2);
    }

    for (const auto& track: *vtx_tracks)
    {
      const float vtx_p = track.GetField<float>(rp);
      const float vtx_dEdx = track.GetField<float>(dE_over_dx);
      // cout << "p: " << vtx_p << '\t';
      // cout << "dEdx: " << vtx_dEdx << '\n';

      dEdx_p.Fill(vtx_p, vtx_dEdx);
    }
      
    
    for (const auto& track: *trd_tracks)
    {
      const float _trd_p = track.GetField<float>(trd_p);
      const float _trd_eloss0 = track.GetField<float>(trd_eloss0);
      const float _trd_eloss1 = track.GetField<float>(trd_eloss1);
      const float _trd_eloss2 = track.GetField<float>(trd_eloss2);
      const float _trd_eloss3 = track.GetField<float>(trd_eloss3);
      // cout << "p: " << _trd_p << '\t';
      // cout << "eloss0: " << _trd_eloss0 << '\n';
      // cout << "eloss0: " << _trd_eloss1 << '\n';
      // cout << "eloss0: " << _trd_eloss2 << '\n';
      // cout << "eloss0: " << _trd_eloss3 << '\n';

      eloss0_p.Fill(_trd_p, _trd_eloss0);
      eloss1_p.Fill(_trd_p, _trd_eloss1);
      eloss2_p.Fill(_trd_p, _trd_eloss2);
      eloss3_p.Fill(_trd_p, _trd_eloss3);
    }
   }


   hc_qp_mass2.Write();
   dEdx_p.Write();
   eloss0_p.Write();
   eloss1_p.Write();
   eloss2_p.Write();
   eloss3_p.Write();
   p_mass2.Write();

  fileOut->Close();
}