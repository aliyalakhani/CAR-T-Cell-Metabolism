function mol=CARTcell_ccm(free_net,free_xch,inp)

mol.EMU=1;

persistent inC
if isempty(inC)
	inC.CO2_IN_1=emufy(inp.CO2_IN__C,[1]);
	inC.AC_IN_2=emufy(inp.AC_IN__C,[2]);
	inC.OAA_IN_2=emufy(inp.OAA_IN__C,[2]);
	inC.AC_IN_1=emufy(inp.AC_IN__C,[1]);
	inC.OAA_IN_4=emufy(inp.OAA_IN__C,[4]);
	inC.OAA_IN_3=emufy(inp.OAA_IN__C,[3]);
	inC.Glu_IN_3=emufy(inp.Glu_IN__C,[3]);
	inC.OAA_IN_1=emufy(inp.OAA_IN__C,[1]);
	inC.Glu_IN_1=emufy(inp.Glu_IN__C,[1]);
	inC.GLC_IN_3=emufy(inp.GLC_IN__C,[3]);
	inC.GLC_IN_2=emufy(inp.GLC_IN__C,[2]);
	inC.Gln_IN_3=emufy(inp.Gln_IN__C,[3]);
	inC.Glu_IN_2=emufy(inp.Glu_IN__C,[2]);
	inC.GLC_IN_1=emufy(inp.GLC_IN__C,[1]);
	inC.Gln_IN_1=emufy(inp.Gln_IN__C,[1]);
	inC.Glu_IN_4=emufy(inp.Glu_IN__C,[4]);
	inC.Glu_IN_5=emufy(inp.Glu_IN__C,[5]);
	inC.Gln_IN_2=emufy(inp.Gln_IN__C,[2]);
	inC.GLC_IN_4=emufy(inp.GLC_IN__C,[4]);
	inC.Gln_IN_4=emufy(inp.Gln_IN__C,[4]);
	inC.Gln_IN_5=emufy(inp.Gln_IN__C,[5]);
	inC.GLC_IN_5=emufy(inp.GLC_IN__C,[5]);
	inC.GLC_IN_6=emufy(inp.GLC_IN__C,[6]);
	inC.AC_IN_1_2=emufy(inp.AC_IN__C,[1 2]);
	inC.OAA_IN_2_3=emufy(inp.OAA_IN__C,[2 3]);
	inC.Glu_IN_2_3=emufy(inp.Glu_IN__C,[2 3]);
	inC.GLC_IN_1_2=emufy(inp.GLC_IN__C,[1 2]);
	inC.Gln_IN_2_3=emufy(inp.Gln_IN__C,[2 3]);
	inC.Glu_IN_4_5=emufy(inp.Glu_IN__C,[4 5]);
	inC.GLC_IN_2_3=emufy(inp.GLC_IN__C,[2 3]);
	inC.Gln_IN_4_5=emufy(inp.Gln_IN__C,[4 5]);
	inC.GLC_IN_4_5=emufy(inp.GLC_IN__C,[4 5]);
	inC.OAA_IN_1_2=emufy(inp.OAA_IN__C,[1 2]);
	inC.OAA_IN_3_4=emufy(inp.OAA_IN__C,[3 4]);
	inC.GLC_IN_5_6=emufy(inp.GLC_IN__C,[5 6]);
	inC.Glu_IN_3_4=emufy(inp.Glu_IN__C,[3 4]);
	inC.Gln_IN_3_4=emufy(inp.Gln_IN__C,[3 4]);
	inC.Glu_IN_1_2=emufy(inp.Glu_IN__C,[1 2]);
	inC.Gln_IN_1_2=emufy(inp.Gln_IN__C,[1 2]);
	inC.OAA_IN_2_3_4=emufy(inp.OAA_IN__C,[2 3 4]);
	inC.OAA_IN_1_2_3=emufy(inp.OAA_IN__C,[1 2 3]);
	inC.Glu_IN_1_2_3=emufy(inp.Glu_IN__C,[1 2 3]);
	inC.GLC_IN_1_2_3=emufy(inp.GLC_IN__C,[1 2 3]);
	inC.Gln_IN_1_2_3=emufy(inp.Gln_IN__C,[1 2 3]);
	inC.GLC_IN_4_5_6=emufy(inp.GLC_IN__C,[4 5 6]);
	inC.Glu_IN_2_3_4=emufy(inp.Glu_IN__C,[2 3 4]);
	inC.Glu_IN_3_4_5=emufy(inp.Glu_IN__C,[3 4 5]);
	inC.Gln_IN_2_3_4=emufy(inp.Gln_IN__C,[2 3 4]);
	inC.Gln_IN_3_4_5=emufy(inp.Gln_IN__C,[3 4 5]);
	inC.OAA_IN_1_2_3_4=emufy(inp.OAA_IN__C,[1 2 3 4]);
	inC.Glu_IN_2_3_4_5=emufy(inp.Glu_IN__C,[2 3 4 5]);
	inC.GLC_IN_3_4_5_6=emufy(inp.GLC_IN__C,[3 4 5 6]);
	inC.Gln_IN_2_3_4_5=emufy(inp.Gln_IN__C,[2 3 4 5]);
	inC.Glu_IN_1_2_3_4_5=emufy(inp.Glu_IN__C,[1 2 3 4 5]);
	inC.Gln_IN_1_2_3_4_5=emufy(inp.Gln_IN__C,[1 2 3 4 5]);
	inC.GLC_IN_2_3_4_5_6=emufy(inp.GLC_IN__C,[2 3 4 5 6]);
	inC.GLC_IN_1_2_3_4_5_6=emufy(inp.GLC_IN__C,[1 2 3 4 5 6]);
end
mol_=CARTcell_ccm_C(free_net,free_xch,inC);
mol.CO2__C=mol_.CO2;
mol.GLC__C=mol_.GLC;
mol.G6P__C=mol_.G6P;
mol.F6P__C=mol_.F6P;
mol.FBP__C=mol_.FBP;
mol.GAP__C=mol_.GAP;
mol.DHAP__C=mol_.DHAP;
mol.BPG__C=mol_.BPG;
mol.PGA__C=mol_.PGA;
mol.PEP__C=mol_.PEP;
mol.PYR__C=mol_.PYR;
mol.LAC__C=mol_.LAC;
mol.OAA__C=mol_.OAA;
mol.MAL__C=mol_.MAL;
mol.m6PG__C=mol_.m6PG;
mol.Ru5P__C=mol_.Ru5P;
mol.R5P__C=mol_.R5P;
mol.X5P__C=mol_.X5P;
mol.S7P__C=mol_.S7P;
mol.E4P__C=mol_.E4P;
mol.AcCoA__C=mol_.AcCoA;
mol.CitICit__C=mol_.CitICit;
mol.OGA__C=mol_.OGA;
mol.SuccCoA__C=mol_.SuccCoA;
mol.Succ__C=mol_.Succ;
mol.Fum__C=mol_.Fum;
mol.Glu__C=mol_.Glu;
mol.Ala__C=mol_.Ala;
mol.Gln__C=mol_.Gln;

persistent inC2
if isempty(inC2)
	inC2.CO2_IN_1=emufy(inp.CO2_IN__C2,[1]);
	inC2.AC_IN_2=emufy(inp.AC_IN__C2,[2]);
	inC2.OAA_IN_2=emufy(inp.OAA_IN__C2,[2]);
	inC2.AC_IN_1=emufy(inp.AC_IN__C2,[1]);
	inC2.OAA_IN_4=emufy(inp.OAA_IN__C2,[4]);
	inC2.OAA_IN_3=emufy(inp.OAA_IN__C2,[3]);
	inC2.Glu_IN_3=emufy(inp.Glu_IN__C2,[3]);
	inC2.OAA_IN_1=emufy(inp.OAA_IN__C2,[1]);
	inC2.Glu_IN_1=emufy(inp.Glu_IN__C2,[1]);
	inC2.GLC_IN_3=emufy(inp.GLC_IN__C2,[3]);
	inC2.GLC_IN_2=emufy(inp.GLC_IN__C2,[2]);
	inC2.Gln_IN_3=emufy(inp.Gln_IN__C2,[3]);
	inC2.Glu_IN_2=emufy(inp.Glu_IN__C2,[2]);
	inC2.GLC_IN_1=emufy(inp.GLC_IN__C2,[1]);
	inC2.Gln_IN_1=emufy(inp.Gln_IN__C2,[1]);
	inC2.Glu_IN_4=emufy(inp.Glu_IN__C2,[4]);
	inC2.Glu_IN_5=emufy(inp.Glu_IN__C2,[5]);
	inC2.Gln_IN_2=emufy(inp.Gln_IN__C2,[2]);
	inC2.GLC_IN_4=emufy(inp.GLC_IN__C2,[4]);
	inC2.Gln_IN_4=emufy(inp.Gln_IN__C2,[4]);
	inC2.Gln_IN_5=emufy(inp.Gln_IN__C2,[5]);
	inC2.GLC_IN_5=emufy(inp.GLC_IN__C2,[5]);
	inC2.GLC_IN_6=emufy(inp.GLC_IN__C2,[6]);
	inC2.AC_IN_1_2=emufy(inp.AC_IN__C2,[1 2]);
	inC2.OAA_IN_2_3=emufy(inp.OAA_IN__C2,[2 3]);
	inC2.Glu_IN_2_3=emufy(inp.Glu_IN__C2,[2 3]);
	inC2.GLC_IN_1_2=emufy(inp.GLC_IN__C2,[1 2]);
	inC2.Gln_IN_2_3=emufy(inp.Gln_IN__C2,[2 3]);
	inC2.Glu_IN_4_5=emufy(inp.Glu_IN__C2,[4 5]);
	inC2.GLC_IN_2_3=emufy(inp.GLC_IN__C2,[2 3]);
	inC2.Gln_IN_4_5=emufy(inp.Gln_IN__C2,[4 5]);
	inC2.GLC_IN_4_5=emufy(inp.GLC_IN__C2,[4 5]);
	inC2.OAA_IN_1_2=emufy(inp.OAA_IN__C2,[1 2]);
	inC2.OAA_IN_3_4=emufy(inp.OAA_IN__C2,[3 4]);
	inC2.GLC_IN_5_6=emufy(inp.GLC_IN__C2,[5 6]);
	inC2.Glu_IN_3_4=emufy(inp.Glu_IN__C2,[3 4]);
	inC2.Gln_IN_3_4=emufy(inp.Gln_IN__C2,[3 4]);
	inC2.Glu_IN_1_2=emufy(inp.Glu_IN__C2,[1 2]);
	inC2.Gln_IN_1_2=emufy(inp.Gln_IN__C2,[1 2]);
	inC2.OAA_IN_2_3_4=emufy(inp.OAA_IN__C2,[2 3 4]);
	inC2.OAA_IN_1_2_3=emufy(inp.OAA_IN__C2,[1 2 3]);
	inC2.Glu_IN_1_2_3=emufy(inp.Glu_IN__C2,[1 2 3]);
	inC2.GLC_IN_1_2_3=emufy(inp.GLC_IN__C2,[1 2 3]);
	inC2.Gln_IN_1_2_3=emufy(inp.Gln_IN__C2,[1 2 3]);
	inC2.GLC_IN_4_5_6=emufy(inp.GLC_IN__C2,[4 5 6]);
	inC2.Glu_IN_2_3_4=emufy(inp.Glu_IN__C2,[2 3 4]);
	inC2.Glu_IN_3_4_5=emufy(inp.Glu_IN__C2,[3 4 5]);
	inC2.Gln_IN_2_3_4=emufy(inp.Gln_IN__C2,[2 3 4]);
	inC2.Gln_IN_3_4_5=emufy(inp.Gln_IN__C2,[3 4 5]);
	inC2.OAA_IN_1_2_3_4=emufy(inp.OAA_IN__C2,[1 2 3 4]);
	inC2.Glu_IN_2_3_4_5=emufy(inp.Glu_IN__C2,[2 3 4 5]);
	inC2.GLC_IN_3_4_5_6=emufy(inp.GLC_IN__C2,[3 4 5 6]);
	inC2.Gln_IN_2_3_4_5=emufy(inp.Gln_IN__C2,[2 3 4 5]);
	inC2.Glu_IN_1_2_3_4_5=emufy(inp.Glu_IN__C2,[1 2 3 4 5]);
	inC2.Gln_IN_1_2_3_4_5=emufy(inp.Gln_IN__C2,[1 2 3 4 5]);
	inC2.GLC_IN_2_3_4_5_6=emufy(inp.GLC_IN__C2,[2 3 4 5 6]);
	inC2.GLC_IN_1_2_3_4_5_6=emufy(inp.GLC_IN__C2,[1 2 3 4 5 6]);
end
mol_=CARTcell_ccm_C(free_net,free_xch,inC2);
mol.CO2__C2=mol_.CO2;
mol.GLC__C2=mol_.GLC;
mol.G6P__C2=mol_.G6P;
mol.F6P__C2=mol_.F6P;
mol.FBP__C2=mol_.FBP;
mol.GAP__C2=mol_.GAP;
mol.DHAP__C2=mol_.DHAP;
mol.BPG__C2=mol_.BPG;
mol.PGA__C2=mol_.PGA;
mol.PEP__C2=mol_.PEP;
mol.PYR__C2=mol_.PYR;
mol.LAC__C2=mol_.LAC;
mol.OAA__C2=mol_.OAA;
mol.MAL__C2=mol_.MAL;
mol.m6PG__C2=mol_.m6PG;
mol.Ru5P__C2=mol_.Ru5P;
mol.R5P__C2=mol_.R5P;
mol.X5P__C2=mol_.X5P;
mol.S7P__C2=mol_.S7P;
mol.E4P__C2=mol_.E4P;
mol.AcCoA__C2=mol_.AcCoA;
mol.CitICit__C2=mol_.CitICit;
mol.OGA__C2=mol_.OGA;
mol.SuccCoA__C2=mol_.SuccCoA;
mol.Succ__C2=mol_.Succ;
mol.Fum__C2=mol_.Fum;
mol.Glu__C2=mol_.Glu;
mol.Ala__C2=mol_.Ala;
mol.Gln__C2=mol_.Gln;

persistent inC3
if isempty(inC3)
	inC3.CO2_IN_1=emufy(inp.CO2_IN__C3,[1]);
	inC3.AC_IN_2=emufy(inp.AC_IN__C3,[2]);
	inC3.OAA_IN_2=emufy(inp.OAA_IN__C3,[2]);
	inC3.AC_IN_1=emufy(inp.AC_IN__C3,[1]);
	inC3.OAA_IN_4=emufy(inp.OAA_IN__C3,[4]);
	inC3.OAA_IN_3=emufy(inp.OAA_IN__C3,[3]);
	inC3.Glu_IN_3=emufy(inp.Glu_IN__C3,[3]);
	inC3.OAA_IN_1=emufy(inp.OAA_IN__C3,[1]);
	inC3.Glu_IN_1=emufy(inp.Glu_IN__C3,[1]);
	inC3.GLC_IN_3=emufy(inp.GLC_IN__C3,[3]);
	inC3.GLC_IN_2=emufy(inp.GLC_IN__C3,[2]);
	inC3.Gln_IN_3=emufy(inp.Gln_IN__C3,[3]);
	inC3.Glu_IN_2=emufy(inp.Glu_IN__C3,[2]);
	inC3.GLC_IN_1=emufy(inp.GLC_IN__C3,[1]);
	inC3.Gln_IN_1=emufy(inp.Gln_IN__C3,[1]);
	inC3.Glu_IN_4=emufy(inp.Glu_IN__C3,[4]);
	inC3.Glu_IN_5=emufy(inp.Glu_IN__C3,[5]);
	inC3.Gln_IN_2=emufy(inp.Gln_IN__C3,[2]);
	inC3.GLC_IN_4=emufy(inp.GLC_IN__C3,[4]);
	inC3.Gln_IN_4=emufy(inp.Gln_IN__C3,[4]);
	inC3.Gln_IN_5=emufy(inp.Gln_IN__C3,[5]);
	inC3.GLC_IN_5=emufy(inp.GLC_IN__C3,[5]);
	inC3.GLC_IN_6=emufy(inp.GLC_IN__C3,[6]);
	inC3.AC_IN_1_2=emufy(inp.AC_IN__C3,[1 2]);
	inC3.OAA_IN_2_3=emufy(inp.OAA_IN__C3,[2 3]);
	inC3.Glu_IN_2_3=emufy(inp.Glu_IN__C3,[2 3]);
	inC3.GLC_IN_1_2=emufy(inp.GLC_IN__C3,[1 2]);
	inC3.Gln_IN_2_3=emufy(inp.Gln_IN__C3,[2 3]);
	inC3.Glu_IN_4_5=emufy(inp.Glu_IN__C3,[4 5]);
	inC3.GLC_IN_2_3=emufy(inp.GLC_IN__C3,[2 3]);
	inC3.Gln_IN_4_5=emufy(inp.Gln_IN__C3,[4 5]);
	inC3.GLC_IN_4_5=emufy(inp.GLC_IN__C3,[4 5]);
	inC3.OAA_IN_1_2=emufy(inp.OAA_IN__C3,[1 2]);
	inC3.OAA_IN_3_4=emufy(inp.OAA_IN__C3,[3 4]);
	inC3.GLC_IN_5_6=emufy(inp.GLC_IN__C3,[5 6]);
	inC3.Glu_IN_3_4=emufy(inp.Glu_IN__C3,[3 4]);
	inC3.Gln_IN_3_4=emufy(inp.Gln_IN__C3,[3 4]);
	inC3.Glu_IN_1_2=emufy(inp.Glu_IN__C3,[1 2]);
	inC3.Gln_IN_1_2=emufy(inp.Gln_IN__C3,[1 2]);
	inC3.OAA_IN_2_3_4=emufy(inp.OAA_IN__C3,[2 3 4]);
	inC3.OAA_IN_1_2_3=emufy(inp.OAA_IN__C3,[1 2 3]);
	inC3.Glu_IN_1_2_3=emufy(inp.Glu_IN__C3,[1 2 3]);
	inC3.GLC_IN_1_2_3=emufy(inp.GLC_IN__C3,[1 2 3]);
	inC3.Gln_IN_1_2_3=emufy(inp.Gln_IN__C3,[1 2 3]);
	inC3.GLC_IN_4_5_6=emufy(inp.GLC_IN__C3,[4 5 6]);
	inC3.Glu_IN_2_3_4=emufy(inp.Glu_IN__C3,[2 3 4]);
	inC3.Glu_IN_3_4_5=emufy(inp.Glu_IN__C3,[3 4 5]);
	inC3.Gln_IN_2_3_4=emufy(inp.Gln_IN__C3,[2 3 4]);
	inC3.Gln_IN_3_4_5=emufy(inp.Gln_IN__C3,[3 4 5]);
	inC3.OAA_IN_1_2_3_4=emufy(inp.OAA_IN__C3,[1 2 3 4]);
	inC3.Glu_IN_2_3_4_5=emufy(inp.Glu_IN__C3,[2 3 4 5]);
	inC3.GLC_IN_3_4_5_6=emufy(inp.GLC_IN__C3,[3 4 5 6]);
	inC3.Gln_IN_2_3_4_5=emufy(inp.Gln_IN__C3,[2 3 4 5]);
	inC3.Glu_IN_1_2_3_4_5=emufy(inp.Glu_IN__C3,[1 2 3 4 5]);
	inC3.Gln_IN_1_2_3_4_5=emufy(inp.Gln_IN__C3,[1 2 3 4 5]);
	inC3.GLC_IN_2_3_4_5_6=emufy(inp.GLC_IN__C3,[2 3 4 5 6]);
	inC3.GLC_IN_1_2_3_4_5_6=emufy(inp.GLC_IN__C3,[1 2 3 4 5 6]);
end
mol_=CARTcell_ccm_C(free_net,free_xch,inC3);
mol.CO2__C3=mol_.CO2;
mol.GLC__C3=mol_.GLC;
mol.G6P__C3=mol_.G6P;
mol.F6P__C3=mol_.F6P;
mol.FBP__C3=mol_.FBP;
mol.GAP__C3=mol_.GAP;
mol.DHAP__C3=mol_.DHAP;
mol.BPG__C3=mol_.BPG;
mol.PGA__C3=mol_.PGA;
mol.PEP__C3=mol_.PEP;
mol.PYR__C3=mol_.PYR;
mol.LAC__C3=mol_.LAC;
mol.OAA__C3=mol_.OAA;
mol.MAL__C3=mol_.MAL;
mol.m6PG__C3=mol_.m6PG;
mol.Ru5P__C3=mol_.Ru5P;
mol.R5P__C3=mol_.R5P;
mol.X5P__C3=mol_.X5P;
mol.S7P__C3=mol_.S7P;
mol.E4P__C3=mol_.E4P;
mol.AcCoA__C3=mol_.AcCoA;
mol.CitICit__C3=mol_.CitICit;
mol.OGA__C3=mol_.OGA;
mol.SuccCoA__C3=mol_.SuccCoA;
mol.Succ__C3=mol_.Succ;
mol.Fum__C3=mol_.Fum;
mol.Glu__C3=mol_.Glu;
mol.Ala__C3=mol_.Ala;
mol.Gln__C3=mol_.Gln;