function C = case_list()
% CASE_LIST  Simulation cases grouped by the figure or table that uses them.
% Each field holds a cell array of case names (see case_params.m for the parameters).
fg = {}; for D = [2 5 10], for T = [2 5 10], fg{end+1} = sprintf('FG_D%d_T%d_F2', D, T); end, end %#ok<AGROW>
for D = [2 5 10], for F = [3 6], fg{end+1} = sprintf('FG_D%d_T5_F%d', D, F); end, end %#ok<AGROW>
C.Fig1to3_S1_S2 = fg;
C.Fig4 = {'R_F4S_s1', 'R_F4M_s1', 'R_F4C_s1'};   % seed 1, shared with Fig. 5
C.Fig5 = {'R_F5C_s1', 'R_F5B_s1'};               % seed 1, shared with Fig. 4
C.Fig6 = {'C6_base', 'F6B'};
ob = {}; for ph = [0 15 30 45 60], for dl = [0 5 10 15], ob{end+1} = sprintf('OB_a60_p%d_t%d', ph, dl); end, end %#ok<AGROW>
C.Fig7 = ob;
ob = {}; for ph = [0 30 60], for dl = [0 10], ob{end+1} = sprintf('OB_a40_p%d_t%d', ph, dl); end, end %#ok<AGROW>
C.FigS7 = ob;
om = {}; for a = [40 60], for dl = [0 5 10 15], for of = [0 6 12], om{end+1} = sprintf('OM_a%d_t%d_o%d', a, dl, of); end, end, end %#ok<AGROW>
C.Fig7_text_tip_to_tip = om;   % second set of runs described in the Fig. 7 subsection
C.FigS3_convergence = {'CS_base','CS_dt5','CS_dt1','CS_dx05','CS_r30','CS_r45','CS_L512', ...
    'CC_base','CC_dt5','CC_dt1','CC_dx05','N6C_dx1','N6C_dx05','N6B_dx1','N6B_dx05','C6_rot','CSZ_dx1','CSZ_dx05'};
C.FigS4_interface_width = {'WS_75','WS_125','WC_75','WC_125'};
sw = {}; for D = [2 5 10], for T = [2 5 10], sw{end+1} = sprintf('SW_D%d_T%d_F2', D, T); end, end %#ok<AGROW>
for D = [2 5 10], for F = [3 6], sw{end+1} = sprintf('SW_D%d_T5_F%d', D, F); end, end %#ok<AGROW>
C.FigS5_descriptors = sw;
C.FigS5_thin_interface = {'WSW_D2_T5_F2','WSW_D5_T5_F2','WSW_D10_T5_F2','WSW_D5_T2_F2','WSW_D5_T10_F2','WSW_D10_T5_F3','WSW_D10_T5_F6'};
C.FigS6_bicrystals = {'BX_0','BX_20','BX_40','BX_60'};
C.FigS8_MoS2 = {'MT_35','MT_20','MT_10','MT_03','MG_35','MG_20','MG_10','MG_03'};
r = {}; for cf = {'F4S','F4M','F4C','F5C','F5B'}, for s = 1:5, r{end+1} = sprintf('R_%s_s%d', cf{1}, s); end, end %#ok<AGROW>
C.TableS1_realisations = r;
C.symmetry_test = {'SYM_F4C_sym','SYM_F4S_sym','SYM_F6_sym','SYM_F5C_sym3','SYM_F5C_sym6','SYM_F5B_sym3','SYM_F5B_sym6'};
C.letter_unstable_dx025 = {'CS_dx025'};   % dx = 0.25 with dt = 1e-3: the theta equation is unstable (response letter, Reviewer 2 comment 6)
end
