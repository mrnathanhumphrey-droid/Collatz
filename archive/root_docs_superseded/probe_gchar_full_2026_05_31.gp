\p 60
default(parisize, "4G");
print("=== gchar full scan v2 ===");
bnf = bnfinit(x^2 + 1);
K = bnf.nf;

p2_dec = idealprimedec(K, 2); p5_dec = idealprimedec(K, 5); p13_dec = idealprimedec(K, 13); p17_dec = idealprimedec(K, 17);
p2 = p2_dec[1]; p5a = p5_dec[1]; p5b = p5_dec[2]; p13a = p13_dec[1]; p13b = p13_dec[2]; p17a = p17_dec[1]; p17b = p17_dec[2];

conds_id = [p5a, p5b, p13a, p13b, p17a, p17b, idealmul(K,p5a,p5b), idealmul(K,p17a,p17b), idealmul(K,p5a,p17a), idealmul(K,p5a,p17b), idealmul(K,p5b,p17a), idealmul(K,p5b,p17b), idealmul(K,p2,p17a), idealmul(K,p2,p17b)];
conds_lab = ["p5a", "p5b", "p13a", "p13b", "p17a", "p17b", "(5)", "(17)", "p5a*p17a", "p5a*p17b", "p5b*p17a", "p5b*p17b", "p2*p17a", "p2*p17b"];
inf_types = [[1, 0], [2, 0], [1, 1], [0, 1], [0, 2], [-1, 0], [1, -1], [2, 1], [1, 2], [3, 0], [0, 3]];
inf_labs  = ["(1,0)", "(2,0)", "(1,1)", "(0,1)", "(0,2)", "(-1,0)", "(1,-1)", "(2,1)", "(1,2)", "(3,0)", "(0,3)"];

safeval(v) = if(type(v) == "t_COMPLEX" || type(v) == "t_REAL" || type(v) == "t_INT" || type(v) == "t_FRAC", v, if(type(v) == "t_SER", polcoef(v, 0), 0));

outf = "C:/Collatz/hecke_gchar_lvalues_2026_05_31.csv";
write(outf, "label;inftype;chi;L1_re;L1_im;L2_re;L2_im");
nlvals = 0;
nerrs = 0;

for(i=1, #conds_id, myid = conds_id[i]; lab = conds_lab[i]; gc = 0; iferr(gc = gcharinit(bnf, myid), E, print("gcharinit fail at ", lab); next); cyc = gc.cyc; print(""); print("--- ", lab, " --- cyc=", cyc); for(t=1, #inf_types, itype = inf_types[t]; ilab = inf_labs[t]; base_list = 0; iferr(base_list = gcharalgebraic(gc, [itype]), E, next); if(#base_list == 0, next); base = base_list[1]; n_twists = if(#cyc>=1 && cyc[1] > 0, cyc[1], 1); for(j=0, n_twists-1, chi = base; if(n_twists > 1, chi[1] = chi[1] + j); Lobj = 0; iferr(Lobj = lfuncreate([gc, chi]), E, nerrs = nerrs+1; next); L1v = 0; L2v = 0; iferr(L1v = safeval(lfun(Lobj, 1)), E, L1v = 0); iferr(L2v = safeval(lfun(Lobj, 2)), E, L2v = 0); iferr(write(outf, Strprintf("%s;%s;%s;%.50f;%.50f;%.50f;%.50f", lab, ilab, chi, real(L1v), imag(L1v), real(L2v), imag(L2v))), E, nerrs = nerrs+1; next); nlvals = nlvals + 1)); print("  done ", lab, ", running total nlvals=", nlvals));
print("");
print("Total L-values written: ", nlvals);
print("Errors: ", nerrs);
print("Output: ", outf);
quit;
