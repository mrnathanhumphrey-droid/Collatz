\p 60
default(parisize, "2G");
print("=== gchar probe ===");
bnf = bnfinit(x^2 + 1);
K = bnf.nf;
p17a = idealprimedec(K, 17)[1];
print("Trying gcharinit at p17a (norm ", idealnorm(K, p17a), ")");
gc = gcharinit(bnf, p17a);
print("gc.cyc = ", gc.cyc);
print("gc.mod = ", gc.mod);
print("");
\\ Algebraic chars at various infinity types
for(k1=0, 4, for(k2=0, 4, algs = gcharalgebraic(gc, [k1, k2]); print("inf=[",k1,",",k2,"]: ", #algs, " algebraic chars"); if(#algs > 0 && #algs <= 4, for(j=1, #algs, print("  ", algs[j])))));
print("");
\\ Test L(s) at s=2 for one char
algs10 = gcharalgebraic(gc, [1, 0]);
algs01 = gcharalgebraic(gc, [0, 1]);
algs11 = gcharalgebraic(gc, [1, 1]);
print("inf (1,0): #=", #algs10);
print("inf (0,1): #=", #algs01);
print("inf (1,1): #=", #algs11);
if(#algs10 > 0, chi = algs10[1]; print("Sample chi=(1,0): ", chi); Lobj = lfuncreate([gc, chi]); L1val = lfun(Lobj, 1); L2val = lfun(Lobj, 2); L0val = lfun(Lobj, 0); print("  L(0) = ", L0val); print("  L(1) = ", L1val); print("  L(2) = ", L2val));
if(#algs11 > 0, chi = algs11[1]; print("Sample chi=(1,1): ", chi); Lobj = lfuncreate([gc, chi]); L1val = lfun(Lobj, 1); L2val = lfun(Lobj, 2); print("  L(1) = ", L1val); print("  L(2) = ", L2val));
print("=== Done ===");
quit;
