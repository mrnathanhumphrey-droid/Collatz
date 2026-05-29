import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=json.load(open("experiments_output/probe_syrac_decay_fit_2026_05_28.json"))
rows=d["rows"]; ns=np.array([r["n"] for r in rows]); mx=np.array([r["max"] for r in rows])
mean=np.array([r["mean"] for r in rows])
fits=d["fits"]
A=fits["4"]["power"]["A"]; rho=fits["4"]["geom"]["rho"]
a,c,b=fits["4"]["stretched"]["a_c_beta"]
xs=np.linspace(2,16,200)
# anchor power & geom at n=4
pw=mx[3]*(xs/4.0)**(-A); gm=mx[3]*rho**(xs-4); st=np.exp(a-c*xs**b)

fig,ax=plt.subplots(1,2,figsize=(15,6))
ax[0].semilogy(ns,mx,"o-",color="C0",label="max_xi |S_chi(n)(xi)|  (worst-case)")
ax[0].semilogy(ns,mean,"s-",color="C1",label="mean_xi |..|  (typical)")
ax[0].semilogy(xs,pw,":",color="C2",label=f"power n^-{A:.2f}")
ax[0].semilogy(xs,gm,"-.",color="C3",label=f"geometric {rho:.3f}^n")
ax[0].semilogy(xs,st,"--",color="C4",lw=2,label=f"stretched exp(-{c:.2f} n^{b:.2f})")
ax[0].set_xlabel("n"); ax[0].set_ylabel("|char fn|"); ax[0].grid(alpha=.3)
ax[0].set_title("Tao Prop 1.17 decay through n=16 (offset = stationary, exact)")
ax[0].legend(fontsize=9)

ratio=[mx[i]/mx[i-1] for i in range(1,len(ns))]
locA=[-np.log(mx[i]/mx[i-1])/np.log(ns[i]/ns[i-1]) for i in range(1,len(ns))]
ax[1].plot(ns[1:],ratio,"o-",color="C3",label="local ratio m_n/m_{n-1}")
ax[1].axhline(1.0,ls=":",color="gray")
ax[1].set_ylabel("local geometric ratio",color="C3"); ax[1].tick_params(axis="y",labelcolor="C3")
ax[1].set_ylim(0.6,1.02)
ax2=ax[1].twinx(); ax2.plot(ns[1:],locA,"s--",color="C2",label="local power A")
ax2.set_ylabel("local power-law A",color="C2"); ax2.tick_params(axis="y",labelcolor="C2")
ax[1].set_xlabel("n"); ax[1].grid(alpha=.3)
ax[1].set_title("ratio climbs to 1 (subgeometric); A creeps up (superpolynomial)")
fig.tight_layout(); fig.savefig("probe_syrac_charfn_decay_2026_05_28.png",dpi=110)
print("updated probe_syrac_charfn_decay_2026_05_28.png")
