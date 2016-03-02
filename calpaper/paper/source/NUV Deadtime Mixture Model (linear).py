import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gPhoton.dbasetools as dt

import pprint, pickle

scl = 1.4

inpath = '.'
outpath = '.'

pkl_file = open('{path}/stimdata.pkl'.format(path=inpath), 'rb')
dtdata = pickle.load(pkl_file)
pkl_file.close()
print dtdata.keys()

band = 'NUV'
ix = np.where(np.isfinite(dtdata[band]['globalcount']) &
              np.isfinite(dtdata[band]['stimcount']) &
              (dtdata[band]['exptime']>0))
rawexpt = dtdata[band]['t1'][ix]-dtdata[band]['t0'][ix]
x = dtdata[band]['globalcount'][ix]/rawexpt
y = dtdata[band]['stimcount'][ix]/rawexpt
xerr = np.sqrt(dtdata[band]['globalcount'][ix])/dtdata[band]['exptime'][ix]
yerr = np.sqrt(dtdata[band]['stimcount'][ix])/dtdata[band]['exptime'][ix]
print x.min(),x.max(),xerr.min(),xerr.max()
print y.min(),y.max(),yerr.min(),yerr.max()

plt.figure()
plt.title('{b} Stim vs. Global Count Rates'.format(b=band))
plt.xlabel('Global Count Rate (ct/s)')
plt.ylabel('Stim Count Rate (ct/s)')
plt.plot(x,y,'.',alpha=0.1,label='{b}'.format(b=band))

plt.figure()
plt.title('{b} Stim vs. Global Count Rates Errors'.format(b=band))
plt.xlabel('Global Count Rate Error (ct/s)')
plt.ylabel('Stim Count Rate Error (ct/s)')
plt.plot(xerr,yerr,'.',alpha=0.1,label='{b}'.format(b=band))

plt.show()

# Build a linear mixture model. Run MCMC on it.
# Based on [example by DFM].
import emcee

# Define the probabilistic model...
def lnprior(p,bounds):
    # We'll just put reasonable uniform priors on all the parameters.
    if not all([b[0] < v < b[1] for v, b in zip(p, bounds)]):
        return -np.inf
    return 0

# The "foreground" linear likelihood:
def lnlike_fg(p,x,y,yerr):
    m, b, _, M, lnV = p
    model = m * x + b
    return -0.5 * (((model - y) / yerr) ** 2 + 2 * np.log(yerr))

# The "background" outlier likelihood:
def lnlike_bg(p,x,y,yerr):
    _, _, Q, M, lnV = p
    var = np.exp(lnV) + yerr**2
    return -0.5 * ((M - y) ** 2 / var + np.log(var))

# Full probabilistic model.
def lnprob(p,x,y,yerr,bounds):
    m, b, Q, M, lnV = p
    lp = lnprior(p,bounds)
    if not np.isfinite(lp):
        return -np.inf, None
    ll_fg = lnlike_fg(p,x,y,yerr)
    arg1 = ll_fg + np.log(Q)
    ll_bg = lnlike_bg(p,x,y,yerr)
    arg2 = ll_bg + np.log(1.0 - Q)
    ll = np.sum(np.logaddexp(arg1, arg2))
    return lp + ll, (arg1, arg2)

# Initialize the walkers at a reasonable location.
ndim, nwalkers = 5, 32
params = ['m', 'b', 'Q', 'M', 'lnV']
bounds = [(-.001,0), (72, 81), (0, 1), (0, 100), (0, 10)]
p0 = np.array([-0.0005, 79, 0.5, 50, 5])

p0 = [p0 + 1e-5 * np.random.randn(ndim) for k in range(nwalkers)]

# Set up the sampler.
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr, bounds))

# Run a burn-in chain and save the final location.
pos, _, _, _ = sampler.run_mcmc(p0, 500)

# Run the production chain.
sampler.reset()
sampler.run_mcmc(pos, 1000);


# In[96]:

fig, axes = plt.subplots(len(params), 1, sharex=True, figsize=(8, len(params)*3))
for i,p in enumerate(params):
    axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
    axes[i].yaxis.set_major_locator(MaxNLocator(5))
    axes[i].set_ylabel("${p}$".format(p=p))


# In[97]:

import triangle
#labels = ["$m$", "$b$", "$Q$", "$M$", "$\ln V$"]
#truths = true_params + [true_frac, true_outliers[0], np.log(true_outliers[1])]
triangle.corner(sampler.flatchain, bins=50, extents=bounds, labels=params);


# In[158]:

samples = sampler.flatchain[:, :2]
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
print("""MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))


# In[99]:

gcr = np.arange(1,100000,1000)
tec2fdead=5.52e-6
feeclkratio=0.966
def dt_old(gcr):
    return tec2fdead*(gcr)/feeclkratio
def dt_new(gcr):
    return 1-(m_mcmc[0] * gcr + b_mcmc[0])/b_mcmc[0]

gcr = 50000
print 1-dt_old(gcr), 1-dt_new(gcr)


# In[100]:

norm = 0.0
post_prob = np.zeros(len(x))
for i in range(sampler.chain.shape[1]):
    for j in range(sampler.chain.shape[0]):
        ll_fg, ll_bg = sampler.blobs[i][j]
        post_prob += np.exp(ll_fg - np.logaddexp(ll_fg, ll_bg))
        norm += 1
post_prob /= norm


# In[154]:

plt.figure(figsize=(8*1,4*1))

# Plot the predition.
x0 = np.linspace(0,80000, 1000)
A = np.vander(x0, 2)
lines = np.dot(sampler.flatchain[:, :2], A.T)
quantiles = np.percentile(lines, [16, 84], axis=0)
plt.fill_between(x0, quantiles[0], quantiles[1], color="#8d44ad", alpha=0.5)

# Plot the "bad" points.
ix_bg = np.where(post_prob<.5)
bg_ct = len(ix_bg[0])
plt.errorbar(x[ix_bg], y[ix_bg], yerr=yerr[ix_bg], fmt=",k", marker='.',alpha=0.1)
# Plot the "good" points.
ix_fg = np.where(post_prob>=.5)
fg_ct = len(ix_fg[0])
plt.errorbar(x[ix_fg], y[ix_fg], yerr=yerr[ix_fg], fmt=",k", marker='o', ms=0,
    capsize=0, lw=1, zorder=999, alpha=0.1)
plt.text(15000, 45,
    r'$scr={m}^{{{mp}}}_{{{mm}}}gcr+{b}^{{{bp}}}_{{{bm}}}$'.format(
        m=round(m_mcmc[0],6),b=round(b_mcmc[0],2),
        mp='+{v}'.format(v=round(m_mcmc[1],6)),
        mm='-{v}'.format(v=round(m_mcmc[2],6)),
        bp='+{v}'.format(v=round(b_mcmc[1],2)),
        bm='-{v}'.format(v=round(b_mcmc[2],2))), fontsize=18)

plt.title('{b} Stim vs. Global Countrate (n={n})'.format(
    b=band,n=fg_ct+bg_ct),fontsize=16)
plt.xlabel("Global Countrate (ct/s)",fontsize=14)
plt.ylabel("Stim Countrate (ct/s)",fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(40, 75)
plt.xlim(10000, 70000)
plt.savefig('{path}/Stim_v_GCR_{b}.pdf'.format(path=outpath,b=band),
    format='pdf',dpi=1000)

print("""MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

print 'Tossing ~{p}% of data. ({n} of {m})'.format(
    p=100*bg_ct/(fg_ct+bg_ct),n=bg_ct,m=fg_ct+bg_ct)
