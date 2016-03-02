
"""Read and repackage the stim data that has been generated within the database
for a large number of eclipses.


There is a lot of redundant code between the FUV and NUV mixture models.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gPhoton.dbasetools as dt
import gPhoton.gQuery as gq
import pprint, pickle
import triangle
import emcee

inpath = '.'
outpath = '.'

#------------

"""Read and repackage the stim data that has been generated within the database for a large number of eclipses."""
bands = ['FUV','NUV']
querydata = pd.read_csv('{path}/stimquery.csv'.format(path=inpath))

dtdata = {}
for band in bands:
    stimdata = pd.read_csv('{b}Stim.csv'.format(b=band))
    n = len(stimdata)
    dtdata[band] = {'t0':np.zeros(n),'t1':np.zeros(n),
                    'stimcount':np.zeros(n),'globalcount':np.zeros(n),
                    'exptime':np.zeros(n), 'deadtime':np.zeros(n)}
    for i,data in enumerate(np.array(stimdata)):
        ix = np.where(querydata['ix']==data[0])
        t0 = np.array(querydata['t0'][ix[0]])[0]/1000.
        t1 = np.array(querydata['t1'][ix[0]])[0]/1000.
        print i,data[0],band,t0,t1
        exptime = t1-t0-dt.compute_shutter(band,[t0,t1])
        stimcount = data[3]
        globalcount = gq.getValue(gq.globalcounts(band,t0,t1))
        deadtime = gq.getValue(gq.deadtime(band,t0,t1))
        dtdata[band]['t0'][i] = t0
        dtdata[band]['t1'][i] = t1
        dtdata[band]['stimcount'][i] = stimcount
        dtdata[band]['globalcount'][i] = globalcount
        dtdata[band]['exptime'][i] = exptime
        dtdata[band]['deadtime'][i] = deadtime

# Note that this actually writes to 'inpath' because it's read in below
output = open('{path}/stimdata.pkl'.format(path=inpath), 'wb')
# Pickle dictionary using protocol 0.
pickle.dump(dtdata, output)
output.close()

#------------------------------------
# FUV linear mixture model

scl = 1.4

pkl_file = open('{path}/stimdata.pkl'.format(path=inpath), 'rb')
dtdata = pickle.load(pkl_file)
pkl_file.close()
band = 'FUV'
print dtdata[band].keys()

ix = np.where(np.isfinite(dtdata[band]['globalcount']) &
    np.isfinite(dtdata[band]['stimcount']) & (dtdata[band]['exptime']>0) &
    (dtdata[band]['stimcount']/dtdata[band]['exptime']<140))
rawexpt = dtdata[band]['t1'][ix]-dtdata[band]['t0'][ix]
x = dtdata[band]['globalcount'][ix]/rawexpt
y = dtdata[band]['stimcount'][ix]/rawexpt
xerr = np.sqrt(dtdata[band]['globalcount'][ix])/dtdata[band]['exptime'][ix]
yerr = np.sqrt(dtdata[band]['stimcount'][ix])/dtdata[band]['exptime'][ix]
print x.min(),x.max(),xerr.min(),xerr.max()
print y.min(),y.max(),yerr.min(),yerr.max()

# Plots the raw data. Useful for diagnostics.
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates'.format(b=band))
# plt.xlabel('Global Count Rate (ct/s)')
# plt.ylabel('Stim Count Rate (ct/s)')
# plt.plot(x,y,'.',alpha=0.25,label='{b}'.format(b=band))
#
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates Errors'.format(b=band))
# plt.xlabel('Global Count Rate Error (ct/s)')
# plt.ylabel('Stim Count Rate Error (ct/s)')
# plt.plot(xerr,yerr,'.',alpha=0.25,label='{b}'.format(b=band))
#
# plt.show()

# Build a linear mixture model. Run MCMC on it.
# Based on [example by DFM].
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
bounds = [(-.001,0), (75, 81), (0, 1), (0, 100), (0, 10)]
p0 = np.array([-0.0005, 79, 0.5, 50, 5])

p0 = [p0 + 1e-5 * np.random.randn(ndim) for k in range(nwalkers)]

# Set up the sampler.
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
    args=(x, y, yerr, bounds))

# Run a burn-in chain and save the final location.
pos, _, _, _ = sampler.run_mcmc(p0, 1000)

# Run the production chain.
sampler.reset()
sampler.run_mcmc(pos, 1000);

# Plots the walkers. Useful for diagnostics, particular of the priors.
# fig, axes = plt.subplots(len(params), 1, sharex=True,
#     figsize=(8, len(params)*3))
# for i,p in enumerate(params):
#     axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
#     axes[i].yaxis.set_major_locator(MaxNLocator(5))
#     axes[i].set_ylabel("${p}$".format(p=p))

# Generates a corner plot (of all variables against each other).
# Useful for diagnostics.
# triangle.corner(sampler.flatchain, bins=50, extents=bounds, labels=params);

samples = sampler.flatchain[:, :2]
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
print("""FUV - MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

# Setup reference variables for emperical deadtime correction
gcr = np.arange(1,100000,1000)
tec2fdead=5.52e-6
feeclkratio=0.966
dt_old = tec2fdead*(gcr)/feeclkratio
dt_new = 1-(m_mcmc[0] * gcr + b_mcmc[0])/b_mcmc[0]

norm = 0.0
post_prob = np.zeros(len(x))
for i in range(sampler.chain.shape[1]):
    for j in range(sampler.chain.shape[0]):
        ll_fg, ll_bg = sampler.blobs[i][j]
        post_prob += np.exp(ll_fg - np.logaddexp(ll_fg, ll_bg))
        norm += 1
post_prob /= norm

# Plot the prediction.
plt.figure(figsize=(8*1,4*1))
x0 = np.linspace(0,18000, 1000)
A = np.vander(x0, 2)
lines = np.dot(sampler.flatchain[:, :2], A.T)
quantiles = np.percentile(lines, [16, 84], axis=0)
plt.fill_between(x0, quantiles[0], quantiles[1], color="#8d44ad", alpha=0.5)

# Plot the "foreground" points.
ix_bg = np.where(post_prob<.5)
bg_ct = len(ix_bg[0])
plt.errorbar(x[ix_bg], y[ix_bg], yerr=yerr[ix_bg], fmt=",k",
    marker='.',alpha=0.1)

# Plot the "noise" points.
ix_fg = np.where(post_prob>=.5)
fg_ct = len(ix_fg[0])
plt.errorbar(x[ix_fg], y[ix_fg], yerr=yerr[ix_fg], fmt=",k", marker='o', ms=0,
    capsize=0, lw=1, zorder=999, alpha=0.1)
plt.text(5000, 76,
    r'$scr={m}^{{{mp}}}_{{{mm}}}gcr+{b}^{{{bp}}}_{{{bm}}}$'.format(
        m=round(m_mcmc[0],6),b=round(b_mcmc[0],2),
        mp='+{v}'.format(v=round(m_mcmc[1],6)),
        mm='-{v}'.format(v=round(m_mcmc[2],6)),
        bp='+{v}'.format(v=round(b_mcmc[1],2)),
        bm='-{v}'.format(v=round(b_mcmc[2],2))), fontsize=18)

plt.title('{b} Stim vs. Global Countrate (n={n})'.format(
    b=band,n=bg_ct+fg_ct),fontsize=16)
plt.xlabel("Global Countrate (ct/s)",fontsize=14)
plt.ylabel("Stim Countrate (ct/s)",fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(68, 78)
plt.xlim(0, 17000)
plt.savefig('{path}/Stim_v_GCR_{b}.pdf'.format(path=outpath,b=band),
    format='pdf',dpi=1000)
plt.close()

# Print a summary
print("""FUV - MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

print 'Tossing ~{p}% of data. ({n} of {m})'.format(
    p=100*bg_ct/(fg_ct+bg_ct),n=bg_ct,m=bg_ct+fg_ct)

#-------------------------------
# NUV linear mixture model
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

# Plot the data. Useful for diagnostics.
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates'.format(b=band))
# plt.xlabel('Global Count Rate (ct/s)')
# plt.ylabel('Stim Count Rate (ct/s)')
# plt.plot(x,y,'.',alpha=0.1,label='{b}'.format(b=band))
#
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates Errors'.format(b=band))
# plt.xlabel('Global Count Rate Error (ct/s)')
# plt.ylabel('Stim Count Rate Error (ct/s)')
# plt.plot(xerr,yerr,'.',alpha=0.1,label='{b}'.format(b=band))
#
# plt.show()

# Build a linear mixture model. Run MCMC on it.
# Based on [example by DFM].
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

# Plot the walkers. Useful for diagnostics.
# fig, axes = plt.subplots(len(params), 1, sharex=True, figsize=(8, len(params)*3))
# for i,p in enumerate(params):
#     axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
#     axes[i].yaxis.set_major_locator(MaxNLocator(5))
#     axes[i].set_ylabel("${p}$".format(p=p))

# Make "triangle" plots of every variable against every other. For diagnostics.
# triangle.corner(sampler.flatchain, bins=50, extents=bounds, labels=params);

samples = sampler.flatchain[:, :2]
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
print("""MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

gcr = np.arange(1,100000,1000)
tec2fdead=5.52e-6
feeclkratio=0.966

norm = 0.0
post_prob = np.zeros(len(x))
for i in range(sampler.chain.shape[1]):
    for j in range(sampler.chain.shape[0]):
        ll_fg, ll_bg = sampler.blobs[i][j]
        post_prob += np.exp(ll_fg - np.logaddexp(ll_fg, ll_bg))
        norm += 1
post_prob /= norm

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

# Report some results.
print("""MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

print 'Tossing ~{p}% of data. ({n} of {m})'.format(
    p=100*bg_ct/(fg_ct+bg_ct),n=bg_ct,m=fg_ct+bg_ct)
