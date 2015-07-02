import dadi
import numpy
import pylab
import demographic_models

from numpy import array

data = dadi.Spectrum.from_file('../data/dadiSpectrum_CanMorIbp.txt')

ns = data.sample_sizes
pts_l = [40, 50, 60]
params = array([1, 1, 1, 1, 1, 1, 1.5, 2])
upper_bound = [100, 100, 100, 100, 100, 100, 100, 100]
lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2]

func = demographic_models.simple_canary
func_ex = dadi.Numerics.make_extrap_log_func(func)

s_model = func_ex(params, ns, pts_l)
ll_s_model = dadi.Inference.ll_multinom(s_model, data)
print 'Model log-likelihood:', ll_s_model
s_theta = dadi.Inference.optimal_sfs_scaling(s_model, data)

p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=1000)
print 'Optimized parameters', repr(popt)
s_model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(s_model, data)
print 'Optimized log-likelihood:', ll_opt

pylab.figure(figsize=(9, 9), dpi=300)
dadi.Plotting.plot_3d_comp_multinom(s_model, data, vmin=1, resid_range=3,
                                    pop_ids =('Ca','Mo', 'Ib'))
pylab.savefig('../results/simple_unfolded.png',dpi=300)
