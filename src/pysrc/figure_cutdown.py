#!/usr/bin/env python3

#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import plotting_style as ps

from point_of_transmission import analytic_p_inocs, neutralize_prob_from_z

def plot_cutdown(
        axis=None,
        f_mut=3e-5,
        min_ratio=1,
        max_ratio=200,
        fineness=100,
        final_bottlenecks=[10, 1],
        cmaps=[plt.cm.Reds, plt.cm.Blues],
        color_low=0.4,
        color_high=0.8,
        drift_cmap=plt.cm.Greys,
        mut_neut_probs=[0.95, 0.99],
        z_wt=0.95,
        linestyle='dashed',
        dash_styles=[(1, 0), (1,1)],
        max_exact = 50,
        inoculum_model = "poisson",
        linealpha=1,
        **kwargs):
    
    if axis is None:
        fig, axis = plt.subplots()

    vb_ratios = np.linspace(min_ratio,
                            max_ratio,
                            fineness)

    n_mps = len(mut_neut_probs)
    colors = np.linspace(color_low, color_high, n_mps)

    for bn, cmap, dashes, in zip(final_bottlenecks, cmaps, dash_styles):
        drifts = [analytic_p_inocs(bn,
                                   f_mut,
                                   p_loss = 0,
                                   mut_neut_prob = 0,
                                   wt_neut_prob = 0,
                                   inoculum_size = bn * vb_ratio,
                                   verbose = False,
                                   exact = (bn * vb_ratio) < max_exact)
                  for vb_ratio in vb_ratios]
        axis.plot(vb_ratios,
                  drifts,
                  color=cmap(0.95),
                  linestyle='dashed',
                  **kwargs)

        for line_id, mnp in enumerate(mut_neut_probs):
                
            probs = [analytic_p_inocs(bn,
                                      f_mut,
                                      p_loss = 0,
                                      mut_neut_prob = mnp,
                                      wt_neut_prob = neutralize_prob_from_z(
                                          z_wt,
                                          bn * vb_ratio,
                                          inoculum_model),
                                      inoculum_size = bn * vb_ratio,
                                      exact = (bn * vb_ratio) < max_exact,
                                      verbose = False)
                     for vb_ratio in vb_ratios]
            axis.plot(vb_ratios,
                      probs,
                      color = cmap(colors[line_id]),
                      label = "$b={}$, $\kappa_m={}$".format(
                          bn, mnp),
                      alpha = linealpha,
                      linestyle = linestyle,
                      dashes = dashes,
                      **kwargs)

    axis.set_yscale('log')
    axis.set_xlim(left=1)
    axis.legend(ncol = 1,
                handlelength = 1)
    axis.set_xlabel('ratio of mucus bottleneck\n'
                    'to final bottleneck ($v / b$)')
    axis.set_ylabel('new variant\nsurvival probability')
