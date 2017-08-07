import seaborn as sbn
import numpy as np


def generate_genome_plot(ax, data, xlims=(0, 9717), top_axis=False):
    """Plots the Percent Cleaved vs the HXB2 start position for each gRNA."""
    
    data.plot(ax=ax, x = 'Start', y = 'Percent cleaved', 
          kind='scatter', color = 'k', s = 30)

    ax.set_xlim(*xlims)
    ax.set_xlabel('')
    ax.set_xticks([])
    ax.set_ylim(0, 100)
    ax.set_ylabel('Predicted percent cleaved', fontsize=16)

    if top_axis:
        tw_ax = ax.twiny()
        tw_ax.set_xlim(*xlims)
        tw_ax.set_xlabel('Position in HXB2', fontsize=16)
        sbn.despine(ax=tw_ax, top=False)
        tw_ax.spines['top'].set_visible(False)
    else:
        ax.set_xlabel('Position in HXB2', fontsize=16)

    sbn.despine(ax=ax)
    
def generate_entropy_plot(ax, data):
    """Plot the Percent cleaved vs the 20-mer entropy for each gRNA."""
    
    data.plot(ax=ax, x = 'Entropy (bits)', y = 'Percent cleaved', 
          kind='scatter', color = 'k', s = 30)
    sbn.regplot(data = data, x = 'Entropy (bits)', 
                y = 'Percent cleaved', 
                ax=ax, color = 'k')

    ax.set_xlabel('Entropy (bits)', fontsize=16)
    ax.set_ylabel('Predicted percent cleaved', fontsize=16)
    ax.set_ylim(0, 100)
    ax.set_xlim(0, 8)
    ax.bar(left = 0.5, width = 4.5,
               bottom = -1, height = 10, linewidth = 2,
               edgecolor='r', facecolor='None')

    sbn.despine(ax=ax)
    
def generate_heatmap(ax, freqs, gRNA, offset, for_strand=True, vmin=0, vmax=12):
    """Generate the frequency heatmap of the alignment and annotate the target sequence."""
    
    # Use seaborn to create the heatmap
    sbn.heatmap(-np.log10(freqs), ax=ax, 
                vmin = vmin, vmax = vmax,
                cbar_kws={'label': '-log(Frequency)'}, cmap='copper_r')
    ax.set_xticklabels(['%s\n%i' % (l, p) for p, l in enumerate(gRNA, offset)], rotation=0)

    # Create green boxes around the target sequence
    for n, l in enumerate(gRNA):
        bottom = (a for a, b in enumerate('TGCA.N') if b == l).next()
        ax.bar(left = n,
               width = 1,
               height = 1,
               bottom = bottom,
               linewidth = 2,
               edgecolor = 'r', 
               facecolor='None')
    [t.set_rotation(0) for t in ax.get_yticklabels()]

    # Add the penalties to the top of the figure
    pen_ax = ax.twiny()
    pen_ax.set_xlim(0, 23)
    pen_ax.set_xticks(np.arange(23)+0.5)
    penalties = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0,
                     0.389, 0.079, 0.445, 0.508, 0.613,
                     0.851, 0.732, 0.828, 0.615, 0.804,
                     0.685, 0.583]
    
    # Deal with orientation
    if for_strand:
        labels = penalties+list('***')
    else:
        labels = list('***') + penalties[::-1]
        
    pen_ax.set_xticklabels(labels, rotation = 90)
    

def generate_hiv_orfs(ax):
    """Create a diagram of the HXB2 ORFs."""
    
    ax.set_yticks([])
    ax.set_xticks([])

    # List of ORF positions from https://www.hiv.lanl.gov/components/sequence/HIV/search/help.html#region
    points = [("5' LTR", 1, 634, 1),
              ("Gag", 790, 2292, 1),
              ("Pol", 2085, 5096, 3),
              ("Vif", 5041, 5619, 1),
              ("Vpr", 5559, 5850, 3),
              ("Tat1", 5831, 6045, 2),
              ("Tat2", 8379, 8469, 1), 
              ("Rev1", 5970, 6045, 3), 
              ("Rev2", 8379, 8653, 2),
              ("Vpu", 6062, 6310, 2), 
              ("Env", 6225, 8795, 3), 
              ("Nef", 8797, 9417, 1), 
              ("3' LTR", 9086,  9719, 2)]

    # Create boxes for each ORF along with 
    for name, start, stop, frame in points:
        ax.bar(bottom = frame+0.1, height = 0.8,
               left = start, width = stop-start,
               edgecolor = 'k', facecolor='k')
        if stop-start > 500:
            ax.annotate(name, xy = (start + (stop-start)/2, frame+0.5, ),
                        ha = 'center', va='center', color='w', fontsize=14)

    a = 0.5
    
    # Add arrows joining Tat and Rev exons
    ax.annotate("Tat", xy = (5900, 2.5), xycoords = 'data', xytext = (7000, 1.5),
                fontsize=14, color='k',ha = 'center', va='center',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle="arc",
                                lw=2, alpha=a
                                ))
    ax.annotate("Tat", xy = (8379, 1.5), xycoords = 'data', xytext = (7000, 1.5),
                fontsize=14, color='k',ha = 'center', va='center',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle="arc",
                                lw=2, alpha=a
                                ))

    ax.annotate("Rev", xy = (6000, 3.5), xycoords = 'data', xytext = (7000, 2.5),
                fontsize=14, color='k',ha = 'center', va='center',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle="arc",
                                lw=2, alpha=a
                                ))
    ax.annotate("Rev", xy = (8379, 2.5), xycoords = 'data', xytext = (7000, 2.5),
                fontsize=14, color='k',ha = 'center', va='center',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle="arc",
                                lw=2, alpha=a
                                ))

    # Add annotations for shorter genes
    ax.annotate('Vpr', xy = (5559, 3.5), xycoords = 'data', fontsize=14,
                ha='right', va='center')
    ax.annotate('Vpr', xy = (5559, 3.5), xycoords = 'data', fontsize=14,
                ha='right', va='center')

    ax.annotate('Vpu', xy = (6200, 2), xycoords = 'data', fontsize=14,
                ha='center', va='bottom')
    ax.annotate('Vpu', xy = (6200, 2), xycoords = 'data', fontsize=14,
                ha='center', va='bottom')

    sbn.despine(ax=ax, left=True, bottom=True)    

    ax.invert_yaxis()
    
