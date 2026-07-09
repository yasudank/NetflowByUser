from collections import defaultdict
import matplotlib.pyplot as plt

def sort_key(c):
    if c.startswith('sci_P'):
        return (0, int(c.split('_P')[1]))
    elif c == 'cal':
        return (1, 0)
    else:
        return (2, 0)

def plot_sky_distribution(res, tgt, plot_file_path):
    """Plot the sky distribution of targeted and un-assigned high priority objects."""
    print("\nPlotting sky distribution...")
    
    # Group targeted and un-targeted
    plot_data = defaultdict(lambda: {'ra': [], 'dec': []})
    unassigned_high_prio = {'ra': [], 'dec': []}

    assigned_tidxs = set()
    for vis in res:
        for tidx in vis.keys():
            assigned_tidxs.add(tidx)

    for tidx, t in enumerate(tgt):
        cls = t.targetclass
        if tidx in assigned_tidxs:
            plot_data[cls]['ra'].append(t.ra)
            plot_data[cls]['dec'].append(t.dec)
        elif cls in ['sci_P0', 'sci_P1', 'sci_P2']:
            unassigned_high_prio['ra'].append(t.ra)
            unassigned_high_prio['dec'].append(t.dec)

    classes = sorted(list(plot_data.keys()), key=sort_key)

    # Define a color map for priorities
    colors = {
        'sci_P0': 'red',
        'sci_P1': 'orange',
        'sci_P2': 'gold',
        'sci_P3': 'yellowgreen',
        'sci_P4': 'green',
        'sci_P6': 'cyan',
        'sci_P7': 'dodgerblue',
        'sci_P9': 'navy',
        'cal': 'magenta',
        'sky': 'gray'
    }

    plt.figure(figsize=(12, 10))

    # 1. Plot Unassigned High Priority targets first (background)
    if unassigned_high_prio['ra']:
        plt.scatter(unassigned_high_prio['ra'], unassigned_high_prio['dec'], 
                    label=f"Unassigned P0, P1, P2 ({len(unassigned_high_prio['ra'])})", 
                    color='black', marker='x', alpha=0.5, s=20)

    # 2. Plot lower priorities first (sky, cal, P9) so high priorities are on top
    for cls in reversed(classes):
        color = colors.get(cls, 'black')
        alpha = 0.2 if cls in ['sky', 'sci_P9'] else 0.8
        
        if cls == 'sky':
            marker_size = 5
        elif cls == 'cal':
            marker_size = 10
        elif cls.startswith('sci_P'):
            pri = int(cls.split('_P')[1])
            # Highest priority (P0) gets largest marker (50), lowest (P9) gets smallest (5)
            marker_size = max(5, 50 - pri * 5)
        else:
            marker_size = 15
        
        plt.scatter(plot_data[cls]['ra'], plot_data[cls]['dec'], 
                    label=f"Assigned {cls} ({len(plot_data[cls]['ra'])})", 
                    color=color, alpha=alpha, s=marker_size, edgecolors='none')

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title('Sky Distribution of Targeted Objects')

    # Adjust legend position
    plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), title="Target Classes")
    plt.grid(True, linestyle='--', alpha=0.5)

    # Invert x-axis to represent RA properly (east to the left)
    plt.gca().invert_xaxis()

    plt.tight_layout()
    plt.savefig(plot_file_path, dpi=200, bbox_inches='tight')
    print(f"Plot saved to {plot_file_path}")
