import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from postdoc.utils import codify, split_by_mask


def setup_figure():
    widths = 1, 4
    fig = plt.figure(figsize=(8, 10))
    gs = fig.add_gridspec(3, 2, width_ratios=widths)
    ax_mean_norm = fig.add_subplot(gs[0, 0])
    ax_mean = fig.add_subplot(gs[1, 0])
    ax_mean_log = fig.add_subplot(gs[2, 0])
    ax_sec_norm = fig.add_subplot(gs[0, 1], sharey=ax_mean_norm)
    ax_sec = fig.add_subplot(gs[1, 1], sharey=ax_mean)
    ax_sec_log = fig.add_subplot(gs[2, 1], sharey=ax_mean_log)
    
    return fig, (ax_mean_norm, ax_mean, ax_mean_log), (ax_sec_norm, ax_sec, ax_sec_log)


def prepare_plot_means(df_int, intensity):
    intensity_norm = f'{intensity}_peak_norm'
    return (df_int
    .groupby(['barcode', 'stage', 'plot_barcode'])
    [intensity].apply(lambda x: x.mean(skipna=False))
    .reset_index()
    .query('~(stage == "SEC" & plot_barcode == False)')
    .assign(**{f'{intensity}_max': lambda x: x.groupby('barcode')[intensity].transform('max')})
    .assign(**{intensity_norm: lambda x: x.eval(f'{intensity} / {intensity}_max')})
    .pipe(codify, barcode_cat='barcode')
    )


def prepare_plot_sec(df_int, intensity):
    intensity_norm = f'{intensity}_peak_norm'
    return (df_int
    .query('stage == "SEC"')
    .assign(**{f'{intensity}_max': lambda x: x.groupby('barcode')[intensity].transform('max')})
    .assign(**{intensity_norm: lambda x: x.eval(f'{intensity} / {intensity}_max')})
    .pipe(codify, barcode_cat='barcode')
    )


def plot(df_int, stages, fraction_centers, intensity):
    def plot_markers(ax, yvar):
        for _, row in df_plot_means.iterrows():
            color = colors[row['barcode']]
            fill = color if row['plot_barcode'] else 'None'
            x = stages.index(row['stage'])
            ax.scatter(x=x, y=row[yvar], 
                    edgecolors=color, facecolors=fill)
            
    def plot_lines(ax, yvar):
        for barcode, row in df_plot_means_wide[yvar][stages].iterrows():
            ax.plot(row, color=colors[barcode], label=barcode)

    intensity_norm = f'{intensity}_peak_norm'

    # stage plots (SEC stage is averaged)
    df_plot_means = prepare_plot_means(df_int, intensity)

    values = [intensity, intensity_norm]
    # return (df_plot_means
    #  .pivot_table(index='barcode', columns='stage', values=values))
    df_plot_means_wide = (df_plot_means
     .pivot_table(index='barcode', columns='stage', values=values)
     .stack(level=0).reindex(columns=stages).unstack(level=1).reorder_levels([1, 0], axis=1)
    )

    barcodes = df_plot_means['barcode'].drop_duplicates()
    palette = sns.color_palette('Set1', n_colors=len(barcodes))
    colors = {b: c for b,c in zip(barcodes, palette)}

    fig, axs_mean, axs_sec = setup_figure()
    ax_mean_norm, ax_mean, ax_mean_log = axs_mean
    ax_sec_norm, ax_sec, ax_sec_log = axs_sec

    ax_mean.set_ylabel(intensity)
    ax_mean_log.set_ylabel(intensity)
    ax_mean_norm.set_ylabel(f'{intensity}_peak_norm')

    plot_lines(ax_mean, intensity)
    plot_lines(ax_mean_log, intensity)
    plot_lines(ax_mean_norm, intensity_norm)

    plot_markers(ax_mean, intensity)
    plot_markers(ax_mean_log, intensity)
    plot_markers(ax_mean_norm, intensity_norm)

    ax_mean_log.set_yscale('log')
    for ax in axs_mean[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in axs_mean:
        ax.set_xlim([-0.25, 0.25 + len(stages) - 1])
    
    axs_mean[-1].set_xticks(np.arange(len(stages)))
    axs_mean[-1].set_xticklabels(stages, rotation=30)

    def plot_lines(ax, yvar):
        for barcode, row in df_sec_means_wide[yvar].iterrows():
            ax.plot(row, color=colors[barcode], label=barcode)

    # SEC plots
    df_plot_sec = prepare_plot_sec(df_int, intensity)
    df_plot_sec_wide = (df_plot_sec
     .pivot_table(index='barcode', columns='fraction_center', 
                  values=[intensity, intensity_norm, 'plot_barcode'])
     .stack(level=0).reindex(columns=fraction_centers).unstack(level=1).reorder_levels([1, 0], axis=1)
    )
    def plot_lines(ax, yvar):
        for barcode, row in df_plot_sec_wide.iterrows():
            traces = split_by_mask(row[yvar], row['plot_barcode'])
            for trace in traces:
                ax.plot(trace, color=colors[barcode], label=barcode)

    def plot_markers(ax, yvar):
        for _, row in df_plot_sec.iterrows():
            color = colors[row['barcode']]
            fill = color if row['plot_barcode'] else 'None'
            x = row['fraction_center']
            ax.scatter(x=x, y=row[yvar], marker='.',
                    edgecolors=color, facecolors=fill)

    plot_lines(ax_sec_norm, intensity_norm)
    plot_lines(ax_sec, intensity)
    plot_lines(ax_sec_log, intensity)

    plot_markers(ax_sec_norm, intensity_norm)
    plot_markers(ax_sec, intensity)
    plot_markers(ax_sec_log, intensity)

    start, end = df_plot_sec['fraction_center'].describe()[['min', 'max']]
    for ax in axs_sec:
        plt.setp(ax.get_yticklabels(), visible=False) 
        ax.set_xlim([start - 0.5, end + 0.5])
        ax.set_xticks(np.arange(int(start), int(end) + 1))
        ax.grid(axis='x')
        
    axs_sec[0].tick_params(top=True, labeltop=True)
    axs_sec[0].xaxis.set_label_position('top')
    
    axs_mean[0].tick_params(top=True, labeltop=True, rotation=30)
    for ax in (axs_sec[0], axs_sec[-1]):
        ax.set_xlabel('SEC fraction center (mL)')
    
    ax_mean.set_ylim([0, 1.1 * df_int['area_ms1'].max()])
    ax_mean.ticklabel_format(useOffset=False, style='plain', axis='y')
    ax_mean.set_yticks(ax_sec.get_yticks()) # avoid mpl warning
    ax_mean.set_yticklabels([f'{y:.1e}' for y in ax_mean.get_yticks()])
    

    fig.tight_layout()

    return fig