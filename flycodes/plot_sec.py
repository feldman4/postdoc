import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from tqdm.auto import tqdm

from postdoc.utils import codify, split_by_mask


def setup_figure():
    widths = 1, 3.75, 1
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(3, 3, width_ratios=widths)
    ax_mean_norm = fig.add_subplot(gs[0, 0])
    ax_mean = fig.add_subplot(gs[1, 0])
    ax_mean_log = fig.add_subplot(gs[2, 0])
    ax_sec_norm = fig.add_subplot(gs[0, 1], sharey=ax_mean_norm)
    ax_sec = fig.add_subplot(gs[1, 1])
    ax_sec_log = fig.add_subplot(gs[2, 1], sharey=ax_mean_log)
    ax_not_used = fig.add_subplot(gs[0, 2])
    # ax_not_used.set_visible(False)
    ax_not_used.axis('off')
    return fig, (ax_mean_norm, ax_mean, ax_mean_log), (ax_sec_norm, ax_sec, ax_sec_log)


def prepare_plot_means(df_int, intensity):
    intensity_norm = f'{intensity}_peak_norm'
    return (df_int
    .groupby(['barcode', 'stage', 'plot_barcode'])
    [intensity].apply(lambda x: x.mean(skipna=False))
    .reset_index()
    .query('~(stage == "SEC" & plot_barcode == False)')
    .pipe(add_intensity_metrics, intensity)
    .pipe(codify, barcode_cat='barcode')
    )


def prepare_plot_sec(df_int, intensity):
    return (df_int
    .query('stage == "SEC"')
    .pipe(add_intensity_metrics, intensity)
    .pipe(codify, barcode_cat='barcode')
    )


def add_intensity_metrics(df_int, intensity):
    return (df_int
    .assign(**{f'{intensity}_max': lambda x: x.groupby('barcode')[intensity].transform('max')})
    .assign(**{f'{intensity}_sum': lambda x: x.groupby('barcode')[intensity].transform('sum')})
    .assign(**{f'{intensity}_peak_norm': lambda x: x.eval(f'{intensity} / {intensity}_max')})
    .assign(**{f'{intensity}_L1_norm': lambda x: x.eval(f'{intensity} / {intensity}_sum')})
    )



def plot(df_int, stages, fraction_centers, intensity):
    def plot_markers(ax, yvar):
        for _, row in df_plot_means.iterrows():
            color = colors.get(row['barcode'], 'gray')
            fill = color if row['plot_barcode'] else 'None'
            x = stages.index(row['stage'])
            ax.scatter(x=x, y=row[yvar], 
                    edgecolors=color, facecolors=fill)
            
    def plot_lines(ax, yvar):
        for barcode, row in df_plot_means_wide[yvar][stages].iterrows():
            ax.plot(row, color=colors[barcode], label=barcode)

    intensity_peak_norm = f'{intensity}_peak_norm'
    intensity_L1_norm = f'{intensity}_L1_norm'
    
    df_plot_sec = prepare_plot_sec(df_int, intensity)    
    df_plot_means = prepare_plot_means(df_int, intensity)
    barcodes = list(df_int['barcode'].drop_duplicates())

    # stage plots (SEC stage is averaged)
    values = [intensity, intensity_peak_norm, intensity_L1_norm]
    df_plot_means_wide = (df_plot_means
     .pivot_table(index='barcode', columns='stage', values=values)
     .stack(level=0).reindex(columns=stages).unstack(level=1).reorder_levels([1, 0], axis=1)
    )

    palette = sns.color_palette('Set1', n_colors=len(barcodes))
    colors = {b: c for b,c in zip(barcodes, palette)}

    fig, axs_mean, axs_sec = setup_figure()
    ax_mean_norm, ax_mean, ax_mean_log = axs_mean
    ax_sec_norm, ax_sec, ax_sec_log = axs_sec

    ax_mean.set_ylabel(intensity_L1_norm)
    ax_mean_log.set_ylabel(intensity)
    ax_mean_norm.set_ylabel(intensity_peak_norm)

    plot_lines(ax_mean_norm, intensity_peak_norm)
    plot_lines(ax_mean, intensity_L1_norm)
    plot_lines(ax_mean_log, intensity)
    
    plot_markers(ax_mean_norm, intensity_peak_norm)
    plot_markers(ax_mean, intensity_L1_norm)
    plot_markers(ax_mean_log, intensity)
    
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
    
    df_plot_sec_wide = (df_plot_sec
     .pivot_table(index='barcode', columns='fraction_center', 
                  values=values + ['plot_barcode'])
     .stack(level=0).reindex(columns=fraction_centers).unstack(level=1).reorder_levels([1, 0], axis=1)
    )
    def plot_lines(ax, yvar):
        for barcode, row in df_plot_sec_wide.iterrows():
            traces = split_by_mask(row[yvar], row['plot_barcode'].fillna(False))
            for trace in traces:
                ax.plot(trace, color=colors[barcode], label=barcode)

    def plot_markers(ax, yvar):
        for _, row in df_plot_sec.iterrows():
            color = colors[row['barcode']]
            fill = color if row['plot_barcode'] else 'None'
            x = row['fraction_center']
            ax.scatter(x=x, y=row[yvar], marker='.',
                    edgecolors=color, facecolors=fill)

    plot_lines(ax_sec_norm, intensity_peak_norm)
    plot_lines(ax_sec, intensity_L1_norm)
    plot_lines(ax_sec_log, intensity)

    plot_markers(ax_sec_norm, intensity_peak_norm)
    plot_markers(ax_sec, intensity_L1_norm)
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
    
    ymax = (df_int
     .query('stage == "SEC"')
     .pipe(add_intensity_metrics, intensity)
     .query('plot_barcode')
     [intensity_L1_norm].max()
    )

    ax_sec.set_ylim([0, 1.1 * ymax])
    # ax_mean.ticklabel_format(useOffset=False, style='plain', axis='y')
    # ax_mean.set_yticks(ax_sec.get_yticks()) # avoid mpl warning
    # ax_mean.set_yticklabels([f'{y:.1e}' for y in ax_mean.get_yticks()])
    
    # legend with design name in title and barcode colors
    handles, labels = ax_mean.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', frameon=False,
               bbox_to_anchor=[0.82, 0, 1, .95])

    fig.tight_layout()

    return fig


def plot_ratio_vs_barcodes(df_barcode_metrics, min_barcodes=20, xlim=(-2, 1), wrap=4):
    """Stratify barcodes by amino acid counts and plot soluble/insoluble ratio cumulative
    distributions within categories.
    """
    bins = np.linspace(xlim[0], xlim[1] + 0.1, 300)

    with sns.plotting_context('notebook', rc={'axes.titlesize': 20}):
        ratio_label = 'log(soluble / insoluble)'
        df_plot = (df_barcode_metrics
         .assign(ratio=lambda x: x.eval('log10(soluble / insoluble)'))
         .pipe(add_aa_counts)
         .filter(regex='[^K]_count$|^ratio$')
         .rename(columns={'hydro_count': 'YVL_count'})
         .set_index('ratio').stack().reset_index()
         .rename(columns={'level_1': 'variable', 0: 'count'})
         .loc[lambda x: x.groupby(['variable', 'count'])['ratio']
                         .transform('size') > min_barcodes]
         .sort_values(['variable', 'count'])
         .assign(variable=lambda x: x['variable'].str.split('_').str[0])
        )

        fg = (df_plot
         .pipe(sns.FacetGrid, hue='count', col='variable', col_wrap=4, palette='nipy_spectral', height=2)
         .map(plt.hist, 'ratio', bins=bins, cumulative=True, density=True, histtype='step', lw=2)
         .add_legend()
         .set_titles(col_template='{col_name}')
        )

        fg.axes.flat[0].set_xlim(xlim)
        [ax.set_xlabel('') for ax in fg.axes.flat[:]]

        fg.fig.text(0.5, 0.04, 'log10(soluble / insoluble)', ha='center', fontsize=18)
        fg.fig.text(0, 0.5, 'cumulative barcode fraction', va='center', rotation='vertical', fontsize=18)

    return fg, df_plot


def add_aa_counts(df, col='barcode'):
    df = df.copy()
    for aa in sorted(set(''.join(df['barcode']))):
        df[f'{aa}_count'] = df['barcode'].str.count(aa)
    return df


def plot_validation_overlays(df_summary, df_uv_data, traces, baseline_range=(7, 21)):

    it = (df_summary
     .dropna(subset=['num_ms_barcodes'])
     [['design_name', 'dataset', 'export_name', 'description']].values
    )

    # for baseline normalization
    v0, v1 = baseline_range

    def minmax(x):
        return (x - x.min()) / (x.max() - x.min())

    def baseline(x):
        return x - x.min()

    def normalize_to_X0(df, X0):
        """First minmax normalize to remove baseline, then divide by integration over volume range of X0.
        """
        df = df.copy()
        df['amplitude_0'] = minmax(df['amplitude'])
        df_ = df.query('@v0_ < volume < @v1_')
        v0_, v1_ = X0.index.min(), X0.index.max()
        integrated = df['amplitude_0'].sum()
        return (df['amplitude_0'] / integrated) * (df_.shape[0] / X0.shape[0])

    for design_name, dataset, export_name, description in tqdm(list(it)):
        X0 = traces[(dataset, 'barcodes')].loc[design_name].T
        X0 = X0.iloc[:, :13]
        X0 = X0 / X0.sum()
        y1 = X0.max().mean() * 1.3

        fig = plt.figure(figsize=(8, 7))

        widths = 0.9, 0.1
        gs = fig.add_gridspec(2, 2, width_ratios=widths)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)

        it_ = (df_uv_data
         .query('@v0 < volume < @v1')
         .query('export_name == @export_name').groupby('channel'))
        for channel, df in it_:
            y = baseline(df['amplitude'])
            y = (y / y.max()) * y1
            ax0.plot(df['volume'], y, label=channel, ls='dashed', lw=4, alpha=0.3, zorder=-10)
            ax1.plot(df['volume'], (df['amplitude']), label=channel)

        X0.plot(ax=ax0, marker='.')
        try:
            X1 = traces[(dataset, 'consensus')].loc[design_name].T
            (X1 + 0.005).plot(ax=ax0, label='consensus', color='black')
        except KeyError:
            pass
        ax0.set_ylim([0, y1])

        title = ': '.join(df.iloc[0][['SystemName', 'MethodStartTime']])
        ax1.set_title(title)
        ax1.set_xlabel('volume')
        ax1.set_ylabel('raw amplitude')
        ax0.set_ylabel('normalized within pool range')
        ax0.set_title(f'{export_name}: {dataset}')

        ax0.legend().remove()
        ax0.legend(bbox_to_anchor=(1, 0, 1, 1), loc='upper left')

        f = f'{export_name}_{dataset}.png'
        fig.savefig(f)
        plt.close(fig)