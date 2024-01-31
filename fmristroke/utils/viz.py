from uuid import uuid4

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec as mgs
from nilearn import plotting
from nilearn.image import index_img, threshold_img
from niworkflows.viz.utils import extract_svg
from svgutils.transform import fromstring


def plot_multicomponents(
    tseries,
    anat_nii,
    stat_map_nii,
    lesion_nii,
):
    try:
        ts = (
            pd.read_csv(tseries, sep="\t") if isinstance(tseries, str) else tseries
        )
    except:
        return plt.figure()

    n_components = ts.shape[1]
    anat = nib.load(anat_nii)
    stat_map = nib.load(stat_map_nii)
    lesion = nib.load(lesion_nii)

    fig = plt.figure(figsize=(10, 5 * n_components))
    gs = mgs.GridSpec(nrows=n_components, ncols=1, wspace=0.0, hspace=0.05)
    for j, IC in enumerate(ts):
        plot_components(
            anat,
            index_img(stat_map, int(IC.split("_")[-1])),
            np.array(ts[IC]),
            gs=gs[j],
            contour=lesion,
            name=f"IC {IC.split('_')[-1]}",
            color="white",
        )

    return fig


def plot_ts(
    tseries,
    gs_ts,
    gs_dist=None,
    name=None,
    units=None,
    hide_x=True,
    color="b",
    nskip=0,
    cutoff=None,
    ylims=None,
    tr=None,
):
    import seaborn as sns

    # Define TR and number of frames
    notr = False
    if tr is None:
        notr = True
        tr = 1.0
    ntsteps = len(tseries)
    tseries = np.array(tseries)

    # Define nested GridSpec
    gs = mgs.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs_ts, width_ratios=[1, 100], wspace=0.0
    )

    ax_ts = plt.subplot(gs[1])
    ax_ts.grid(False)

    # Set 10 frame markers in X axis
    interval = max((ntsteps // 10, ntsteps // 5, 1))
    xticks = list(range(0, ntsteps)[::interval])
    ax_ts.set_xticks(xticks)

    if not hide_x:
        if notr:
            ax_ts.set_xlabel("time (frame #)")
        else:
            ax_ts.set_xlabel("time (s)")
            labels = tr * np.array(xticks)
            ax_ts.set_xticklabels(["%.02f" % t for t in labels.tolist()])
    else:
        ax_ts.set_xticklabels([])

    if name is not None:
        if units is not None:
            name += " [%s]" % units

        ax_ts.annotate(
            name,
            xy=(0.0, 0.7),
            xytext=(0, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            va="center",
            ha="left",
            color=color,
            size=8,
            bbox={
                "boxstyle": "round",
                "fc": "w",
                "ec": "none",
                "color": "none",
                "lw": 0,
                "alpha": 0.8,
            },
        )

    for side in ["top", "right"]:
        ax_ts.spines[side].set_color("none")
        ax_ts.spines[side].set_visible(False)

    if not hide_x:
        ax_ts.spines["bottom"].set_position(("outward", 20))
        ax_ts.xaxis.set_ticks_position("bottom")
    else:
        ax_ts.spines["bottom"].set_color("none")
        ax_ts.spines["bottom"].set_visible(False)

    # ax_ts.spines["left"].set_position(('outward', 30))
    ax_ts.spines["left"].set_color("none")
    ax_ts.spines["left"].set_visible(False)
    # ax_ts.yaxis.set_ticks_position('left')

    ax_ts.set_yticks([])
    ax_ts.set_yticklabels([])

    nonnan = tseries[~np.isnan(tseries)]
    if nonnan.size > 0:
        # Calculate Y limits
        valrange = nonnan.max() - nonnan.min()
        def_ylims = [
            nonnan.min() - 0.1 * valrange,
            nonnan.max() + 0.1 * valrange,
        ]
        if ylims is not None:
            if ylims[0] is not None:
                def_ylims[0] = min([def_ylims[0], ylims[0]])
            if ylims[1] is not None:
                def_ylims[1] = max([def_ylims[1], ylims[1]])

        # Add space for plot title and mean/SD annotation
        def_ylims[0] -= 0.1 * (def_ylims[1] - def_ylims[0])

        ax_ts.set_ylim(def_ylims)

        # Annotate stats
        maxv = nonnan.max()
        mean = nonnan.mean()
        stdv = nonnan.std()
        p95 = np.percentile(nonnan, 95.0)
    else:
        maxv = 0
        mean = 0
        stdv = 0
        p95 = 0

    if cutoff is None:
        cutoff = []

    for i, thr in enumerate(cutoff):
        ax_ts.plot((0, ntsteps - 1), [thr] * 2, linewidth=0.2, color="dimgray")

        ax_ts.annotate(
            "%.2f" % thr,
            xy=(0, thr),
            xytext=(-1, 0),
            textcoords="offset points",
            va="center",
            ha="right",
            color="dimgray",
            size=3,
        )

    ax_ts.plot(tseries, color=color, linewidth=0.8)
    ax_ts.set_xlim((0, ntsteps - 1))

    if gs_dist is not None:
        ax_dist = plt.subplot(gs_dist)
        sns.displot(tseries, vertical=True, ax=ax_dist)
        ax_dist.set_xlabel("Timesteps")
        ax_dist.set_ylim(ax_ts.get_ylim())
        ax_dist.set_yticklabels([])

        return [ax_ts, ax_dist], gs
    return ax_ts, gs


def plot_components(
    anat_nii,
    stat_map_nii,
    ts,
    gs,
    plot_params=None,
    cuts=None,
    label=None,
    name=None,
    contour=None,
    color="white",
    filled=False,
):
    plot_params = {} if plot_params is None else plot_params

    # nilearn 0.10.0 uses Nifti-specific methods
    anat_nii = nib.Nifti1Image.from_image(anat_nii)

    if contour:
        contour = nib.Nifti1Image.from_image(contour)

    # figure = plt.gcf()

    # Create grid
    grid = mgs.GridSpecFromSubplotSpec(
        nrows=2, ncols=1, subplot_spec=gs, height_ratios=[1] + [0.4]
    )

    # Components spatial map
    ax = plt.subplot(grid[0])

    # Plot each cut axis
    plot_params["title"] = label
    # Generate nilearn figure
    d = plotting.plot_stat_map(
        threshold_img(stat_map_nii, 1),
        display_mode="z",
        colorbar=False,
        cut_coords=6,
        bg_img=anat_nii,
        axes=ax,
        draw_cross=False,
        **plot_params,
    )

    if contour is not None:
        d.add_contours(
            contour, colors=color, filled=filled, levels=[0.5], linewidths=1
        )

    # ts of component
    plot_ts(ts, grid[1], name=name)

    return grid


def plot_lagmaps(
    anat_nii,
    stat_map_nii,
    div_id,
    plot_params=None,
    order=("z", "x", "y"),
    cuts=None,
    label=None,
    contour=None,
    compress="auto",
    vmax=10,
    color="black",
    filled=True,
):
    """
    Plots the foreground and background views
    Default order is: axial, coronal, sagittal
    """

    plot_params = {} if plot_params is None else plot_params

    # Use default MNI cuts if none defined
    if cuts is None:
        raise NotImplementedError  # TODO

    # nilearn 0.10.0 uses Nifti-specific methods
    anat_nii = nib.Nifti1Image.from_image(anat_nii)

    out_files = []
    if contour:
        contour = nib.Nifti1Image.from_image(contour)

    # Plot each cut axis
    for i, mode in enumerate(list(order)):
        plot_params["display_mode"] = mode
        plot_params["cut_coords"] = cuts[mode]
        if i == 0:
            plot_params["title"] = label
        else:
            plot_params["title"] = None

        # Generate nilearn figure
        display = plotting.plot_stat_map(
            stat_map_nii, bg_img=anat_nii, vmax=vmax, cmap="jet", **plot_params
        )
        if contour is not None:
            display.add_contours(
                contour,
                colors=color,
                filled=filled,
                levels=[0.5],
                linewidths=0.5,
            )

        svg = extract_svg(display, compress=compress)
        display.close()

        # Find and replace the figure_1 id.
        svg = svg.replace("figure_1", "%s-%s-%s" % (div_id, mode, uuid4()), 1)
        out_files.append(fromstring(svg))

    return out_files


def plot_multilesionconn(
    anat_nii, lesion_nii, stat_map_nii, titles, **plot_params
):
    """
    Plot representing lesion to voxels connectivity
    """
    plot_params = {} if plot_params is None else plot_params

    anat = nib.load(anat_nii)
    anat = nib.Nifti1Image.from_image(anat)

    if lesion_nii is not None:
        lesion_nii = nib.load(lesion_nii)
        lesion_nii = nib.Nifti1Image.from_image(lesion_nii)

    fig = plt.figure(figsize=(10, 5 * len(stat_map_nii)))
    gs = mgs.GridSpec(
        nrows=len(stat_map_nii), ncols=1, wspace=0.0, hspace=0.05
    )

    for i, (name, stat_map) in enumerate(zip(titles, stat_map_nii)):
        stat_map = nib.load(stat_map)
        stat_map = nib.Nifti1Image.from_image(stat_map)

        display = plotting.plot_stat_map(
            stat_map,
            bg_img=anat,
            threshold=0.05,
            cmap="jet",
            title=name,
            axes=plt.subplot(gs[i]),
            **plot_params,
        )

        if lesion_nii is not None:
            display.add_contours(
                lesion_nii,
                colors="black",
                filled=True,
                levels=[0.5],
                linewidths=0.5,
            )
    return fig


def plot_kdeplot(data, title=None):
    """
    Plot representing voxels correlations with lesion
    """

    fig, ax = plt.subplots(1, 1)

    for col in data:
        sns.kdeplot(data[col], shade=True, ax=ax, label=col)
    ax.axvline(x=0, linestyle="dashed", color="k")
    ax.set_title(title)
    ax.set_xlim([-1, 1])
    ax.set_ylabel("Probability density")
    ax.set_xlabel("Connection strength")
    ax.legend(bbox_to_anchor=(1.025, 1), borderaxespad=0)

    return fig

def plot_catplot(x, y, data, xlabel=None, ylabel=None):
    """
    Plot representing quality measure value for each pipeline.

    Args:
        x (str):
            Column name of data corresponding to quality measure.
        y (str):
            Column name of data used to group values.
        data (pd.DataFrame):
            Table of quality measures. It should contain at least one one column
            representing quality measure and one column for grouping variable.
            Usually each row corresponds to one pipeline.
        xlabel (str, optional):
            Custom x-axis label.
        ylabel (str, optional):
            Custom y-axis label.

    Returns:
        generated plot.
    """

    fig = sns.catplot(x=x, y=y, kind="bar", data=data)
    if xlabel:
        fig.ax.set_xlabel(xlabel)
    if ylabel:
        fig.ax.set_ylabel(ylabel)

    sns.despine(top=False, right=False)
    fig.ax.axvline(x=0, color="k")

    return fig


def make_motion_plot(
    group_conf_summary: pd.DataFrame,
    x=None,
    y= None,):
    """
    Generates plot presenting number of subjects excluded with high motion
    according specified thresholds.
    Args:
        group_conf_summary (DataFrame): dataframe, output from GroupConfounds
        output_path (str): output path where plot is saved

    Returns:
        created figure, ready to save
    """

    sns.set_style("ticks")
    colors = ['#00a074', '#fe6863']
    palette = sns.set_palette(colors, 2)

    fig, axes = plt.subplots(1, 1, figsize=(16, 7))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)


    p = sns.swarmplot(
        y=y,
        x=x,
        data=group_conf_summary,
        alpha=0.8,
        s=10,
        palette=palette,
        ax=axes,
    )

    p = sns.boxplot(
        y=y,
        x=x,
        data=group_conf_summary,
        showcaps=False,
        boxprops={'facecolor': 'None'},
        showfliers=False,
        ax=axes,
    )

    p.set(xlabel='')
    p.set(ylabel=y)
    p.tick_params(axis='both', which='both', length=6, width=2.2)

    return fig
