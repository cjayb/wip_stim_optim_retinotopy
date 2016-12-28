# Authors: Christopher Bailey (cjb@cfin.au.dk)
#
# License: BSD (3-clause)

import matplotlib.pyplot as plt

import re
import numpy as np
import subprocess as subp
from collections import OrderedDict
import colorsys
from os.path import join, isfile
from os import environ, makedirs
from matplotlib.pylab import cm
from mne import read_labels_from_annot, read_forward_solution
from mne.io.pick import pick_types, pick_types_forward, pick_info
from mne.forward import (is_fixed_orient, convert_forward_solution)
from mne import Label, BiHemiLabel, SourceEstimate
from mne.simulation import simulate_evoked
from mne.preprocessing.maxwell import _get_coil_scale as mne_get_coil_scale


def hsv2rgb_int(h, s, v):
    return tuple(int(i * 255) for i in colorsys.hsv_to_rgb(h, s, v))


def execute(cmd, verbose=True, fake=False):
    if verbose:
        print(cmd)
    if not fake:
        subp.call(cmd, shell=True)


# radii = [(3., 4.2), (4.2, 6.2), (6.2, 10.2)]
radii = [(3., 4.2), (4.4, 6.4), (6.6, 10.2)]
wedge_start = 12.
wedge_width = 22.
tmp = np.arange(wedge_start, 90. - wedge_start, wedge_width)

theta_width = 2.*np.pi*wedge_width/360.

# This is with the horizontal targets
# theta_starts_deg = np.r_[np.r_[tmp, 90. - wedge_width/2.],
#                                tmp + 90.]
theta_starts_deg = np.r_[tmp, tmp + 90.]

# theta_starts_deg = np.r_[theta_starts_deg, -1.*theta_starts_deg]
theta_starts = 2.*np.pi*theta_starts_deg/360.

code_map = OrderedDict()
for i_th, thetas in enumerate(theta_starts):
    angle_bit = i_th * '0' + '1' + (len(theta_starts) - (i_th + 1)) * '0'
    for i_ecc in range(len(radii)):
        ecc_bit = i_ecc * '0' + '1' + (len(radii) - (i_ecc + 1)) * '0'
        index = i_th * len(radii) + i_ecc
        bin_code = "{:s}{:s}".format(ecc_bit, angle_bit)
        code_map[bin_code] = int(bin_code, 2)

# PLOT STUFF
polax_size = 4  # inches square

global th_plot, r_plot, w_plot, b_plot, val_plot


def polax_bardata_reset():
    global th_plot, r_plot, w_plot, b_plot, val_plot
    th_plot = []  # these are the starting angles, NOT centers
    r_plot = []
    w_plot = []
    b_plot = []
    c_plot = []
    val_plot = []


def polax_bardata_append(hemi_bit, ecc_ind, ang_ind, val):
    '''Modifies X_plot-arrays from outer scope...dangerous...
    '''
    global th_plot, r_plot, w_plot, b_plot, val_plot

    # NB FLIPPED! Left hemisphere gets right hemifield!
    wedge = [2*np.pi - (theta_starts[ang_ind] + theta_width),
             theta_starts[ang_ind]]
    # wedge = [theta_starts[ang_ind],
    #          2*np.pi - (theta_starts[ang_ind] + theta_width)]

    th_plot.append(wedge[hemi_bit])
    w_plot.append(theta_width)
    r_plot.append(np.diff(radii[ecc_ind])[0])
    b_plot.append(radii[ecc_ind][0])
    val_plot.append(val)


def polax_bardata_setcols(ax, cmap=cm.hot, normalizer=1.):
    global th_plot, r_plot, w_plot, b_plot, val_plot
    ax.set_theta_zero_location("N")
    bars = ax.bar(th_plot, r_plot, width=w_plot, bottom=b_plot)
    for val, bar in zip(val_plot, bars):
        cval = polax_get_colour(val, cmap=cmap, normalizer=normalizer)
        bar.set_facecolor(cval)


def polax_get_colour(val, cmap=cm.hot, normalizer=1.):
    cval = val / normalizer if val is not None else 0.0
    return cmap(cval)

# out_dict = dict(radii=radii, wedge_start=wedge_start, wedge_width=wedge_width,
#                 theta_width=theta_width, theta_starts_deg=list(theta_starts_deg),
#                 theta_starts=list(theta_starts))
# with open('polar_plot_data.json', 'wt') as fp:
#     json.dump(out_dict, fp)


def create_LUT():
    LUT = """#$Id: Retinotopy Annotations
    # No. Label Name:           R   G   B   A
    0   None                   255  255  255   255
    """

    c_plot = []
    for ii, th_start in enumerate(theta_starts):
        color = cm.jet((np.abs(th_start) + theta_width/2)/np.pi)
        h, s, v = colorsys.rgb_to_hsv(*color[:-1])
        for jj, rlim in enumerate(radii):

            v = 1 - (jj / len(radii)) / 2
            col = colorsys.hsv_to_rgb(h, s, v)
            c_plot.extend([col, col])

    for ii, (col, bin_code) in enumerate(zip(c_plot[::2], code_map.keys())):
        r, g, b = np.round(np.array(col)*255)
    #     cp = (ii % len(radii)) * 5
        LUT += "{ind:d}\tRM{code:d}      {r:.0f}  {g:.0f}  {b:.0f}  {t:.0f}\n".\
            format(ind=ii + 1, code=code_map[bin_code], r=r, g=g, b=b,
                   t=0)

    with open('RMLUT.txt', 'wt') as fp:
        fp.write(LUT)


def get_ecc_ang_inds(lab_name, n_radii=len(radii),
                     n_angles=len(theta_starts)):
    code_int = int(lab_name.split('-')[0].split('RM')[1])
    bin_code = bin(code_int)[2:]
    bin_code = (n_radii + n_angles - len(bin_code)) * '0' + bin_code  # pad
    ups = [m.start() for m in re.finditer('1', bin_code)]
    ecc_ind = ups[0]
    ang_ind = ups[1] - n_radii
    return ecc_ind, ang_ind


def get_RM_labels(subject, regions=['V1', 'V2', 'V3'], hemis=['lh', 'rh'],
                  subjects_dir=None, verbose=None):
    labels = OrderedDict()
    for reg in regions:
        labels[reg] = OrderedDict()
        for hemi in hemis:
            labels[reg][hemi] = read_labels_from_annot(
                subject, parc='RM.{:s}'.format(reg), subjects_dir=subjects_dir,
                hemi=hemi, verbose=verbose, regexp='RM')
    return labels


def get_mag_scaling_factor(info):
    meg_picks = pick_types(info, meg=True, eeg=False)

    meg_info = pick_info(info, meg_picks)
    meg_picks = pick_types(meg_info, meg=True, eeg=False)
    mag_picks = pick_types(meg_info, meg='mag', eeg=False)
    grad_picks = pick_types(meg_info, meg='grad', eeg=False)

    mag_scale = 'auto'
    coil_scale, mag_scale = mne_get_coil_scale(
        meg_picks, mag_picks, grad_picks, mag_scale, meg_info)

    return(coil_scale, mag_scale)


def prepare_gain(fwd, ch_type='meg', exclude=[], verbose=None):
    """Convert gain matrix to surface orientation

    NB: Assumes fwd is a surface source space

    Parameters
    ----------
    fwd : Forward
        The forward operator, which must be defined on a discrete
        source space.
    ch_type : 'meg' | grad' | 'mag' | 'eeg'
        The type of sensors to use.
    exclude : list of string | str
        List of channels to exclude. If empty do not exclude any (default).
        If 'bads', exclude channels in fwd['info']['bads'].
    verbose : bool, str, int, or None
        If not None, override default verbose level (see :func:`mne.verbose`
        and :ref:`Logging documentation <tut_logging>` for more).

    Returns
    -------
    fwd : Forward
        The prepared forward operator (including gain matrix) in surf_ori
    """
    # check strings
    if ch_type not in ['eeg', 'grad', 'mag', 'meg']:
        raise ValueError("ch_type should be 'eeg', 'meg', mag' or "
                         "'grad (got %s)" % ch_type)

    # check forward
    if is_fixed_orient(fwd, orig=True):
        raise ValueError('fwd should must be computed with free orientation')

    # bit of a hack to allow all MEG channels to be used
    if ch_type == 'meg':
        ch_type = True

    # limit forward (this will make a copy of the data for us)
    # NB because of the copy, the variable fwd is now a local variable
    # -> Python's scope rules mean that the fwd in the calling scope
    # will then remain unchanged!! Wow...
    if ch_type == 'eeg':
        fwd = pick_types_forward(fwd, meg=False, eeg=True, exclude=exclude)
    else:
        fwd = pick_types_forward(fwd, meg=ch_type, eeg=False, exclude=exclude)

    # This operates in-place
    convert_forward_solution(fwd, surf_ori=True, force_fixed=False,
                             copy=False, verbose=False)
    if not fwd['surf_ori'] or is_fixed_orient(fwd):
        raise RuntimeError('Error converting solution, please notify '
                           'mne-python developers')
    return fwd


def pick_gain_vertices(fwd, label):
    """Find the indices to the gain matrix corresponding to label vertices

    Parameters
    ----------
    fwd : Forward
        The forward operator.
    labels : Label | BiHemiLabel
        The patch of cortex in question.
    """
    if fwd['src'][0]['subject_his_id'] != label.subject:
        raise(RuntimeError('Gain and label subjects do not match.'))

    # NB strong assumption (_check!_): hemispheres in gain are stacked L-to-R!
    gain_vertices = fwd['src'][0]['vertno'], fwd['src'][1]['vertno']

    if label.hemi == 'lh':
        lab_verts = [label.vertices, []]
    elif label.hemi == 'rh':
        lab_verts = [[], label.vertices]
    elif label.hemi == 'both':  # BiHemiLabel!
        lab_verts = [label.lh.vertices, label.rh.vertices]

    gain_idx = []
    for hemi_idx, verts in enumerate(lab_verts):
        vert_idx = np.where(np.in1d(gain_vertices[hemi_idx],
                                    verts))[0]
        # NB _check!_ that column stacking is L-R!
        gain_idx.extend(vert_idx + hemi_idx * len(gain_vertices[0]))

    n_locations = len(gain_idx)
    assert(n_locations == len(lab_verts[0]) + len(lab_verts[1]))

    return gain_idx


def patch_sensitivity(fwd, labels):
    """Compute cancellation indices for cortical sources.

    Parameters
    ----------
    fwd : Forward
        The forward operator, which must be in surface orientation.
    labels : Label | list of Label
        Simultaneously active cortical patches.

    Returns
    -------
    sensitivity : dict
        ci -> np.array of cancellation index for each label
        tot_sens -> np.array of vector sum-sensitivity for each label
        w_sens -> np.array of average sensitivity for each label,
                  weighted by the amount cancellation due to anatomy
        sigvecs -> list (of len(labels)) of _total_ signal vectors
    """
    if not fwd['surf_ori']:
        raise RuntimeError('Forward operator must be in surface orientation.')

    if not isinstance(labels, list):
        labels = [labels]

    gain = fwd['sol']['data']
    # # NB strong assumption (_check!_): hemispheres in gain are stacked L-to-R!
    # gain_vertices = fwd['src'][0]['vertno'], fwd['src'][1]['vertno']

    n_sensors, n_dipoles = gain.shape

    alpha = np.empty(len(labels))
    beta = np.empty(len(labels))
    w_sens = np.empty(len(labels))
    mean_sens = np.empty(len(labels))
    ci = np.empty(len(labels))
    sigvec = np.empty((n_sensors, len(labels)), dtype=np.float)

    for nl, lab in enumerate(labels):
        gain_idx = pick_gain_vertices(fwd, lab)
        n_locations = len(gain_idx)

        a_vec = np.zeros(n_sensors)  # alpha: norm of vector sums
        b_vec = np.empty(n_locations)  # beta: sum of norms
        for n, gix in enumerate(gain_idx):
            gz = gain[:, 3 * gix + 2]  # the normal component
            a_vec += gz  # cancellation occurs
            b_vec[n] = np.linalg.norm(gz)

        alpha[nl] = np.linalg.norm(a_vec)
        beta[nl] = np.sum(b_vec)
        mean_sens[nl] = np.mean(b_vec)  # if whole patch active simult.
        sigvec[:, nl] = a_vec

    # These are elementwise operations on numpy arrays
    ci = 1 - alpha/beta
    w_sens = mean_sens * (alpha/beta)  # could have done this two lines up...

    return dict(ci=ci, tot_sens=alpha, w_sens=w_sens, mean_sens=mean_sens,
                sigvec=sigvec)


def find_labels_in_list(lablist, labname):
    '''Return a list of labels matching a name.'''
    labels = [l for l in lablist if labname in l.name]
    if len(labels) == 0:
        raise RuntimeError('No such label name found: {}'.format(labname))
    return labels


def get_2D_connectivity_matrix_value(Cmat, reflab, trglab):
    n_hemi = Cmat.shape[0] // 2  # square

    r_ecc_ind, r_ang_ind = get_ecc_ang_inds(reflab.name)
    ref_hidx = 0 if reflab.hemi == 'lh' else 1

    ecc_ind, ang_ind = get_ecc_ang_inds(trglab.name)
    trg_hidx = 0 if trglab.hemi == 'lh' else 1

    ridx = ref_hidx * n_hemi + len(radii)*r_ang_ind + r_ecc_ind
    cidx = trg_hidx * n_hemi + len(radii)*ang_ind + ecc_ind

    return Cmat[ridx, cidx]


def _stc_from_labels(labels, amp=1e-10, tmin=0.0, tstep=0.001):
    """Create an stc with timepoints = len(labels)
    """
    if not isinstance(labels, list):
        if not isinstance(labels, (Label, BiHemiLabel)):
            raise ValueError('labels must be a list of Labels')
        labels = [labels]

    vertices = []
    all_vertices = [np.array([], dtype=np.int), np.array([], dtype=np.int)]
    for label in labels:
        if label.hemi == 'both':
            vertices.append([label.lh.vertices, label.rh.vertices])
            all_vertices[0] = np.r_[all_vertices[0], vertices[-1][0]]
            all_vertices[1] = np.r_[all_vertices[1], vertices[-1][1]]
        elif label.hemi == 'lh':
            vertices.append([label.vertices, np.array([], dtype=np.int)])
            all_vertices[0] = np.r_[all_vertices[0], vertices[-1][0]]
        elif label.hemi == 'rh':
            vertices.append([np.array([], dtype=np.int), label.vertices])
            all_vertices[1] = np.r_[all_vertices[1], vertices[-1][1]]

    all_vertices = [np.unique(v) for v in all_vertices]
    n_src = sum([len(v) for v in all_vertices])

    data = np.zeros((n_src, len(labels)))
    for idt, verts in enumerate(vertices):
        for hemi in [0, 1]:
            idx = np.where(np.in1d(all_vertices[hemi], verts[hemi]))[0]
            data[idx + hemi * len(all_vertices[0]), idt] = amp

    stc = SourceEstimate(data, vertices=all_vertices, tmin=tmin, tstep=tstep)

    return stc


def plot_stc_topomap(fwd_fixed, stc, info, cov=None, snr=np.inf,
                     axes=None, plot_params=None, title=None):
    if plot_params is None:
        plot_params = dict(times=stc.times, ch_type='mag',
                           outlines='skirt', colorbar=False)
    if axes is None:
        fig, axes = plt.subplots(1, len(stc.times))

    if title is None:
        title = 'Field of label %01d'

    evoked = simulate_evoked(fwd_fixed, stc, info, cov=cov, snr=snr)
    evoked.plot_topomap(time_format=title, axes=axes, **plot_params)


def plot_region_interaction_topomap(lablist, fwd_fixed, info, axes=None,
                                    fig=None):
    if len(lablist) != 2:
        raise RuntimeError('Provide 2 labels as a list')
    lablist.append(lablist[0] + lablist[1])

    if fig is None:
        fig = plt.figure()
    if axes is None or len(axes) != 3:
        axes = [fig.add_subplot(1, len(lablist), ii + 1) for
                ii in range(len(lablist))]

    stc = _stc_from_labels(lablist, tmin=1e-3)
    # for n_lab, label in enumerate(lablist):
        # stc = stc_from_label(label, tmin=(n_lab + 1) * 1e-3, tstep=1e-3)
    plot_stc_topomap(fwd_fixed, stc, info, axes=axes)
