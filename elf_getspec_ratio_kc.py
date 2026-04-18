"""
GOAL: Produce energy-time spectrograms of the energy-flux for electrons and ions
      and their omni-directional flux, perpendicular flux, precipitating flux (i.e.,
      parallel if direction faces north or antiparallel if direction faces south), and
      ratio of precipitating / perpendicular, as well as change in magnetic latitude
      as a function of universal time (UT).

ORIGINAL EXAMPLE PROVIDED BY: Colin Wilkins (colinwilkins@ucla.edu) - 2024-06-12 (YYYY-MM-DD)

AUTHORS: Kai Chen, Drs. Mei-Ching Fok, Suk-Bin Kang, and Cristian Ferradas from NASA GSFC
CONTRIBUTORS: Colin Wilkins, Anton Artemyev, and Dr. Vassilis Angeolopoulos from UCLA's ELFIN team
STARTED: 2024-06-17 (IDL), Python conversion: 2026-04-17

DEPENDENCIES:
    pip install pyspedas pytplot numpy scipy
    geopack is required for tt89 / ttrace2equator:
        pip install geopack
"""

import os
import numpy as np
from scipy.interpolate import interp1d

import pyspedas
from pyspedas import cotrans, tt89, xyz_to_polar, tvector_rotate, tvectot
from pyspedas.utilities.time_double import time_double
from pyspedas.utilities.time_string import time_string

import pytplot
from pytplot import (
    get_data, store_data, options, ylim, zlim,
    tplot, tplot_options, del_data
)

# pyspedas.elfin exposes these two as the primary data loaders
from pyspedas.elfin import epd as elf_load_epd   # replaces elf_load_epd + elf_getspec
from pyspedas.elfin import state as elf_load_state

# ELFIN-specific helper: compute MLT, L-shell, and magnetic latitude under
# dipole configuration — equivalent to IDL's elf_mlt_l_lat.
# If this import path changes across pyspedas versions, adjust accordingly.
try:
    from pyspedas.elfin.epd.elf_mlt_l_lat import elf_mlt_l_lat
except ImportError:
    elf_mlt_l_lat = None  # handled below where it is used

# Optional: shadow/sunlight bar — may not exist in all pyspedas versions
try:
    from pyspedas.elfin.epd.elf_load_sun_shadow_bar import elf_load_sun_shadow_bar
except ImportError:
    elf_load_sun_shadow_bar = None


# =============================================================================
# Configuration
# Look at data_availability from /ela or /elb online to determine these values
# =============================================================================

localdir = '/'                    # directory to store spectrogram plots
probe    = 'a'                    # 'a' for ELFIN-A or 'b' for ELFIN-B (ELFIN-STAR)
myspecies  = 'e'                  # 'e' for electron or 'i' for ion
mydatatype = 'pef'                # 'pef' for electron or 'pif' for ion
direction  = 'South-Descending'   # CHECK el?_epd?_all.csv for 'Descending' or 'Ascending'
tstart = ['2022-07-18/00:00:00']
tend   = ['2022-07-19/00:00:00']

# ---------------------------------------------------------------------------
# Investigation examples of the 2022-07-18 to 2022-07-20 geomagnetic storm
# Uncomment/comment each block as needed (same logic as the IDL script)
# ---------------------------------------------------------------------------

# localdir = '/data/share/elfin/'

# # Both ELFIN-A&B ions and ELFIN-A electron 24-hr plots
# probe      = 'a'
# myspecies  = 'i'
# mydatatype = 'pif'
# tstart = ['2022-07-18/00:00:00', '2022-07-19/00:00:00', '2022-07-20/00:00:00']
# tend   = ['2022-07-19/00:00:00', '2022-07-20/00:00:00', '2022-07-21/00:00:00']

# # ELFIN-B electron 24-hr plots
# probe      = 'b'
# myspecies  = 'e'
# mydatatype = 'pef'
# tstart = ['2022-07-19/00:00:00', '2022-07-20/00:00:00']
# tend   = ['2022-07-20/00:00:00', '2022-07-21/00:00:00']

# # ELFIN-A specific times — electrons and ions share the same pass windows
# probe      = 'a'
# myspecies  = 'i'
# mydatatype = 'pif'
# direction  = 'south'
# tstart = ['2022-07-18/11:29:57', '2022-07-18/12:46:48', '2022-07-18/23:42:26',
#           '2022-07-19/10:21:23', '2022-07-19/11:38:19', '2022-07-19/13:22:20',
#           '2022-07-20/05:57:38']
# tend   = ['2022-07-18/11:35:58', '2022-07-18/13:06:42', '2022-07-18/23:48:28',
#           '2022-07-19/10:27:24', '2022-07-19/11:58:20', '2022-07-19/13:28:20',
#           '2022-07-20/06:17:23']

# # ELFIN-B specific times — electrons (south-pointing)
# probe      = 'b'
# myspecies  = 'e'
# mydatatype = 'pef'
# direction  = 'south'
# tstart = ['2022-07-19/10:06:23', '2022-07-19/13:08:34', '2022-07-20/13:32:32']
# tend   = ['2022-07-19/10:09:13', '2022-07-19/13:13:42', '2022-07-20/13:37:31']

# # ELFIN-B specific times — ions (north-pointing)
# probe      = 'b'
# myspecies  = 'i'
# mydatatype = 'pif'
# direction  = 'north'
# tstart = ['2022-07-18/11:30:26', '2022-07-18/23:43:38',
#           '2022-07-19/10:21:52', '2022-07-19/13:22:20']
# tend   = ['2022-07-18/11:35:50', '2022-07-18/23:48:28',
#           '2022-07-19/10:27:24', '2022-07-19/13:27:05']

# # ELFIN-B specific times — ions (south-pointing)
# probe      = 'b'
# myspecies  = 'i'
# mydatatype = 'pif'
# direction  = 'south'
# tstart = ['2022-07-19/10:06:23', '2022-07-19/13:07:41', '2022-07-20/13:32:32']
# tend   = ['2022-07-19/10:09:24', '2022-07-19/13:13:42', '2022-07-20/13:37:31']


# =============================================================================
# Validate equal number of start and end times (equivalent to IDL stop)
# =============================================================================
if len(tstart) != len(tend):
    raise ValueError("tstart and tend must contain an equal number of time strings.")


# =============================================================================
# Main loop — automates plotting of multiple spectrograms
# =============================================================================
for k in range(len(tstart)):

    # -------------------------------------------------------------------------
    # Colin's contribution from elf_getspec_example_cw.pro
    # -------------------------------------------------------------------------
    time2plot   = [tstart[k], tend[k]]
    timeduration = time_double(tend[k]) - time_double(tstart[k])

    mytype = 'eflux'

    # Load position / state data
    elf_load_state(probe=probe, trange=time2plot)

    # Load EPD data and compute spectra.
    # In PySPEDAS, elf_load_epd() handles both loading *and* the getspec step
    # (equivalent to IDL's elf_load_epd + elf_getspec /get3Dspec).
    # The first call with type_='raw' mirrors the IDL two-pass approach; it is
    # kept here for completeness but the eflux pass below is what drives the plots.
    elf_load_epd(probe=probe, datatype=mydatatype, level='l2',
                 type_='raw', trange=time2plot)
    elf_load_epd(probe=probe, datatype=mydatatype, level='l2',
                 type_=mytype, trange=time2plot)

    # -------------------------------------------------------------------------
    # Get position data and perform coordinate transforms
    # IDL: cotrans with /GEI2GSE, /GSE2GSM, /GSM2SM, /GEI2GEO, /GEO2MAG
    # -------------------------------------------------------------------------
    cotrans('el' + probe + '_pos_gei', 'el' + probe + '_pos_gse',
            coord_in='gei', coord_out='gse')
    cotrans('el' + probe + '_pos_gse', 'el' + probe + '_pos_gsm',
            coord_in='gse', coord_out='gsm')
    # Note: reaching SM from GEI requires the intermediate GSE→GSM step first
    cotrans('el' + probe + '_pos_gsm', 'el' + probe + '_pos_sm',
            coord_in='gsm', coord_out='sm')
    cotrans('el' + probe + '_pos_gei', 'el' + probe + '_pos_geo',
            coord_in='gei', coord_out='geo')
    cotrans('el' + probe + '_pos_geo', 'el' + probe + '_pos_mag',
            coord_in='geo', coord_out='mag')

    # -------------------------------------------------------------------------
    # Calculate IGRF magnetic field
    # -------------------------------------------------------------------------
    datgsm = get_data('el' + probe + '_pos_gsm')

    # Subsample to one point per minute for efficiency (every 60th sample),
    # then interpolate back to the full cadence — mirrors the IDL quick_run path.
    store_data('el' + probe + '_pos_gsm_mins', data={
        'x': datgsm.times[::60],
        'y': datgsm.y[::60, :]
    })
    tt89('el' + probe + '_pos_gsm_mins', igrf_only=True,
         newname='el' + probe + '_bt89_gsm_mins', period=1.)

    # Interpolate minute-cadence IGRF back to the full time array
    gsm_mins = get_data('el' + probe + '_bt89_gsm_mins')
    bt89_interp = np.zeros((len(datgsm.times), 3))
    for comp in range(3):
        f_interp = interp1d(gsm_mins.times, gsm_mins.y[:, comp],
                            bounds_error=False, fill_value='extrapolate')
        bt89_interp[:, comp] = f_interp(datgsm.times)
    store_data('el' + probe + '_bt89_gsm',
               data={'x': datgsm.times, 'y': bt89_interp})
    del_data('*_mins')

    # Transform IGRF to SM coordinates
    state_pos_sm = get_data('el' + probe + '_pos_sm')
    cotrans('el' + probe + '_bt89_gsm', 'el' + probe + '_bt89_sm',
            coord_in='gsm', coord_out='sm')

    # Convert SM position to spherical (co-latitude) coordinates
    xyz_to_polar('el' + probe + '_pos_sm', co_latitude=True)
    pos_sm_th  = get_data('el' + probe + '_pos_sm_th')
    pos_sm_phi = get_data('el' + probe + '_pos_sm_phi')

    csth = np.cos(np.pi * pos_sm_th.y  / 180.)
    csph = np.cos(np.pi * pos_sm_phi.y / 180.)
    snth = np.sin(np.pi * pos_sm_th.y  / 180.)
    snph = np.sin(np.pi * pos_sm_phi.y / 180.)

    # Build the Cartesian → spherical rotation matrix (3 × 3 per time step)
    n = len(pos_sm_th.times)
    rot2rthph = np.zeros((n, 3, 3))
    rot2rthph[:, 0, 0] =  snth * csph
    rot2rthph[:, 0, 1] =  csth * csph
    rot2rthph[:, 0, 2] = -snph
    rot2rthph[:, 1, 0] =  snth * snph
    rot2rthph[:, 1, 1] =  csth * snph
    rot2rthph[:, 1, 2] =  csph
    rot2rthph[:, 2, 0] =  csth
    rot2rthph[:, 2, 1] = -snth
    rot2rthph[:, 2, 2] =  0. * csth
    store_data('rot2rthph', data={'x': pos_sm_th.times, 'y': rot2rthph})
    tvector_rotate('rot2rthph', 'el' + probe + '_bt89_sm',
                   newname='el' + probe + '_bt89_sm_sph')

    # Rotate spherical (r, θ, φ) → NED (North, East, Down)
    rotSMSPH2NED = np.zeros((n, 3, 3))
    rotSMSPH2NED[:, 0, 2] = -1.  # N  = -r_hat
    rotSMSPH2NED[:, 1, 0] = -1.  # E  = -theta_hat
    rotSMSPH2NED[:, 2, 1] =  1.  # D  = +phi_hat
    store_data('rotSMSPH2NED', data={'x': pos_sm_th.times, 'y': rotSMSPH2NED})
    tvector_rotate('rotSMSPH2NED', 'el' + probe + '_bt89_sm_sph',
                   newname='el' + probe + '_bt89_sm_NED')
    tvectot('el' + probe + '_bt89_sm_NED', newname='el' + probe + '_bt89_sm_NEDT')

    options('el' + probe + '_bt89_sm_NED',  'legend_names', ['N', 'E', 'D'])
    options('el' + probe + '_bt89_sm_NEDT', 'legend_names', ['N', 'E', 'D', 'T'])
    options('el' + probe + '_bt89_sm_NED',  'ytitle', 'IGRF')
    options('el' + probe + '_bt89_sm_NED',  'ysubtitle', '[nT]')
    options('el' + probe + '_bt89_sm_NEDT', 'ytitle', 'IGRF')
    options('el' + probe + '_bt89_sm_NEDT', 'ysubtitle', '[nT]')

    # -------------------------------------------------------------------------
    # MLT and LAT — dipole approximation
    # elf_mlt_l_lat computes magnetic local time, L-shell, and magnetic latitude
    # under a simple dipole field configuration.
    # -------------------------------------------------------------------------
    if elf_mlt_l_lat is not None:
        MLT0, L0, lat0 = elf_mlt_l_lat('el' + probe + '_pos_sm')
    else:
        raise ImportError(
            "elf_mlt_l_lat not found. Check your pyspedas installation "
            "or update the import path at the top of this script."
        )

    elfin_pos = get_data('el' + probe + '_pos_sm')
    store_data('el' + probe + '_MLT_dip',  data={'x': elfin_pos.times, 'y': MLT0})
    store_data('el' + probe + '_L_dip',    data={'x': elfin_pos.times, 'y': L0})
    store_data('el' + probe + '_MLAT_dip', data={'x': elfin_pos.times,
                                                  'y': lat0 * 180. / np.pi})
    options('el' + probe + '_MLT_dip',  'ytitle', 'dip')
    options('el' + probe + '_L_dip',    'ytitle', 'dip')
    options('el' + probe + '_MLAT_dip', 'ytitle', 'dip')

    # Median satellite altitude above Earth's surface (km)
    alt = np.median(np.sqrt(np.sum(elfin_pos.y ** 2, axis=1))) - 6371.

    # Geographic longitude (GLON)
    dat_geo  = get_data('el' + probe + '_pos_geo')
    phi_geo  = np.degrees(np.arctan2(dat_geo.y[:, 1], dat_geo.y[:, 0]))
    phi_geo[phi_geo < 0] += 360.
    store_data('el' + probe + '_GLON', data={'x': dat_geo.times, 'y': phi_geo})
    options('el' + probe + '_GLON', 'ytitle', 'GLON (east)')

    # -------------------------------------------------------------------------
    # MLT in IGRF — trace field lines to the magnetic equator
    # -------------------------------------------------------------------------
    Rem = 6371.0   # Earth mean radius [km]
    Re  = 6378.0   # Earth equatorial radius [km]

    cotrans('el' + probe + '_pos_gei', 'elx_pos_gse',
            coord_in='gei', coord_out='gse')
    cotrans('elx_pos_gse', 'elx_pos_gsm',
            coord_in='gse', coord_out='gsm')
    datgsm_d = get_data('elx_pos_gsm')

    store_data('elx_pos_gsm_mins', data={
        'x': datgsm_d.times[::60],
        'y': datgsm_d.y[::60, :]
    })
    tt89('elx_pos_gsm_mins', igrf_only=True,
         newname='elx_bigrf_gsm_mins', period=0.1)

    # Trace to the magnetic equator (requires geopack)
    pyspedas.ttrace2equator(
        'elx_pos_gsm_mins',
        external_model='none', internal_model='igrf',
        km=True, in_coord='gsm', out_coord='gsm',
        rlim=100. * Rem
    )

    cotrans('elx_pos_gsm_mins_foot', 'elx_pos_sm_mins_foot',
            coord_in='gsm', coord_out='sm')
    xyz_to_polar('elx_pos_sm_mins_foot', co_latitude=True)

    elx_pos_sm_foot    = get_data('elx_pos_sm_mins_foot')
    foot_mag           = get_data('elx_pos_sm_mins_foot_mag')
    foot_th            = get_data('elx_pos_sm_mins_foot_th')
    elx_bigrf_gsm_mins = get_data('elx_bigrf_gsm_mins')
    elx_pos_gsm_mins_d = get_data('elx_pos_gsm_mins')

    # L-shell in IGRF
    Ligrf = (foot_mag.y / Rem) / np.sin(foot_th.y * np.pi / 180.) ** 2

    # Magnetic hemisphere from sign of radial IGRF component
    Br_tmp   = np.sum(elx_bigrf_gsm_mins.y * elx_pos_gsm_mins_d.y, axis=1)
    hemisphere = np.sign(-Br_tmp)

    # MLAT at the satellite footpoint
    r_ift_dip = 1. + 100. / Rem
    MLAT_mins = (180. / np.pi) * np.arccos(
        np.sqrt(Rem * r_ift_dip / foot_mag.y) *
        np.sin(foot_th.y * np.pi / 180.)
    ) * hemisphere

    # Interpolate minute-cadence quantities to full time array
    f_mlat = interp1d(elx_pos_sm_foot.times, MLAT_mins,
                      bounds_error=False, fill_value='extrapolate')
    store_data('el' + probe + '_MLAT_igrf',
               data={'x': datgsm_d.times, 'y': f_mlat(datgsm_d.times)})

    L_mins = np.sqrt(np.sum(elx_pos_sm_foot.y ** 2, axis=1)) / Re
    f_L = interp1d(elx_pos_sm_foot.times, L_mins,
                   bounds_error=False, fill_value='extrapolate')
    store_data('el' + probe + '_L_igrf',
               data={'x': datgsm_d.times, 'y': f_L(datgsm_d.times)})

    MLT0_igrf, _, _ = elf_mlt_l_lat('elx_pos_sm_mins_foot')
    sm_mins_d = get_data('elx_pos_sm_mins_foot')
    f_mlt = interp1d(sm_mins_d.times, MLT0_igrf,
                     bounds_error=False, fill_value='extrapolate')
    store_data('el' + probe + '_MLT_igrf',
               data={'x': datgsm_d.times, 'y': f_mlt(datgsm_d.times)})

    # Clean up intermediate tplot variables
    del_data(['*_mins', 'elx_*', 'Ligrf', 'MLAT', 'rot2rthph', 'rotSMSPH2NED'])

    options('el' + probe + '_L_igrf',    'ytitle', 'L-igrf')
    options('el' + probe + '_MLAT_igrf', 'ytitle', 'MLAT-igrf')
    options('el' + probe + '_MLT_igrf',  'ytitle', 'MLT-igrf')

    # Adjust character/label size depending on plot type
    char_size = 1.0 if timeduration == 86400 else 1.25
    for var in ['el' + probe + '_MLT_dip',  'el' + probe + '_L_dip',
                'el' + probe + '_MLAT_dip', 'el' + probe + '_L_igrf',
                'el' + probe + '_MLAT_igrf','el' + probe + '_MLT_igrf']:
        options(var, 'char_size', char_size)

    # -------------------------------------------------------------------------
    # Sunlight / shadow bar
    # -------------------------------------------------------------------------
    if elf_load_sun_shadow_bar is not None:
        elf_load_sun_shadow_bar(tplotname='el' + probe + '_pos_gse')
        options('shadow_bar',
                thick=5.5, xstyle=4, ystyle=4, yrange=[-0.1, 0.1],
                ytitle='', panel_size=0.1, ztitle='')
        options('sun_bar',
                thick=5.5, xstyle=4, ystyle=4, yrange=[-0.1, 0.1],
                ytitle='', panel_size=0.1, ztitle='')
        store_data('sunlight_bar', data=['sun_bar', 'shadow_bar'])
        options('sunlight_bar', 'panel_size', 0.1)
        options('sunlight_bar', 'yrange', [-0.1, 0.1])
    else:
        # Placeholder: create an empty bar so downstream tplot calls still work
        store_data('sunlight_bar', data={
            'x': datgsm_d.times,
            'y': np.zeros(len(datgsm_d.times))
        })
        options('sunlight_bar', 'panel_size', 0.1)

    # =========================================================================
    # Build spectrogram tplot variables
    # =========================================================================
    prefix    = 'el' + probe + '_' + mydatatype + '_en_spec2plot_'
    omni_plot = prefix + 'omni'
    anti_plot = prefix + 'anti'
    para_plot = prefix + 'para'
    perp_plot = prefix + 'perp'

    anti_data = get_data(anti_plot)
    para_data = get_data(para_plot)
    perp_data = get_data(perp_plot)

    # Sanity check: all arrays must have the same shape
    if not (anti_data.y.shape == para_data.y.shape == perp_data.y.shape):
        raise ValueError(
            f"Mismatched spectrogram array sizes for interval {tstart[k]}–{tend[k]}. "
            "Consider interpolating to a common time grid."
        )

    if timeduration == 86400:
        # -----------------------------------------------------------------
        # 24-hr plot: precipitating flux = anti-parallel + parallel combined
        # -----------------------------------------------------------------
        antiandpara = anti_data.y + para_data.y

        store_data(prefix + 'anti_and_para', data={
            'x': para_data.times, 'y': antiandpara, 'v': para_data.v
        })
        anti_and_para_plot = prefix + 'anti_and_para'
        anti_and_para_data = get_data(anti_and_para_plot)

        options(anti_and_para_plot, 'spec', 1)
        ylim(anti_and_para_plot, 50, 5e3, 1)
        zlim(anti_and_para_plot, 1e4, 1e9, 1)

        # Ratio = (anti + para) / perp  — NaN where both are zero
        ratio_data = np.full_like(antiandpara, np.nan)
        mask = (antiandpara != 0) | (perp_data.y != 0)
        ratio_data[mask] = antiandpara[mask] / perp_data.y[mask]

        store_data(prefix + 'anti+para_ovr_perd', data={
            'x': anti_and_para_data.times,
            'y': ratio_data,
            'v': anti_and_para_data.v
        })
        ratio_plot = prefix + 'anti+para_ovr_perd'

    else:
        # -----------------------------------------------------------------
        # Time-specific plot: precipitating flux depends on pointing direction
        # parallel (para) if north-pointing, anti-parallel (anti) if south-pointing
        # -----------------------------------------------------------------
        if 'north' in direction.lower():
            prec_y, prec_times, prec_v = para_data.y, para_data.times, para_data.v
        else:
            prec_y, prec_times, prec_v = anti_data.y, anti_data.times, anti_data.v

        store_data(prefix + 'prec', data={'x': prec_times, 'y': prec_y, 'v': prec_v})
        prec_plot = prefix + 'prec'
        prec_data = get_data(prec_plot)

        options(prec_plot, 'spec', 1)
        ylim(prec_plot, 50, 5e3, 1)
        zlim(prec_plot, 1e4, 1e9, 1)

        # Ratio = prec / perp  — NaN where both are zero
        ratio_data = np.full_like(prec_data.y, np.nan)
        mask = (prec_data.y != 0) | (perp_data.y != 0)
        ratio_data[mask] = prec_data.y[mask] / perp_data.y[mask]

        store_data(prefix + 'prec_ovr_perd', data={
            'x': prec_data.times, 'y': ratio_data, 'v': prec_data.v
        })
        ratio_plot = prefix + 'prec_ovr_perd'

    # Shared axes for ratio panel
    # EPD-E: 50 keV – 4.5 MeV  |  EPD-I: 50 keV – 300 keV
    options(ratio_plot, 'spec', 1)
    zlim(ratio_plot, 0.02, 2, 1)
    ylim(ratio_plot, 50, 5e3, 1)

    # Z-axis titles
    options(omni_plot, 'ztitle', 'nflux')
    options(perp_plot, 'ztitle', 'nflux')
    if timeduration == 86400:
        options(anti_and_para_plot, 'ztitle', 'eflux')
    else:
        options(prec_plot, 'ztitle', 'eflux')
    options(ratio_plot, 'ztitle', 'ratio')

    # =========================================================================
    # Plot and save PNG
    # =========================================================================
    tplot_options('xmargin', [18, 11])   # left/right margins in character widths
    tplot_options('ymargin', [4, 4])     # bottom/top margins in line units

    date = tstart[k][:10]   # 'YYYY-MM-DD'

    if timeduration == 86400:
        mytitle = (
            f'PRELIMINARY ELFIN-{probe.upper()} EPD-{myspecies.upper()}, '
            f'alt={int(alt)}km, {date}, 24hrs'
        )
    else:
        t0_str = tstart[k][11:16]   # 'HH:MM'
        t1_str = tend[k][11:16]
        mytitle = (
            f'PRELIMINARY ELFIN-{probe.upper()} EPD-{myspecies.upper()}, '
            f'alt={int(alt)}km, {date}, {t0_str} to {t1_str}, {direction}'
        )

    # Bottom-of-panel position labels (GLON, MLAT, MLT, L in both IGRF and dipole)
    varstring = [
        f'el{probe}_GLON',
        f'el{probe}_MLAT_igrf[el{probe}_MLAT_dip]',
        f'el{probe}_MLT_igrf[el{probe}_MLT_dip]',
        f'el{probe}_L_igrf[el{probe}_L_dip]',
    ]

    if timeduration == 86400:
        vars_to_plot = [omni_plot, perp_plot, anti_and_para_plot, ratio_plot,
                        'sunlight_bar', 'el' + probe + '_MLAT_igrf']
    else:
        vars_to_plot = [omni_plot, perp_plot, prec_plot, ratio_plot,
                        'sunlight_bar', 'el' + probe + '_MLAT_igrf']

    # Create output directory: localdir/el{probe}/overplots/YYYY/MM/DD/
    dir_path = os.path.join(localdir, 'el' + probe, 'overplots',
                            tstart[k][:4], tstart[k][5:7], tstart[k][8:10])
    os.makedirs(dir_path, exist_ok=True)

    if timeduration == 86400:
        filename = f'el{probe}_epd{myspecies}_scizone_specplots_{date}'
    else:
        time_tag = tstart[k][11:16].replace(':', '')
        filename = f'el{probe}_epd{myspecies}_scizone_specplots_{date}_{time_tag}'

    tplot(vars_to_plot,
          title=mytitle,
          var_label=varstring,
          save_png=os.path.join(dir_path, filename))

    # =========================================================================
    # Write SM position data to ASCII file
    # Columns: time (ISO string), X_SM, Y_SM, Z_SM  [all in km / Earth radii]
    # Provided for comparison with the CIMI model (Dr. Mei-Ching Fok et al.)
    # =========================================================================
    pos_sm_data = get_data('el' + probe + '_pos_sm')

    if timeduration == 86400:
        dat_filename = f'pos_sm_{date}.dat'
    else:
        time_tag = tstart[k][11:16].replace(':', '')
        dat_filename = f'pos_sm_{date}_{time_tag}.dat'

    dat_filepath = os.path.join(dir_path, dat_filename)

    with open(dat_filepath, 'w') as fh:
        for count in range(len(pos_sm_data.times)):
            t_str = time_string(pos_sm_data.times[count])
            x = pos_sm_data.y[count, 0]
            y = pos_sm_data.y[count, 1]
            z = pos_sm_data.y[count, 2]
            fh.write(f'{t_str}  {x:.6f}  {y:.6f}  {z:.6f}\n')
