; GOAL: Produce energy-time spectrograms of the energy-flux for electrons and ions 
; 		and their omni-directional flux, perpendicular flux, precipitating flux (i.e.,
; 		parallel if north-descending or antiparallel if south-ascending), and 
; 		ratio of precipitating / perpendicular, as well as change in magnetic latitude
; 		as a function of universial time (UT).
; 
; ORIGINAL BY: Colin Wilkins (colinwilkins@ucla.edu) - 2024-06-12 (YYYY-MM-DD)
; 
; AUTHORS: Kai Chen, Drs. Mei-Ching Fok, Suk-Bin Kang, and Cristian Ferradas from NASA GSFC
; CONTRIBUTORS: Colin Wilkins, Anton Artemyev, and Dr. Vassilis Angeolopoulos from UCLA's ELFIN team
; STARTED: 2024-06-17
; ENDED: TBD

pro elf_getspec_ratio_kc

	elf_init
	;!elf.local_data_dir='D:\docs\data\elfin_data\'

	; the following are the presets to investigate available data for the geomagnetic storm between 2022-07-18 and 2022-07-20
	yourdir = '/data/share/elfin/'		; directory to store spectrogram plots

	; both ELFIN-A&B ions and ELFIN-A electron
	;probe='a'			; 'a' for ELFIN-A probe or 'b' for ELFIN-B (ELFIN-STAR) probe
	;myspecies='i'		; 'i' for ion or 'e' for electron
	;mydatatype='pif'		; 'pif' for ion or 'pef' for electron ; be sure to also change when calling elf_load_epd and elf_getspec
	;;direction is not needed
	;direction = 'south'
	;tstart=['2022-07-18/11:29:57']
	;tend=['2022-07-18/11:35:58']

	; only for ELFIN-B electron
	probe='b'
	myspecies='e'
	mydatatype='pef'
	;direction is not needed
	tstart=['2022-07-19/00:00:00','2022-07-20/00:00:00']
	tend=['2022-07-20/00:00:00','2022-07-21/00:00:00']

	; ELFIN-A specific times (electron and ion are the same) (got times from ela_epdi_all.csv and ela_epde_all.csv) (all of them are south)
	;probe='a'
	;myspecies='i'
	;mydatatype='pif'		; be sure to also change when calling elf_load_epd and elf_getspec
	;direction='south'
	;tstart=['2022-07-18/11:29:57', $
	;		'2022-07-18/12:46:48', $
	;		'2022-07-18/23:42:26', $
	;		'2022-07-19/10:21:23', $
	;		'2022-07-19/11:38:19', $
	;		'2022-07-19/13:22:20', $
	;		'2022-07-20/05:57:38']
	;tend = ['2022-07-18/11:35:58', $
	;		'2022-07-18/13:06:42', $
	;		'2022-07-18/23:48:28', $
	;		'2022-07-19/10:27:24', $
	;		'2022-07-19/11:58:20', $
	;		'2022-07-19/13:28:20', $
	;		'2022-07-20/06:17:23']

	; ELFIN-B specific times for electrons (all of them are south)
	;probe='b'
	;myspecies='e'
	;mydatatype='pef'		; be sure to also change when calling elf_load_epd and elf_getspec
	;direction='south'
	;tstart=['2022-07-19/10:06:23', $
	;		'2022-07-19/13:08:34', $
	;		'2022-07-20/13:32:32']
	;tend = ['2022-07-19/10:09:13', $
	;		'2022-07-19/13:13:42', $
	;		'2022-07-20/13:37:31']

	; ELFIN-B specific times for ions (north ones)
	;probe='b'
	;myspecies='i'
	;mydatatype='pif'		; be sure to also change when calling elf_load_epd and elf_getspec
	;direction='north'
	;tstart=['2022-07-18/11:30:26', $
	;		'2022-07-18/23:43:38', $
	;		'2022-07-19/10:21:52', $
	;		'2022-07-19/13:22:20']
	;tend = ['2022-07-18/11:35:50', $
	;		'2022-07-18/23:48:28', $
	;		'2022-07-19/10:27:24', $
	;		'2022-07-19/13:27:05']

	; ELFIN-B specific times for ions (south ones)
	;probe='b'
	;myspecies='i'
	;mydatatype='pif'		; be sure to also change when calling elf_load_epd and elf_getspec
	;direction='south'
	;tstart=['2022-07-19/10:06:23', $
	;		'2022-07-19/13:07:41', $
	;		'2022-07-20/13:32:32']
	;tend = ['2022-07-19/10:09:24', $
	;		'2022-07-19/13:13:42', $
	;		'2022-07-20/13:37:31']

	tstart_size = size(tstart)
	tend_size = size(tend)
	IF(tstart_size[3] NE tend_size[3]) THEN BEGIN
		stop			; stop if tstart and tend are not the same length
	ENDIF

	FOR k=0, tstart_size[3]-1 DO BEGIN

		; Contribution by CW START
		; Contribution by CW START
		time2plot=[tstart[k],tend[k]]
		timeduration=time_double(tend[k])-time_double(tstart[k])
		timespan,tstart[k],timeduration,/seconds

		mytype='eflux'

		elf_load_state,probe=probe,trange=time2plot

		; loads pef/pif data
		elf_load_epd,probe=probe,datatype='pef',type = 'raw',trange=time2plot
		elf_getspec,probe=probe,datatype='pef',type = 'raw',fullspin=fullspin,/get3Dspec ; this simply acts on the data loaded previously

		; reload data in eflux units now (two calls: first reload data, then recompute spectra)
		elf_load_epd,probe=probe,datatype='pef',type = mytype,trange=time2plot
		elf_getspec,probe=probe,datatype='pef',type = mytype,fullspin=fullspin,/get3Dspec ; this simply acts on the data loaded previously
		; Contribution by CW END
		; Contribution by CW END

		; COPY PASTA START
		; COPY PASTA START
		; COPY PASTA FROM spedas/projects/elfin/plots/epde_plot_overviews.pro
		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  ; Get position data
		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  elf_load_state, probes=probe, no_download=no_download
		  get_data, 'el'+probe+'_pos_gei', data=dat_gei
		  cotrans,'el'+probe+'_pos_gei','el'+probe+'_pos_gse',/GEI2GSE
		  cotrans,'el'+probe+'_pos_gse','el'+probe+'_pos_gsm',/GSE2GSM	; note from CW - to get GSM, you need to first transform to GSE
		  cotrans,'el'+probe+'_pos_gsm','el'+probe+'_pos_sm',/GSM2SM 	; in SM
		  cotrans,'el'+probe+'_pos_gei','el'+probe+'_pos_geo',/GEI2GEO
		  cotrans,'el'+probe+'_pos_geo','el'+probe+'_pos_mag',/GEO2MAG

		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  ; Calculate IGRF
		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  threeones=[1,1,1]
		  ; quick_run -> do only every 60th point (i.e. per minute)
		  if keyword_set(quick_run) then begin
			get_data, 'el'+probe+'_pos_gsm', data=datgsm, dlimits=dl, limits=l
			store_data, 'el'+probe+'_pos_gsm_mins', data={x: datgsm.x[0:*:60], y: datgsm.y[0:*:60,*]}, dlimits=dl, limits=l
			tt89,'el'+probe+'_pos_gsm_mins',/igrf_only,newname='el'+probe+'_bt89_gsm_mins',period=1.
			; interpolate the minute-by-minute data back to the full array
			get_data,'el'+probe+'_bt89_gsm_mins',data=gsm_mins, dlimits=dl, limits=l
			store_data,'el'+probe+'_bt89_gsm',data={x: datgsm.x, y: interp(gsm_mins.y[*,*], gsm_mins.x, datgsm.x)},dlimits=dl, limits=l
			; clean up the temporary data
			del_data, '*_mins'
		  endif else begin
			tt89,'el'+probe+'_pos_gsm',/igrf_only,newname='el'+probe+'_bt89_gsm',period=1.
		  endelse

		  get_data, 'el'+probe+'_pos_sm', data=state_pos_sm, dlimits=dl, limits=l
		  ; calculate IGRF in nT
		  cotrans,'el'+probe+'_bt89_gsm','el'+probe+'_bt89_sm',/GSM2SM ; Bfield in SM coords as well
		  xyz_to_polar,'el'+probe+'_pos_sm',/co_latitude
		  get_data,'el'+probe+'_pos_sm_th',data=pos_sm_th;,dlim=myposdlim,lim=myposlim
		  get_data,'el'+probe+'_pos_sm_phi',data=pos_sm_phi
		  csth=cos(!PI*pos_sm_th.y/180.)
		  csph=cos(!PI*pos_sm_phi.y/180.)
		  snth=sin(!PI*pos_sm_th.y/180.)
		  snph=sin(!PI*pos_sm_phi.y/180.)
		  rot2rthph=[[[snth*csph],[csth*csph],[-snph]],[[snth*snph],[csth*snph],[csph]],[[csth],[-snth],[0.*csth]]]
		  store_data,'rot2rthph',data={x:pos_sm_th.x,y:rot2rthph},dlimits=dl, limits=l ;dlim=myposdlim,lim=myposlim
		  tvector_rotate,'rot2rthph','el'+probe+'_bt89_sm',newname='el'+probe+'_bt89_sm_sph'
		  rotSMSPH2NED=[[[snth*0.],[snth*0.],[snth*0.-1.]],[[snth*0.-1.],[snth*0.],[snth*0.]],[[snth*0.],[snth*0.+1.],[snth*0.]]]
		  store_data,'rotSMSPH2NED',data={x:pos_sm_th.x,y:rotSMSPH2NED},dlimits=dl, limits=l;dlim=myposdlim,lim=myposlim
		  tvector_rotate,'rotSMSPH2NED','el'+probe+'_bt89_sm_sph',newname='el'+probe+'_bt89_sm_NED' ; North (-Spherical_theta), East (Spherical_phi), Down (-Spherical_r)
		  tvectot,'el'+probe+'_bt89_sm_NED',newname='el'+probe+'_bt89_sm_NEDT'
		  get_data, 'el'+probe+'_bt89_sm_NEDT', data=d, dlimits=dl, limits=l
		  dl.labels=['N','E','D','T']
		  dl.colors=[60,155,254,1]
		  store_data, 'el'+probe+'_bt89_sm_NEDT', data=d, dlimits=dl, limits=l
		  options,'el'+probe+'_bt89_sm_N*','ytitle', 'IGRF'
		  options,'el'+probe+'_bt89_sm_N*','ysubtitle','[nT]'
		  options,'el'+probe+'_bt89_sm_NED','labels',['N','E','D']
		  options,'el'+probe+'_bt89_sm_NEDT','labels',['N','E','D','T']
		  options,'el'+probe+'_bt89_sm_NEDT','colors',[60,155,254,1]
		  options,'el'+probe+'_bt89_sm_N*','databar',0.

		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  ; Get MLT amd LAT (dipole)
		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  elf_mlt_l_lat,'el'+probe+'_pos_sm',MLT0=MLT0,L0=L0,lat0=lat0 ;;subroutine to calculate mlt,l,mlat under dipole configuration
		  get_data, 'el'+probe+'_pos_sm', data=elfin_pos
		  store_data,'el'+probe+'_MLT_dip',data={x:elfin_pos.x,y:MLT0}
		  store_data,'el'+probe+'_L_dip',data={x:elfin_pos.x,y:L0}
		  store_data,'el'+probe+'_MLAT_dip',data={x:elfin_pos.x,y:lat0*180./!pi}
		  options,'el'+probe+'_MLT_dip',ytitle='dip'
		  options,'el'+probe+'_L_dip',ytitle='dip'
		  options,'el'+probe+'_MLAT_dip',ytitle='dip'
		  ; charsize options was here
		  alt = median(sqrt(elfin_pos.y[*,0]^2 + elfin_pos.y[*,1]^2 + elfin_pos.y[*,2]^2))-6371.

		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  ; GLON
		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  GLON=lat0
		  get_data, 'el'+probe+'_pos_geo', data=dat_geo
		  cart_to_sphere,dat_geo.y[*,0],dat_geo.y[*,1],dat_geo.y[*,2],r_geo,theta_geo,phi_geo
		  gidx=where(phi_geo LT 0, ncnt)
		  if ncnt GT 0 then phi_geo[gidx]=360.+phi_geo[gidx]
		  store_data,'el'+probe+'_GLON',data={x:dat_geo.x,y:phi_geo}
		  options,'el'+probe+'_GLON',ytitle='GLON (east)'
		  ; charsize options was here
		  re=6378.

		  ;;;;;;;;;;;;;;;;;;;;
		  ; MLT IGRF
		  ;;;;;;;;;;;;;;;;;;;;
		  pival=!PI
		  Rem=6371.0 ; Earth mean radius in km
		  cotrans,'el'+probe+'_pos_gei','elx_pos_gse',/GEI2GSE
		  cotrans,'elx_pos_gse','elx_pos_gsm',/GSE2GSM
		  get_data, 'elx_pos_gsm',data=datgsm, dlimits=datgsmdl, limits=datgsml
		  store_data, 'elx_pos_gsm_mins', data={x: datgsm.x[0:*:60], y: datgsm.y[0:*:60,*]}, dlimits=datgsmdl, limits=datgsml
		  tt89,'elx_pos_gsm_mins',/igrf_only,newname='elx_bigrf_gsm_mins',period=0.1; gets IGRF field at ELF location
		  ; find igrf coordinates for satellite, same as for footpoint: Ligrf, MLATigrf, MLTigrf
		  ttrace2equator,'elx_pos_gsm_mins',external_model='none',internal_model='igrf',/km,in_coord='gsm',out_coord='gsm',rlim=100.*Rem ; native is gsm
		  cotrans,'elx_pos_gsm_mins_foot','elx_pos_sm_mins_foot',/GSM2SM ; now in SM
		  get_data,'elx_pos_sm_mins_foot',data=elx_pos_sm_foot
		  xyz_to_polar,'elx_pos_sm_mins_foot',/co_latitude ; get position in rthphi (polar) coords
		  calc," 'Ligrf'=('elx_pos_sm_mins_foot_mag'/Rem)/(sin('elx_pos_sm_mins_foot_th'*pival/180.))^2 " ; uses 1Rem (mean E-radius, the units of L) NOT 1Rem+100km!
		  tdotp,'elx_bigrf_gsm_mins','elx_pos_gsm_mins',newname='elx_br_tmp'
		  get_data,'elx_br_tmp',data=Br_tmp
		  hemisphere=sign(-Br_tmp.y)
		  r_ift_dip = (1.+100./Rem)
		  calc," 'MLAT' = (180./pival)*arccos(sqrt(Rem*r_ift_dip/'elx_pos_sm_mins_foot_mag')*sin('elx_pos_sm_mins_foot_th'*pival/180.))*hemisphere " ; at footpoint
		  ; interpolate the minute-by-minute data back to the full array
		  get_data,'MLAT',data=MLAT_mins  
		  store_data,'el'+probe+'_MLAT_igrf',data={x: datgsm.x, y: interp(MLAT_mins.y, MLAT_mins.x, datgsm.x)}

		  ;;trace to equator to get L, MLAT, and MLT in IGRF
		  get_data,'elx_pos_gsm_mins_foot',data=elx_pos_eq
		  L1=sqrt(total(elx_pos_eq.y^2.0,2,/nan))/Re
		  store_data,'el'+probe+'_L_igrf',data={x: datgsm.x, y: interp(L1, elx_pos_eq.x, datgsm.x)}
		  
		  elf_mlt_l_lat,'elx_pos_sm_mins_foot',MLT0=MLT0,L0=L0,lat0=lat0
		  get_data, 'elx_pos_sm_mins_foot', data=sm_mins
		  store_data,'el'+probe+'_MLT_igrf',data={x: datgsm.x, y: interp(MLT0, sm_mins.x, datgsm.x)}
		  del_data, '*_mins'
		  ; charsize option was here
		  options,'el'+probe+'_L_igrf',ytitle='L-igrf'
		  ; charsize option was here
		  options,'el'+probe+'_MLAT_igrf',ytitle='MLAT-igrf'
		  ; charsize option was here
		  options,'el'+probe+'_MLT_igrf',ytitle='MLT-igrf'
		  ; charsize option was here

		  ;	EDIT HERE
		  ; EDIT HERE
		  ; adjusting charsize based on 24hr plot or time-specific plot - KC
		  IF(timeduration EQ 86400) THEN BEGIN
			;day plot
			options,'el'+probe+'_MLT_dip',charsize=1.0
			options,'el'+probe+'_L_dip',charsize=1.0
			options,'el'+probe+'_MLAT_dip',charsize=1.0
			options,'el'+probe+'_L_igrf',charsize=1.0
			options,'el'+probe+'_MLAT_igrf',charsize=1.0
			options,'el'+probe+'_MLT_igrf',charsize=1.0
		  ENDIF ELSE BEGIN
			;time-specific plot
			options,'el'+probe+'_MLT_dip',charsize=1.25
			options,'el'+probe+'_L_dip',charsize=1.25
			options,'el'+probe+'_MLAT_dip',charsize=1.25
			options,'el'+probe+'_L_igrf',charsize=1.25
			options,'el'+probe+'_MLAT_igrf',charsize=1.25
			options,'el'+probe+'_MLT_igrf',charsize=1.25
		  ENDELSE
		  ; EDIT END
		  ; EDIT END

		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  ; ... shadow/sunlight bar 0 (shadow) or 1 (sunlight)
		  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		  elf_load_sun_shadow_bar, tplotname='el'+probe+'_pos_gse', no_download=no_download
		  options,'shadow_bar',thick=5.5,xstyle=4,ystyle=4,yrange=[-0.1,0.1],ytitle='',$
			ticklen=0,panel_size=0.1, charsize=2., ztitle=''
		  options,'sun_bar',thick=5.5,xstyle=4,ystyle=4,yrange=[-0.1,0.1],ytitle='',$
			ticklen=0,panel_size=0.1,colors=195, charsize=2., ztitle=''

		  ; create one bar for both sun and shadow
		  store_data, 'sunlight_bar', data=['sun_bar','shadow_bar']
		  options, 'sunlight_bar', panel_size=0.1
		  options, 'sunlight_bar',ticklen=0
		  options, 'sunlight_bar', 'ystyle',4
		  options, 'sunlight_bar', 'xstyle',4
		  options, 'sunlight_bar', yrange=[-0.1,0.1]
		; COPY PASTA END
		; COPY PASTA END


		; CREATING TPLOTS

		; get data from anti, para, and perp tplots
		prefix = 'el' + probe + '_' + mydatatype + '_en_spec2plot_'
		omni_plot = prefix + 'omni'
		anti_plot = prefix + 'anti'
		para_plot = prefix + 'para'
		perp_plot = prefix + 'perp'
		get_data, anti_plot, data = anti_data		; anti_plot is the tplot while anti_data is the data
		get_data, para_plot, data = para_data
		get_data, perp_plot, data = perp_data

		; check that dimensions of all ?_data.y arrays are the same and store them
		anti_size = size(anti_data.y)
		para_size = size(para_data.y)
		perp_size = size(perp_data.y)
		IF((anti_size[4] NE para_size[4]) || (para_size[4] NE perp_size[4])) THEN BEGIN		; ?_size[4] = total num of elements
			stop		; stop for now, maybe interpolate in the future
		ENDIF

		IF(timeduration EQ 86400) THEN BEGIN		; if 24hr plot
			; create anti_and_para tplot
			antiandpara = fltarr(para_size[1], para_size[2])
			antiandpara[*] = !Values.F_NAN
			FOR i=0, para_size[1]-1 DO BEGIN
				FOR j=0, para_size[2]-1 DO BEGIN
					antiandpara[i,j] = anti_data.y[i,j] + para_data.y[i,j]	; add y values
				ENDFOR
			ENDFOR
			
			store_data, prefix + 'anti_and_para', data = {x:para_data.x, y:antiandpara, v:para_data.v}
			anti_and_para_plot = prefix + 'anti_and_para'
			get_data, anti_and_para_plot, data = anti_and_para_data

			options, anti_and_para_plot, 'spec', 1
			ylim, anti_and_para_plot, 50, 5e3, 1
			zlim, anti_and_para_plot, 1e4, 1e9, 1		; 1e4 means 1 x 10^4
			
			; calculate ratio iff anti_and_para and perp are both non-zero
			ratio_data = fltarr(para_size[1], para_size[2])
			ratio_data[*] = !Values.F_NAN					; gives all values in ratio 2D-array a NAN value
			anti_and_para_size = size(antiandpara)
			FOR i=0, anti_and_para_size[1]-1 DO BEGIN
				FOR j=0, anti_and_para_size[2]-1 DO BEGIN
					IF((antiandpara[i,j] NE 0) || (perp_data.y[i,j] NE 0)) THEN BEGIN	; if not the dividing value is not zero, put it into ratio
						ratio_data[i,j] = antiandpara[i,j] / perp_data.y[i,j]
					ENDIF
			  	ENDFOR
			ENDFOR

			; create new ratio tplot
			store_data, prefix + 'anti+para_ovr_perd', data = {x:anti_and_para_data.x, y:ratio_data, v:anti_and_para_data.v}
			ratio_plot = prefix + 'anti+para_ovr_perd'
		ENDIF ELSE BEGIN							; if time-specific plot
			; create new prec tplot
			IF(STRCMP(direction, 'north') EQ 1) THEN BEGIN
				store_data, prefix + 'prec', data = {x:para_data.x, y:para_data.y, v:para_data.v}	; PARA_DATA (if north descending)
			ENDIF ELSE BEGIN
				store_data, prefix + 'prec', data = {x:anti_data.x, y:anti_data.y, v:anti_data.v}	; ANTI_DATA (if south ascending)
			ENDELSE
			prec_plot = prefix + 'prec'			; used for tplot and get_data
			get_data, prec_plot, data = prec_data

			options, prec_plot, 'spec', 1		; 1 means logarithmic scale, 0 means linear
			ylim, prec_plot, 50, 5e3, 1	
			zlim, prec_plot, 1e4, 1e9, 1		; 1e4 means 1 x 10^4

			; calculate ratio iff prec and perp are both non-zero
			ratio_data = fltarr(para_size[1], para_size[2])	; doesn't matter, para_size and anti_size are the same
			ratio_data[*] = !Values.F_NAN					; gives all values in ratio 2D-array a NAN value
			prec_size = size(prec_data.y)
			FOR i=0, prec_size[1]-1 DO BEGIN
				FOR j=0, prec_size[2]-1 DO BEGIN
					IF((prec_data.y[i,j] NE 0) || (perp_data.y[i,j] NE 0)) THEN BEGIN	; if not the dividing value is not zero, put it into ratio
						ratio_data[i,j] = prec_data.y[i,j] / perp_data.y[i,j]
					ENDIF
			  	ENDFOR
			ENDFOR

			; create new ratio tplot
			store_data, prefix + 'prec_ovr_perd', data = {x:prec_data.x, y:ratio_data, v:prec_data.v}
			ratio_plot = prefix + 'prec_ovr_perd'
		ENDELSE

		; scaling z and y axis for ratio_plot
		options, ratio_plot, 'spec', 1		; unrelated but fyi,
		zlim, ratio_plot, 0.02, 2, 1		; EPD-E ranges from 50keV to 4.5MeV
		ylim, ratio_plot, 50, 5e3, 1		; EPD-I ranges from 50keV to 300keV.

		; ztitle options
		options, omni_plot, 'ztitle','nflux'
		options, perp_plot, 'ztitle','nflux'
		IF(timeduration EQ 86400) THEN BEGIN
			options, anti_and_para_plot, 'ztitle','nflux'
		ENDIF ELSE BEGIN
			options, prec_plot, 'ztitle','nflux'
		ENDELSE
		options, ratio_plot, ztitle='ratio'


		; DISPLAYING TPLOTS AND PNGing THEM

		; tplot_options affect global settings
		tplot_options, 'xmargin', [18,11]		; set left/right marjins to 18 and 11 characters, respectively
		tplot_options, 'ymargin', [4, 4]		; set bottom/top margins to 4 and 4 lines, respectively

		date = tstart[k].Substring(0,9)
		IF(timeduration EQ 86400) THEN BEGIN
			;day plot title
			mytitle='PRELIMINARY ELFIN-'+strupcase(probe)+' EPD-'+strupcase(myspecies)+', alt='+strmid(strtrim(alt,1),0,3)+'km, ' + date + ', 24hrs'
		ENDIF ELSE BEGIN
			;time-specific title
			mytitle='PRELIMINARY ELFIN-'+strupcase(probe)+' EPD-'+strupcase(myspecies)+', alt='+strmid(strtrim(alt,1),0,3)+'km, ' + date + ', ' +  tstart[k].Substring(11,15)+' to '+tend[k].Substring(11,15) + ', '
			IF(STRCMP(direction, 'north') EQ 1) THEN BEGIN
				mytitle = mytitle + 'North Descending'
			ENDIF ELSE BEGIN
				mytitle = mytitle + 'South Ascending'
			ENDELSE
		ENDELSE

		if strlowcase(probe) eq 'a' then  $
			varstring=['ela_GLON','ela_MLAT_igrf[ela_MLAT_dip]', 'ela_MLT_igrf[ela_MLT_dip]', 'ela_L_igrf[ela_L_dip]'] else $
			varstring=['elb_GLON','elb_MLAT_igrf[elb_MLAT_dip]', 'elb_MLT_igrf[elb_MLT_dip]', 'elb_L_igrf[elb_L_dip]']

		IF(timeduration EQ 86400) THEN BEGIN
			tplot, [omni_plot, $	
				perp_plot, $
				anti_and_para_plot, $
				ratio_plot, $
				'sunlight_bar', $
				'el'+probe+'_MLAT_igrf'], $
				title = mytitle, $
				var_label = varstring
		ENDIF ELSE BEGIN
			tplot, [omni_plot, $	
				perp_plot, $
				prec_plot, $
				ratio_plot, $
				'sunlight_bar', $
				'el'+probe+'_MLAT_igrf'], $
				title = mytitle, $
				var_label = varstring
		ENDELSE

		; check if directory exists
		dir = yourdir+prefix.Substring(0,2)+'/overplots/'+tstart[k].Substring(0,3)+'/'	; e.g., ...elfin/ela/overplots/2022/
		result = FILE_TEST(dir, /DIRECTORY)	; make year directory if it doesn't exist

		IF(result EQ 0) THEN BEGIN			; 1 means directory exists, 0 means it does not exist
			FILE_MKDIR, dir					; makes directory
		ENDIF
		dir = dir+tstart[k].Substring(5,6)+'/'
		result = FILE_TEST(dir, /DIRECTORY)
		IF(result EQ 0) THEN BEGIN
			FILE_MKDIR, dir					; make month directory if it doesn't exist
		ENDIF
		dir = dir+tstart[k].Substring(8,9)+'/'
		result = FILE_TEST(dir, /DIRECTORY)
		IF(result EQ 0) THEN BEGIN
			FILE_MKDIR, dir					; make day directory if it doesn't exist
		ENDIF

		;day plot or time-specific plot
		IF(timeduration EQ 86400) THEN BEGIN
			;day plot
			makepng, dir + 'el' + probe + '_epd' + myspecies + '_scizone_specplots_' + date
		ENDIF ELSE BEGIN
			;time-specific plot
			time_specific = tstart[k].Substring(11,15)
			makepng, dir + 'el' + probe + '_epd' + myspecies + '_scizone_specplots_' + date + '_' + time_specific
		ENDELSE

	; Writing SM position data onto ASCII file
	get_data, 'el'+probe+'_pos_sm', data=state_pos_sm
	IF(timeduration EQ 86400) THEN BEGIN
		;day data file
		OPENW, 1, FILEPATH('pos_sm_'+date+'.dat', ROOT_DIR=dir)
	ENDIF ELSE BEGIN
		;time-specific data file
		OPENW, 1, FILEPATH('pos_sm_'+date+'_'+tstart[k].Substring(11,15)+'.dat', ROOT_DIR=dir)
	ENDELSE
	
	limit = size(state_pos_sm.x)	; size yields [dimensions,num_of_rows,type_code,num_of_elements]
									; x = time in seconds, y[*,0] = x-coords, y[*,1] = y-coords, y[*,2] = z-coords
	FOR count=0, limit[3]-1 DO BEGIN
		; columns are time as a float, x, y, and z SM coordinates
		state_pos_sm_time = time_string(state_pos_sm.x[count])
		PRINTF, 1, state_pos_sm_time, state_pos_sm.y[count,0], state_pos_sm.y[count,1], state_pos_sm.y[count,2]
	ENDFOR

	CLOSE, 1

	ENDFOR
end
