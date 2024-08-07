pro elf_getspec_example_cw
; A quick example to load ELFIN EPD-E spectra
; Colin Wilkins (colinwilkins@ucla.edu) - 2024/6/12

elf_init
;!elf.local_data_dir='D:\docs\data\elfin_data\'
sclet='a'
;tstart='2020-08-04/12:40:00'
;tend='2020-08-04/12:47:00'
tstart='2021-11-04/00:00:00'
tend='2021-11-05/00:00:00'
time2plot=[tstart,tend]
timeduration=time_double(tend)-time_double(tstart)
timespan,tstart,timeduration,/seconds

mytype='eflux'
species = 'pef'

elf_load_state,probe=sclet,trange=time2plot
;
elf_load_epd,probe=sclet,datatype=species,type = 'raw',trange=time2plot
elf_getspec,probe=sclet,datatype=species,type = 'raw',fullspin=fullspin,/get3Dspec ; this simply acts on the data loaded previously

; reload data in eflux units now (two calls: first reload data, then recompute spectra)
elf_load_epd,probe=sclet,datatype=species,type = mytype,trange=time2plot
elf_getspec,probe=sclet,datatype=species,type = mytype,fullspin=fullspin,/get3Dspec ; this simply acts on the data loaded previously

tplot, ['ela_pef_en_spec2plot_omni','ela_pef_en_spec2plot_perp','ela_pef_en_spec2plot_para','ela_pef_en_spec2plot_anti']
makepng,'elfin_epde_example_scizone_specplots'
;makepng,'/data/share/elfin/plots/elfin_epde_example_scizone_specplots'
stop
end
