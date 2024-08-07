# How to use
This IDL routine produces energy-time spectrograms, specifically a ratio of precipitating over perpendicular, from the ELFIN mission.

Be sure you have the [SPEDAS framework installed]([url](http://spedas.org/wiki/index.php?title=Downloads_and_Installation)) before attempting to run this routine.

Once you have the SPEDAS folder installed, download [elf_getspec_ratio_kc.pro]([url](https://github.com/kchen3490/elf-getspec-ratio/blob/main/elf_getspec_ratio_kc.pro)) and open it with your favorite IDE or text editor.

At the top, you will see some variables to hardcode in. They are as such:
  localdir = '/'						    ; directory to store spectrogram plots
	probe = 'a'								    ; 'a' for ELFIN-A or 'b' for ELFIN-B, aka ELFIN-STAR
	myspecies = 'e'							  ; 'e' for electron or 'i' for ion
	datatype = 'pef'						  ; 'pef' for electron or 'pif' for ion (change this also on lines 133-134, 138-139)
	direction = 'south'						; 'north' for north-descending or 'south' for south-ascending
	tstart=['2022-07-18/00:00:00']; time to start investigation in format 'YYYY-MM-DD/hh:mm:ss'
	tend = ['2022-07-19/00:00:00']; time to end investigation in format 'YYYY-MM-DD/hh:mm:ss'

To avoid running into errors (because some times will not have the data you wish to observe), download the [data_availability files]([url](https://data.elfin.ucla.edu/ela/data_availability/)) in directories '/ela' or '/elb' (depending on the probe) to understand which specific times have data as well as the direction of each specific time.

# Limitations
In addition to the limitations provided in the data_availability files, understand the ELFIN mission's 3U cubesats could only capture so much data, so its primary objective was to investigate electron losses and not as much ion losses during its 4-year orbit. As such, ion data should only be investigated from 2022-06-29 to 2022-09-15 for ELFIN-A and from 2022-07-06 to 2022-09-30 for ELFIN-B ([source]([url](https://elfin.igpp.ucla.edu/data-notes))).
