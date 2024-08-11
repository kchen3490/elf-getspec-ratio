# elf-getspec-ratio routine

### How to use 🤔
This single-use utility program/IDL routine produces energy-time spectrograms, specifically a ratio of precipitating over perpendicular, from the ELFIN mission.

Be sure you have the [SPEDAS framework](https://spedas.org/wiki/index.php?title=Downloads_and_Installation) installed before attempting to run this routine.

Once you have the SPEDAS folder installed, download **[elf_getspec_ratio_kc.pro](elf_getspec_ratio_kc.pro)** and open it with your favorite IDE or text editor.

At the top, you will see some variables to hardcode in. They are as such:
![image](https://github.com/user-attachments/assets/a08853ad-b7fd-48ee-8fdd-49ea6169246a)

To avoid running into errors (because some times will not have the data you wish to observe), download the [data_availability files](https://data.elfin.ucla.edu/ela/data_availability/) in directories '/ela' or '/elb' (depending on the probe) to understand which specific times have data as well as the direction of each specific time.

### Limitations 💀
In addition to the limitations provided in the data_availability files, understand the ELFIN mission's 3U cubesats could only capture so much data, so its primary objective was to investigate electron losses and not as much ion losses during its 4-year orbit. As such, ion data should only be investigated from **2022-06-29 to 2022-09-15 for ELFIN-A and from 2022-07-06 to 2022-09-30 for ELFIN-B** ([source](https://elfin.igpp.ucla.edu/data-notes)).
