## 4. Delomas et al. method
Get mtype2 installed: follow instructions on Delomas github: for university servers, you will likely have to download it you your files and call it from there.
download microhapWrap.py from Delomas github

If using full set of loci (341): wolfAmpRef_5.fa, with associated files: wolfAmpRef_5.1.bt2, wolfAmpRef_5.2.bt2, wolfAmpRef_5.3.bt2, wolfAmpRef_5.4.bt2, wolfAmpRef_5.rev.1.bt2, wolfAmpRef_5.rev.2.bt2
If using fecal set of loci (200): wolfAmpRef_200.fa, with associated files: wolfAmpRef_200.1.bt2, wolfAmpRef_200.2.bt2, wolfAmpRef_200.3.bt2, wolfAmpRef_200.4.bt2, wolfAmpRef_200.rev.1.bt2, wolfAmpRef_200.rev.2.bt2

Create homebrew genome: you can generate in a text editer: two lines per locus. first line has name of locus in format: >Clu_CM00000009_###; second line is the locus sequence: for this panel it if 79 basepairs
Save it as .fa

Generate associated files:
>module load bowtie2
>bowtie2-build wolfAmpRef_5.fa wolfAmpRef_5

Second argument is the prefix that will be used to make the .bt2 files

Other files needed: Clu200_pos_4.txt, primers_with_adapters.txt, Clu341_pres_abs.txt

For fecal: you should change
