R=-2500/2500/0.05/0.6
PS=COR_CO.CSB.00.LHZ_CO.CSB.01.LHZ.SAC.ps
PDF=COR_CO.CSB.00.LHZ_CO.CSB.01.LHZ.SAC.pdf

gmt gmtset MAP_TICK_LENGTH_PRIMARY 0.1c
gmt gmtset MAP_TICK_LENGTH_SECONDARY 0.05c
gmt gmtset MAP_LABEL_OFFSET 12P
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 9p
gmt gmtset FONT_LABEL 18p,Times-Bold
gmt gmtset FONT_ANNOT_PRIMARY 15p,5

awk '{print $1,$2,$3/6.2439}' out.txt >tmp.file
gmt makecpt -Cjet.cpt -T0/1.0/0.1 -Z >tmp1.cpt
gmt surface tmp.file -R$R -I1.66667/0.00275 -Gtmp.grd
gmt grd2cpt tmp.grd -Cjet>tmp.cpt

gmt psxy -R$R -JX8i/4i -K -T > $PS
gmt grdimage tmp.grd -R -J -K -O -Ctmp.cpt -B500:"Time(sec)":/0.055:"Frequency(Hz)":WSen -Xa0.5i -Ya1i >> $PS
echo COR_CO.CSB.00.LHZ_CO.CSB.01.LHZ.SAC -2500 1 | gmt pssac -R-2500/2500/0/2 -JX8i/1i -K -O -B500/0:"Amplitude":WSen -C-2500/2500 -M0.8 -W0.8p -Xa0.5i -Ya5.5i >> $PS
gmt psscale -Ctmp1.cpt -D8.5i/3i/10.15/0.5 -Ba0.1g0:"Normalized amplitude power": -K -O >> $PS
echo 0 1 COR_CO.CSB.00.LHZ_CO.CSB.01.LHZ.SAC | gmt pstext -R -J -K -F+f15,5 -O -Xa0.5i -Ya6.2i>> $PS
gmt psxy -R -J -O -T >> $PS

gmt psconvert -Tt -A -P -E300 $PS
ps2pdf $PS $PDF
rm tmp* out.txt.* gmt.* $PS
evince $PDF
