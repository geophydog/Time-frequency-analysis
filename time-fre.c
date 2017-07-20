#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "sacio.h"

int near_pow2(int n) {
    int m;
    float f;
    f = log((float)n) / log(2.) + 1.;
    m = (int)pow(2.,(int)f);
    return m;
}

int main(int argc, char *argv[]) {
    int i, k, npts, seg_npts, seg_num, begin_index;
    float seg_len, *data, t1, t2, time = 0., fre, amp, high_f, peak = 0.;
    SACHEAD hd;
    fftw_complex *in, *out;
    fftw_plan p;
    FILE *ft, *fp;

    if( argc != 7 ) {
        fprintf(stderr,"Usage: time-fre <sacfile> <t1> <t2> <segmentation_length> <high_frequency_limit> <output_result_file_name>\n");
        fprintf(stderr,"       return power-spectral in time-frequency domains,which will be saved output_result_file_name\n");
        fprintf(stderr,"       <sacfile>                 The name of inputing SAC format file;\n");
        fprintf(stderr,"       <t1>                      The beginning time of time-fre analysis;\n");
        fprintf(stderr,"       <t2>                      The ending time of time-fre analysis;\n");
        fprintf(stderr,"       <segmentation_length>     The short FFT time length;\n");
        fprintf(stderr,"       <high_frerquency_limit>   The maximum of frequency of outputing power spectral;\n");
        fprintf(stderr,"       <output_result_file_name> The name of outputing results;\n");
        fprintf(stderr,"         Power spectral plot by will be saved in file \"plot.sh\",just run with command \"sh plot.sh\".\n");
        fprintf(stderr,"         Attention!!! Executing \"sh plot.sh\" requires GMT(the Generic Mapping Tools).\n");
        exit(1);
    }

    t1 = atof(argv[2]);
    t2 = atof(argv[3]);
    seg_len= atof(argv[4]);
    high_f = atof(argv[5]);
    ft = fopen(argv[6],"w");
    fp = fopen("plot.sh","w");

    data = read_sac(argv[1],&hd);
    seg_num = (int)((t2-t1)/seg_len);
    seg_npts = (int)(seg_len/hd.delta);
    begin_index = (int)(t1/hd.delta);
    npts = near_pow2(seg_npts);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);

    while ( k < seg_num ) {
        time = t1 + k*seg_len+seg_len/2.;

        for ( i = 0; i < npts; i ++ ) {
            if( i < seg_npts ) in[i][0] = data[begin_index+i+k*seg_npts];
            else in[i][0] = 0.;
            in[i][1] = 0.;
        }

        p = fftw_plan_dft_1d(npts, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);

        for(i = 0; i < seg_npts; i ++) {
            amp = sqrt(pow(out[i][0],2.)+pow(out[i][1],2.));
            fre = i/hd.delta/npts;
            if ( peak < amp ) peak = amp;
            fprintf(ft,"%6.4f %6.4f %6.4f\n",time, fre, amp);
        }

        k += 1;
        fftw_destroy_plan(p);
    }
    fftw_free(in); fftw_free(out);

    fprintf(fp,"gmtset LABEL_FONT Times-Bold\n");
    fprintf(fp,"gmtset LABEL_FONT_SIZE 15p\n");
    fprintf(fp,"awk '{print $1,$2,$3/%f}' %s >tmp.file\n", peak, argv[6]);
    fprintf(fp,"rm %s\n", argv[6]);
    fprintf(fp,"mv tmp.file %s\n", argv[6]);
    fprintf(fp,"makecpt -Cjet.cpt -T0/1.0/0.1 -Z >tmp1.cpt\n");
    fprintf(fp,"psxy -R%f/%.1f/0/%f -JX8i/5i -K -T >plot.ps\n", t1, t2, high_f);
    fprintf(fp,"surface %s -R -I%d/%f -G%s.grd\n" ,argv[6], (int)(seg_len/2), high_f/100., argv[6]);
    fprintf(fp,"grd2cpt %s.grd -Cjet>tmp.cpt\n",argv[6]);
    fprintf(fp,"grdimage %s.grd -R -J -K -O -Ctmp.cpt -B%d:\"Time(sec)\":/%f:\"Frequency(Hz)\":WSen>>plot.ps\n", argv[6], (int)((t2-t1)/10.), high_f/10.);
    fprintf(fp,"echo %s 0 1 | pssac2 -R%f/%f/0/2 -JX8i/1i -K -O -B%d/0:\"Amplitude\":WSen -C%f/%f -M0.8 -W0.8p -Xa0i -Ya5.4i>>plot.ps\n", argv[1], t1, t2,(int)((t2-t1)/10), t1, t2);
    fprintf(fp,"psscale -Ctmp1.cpt -D8.5i/2.5i/12.55/0.8 -Ba0.1g0:\"Normalized amplitude power\": -K -O >>plot.ps\n");
    fprintf(fp,"echo %.1f 1 20 0 Times-Bold 0 %s | pstext -R -J -K -O -Xa-1i -Ya6.2i>>plot.ps\n", (t1+t2)/2., argv[1]);
    fprintf(fp,"psxy -R -J -O -T >>plot.ps\n");
    fprintf(fp,"ps2pdf plot.ps plot.pdf\n");
    fprintf(fp,"rm tmp*.cpt\n");
    fprintf(fp,"evince plot.pdf\n");
    fclose(fp);

    free(data);
    fclose(ft);
    return 0;
}

