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
    int i, k, npts, seg_npts, seg_num;
    float seg_len, *data, time = 0., fre, amp, high_f, peak = 0.;
    SACHEAD hd;
    fftw_complex *in, *out;
    fftw_plan p;
    FILE *ft, *fp;

    if( argc != 5 ) {
        fprintf(stderr,"Usage: time-fre <sacfile> <segmentation_length> <high_frequency_limit> <output_result_file_name>\n");
        fprintf(stderr,"       return power-spectral in time-frequency domains,which will be saved output_result_file_name\n");
        fprintf(stderr,"       [sacfile] the name of inputing SAC format file;\n");
        fprintf(stderr,"       [segmentation_length] the short FFT time length;\n");
        fprintf(stderr,"       [high_frerquency_limit] the maximum of frequency of outputing power spectral;\n");
        fprintf(stderr,"       [output_result_file_name] the name of outputing results;\n");
        fprintf(stderr,"         Power spectral plot by will be saved in file \"plot.sh\",just run with command \"sh plot.sh\".\n");
        fprintf(stderr,"         Attention!!! Executing \"sh plot.sh\" requires GMT(the Generic Mapping Tools).\n");
        exit(1);
    }

    seg_len= atof(argv[2]);
    high_f = atof(argv[3]);
    ft = fopen(argv[4],"w");
    fp = fopen("plot.sh","w");

    data = read_sac(argv[1],&hd);
    seg_num = (int)(hd.e/seg_len);
    seg_npts = (int)(seg_len/hd.delta);
    npts = near_pow2(seg_npts);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);

    while ( k < seg_num ) {
        time = k*seg_len+seg_len/2.;

        for ( i = 0; i < npts; i ++ ) {
            if( i < seg_npts ) in[i][0] = data[i+k*seg_npts];
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
    fprintf(fp,"awk '{print $1,$2,$3/%f}' %s >tmp.file\n", peak, argv[4]);
    fprintf(fp,"rm %s\n", argv[4]);
    fprintf(fp,"mv tmp.file %s\n", argv[4]);
    fprintf(fp,"makecpt -Cjet.cpt -T0/1.0/0.1 -Z >tmp1.cpt\n");
    fprintf(fp,"psxy -R0/%.1f/0/%f -JX8i/5i -K -T >plot.ps\n", hd.e, high_f);
    fprintf(fp,"surface %s -R -I%d/%f -G%s.grd\n" ,argv[4], (int)(seg_len/2), high_f/100., argv[4]);
    fprintf(fp,"grd2cpt %s.grd -Cjet>tmp.cpt\n",argv[4]);
    fprintf(fp,"grdimage %s.grd -R -J -K -O -Ctmp.cpt -B%d:\"Time(sec)\":/%f:\"Frequency(Hz)\":WSen>>plot.ps\n", argv[4], (int)(hd.e/10.), high_f/10.);
    fprintf(fp,"echo %s 0 1 | pssac2 -R0/%d/0/2 -JX8i/1i -K -O -B%d/1:\"amp\":WSen -C0/%d -M0.8 -W0.8p -Xa0i -Ya5.4i>>plot.ps\n", argv[1], (int)hd.e, (int)(hd.e/10.), (int)hd.e);
    fprintf(fp,"psscale -Ctmp1.cpt -D8.5i/2.5i/12.55/0.8 -Ba0.1g0:\"Amplitude power\": -K -O >>plot.ps\n");
    fprintf(fp,"echo %.1f 1 20 0 Times-Bold 0 %s | pstext -R -J -K -O -Xa-1i -Ya6.2i>>plot.ps\n", hd.e/2, argv[1]);
    fprintf(fp,"psxy -R -J -O -T >>plot.ps\n");
    fprintf(fp,"ps2pdf plot.ps plot.pdf\n");
    fprintf(fp,"rm tmp*.cpt\n");
    fclose(fp);

    free(data);
    fclose(ft);
    return 0;
}

