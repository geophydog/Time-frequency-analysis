#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include "sacio.h"

int near_pow2(int n) {
    int m, i;
    i = (int)(log((float)n) / log(2.) + 1.);
    m = (int)pow(2., i);
    return m;
}

int main(int argc, char *argv[]) {
    int i, k, npts, seg_npts, seg_num, begin_index, f1_i, f2_i;
    float seg_len, *data, t1, t2, time = 0., fre, amp, f1, f2, sp, peak = 0.;
    SACHEAD hd;
    fftw_complex *in, *out;
    fftw_plan p;
    FILE *ft, *fp;
    char str[256] = {""};

    if( argc != 8 ) {
        fprintf(stderr,"Usage: time-fre <sacfile> <t1> <t2> <segmentation_length> <low_frequency> <high_frequency> <output_result_file_name>\n");
        fprintf(stderr,"       return power-spectral in time-frequency domains,which will be saved output_result_file_name\n");
        fprintf(stderr,"       <sacfile>                 The name of inputing SAC format file;\n");
        fprintf(stderr,"       <t1>                      The beginning time of time-fre analysis;\n");
        fprintf(stderr,"       <t2>                      The ending time of time-fre analysis;\n");
        fprintf(stderr,"       <segmentation_length>     The short FFT time length;\n");
        fprintf(stderr,"       <low_frequency>           The low corner frequency;\n");
        fprintf(stderr,"       <high_frerquency>         The high corner frequency;\n");
        fprintf(stderr,"       <output_result_file_name> The name of outputing results;\n");
        fprintf(stderr,"       <<<Power spectral plot by will be saved in file \"plot.sh\",just run with command \"sh plot.sh\".\n");
        fprintf(stderr,"       <<<Attention!!! Executing \"sh plot.sh\" requires GMT(the Generic Mapping Tools).\n");
        exit(1);
    }

    t1 = atof(argv[2]);
    t2 = atof(argv[3]);
    seg_len= atof(argv[4]);
    f1 = atof(argv[5]);
    f2 = atof(argv[6]);
    ft = fopen(argv[7],"w");
    fp = fopen("plot.sh","w");

    data = read_sac(argv[1],&hd);
    sp = 1/hd.delta;
    if ( f2 >= sp/2 ) {
        f2 = sp/2*0.99;
        fprintf(stderr, "ATTENTION!!! %g(your high_fre) >= %g(sp)\n", f2, sp);
    }

    seg_num = (int)((t2-t1)/seg_len);
    seg_npts = (int)(seg_len/hd.delta);
    begin_index = (int)((t1-hd.b)/hd.delta);
    npts = near_pow2(seg_npts);
    f1_i = (int) (f1*npts/sp);
    f2_i = (int) (f2*npts/sp);
    printf("fft npts: %d low_f: %g  high_f: %g  high_f_i: %d   low_f_i: %d  t1: %g  t2: %g seg_len: %g\n", \
            npts, f1, f2, f1_i, f2_i, t1, t2, seg_len);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);

    while ( k <= seg_num ) {
        time = t1 + k*seg_len + seg_len/2.;

        for ( i = 0; i < npts; i ++ ) {
            if( i < seg_npts ) in[i][0] = data[begin_index+i+k*seg_npts];
            else in[i][0] = 0.;
            in[i][1] = 0.;
        }

        p = fftw_plan_dft_1d(npts, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);

        for(i = f1_i; i <= f2_i; i ++) {
            fre = i*sp/npts;
            //printf("%d  %g\n", i, fre);
            amp = sqrt(pow(out[i][0],2.)+pow(out[i][1],2.));
            if ( peak < amp ) peak = amp;
            fprintf(ft,"%g %g %g\n",time, fre, amp);
        }

        k += 1;
        fftw_destroy_plan(p);
    }
    fftw_free(in); fftw_free(out);

//----------------------------------------------------shell plot script-------------------------------------------------------------
    fprintf(fp,"gmtset LABEL_FONT Times-Bold\n");
    fprintf(fp,"gmtset LABEL_FONT_SIZE 20\n");
    fprintf(fp,"gmtset ANNOT_FONT_PRIMARY 5\ngmtset ANNOT_FONT_SIZE_PRIMARY 15\n");
    fprintf(fp,"awk '{print $1,$2,$3/%g}' %s >tmp.file\n", peak, argv[7]);
    fprintf(fp,"makecpt -Cjet.cpt -T0/1.0/0.1 -Z >tmp1.cpt\n");
    fprintf(fp,"psxy -R%g/%g/%g/%g -JX8i/5i -K -T >%s.ps\n", t1, t2, f1, f2, argv[1]);
    fprintf(fp,"surface tmp.file -R -I%g/%g -Gtmp.grd\n", (t2-t1)/3333, (f2-f1)/200);
    fprintf(fp,"grd2cpt tmp.grd -Cjet>tmp.cpt\n");
    fprintf(fp,"grdimage tmp.grd -R -J -K -O -Ctmp.cpt -B%d:\"Time(sec)\":/%g:\"Frequency(Hz)\":WSen>>%s.ps\n", \
            (int)((t2-t1)/10.), (f2-f1)/10., argv[1]);
    fprintf(fp,"echo %s 0 1 | pssac2 -R%g/%g/0/2 -JX8i/1i -K -O -B%d/0:\"Amplitude\":WSen -C%g/%g -M0.8 -W0.8p -Xa0i -Ya5.4i>>%s.ps\n", \
            argv[1], t1, t2,(int)((t2-t1)/10), t1, t2, argv[1]);
    fprintf(fp,"psscale -Ctmp1.cpt -D8.5i/2.5i/12.55/0.6 -Ba0.1g0:\"Normalized amplitude power\": -K -O >>%s.ps\n", argv[1]);
    fprintf(fp,"echo %g 1 20 0 Times-Bold 0 %s | pstext -R -J -K -O -Xa-1i -Ya6.2i>>%s.ps\n", (t1+t2)/2., argv[1], argv[1]);
    fprintf(fp,"psxy -R -J -O -T >>%s.ps\n", argv[1]);
    fprintf(fp,"ps2raster -Tt -A -P -E300 %s.ps\n", argv[1]);
    fprintf(fp,"ps2pdf %s.ps %s.pdf\n", argv[1], argv[1]);
    fprintf(fp,"rm tmp*\n");
    fprintf(fp,"evince %s.pdf\n", argv[1]);
    fclose(fp);

    free(data);
    fclose(ft);
    return 0;
}
