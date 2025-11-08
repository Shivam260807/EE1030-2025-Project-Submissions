#include"../c_libs/main.h"
  

int main() {
    char input[100], output[100];
    printf("Enter input File name: \n");
    scanf("%s", input);
    printf("Enter output File name: \n");
    scanf("%s", output);
    int k;
    printf("Enter k\n");
    scanf("%d", &k);

    char pgm_input[128], pgm_output[128];
    strcpy(pgm_input, input);
    strcpy(pgm_output, output);

    if (strstr(input, ".jpg") || strstr(input, ".jpeg")) {
        strcpy(pgm_input, "temp_input.pgm");
        jpg_to_pgm(input, pgm_input);
    }

    int h, w, maxval;
    char fmt[3];
    if (pgm_in(pgm_input, &w, &h, &maxval, fmt)) {
        printf("Error reading input file\n");
        return 1;
    }

    if (k > h || k > w) k = (h < w) ? h : w;

    double *img = malloc(sizeof(double) * h * w);
    double *imgt = malloc(sizeof(double) * w * h);
    double *a = malloc(sizeof(double) * w * w);
    pgm_read(pgm_input, w, h, maxval, img, fmt);

    trans(h, w, img, imgt);
    mul(w, h, w, imgt, img, a);

    double *V = calloc(w * k, sizeof(double));
    double *U = calloc(h * k, sizeof(double));
    double *Sig = calloc(k * k, sizeof(double));

    for (int j = 0; j < k; j++) {
        double *v = malloc(sizeof(double) * w);
        for (int i = 0; i < w; i++) v[i] = 1;
        double err = 1;
        int iter = 0;
        while (fabs(err) > 1e-8 && iter < 100) {
            double *v1 = malloc(sizeof(double) * w);
            mul(w, w, 1, a, v, v1);
            double temp = 0;
            for (int i = 0; i < w; i++) temp += v1[i] * v[i];
            err = temp - 1;
            memcpy(v, v1, sizeof(double) * w);
            normalize(w, v);
            free(v1);
            iter++;
        }

        double *z = malloc(sizeof(double) * w);
        mul(w, w, 1, a, v, z);
        double y = 0, vtv = 0;
        for (int i = 0; i < w; i++) {
            y += v[i] * z[i];
            vtv += v[i] * v[i];
        }
        double sig = (vtv != 0 && y > 0) ? sqrt(y / vtv) : 0;
        Sig[j * k + j] = sig;

        double *u = malloc(sizeof(double) * h);
        mul(h, w, 1, img, v, u);
        if (sig > 1e-12)
            for (int i = 0; i < h; i++) u[i] /= sig;

        for (int i = 0; i < w; i++) V[i * k + j] = v[i];
        for (int i = 0; i < h; i++) U[i * k + j] = u[i];

        double *vvt = malloc(sizeof(double) * w * w);
        for (int i = 0; i < w; i++)
            for (int l = 0; l < w; l++)
                vvt[i * w + l] = v[i] * v[l];

        double factor = (vtv != 0) ? (y / vtv) : 0;
        for (int i = 0; i < w; i++)
            for (int l = 0; l < w; l++)
                a[i * w + l] -= factor * vvt[i * w + l];

        free(v); free(z); free(u); free(vvt);
    }

    double *Vt = malloc(sizeof(double) * k * w);
    trans(w, k, V, Vt);
    double *Usig = malloc(sizeof(double) * h * k);
    double *out = malloc(sizeof(double) * h * w);
    mul(h, k, k, U, Sig, Usig);
    mul(h, k, w, Usig, Vt, out);

    strcpy(pgm_output, output);
    if (strstr(output, ".jpg") || strstr(output, ".jpeg"))
        strcpy(pgm_output, "temp_output.pgm");

    if (pgm_out(pgm_output, w, h, maxval, out))
        printf("Error writing output file\n");

    if (strstr(output, ".jpg") || strstr(output, ".jpeg"))
        pgm_to_jpg(pgm_output, output);
        
    printf("Frobenius Norm ||A - A_k||_F = %lf\n", fnor(img, out, w, h));
    double *temp= malloc(sizeof(double)*w*h);
    for(int i=0; i<w*h; i++){
      temp[i]=0;
    }
    double p=k;
    printf("Relative error = %lf\n", fnor(img, out, w, h)/fnor(img, temp, w, h));
    printf("Compression ratio = %lf\n", (h*w)/(p*(h+w+1)));

    free(img); free(imgt); free(a);
    free(V); free(U); free(Sig);
    free(Vt); free(Usig); free(out);
    return 0;
}

