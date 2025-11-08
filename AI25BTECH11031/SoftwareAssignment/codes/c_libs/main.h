#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <jpeglib.h>

void jpg_to_pgm(const char *jpg_file, const char *pgm_file) {
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *infile, *outfile;
    JSAMPARRAY buffer;
    int row_stride;

    if ((infile = fopen(jpg_file, "rb")) == NULL) {
        fprintf(stderr, "Cannot open %s\n", jpg_file);
        return;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);

    row_stride = cinfo.output_width * cinfo.output_components;
    buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);

    outfile = fopen(pgm_file, "wb");
    if (!outfile) {
        fprintf(stderr, "Cannot create %s\n", pgm_file);
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return;
    }

    fprintf(outfile, "P5\n%d %d\n255\n", cinfo.output_width, cinfo.output_height);
    unsigned char *row = malloc(cinfo.output_width);
    if (!row) {
        fprintf(stderr, "Memory error\n");
        fclose(outfile);
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return;
    }

    while (cinfo.output_scanline < cinfo.output_height) {
        jpeg_read_scanlines(&cinfo, buffer, 1);
        if (cinfo.output_components == 3) {
            for (int x = 0; x < cinfo.output_width; x++) {
                int r = buffer[0][x * 3];
                int g = buffer[0][x * 3 + 1];
                int b = buffer[0][x * 3 + 2];
                row[x] = (unsigned char)(0.3 * r + 0.59 * g + 0.11 * b);
            }
        } else if (cinfo.output_components == 1) {
            for (int x = 0; x < cinfo.output_width; x++)
                row[x] = buffer[0][x];
        }
        fwrite(row, 1, cinfo.output_width, outfile);
    }

    free(row);
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    fclose(outfile);
}

void pgm_to_jpg(const char *pgm_file, const char *jpg_file) {
    FILE *infile = fopen(pgm_file, "r");
    if (!infile) {
        fprintf(stderr, "Cannot open %s\n", pgm_file);
        return;
    }

    char magic[3];
    int width, height, maxval;
    fscanf(infile, "%s", magic);
    fscanf(infile, "%d %d %d", &width, &height, &maxval);
    fgetc(infile);

    unsigned char *gray = malloc(width * height);
    fread(gray, 1, width * height, infile);
    fclose(infile);

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile = fopen(jpg_file, "wb");
    if (!outfile) {
        fprintf(stderr, "Cannot create %s\n", jpg_file);
        free(gray);
        return;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 90, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    JSAMPROW row_pointer;
    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer = (JSAMPROW)&gray[cinfo.next_scanline * width];
        jpeg_write_scanlines(&cinfo, &row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    fclose(outfile);
    free(gray);
}

int pgm_in(char *filename, int *width, int *height, int *maxval, char fmt[3]) {
    FILE *f = fopen(filename, "rb");
    if (!f) return 1;
    if (fscanf(f, "%2s", fmt) != 1) { fclose(f); return 1; }
    if (!(strcmp(fmt, "P2") == 0 || strcmp(fmt, "P5") == 0)) { fclose(f); return 1; }

    fscanf(f, "%d", width);
    fscanf(f, "%d", height);
    fscanf(f, "%d", maxval);
    fgetc(f);
    fclose(f);
    return 0;
}

int pgm_read(char *filename, int width, int height, int maxval, double *image, char fmt[3]) {
    FILE *f = fopen(filename, "rb");
    if (!f) return 1;
    if (fscanf(f, "%2s", fmt) != 1) { fclose(f); return 1; }
    int wtmp, htmp, mtmp;
    fscanf(f, "%d", &wtmp);
    fscanf(f, "%d", &htmp);
    fscanf(f, "%d", &mtmp);
    fgetc(f);

    if (strcmp(fmt, "P5") == 0) {
        unsigned char *temp = malloc(width * height);
        fread(temp, 1, width * height, f);
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                image[i * width + j] = (double)temp[i * width + j];
        free(temp);
    } else {
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                fscanf(f, "%lf", &image[i * width + j]);
    }

    fclose(f);
    return 0;
}

int pgm_out(char *filename, int width, int height, int maxval, double *image) {
    FILE *f = fopen(filename, "wb");
    if (!f) return 1;
    fprintf(f, "P5\n%d %d\n%d\n", width, height, maxval);
    unsigned char *temp = malloc(width * height);
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            double val = image[i * width + j];
            if (val < 0) val = 0;
            if (val > maxval) val = maxval;
            temp[i * width + j] = (unsigned char)(val + 0.5);
        }
    fwrite(temp, 1, width * height, f);
    free(temp);
    fclose(f);
    return 0;
}

void mul(int m, int n, int p, double *a, double *b, double *result) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            double sum = 0;
            for (int l = 0; l < n; l++)
                sum += a[i * n + l] * b[l * p + j];
            result[i * p + j] = sum;
        }
    }
}

void trans(int n, int m, double *a, double *b) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            b[j * n + i] = a[i * m + j];
}

double norm(int n, double *v) {
    double s = 0;
    for (int i = 0; i < n; i++)
        s += v[i] * v[i];
    return sqrt(s);
}

void normalize(int n, double *v) {
    double nrm = norm(n, v);
    if (nrm == 0) return;
    for (int i = 0; i < n; i++)
        v[i] /= nrm;
}

double fnor(double *a, double *b, int n, int m) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double d = a[i * m + j] - b[i * m + j];
            sum += d * d;
        }
    }
    return sqrt(sum);
}

