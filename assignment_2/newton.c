#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define lower_bound 0.000001
#define upper_bound 10000000000

// arguments
int threads = -1, lines = -1, degree;

// mutices
pthread_mutex_t result_mutex, done_mutex;

// iteration functions for each degree
typedef void (*complex_function)(double complex *);
complex_function iterate;

void iterate_1(complex *z)
{
    *z = 1.0;
}

void iterate_2(complex *z)
{
    *z = *z * 0.5 + 1.0 / (2.0 * *z);
}

void iterate_3(complex *z)
{
    *z = *z * (2.0 / 3) + 1.0 / (3 * *z * *z);
}

void iterate_4(complex *z)
{
    *z = *z * (3.0 / 4) + 1.0 / (4 * *z * *z * *z);
}

void iterate_5(complex *z)
{
    *z = *z * (4.0 / 5) + 1.0 / (5 * *z * *z * *z * *z);
}

void iterate_6(complex *z)
{
    *z = *z * (5.0 / 6) + 1.0 / (6 * *z * *z * *z * *z * *z);
}

void iterate_7(complex *z)
{
    *z = *z * (6.0 / 7) + 1.0 / (7 * *z * *z * *z * *z * *z * *z);
}

void iterate_8(complex *z)
{
    *z = *z * (7.0 / 8) + 1.0 / (8 * *z * *z * *z * *z * *z * *z * *z);
}

void iterate_9(complex *z)
{
    *z = *z * (8.0 / 9) + 1.0 / (9 * *z * *z * *z * *z * *z * *z * *z * *z);
}

void iterate_10(complex *z)
{
    *z = *z * (9.0 / 10) + 1.0 / (10 * *z * *z * *z * *z * *z * *z * *z * *z * *z);
}

// known roots for each degree
double complex *known_roots;

// result matrices
int **roots;
int **iters;
char *done;

void compute_line(int line)
{
    int *result_roots = malloc(lines * sizeof(int));
    int *result_iters = malloc(lines * sizeof(int));

    double complex z;
    double old_im_z = 2.0 - (4.0 * line / lines);
    double sq_abs_z;

    double distance;
    short int conv;
    short int attr;

    int cx = 0;

    for (double re_z = -2.0; cx < lines; ++cx, re_z += (4.0 / lines))
    {
        z = re_z + old_im_z * I;
        sq_abs_z = creal(z) * creal(z) + cimag(z) * cimag(z);
        for (conv = 0, attr = -1;; ++conv)
        {
            // if real or imag part of z is bigger than upper_bound, assume it will diverge
            if (creal(z) > upper_bound || cimag(z) > upper_bound)
            {
                attr = 0;
                break;
            }
            // if z squared is close to origin, assume it will diverge
            if (sq_abs_z < lower_bound)
            {
                attr = 0;
                break;
            }
            for (size_t ix = 0; ix < degree; ++ix)
            {
                distance = (creal(known_roots[ix]) - creal(z)) *
                               (creal(known_roots[ix]) - creal(z)) +
                           (cimag(known_roots[ix]) - cimag(z)) *
                               (cimag(known_roots[ix]) - cimag(z));
                if (distance < lower_bound)
                {
                    attr = ix + 1;
                    break;
                }
            }
            iterate(&z);
            if (attr != -1)
            {
                break;
            }
        }
        result_roots[cx] = attr;
        result_iters[cx] = fmin(conv, 55);
    }

    pthread_mutex_lock(&result_mutex);
    roots[line] = result_roots;
    iters[line] = result_iters;
    pthread_mutex_unlock(&result_mutex);

    pthread_mutex_lock(&done_mutex);
    done[line] = 1;
    pthread_mutex_unlock(&done_mutex);
}

void *compute_lines(void *restrict arg)
{
    int start_line = ((int *)arg)[0];
    int offset = ((int *)arg)[1];
    free(arg);

    for (int line = start_line; line < lines; line += offset)
    {
        compute_line(line);
    }

    return NULL;
}

void *writer_function()
{
    char root_colors[] = "0 0 8 0 8 0 8 0 0 1 1 5 1 5 1 5 1 1 2 2 2 3 2 2 2 3 2 3 2 2";
    char iter_colors[] = "00 00 00 01 01 01 02 02 02 03 03 03 04 04 04 05 05 05 06 06 06 07 07 07 08 08 08 09 09 09 \
10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 17 17 17 18 18 18 19 19 19 \
20 20 20 21 21 21 22 22 22 23 23 23 24 24 24 25 25 25 26 26 26 27 27 27 28 28 28 29 29 29 \
30 30 30 31 31 31 32 32 32 33 33 33 34 34 34 35 35 35 36 36 36 37 37 37 38 38 38 39 39 39 \
40 40 40 41 41 41 42 42 42 43 43 43 44 44 44 45 45 45 46 46 46 47 47 47 48 48 48 49 49 49 \
50 50 50 51 51 51 52 52 52 53 53 53 54 54 54 55 55 55 ";

    char **root_pixels = malloc(lines * sizeof root_pixels);
    for (int i = 0; i < lines; i++)
        root_pixels[i] = malloc(6 * lines * sizeof(char));

    char **iter_pixels = malloc(lines * sizeof iter_pixels);
    for (int i = 0; i < lines; i++)
        iter_pixels[i] = malloc(9 * lines * sizeof(char));

    FILE *fptr_roots;
    FILE *fptr_iters;

    // initialize ppm by printing headers
    char filenamn_roots[26];
    char filenamn_iters[27];
    snprintf(filenamn_roots, sizeof(char) * 26, "newton_attractors_x%d.ppm", degree);
    snprintf(filenamn_iters, sizeof(char) * 27, "newton_convergence_x%d.ppm", degree);
    fptr_roots = fopen(filenamn_roots, "w");
    fptr_iters = fopen(filenamn_iters, "w");
    fprintf(fptr_roots, "P3\n%d %d\n8\n", lines, lines);
    fprintf(fptr_iters, "P3\n%d %d\n55\n", lines, lines);

    struct timespec sleep_timespec = {0};

    // sleep time one microsecond
    sleep_timespec.tv_nsec = 1 * 100000;

    for (int current_line = 0; current_line < lines;)
    {
        pthread_mutex_lock(&done_mutex);
        if (done[current_line] == 1)
        {
            pthread_mutex_unlock(&done_mutex);

            pthread_mutex_lock(&result_mutex);
            int *local_roots = roots[current_line];
            int *local_iter = iters[current_line];
            pthread_mutex_unlock(&result_mutex);

            for (int column = 0, root_index = 0, iter_index = 0; column < lines; column++, root_index += 6, iter_index += 9)
            {
                memcpy(root_pixels[current_line] + root_index, root_colors + local_roots[column] * 6, 6);
                memcpy(iter_pixels[current_line] + iter_index, iter_colors + local_iter[column] * 9, 9);
            }

            fwrite(root_pixels[current_line], sizeof(char), 6 * lines, fptr_roots);
            fwrite(iter_pixels[current_line], sizeof(char), 9 * lines, fptr_iters);

            fwrite("\n", 1, 1, fptr_roots);
            fwrite("\n", 1, 1, fptr_iters);
            current_line++;
        }
        else
        {
            pthread_mutex_unlock(&done_mutex);
            nanosleep(&sleep_timespec, NULL);
        }
    }

    for (int line = 0; line < lines; line++)
    {
        free(root_pixels[line]);
        free(iter_pixels[line]);
    }
    free(root_pixels);
    free(iter_pixels);

    fclose(fptr_roots);
    fclose(fptr_iters);

    pthread_exit(NULL);
}

int main(int argc, char **argv)
{
    // ########################################################################
    // parse arguments
    // ########################################################################
    extern char *optarg;
    extern int optind;
    int c, err = 0;

    static char usage[] = "\nusage: %s -t threads -l lines d degree\n";

    while ((c = getopt(argc, argv, "t:l:")) != -1)
    {
        switch (c)
        {
        case 't':
            threads = atoi(optarg);
            break;
        case 'l':
            lines = atoi(optarg);
            break;
        case '?':
            err = 1;
            break;
        }
    }

    if (threads == -1)
    {
        printf("error: -t is required!\n");
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }
    if (lines == -1)
    {
        printf("error: -l is required!\n");
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }
    if (err)
    {
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }

    if (optind < argc)
    {
        degree = atoi(argv[optind]);
        // printf("threads: %d\nlines:   %d\ndegree:  %d\n", threads, lines, degree);

        if (threads < 1 || threads > lines)
        {
            printf("error: t should be an integer between 1 to lines, or in other words 1 to %d!\n", lines);
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
        if (lines < 1 || lines > 100000)
        {
            printf("error: l should be an integer between 1 to 100000!\n");
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
        if (degree < 1 || degree > 9)
        {
            printf("error: d should be an integer between 1 to 9!\n");
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
    }
    else
    {
        printf("error: d is required!\n");
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }
    // argument parsing done!

    // ########################################################################
    // computation part below
    // ########################################################################
    known_roots = malloc(degree * sizeof(double complex));

    switch (degree)
    {
    case 1:
        iterate = iterate_1;
        known_roots[0] = 1.0 + 0.0 * I;
        break;
    case 2:
        iterate = iterate_2;
        known_roots[0] = 1.0 + 0.0 * I;
        known_roots[1] = -1.0 + 0.0 * I;
        break;
    case 3:
        iterate = iterate_3;
        known_roots[0] = 1.0 + 0.0 * I;
        known_roots[1] = -0.5 - 0.86603 * I;
        known_roots[2] = -0.5 + 0.86603 * I;
        break;
    case 4:
        iterate = iterate_4;
        known_roots[0] = -1.0 + 0.0 * I;
        known_roots[1] = 1.0 + 0.0 * I;
        known_roots[2] = 0.0 - 1.0 * I;
        known_roots[3] = 0.0 + 1.0 * I;
        break;
    case 5:
        iterate = iterate_5;
        known_roots[0] = 1.0 + 0.0 * I;
        known_roots[1] = -0.80902 - 0.58779 * I;
        known_roots[2] = 0.30902 + 0.95106 * I;
        known_roots[3] = 0.30902 - 0.95106 * I;
        known_roots[4] = -0.80902 + 0.58779 * I;
        break;
    case 6:
        iterate = iterate_6;
        known_roots[0] = -1.0 + 0.0 * I;
        known_roots[1] = 1.0 + 0.0 * I;
        known_roots[2] = -0.5 - 0.86603 * I;
        known_roots[3] = 0.5 + 0.86603 * I;
        known_roots[4] = 0.5 - 0.86603 * I;
        known_roots[5] = -0.5 + 0.86603 * I;
        break;
    case 7:
        iterate = iterate_7;
        known_roots[0] = 1.0 + 0.0 * I;
        known_roots[1] = -0.90097 - 0.43388 * I;
        known_roots[2] = 0.62349 + 0.78183 * I;
        known_roots[3] = -0.22252 - 0.97493 * I;
        known_roots[4] = -0.22252 + 0.97493 * I;
        known_roots[5] = 0.62349 - 0.78183 * I;
        known_roots[6] = -0.90097 + 0.43388 * I;
        break;
    case 8:
        iterate = iterate_8;
        known_roots[0] = -1.0 + 0.0 * I;
        known_roots[1] = 1.0 + 0.0 * I;
        known_roots[2] = 0.0 - 1.0 * I;
        known_roots[3] = 0.0 + 1.0 * I;
        known_roots[4] = -0.70711 - 70711 * I;
        known_roots[5] = 0.70711 + 0.70711 * I;
        known_roots[6] = 0.70711 - 0.70711 * I;
        known_roots[7] = -0.70711 + 0.70711 * I;
        break;
    case 9:
        iterate = iterate_9;
        known_roots[0] = 1.0 + 0.0 * I;
        known_roots[1] = -0.93969 - 0.34202 * I;
        known_roots[2] = 0.76604 + 0.64279 * I;
        known_roots[3] = -0.5 - 0.86603 * I;
        known_roots[4] = 0.17365 + 0.98481 * I;
        known_roots[5] = 0.17365 - 0.98481 * I;
        known_roots[6] = 0.5 + 0.86603 * I;
        known_roots[7] = 0.76604 - 0.64279 * I;
        known_roots[8] = -0.93969 + 0.34202 * I;
        break;
    case 10:
        iterate = iterate_10;
        known_roots[0] = -1.0 + 0.0 * I;
        known_roots[1] = 1.0 + 0.0 * I;
        known_roots[2] = -0.80902 - 0.58779 * I;
        known_roots[3] = 0.80902 + 0.58779 * I;
        known_roots[4] = -0.30902 - 0.95106 * I;
        known_roots[5] = 0.30902 + 0.95106 * I;
        known_roots[6] = 0.30902 - 0.95106 * I;
        known_roots[7] = -0.30902 + 0.95106 * I;
        known_roots[8] = 0.80902 - 0.58779 * I;
        known_roots[9] = -0.80902 + 0.58779 * I;
        break;
    }

    int ret, thread;
    pthread_t pthreads[threads];

    // create result matrices
    iters = calloc(lines, sizeof iters);
    roots = calloc(lines, sizeof roots);
    done = calloc(lines, sizeof done);

    // create threads
    pthread_mutex_init(&result_mutex, NULL);
    pthread_mutex_init(&done_mutex, NULL);

    for (thread = 0; thread < threads; thread++)
    {
        int *arg = malloc(4 * sizeof arg);
        arg[0] = thread;
        arg[1] = threads;

        if ((ret = pthread_create(pthreads + thread, NULL, compute_lines, (void *)arg)))
        {
            printf("Error creating thread: %d\n", ret);
            exit(1);
        }
    }

    pthread_t writer_thread;
    pthread_create(&writer_thread, NULL, writer_function, NULL);

    // join threads
    for (thread = 0; thread < threads; thread++)
    {
        if ((ret = pthread_join(pthreads[thread], NULL)))
        {
            printf("Error joining thread: %d\n", ret);
            exit(1);
        }
    }

    pthread_join(writer_thread, NULL);

    pthread_mutex_destroy(&result_mutex);
    pthread_mutex_destroy(&done_mutex);

    // free memory
    free(known_roots);
    for (int line = 0; line < lines; line++)
    {
        free(roots[line]);
        free(iters[line]);
    }
    free(roots);
    free(iters);
    free(done);

    return 0;
}
