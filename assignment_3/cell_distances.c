#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// arguments
int threads = -1;
int number_of_coords = 0;
int block_size;

short *coords_1;
short *coords_2;

int *distances;

short raw_coord_to_int(const char *p)
{
    short x;
    x = (p[1] - '0') * 10000 + (p[2] - '0') * 1000 + (p[4] - '0') * 100 + (p[5] - '0') * 10 + (p[6] - '0');
    if (*p == '-')
    {
        return -x;
    }
    return x;
}

int main(int argc, char **argv)
{
    // ########################################################################
    // parse arguments
    // ########################################################################
    extern char *optarg;
    extern int optind;
    int c, err = 0;

    static char usage[] = "\nusage: %s -t threads\n";

    while ((c = getopt(argc, argv, "t:")) != -1)
    {
        switch (c)
        {
        case 't':
            threads = atoi(optarg);
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
    if (err)
    {
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }

    if (threads < 1)
    {
        printf("error: t should be an positive integer!\n");
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }
    // argument parsing done!

    // ########################################################################
    // computation part below
    // ########################################################################
    omp_set_num_threads(threads);

    FILE *input_file;

    // input_file = fopen("/home/hpc2019/a3_grading/test_data/cell_e4", "r");
    // input_file = fopen("/home/hpcuser004/assignment_3/test_input", "r");
    input_file = fopen("cells", "r");

    char ch;
    while (!feof(input_file))
    {
        ch = fgetc(input_file);
        if (ch == '\n')
        {
            number_of_coords++;
        }
    }

    block_size = 10000;

    char *raw_coords_1 = malloc(sizeof(char) * block_size * 24);
    char *raw_coords_2 = malloc(sizeof(char) * block_size * 24);
    coords_1 = malloc(sizeof(short) * block_size * 3);
    coords_2 = malloc(sizeof(short) * block_size * 3);

    int i = 0, j = 0;
    int k = 0, l = 0;
    int coord;

    distances = calloc(3465, sizeof distances);

    for (i = 0; i < number_of_coords; i += block_size)
    {
        for (j = i; j < number_of_coords; j += block_size)
        {
            fseek(input_file, i * 24, SEEK_SET);
            fread(raw_coords_1, sizeof(char), 24 * block_size, input_file);

            fseek(input_file, j * 24, SEEK_SET);
            fread(raw_coords_2, sizeof(char), 24 * block_size, input_file);

            for (coord = 0; coord < block_size * 3; coord += 3)
            {
                coords_1[coord] = raw_coord_to_int(raw_coords_1 + 8 * coord);
                coords_1[coord + 1] = raw_coord_to_int(raw_coords_1 + 8 + 8 * coord);
                coords_1[coord + 2] = raw_coord_to_int(raw_coords_1 + 16 + 8 * coord);

                coords_2[coord] = raw_coord_to_int(raw_coords_2 + 8 * coord);
                coords_2[coord + 1] = raw_coord_to_int(raw_coords_2 + 8 + 8 * coord);
                coords_2[coord + 2] = raw_coord_to_int(raw_coords_2 + 16 + 8 * coord);
            }

            {
                #pragma omp parallel for shared(coords_1, coords_2) reduction(+:distances[:3465])
                for (k = 0; k < block_size; k++)
                {
                    if (i == j)
                    {
                        for (l = k + 1; l < block_size; l++)
                        {
                            int index1 = 3 * l;
                            int index2 = 3 * k;

                            int diff_1 = coords_1[index1] - coords_2[index2];
                            int diff_2 = coords_1[index1 + 1] - coords_2[index2 + 1];
                            int diff_3 = coords_1[index1 + 2] - coords_2[index2 + 2];

                            float temp = sqrtf(diff_1 * diff_1 + diff_2 * diff_2 + diff_3 * diff_3);
                            int distance = (int)((temp * 0.1) + 0.5);

                            distances[distance] += 1;
                        }
                    }
                    else
                    {
                        for (l = 0; l < block_size; l++)
                        {
                            int index1 = 3 * l;
                            int index2 = 3 * k;

                            int diff_1 = coords_1[index1] - coords_2[index2];
                            int diff_2 = coords_1[index1 + 1] - coords_2[index2 + 1];
                            int diff_3 = coords_1[index1 + 2] - coords_2[index2 + 2];

                            float temp = sqrtf(diff_1 * diff_1 + diff_2 * diff_2 + diff_3 * diff_3);
                            int distance = (int)((temp * 0.1) + 0.5);

                            distances[distance] += 1;
                        }
                    }
                }
            }
        }
    }

    int distance;
    int number;
    int rest;

    for (int i = 0; i < 3465; i++)
    {
        distance = distances[i];
        if (distance != 0)
        {
            number = i / 100;
            rest = i % 100;
            printf("%02d.%02d %d\n", number, rest, distance);
        }
    }
}
