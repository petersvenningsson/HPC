#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>

MPI_Status status;
MPI_Win win;
MPI_Aint size;

#define master 0

int main(int argc, char **argv)
{
    // ########################################################################
    // init MPI
    // ########################################################################
    int number_of_tasks, task_id;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    // ########################################################################
    // master
    // ########################################################################
    if (task_id == 0)
    {
        // ########################################################################
        // parse arguments
        // ########################################################################

        // arguments
        int iterations = -1;
        double diff_const = -1;

        int width;
        int height;

        int worker;

        extern char *optarg;
        extern int optind;
        int c, err = 0;

        static char usage[] = "\nusage: %s -n iterations -d diff_const\n";

        while ((c = getopt(argc, argv, "n:d:")) != -1)
        {
            switch (c)
            {
            case 'n':
                iterations = atoi(optarg);
                break;
            case 'd':
                diff_const = atof(optarg);
                break;
            case '?':
                err = 1;
                break;
            }
        }

        if (iterations == -1)
        {
            printf("error: -n is required!\n");
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
        if (diff_const == -1)
        {
            printf("error: -d is required!\n");
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
        if (err)
        {
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }

        if (iterations < 1)
        {
            printf("error: n should be an positive integer!\n");
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
        if (diff_const < 0)
        {
            printf("error: d should be an positive float!\n");
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
        // argument parsing done!

        // ########################################################################
        // Read width, height and initial data from file
        // ########################################################################
        FILE *input_file;

        // input_file = fopen("test_input", "r");
        // input_file = fopen("/home/hpc2019/a4_grading/test_data/diffusion_100_100", "r");
        // input_file = fopen("/home/hpc2019/a4_grading/test_data/diffusion_10000_1000", "r");
        // input_file = fopen("/home/hpc2019/a4_grading/test_data/diffusion_100000_100", "r");
        input_file = fopen("diffusion", "r");

        char char_line[40];
        char *char_token;
        const char deliminator[2] = " ";

        // Read width and height
        fgets(char_line, 40, input_file);
        char_token = strtok(char_line, deliminator);
        height = atoi(char_token);
        char_token = strtok(NULL, deliminator);
        width = atoi(char_token);

        int number_of_rows_per_worker, start_row;
        number_of_rows_per_worker = height / number_of_tasks;

        int full_array_length = width * height + 2 * width + 2 * height + 4;

        size = full_array_length * sizeof(double);
        double *current_temperatures;
        double *new_temperatures;
        MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &current_temperatures, &win);
        MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &new_temperatures, &win);

        // Read initial values for temperature map
        int index_1;
        int index_2;
        double initial_value;

        while (fgets(char_line, 50, input_file))
        {
            fgets(char_line, 50, input_file);
            char_token = strtok(char_line, deliminator);
            index_2 = atoi(char_token);
            char_token = strtok(NULL, deliminator);
            index_1 = atoi(char_token);
            char_token = strtok(NULL, deliminator);
            initial_value = atof(char_token);
            current_temperatures[(index_1 + 1) * (width + 2) + index_2 + 1] = initial_value;
        }

        // ########################################################################
        // computation part below
        // ########################################################################
        // TODO: handle rest

        for (worker = 1, start_row = number_of_rows_per_worker; worker < number_of_tasks; worker++, start_row += number_of_rows_per_worker)
        {
            MPI_Send(&iterations, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
            MPI_Send(&diff_const, 1, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD);
            MPI_Send(&width, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
            MPI_Send(&height, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
            MPI_Send(&start_row, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
            MPI_Send(&number_of_rows_per_worker, 1, MPI_INT, worker, 1, MPI_COMM_WORLD);
        }

        for (int iter = 0; iter < iterations; iter++)
        {
            double p1, p2, p3, p4;

            for (int row = 0; row < number_of_rows_per_worker; row++)
            {
                int first_index = (row + 1) * (width + 2) + 1;
                for (int i = first_index; i < width + first_index; i++)
                {
                    p1 = current_temperatures[i - 1];
                    p2 = current_temperatures[i + 1];
                    p3 = current_temperatures[i - width - 2];
                    p4 = current_temperatures[i + width + 2];

                    new_temperatures[i] = diff_const * ((p1 + p2 + p3 + p4) / 4 - current_temperatures[i]) + current_temperatures[i];
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
            double *temp;
            temp = current_temperatures;
            current_temperatures = new_temperatures;
            new_temperatures = temp;
            MPI_Barrier(MPI_COMM_WORLD);
        }

        // for (size_t i = 1; i < (width + 1); i++)
        // {
        //     for (size_t j = 1; j < (height + 1); j++)
        //     {
        //         printf("%3e ", current_temperatures[i * (width + 2) + j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // Display the result to the screen
        long double sum = 0.0;
        for (size_t i = 1; i < (width + 1); i++)
        {
            for (size_t j = 1; j < (height + 1); j++)
            {
                sum += current_temperatures[i * (width + 2) + j];
            }
        }
        long double average = sum / (width * height);
        printf("%Lf\n", average);

        long double diff = 0.0;
        for (size_t i = 1; i < (width + 1); i++)
        {
            for (size_t j = 1; j < (height + 1); j++)
            {
                diff += fabs(average - current_temperatures[i * (width + 2) + j]);
            }
        }
        long double average_diff = diff / (width * height);
        printf("%Lf\n", average_diff);
    }

    // ########################################################################
    // slave
    // ########################################################################
    if (task_id > 0)
    {
        int iters, width, height, start_row, number_of_rows_per_worker;
        double diff_const;

        double *current_temperatures;
        double *new_temperatures;
        int disp_unit;
        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &current_temperatures, &win);
        MPI_Win_shared_query(win, 0, &size, &disp_unit, &current_temperatures);
        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &new_temperatures, &win);
        MPI_Win_shared_query(win, 0, &size, &disp_unit, &new_temperatures);

        MPI_Recv(&iters, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&diff_const, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&width, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&height, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&start_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&number_of_rows_per_worker, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        double p1, p2, p3, p4;

        for (int iter = 0; iter < iters; iter++)
        {
            for (int row = start_row; row < start_row + number_of_rows_per_worker; row++)
            {
                int first_index = (row + 1) * (width + 2) + 1;
                int first_return_index = (row - start_row) * (width + 2) + 1;

                for (int i = first_index, j = first_return_index; i < width + first_index; i++, j++)
                {
                    p1 = current_temperatures[i - 1];
                    p2 = current_temperatures[i + 1];
                    p3 = current_temperatures[i - width - 2];
                    p4 = current_temperatures[i + width + 2];

                    new_temperatures[i] = diff_const * ((p1 + p2 + p3 + p4) / 4 - current_temperatures[i]) + current_temperatures[i];
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            double *temp;
            temp = current_temperatures;
            current_temperatures = new_temperatures;
            new_temperatures = temp;
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Win_free(&win);
    MPI_Finalize();
    return 0;
}
