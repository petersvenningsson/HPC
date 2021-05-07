#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <string.h>
#define MAX_SOURCE_SIZE (0x100000)
#include <math.h>

// arguments
int iterations = -1;
double diff_const = -1;

int width;
int height;
int i;

double *current_temperatures;
double *new_temperatures;

int main(int argc, char **argv)
{
    // ########################################################################
    // parse arguments
    // ########################################################################
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
    input_file = fopen("diffusion", "r");

    char char_line[40];
    char *char_token;
    const char deliminator[2] = " ";

    // Read width and height
    fgets(char_line, 40, input_file);
    char_token = strtok(char_line, deliminator);
    width = atoi(char_token);
    char_token = strtok(NULL, deliminator);
    height = atoi(char_token);

    int array_length = width * height;
    double *current_temperatures = calloc(array_length, sizeof(double));
    double *new_temperatures = calloc(array_length, sizeof(double));

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
        current_temperatures[index_1 * width + index_2] = initial_value;
    }

    // ########################################################################
    // computation part below
    // ########################################################################

    FILE *fp;
    char *source_str;
    size_t source_size;

    fp = fopen("kernel.cl", "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char *)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    // Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, 1, &device_id, &ret_num_devices);

    // Create an OpenCL context
    cl_context context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    // Create memory buffers on the device for each vector
    // TODO: may increase perf with more specific memeory type
    cl_mem cur_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, array_length * sizeof(double), NULL, &ret);
    cl_mem new_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, array_length * sizeof(double), NULL, &ret);

    // Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);

    // Build the program
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);

    // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "heat_diff", &ret);

    // Set the arguments of the kernel
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&cur_mem_obj);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&new_mem_obj);
    ret = clSetKernelArg(kernel, 2, sizeof(double), &diff_const);
    ret = clSetKernelArg(kernel, 3, sizeof(int), &width);
    ret = clSetKernelArg(kernel, 4, sizeof(int), &height);

    // Copy the lists A and B to their respective memory buffers
    // TODO: maybe write/copy new to curr every iteration except the first?
    // DO ONCE
    ret = clEnqueueWriteBuffer(command_queue, cur_mem_obj, CL_TRUE, 0, array_length * sizeof(double), current_temperatures, 0, NULL, NULL);

    // Execute the OpenCL kernel on the list
    size_t global_item_size = array_length; // Process the entire matrix/array
    size_t local_item_size = 64;            // Process in groups of 64

    for (size_t i = 0; i < iterations; i++)
    {
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
        ret = clEnqueueCopyBuffer(command_queue, new_mem_obj, cur_mem_obj, 0, 0, array_length * sizeof(double), 0, NULL, NULL);
    }

    // DO ONCE
    // Read the memory buffer C on the device to the local variable C
    ret = clEnqueueReadBuffer(command_queue, new_mem_obj, CL_TRUE, 0, array_length * sizeof(double), new_temperatures, 0, NULL, NULL);

    // Display the result to the screen
    long double sum = 0.0;
    for (size_t i = 0; i < array_length; i++)
    {
        sum += new_temperatures[i];
    }
    long double average = sum / array_length;
    printf("%Lf\n", average);

    long double diff = 0.0;
    for (size_t i = 0; i < array_length; i++)
    {
        diff += fabs(average - new_temperatures[i]);
    }
    long double average_diff = diff / array_length;
    printf("%Lf\n", average_diff);

    // Clean up
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(cur_mem_obj);
    ret = clReleaseMemObject(new_mem_obj);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);
    free(current_temperatures);
    free(new_temperatures);

    return 0;
}
