__kernel void heat_diff(__global double *cur_temps, __global double *new_temps, double diff_const, int width, int height) {

    // Get the index of the current element
    int i = get_global_id(0);

    double p1;
    double p2;
    double p3;
    double p4;

    // Do the operation
    // p1 = min(i % width, 1) * cur_temps[i - 1];
    // p2 = min(i % width - 1, 1) * cur_temps[i + 1];

    // p3 = min(i % height, 1) * cur_temps[i - width];
    // p4 = min() * cur_temps[i + width];

    if(i % width == 0)
    {
        p1 = 0;
    } else {
        p1 = cur_temps[i - 1];
    }

    if((i + 1) % width == 0)
    {
        p2 = 0;
    } else {
        p2 = cur_temps[i + 1];
    }

    if(i < width)
    {
        p3 = 0;
    } else {
        p3 = cur_temps[i - width];
    }

    if(i >= (width * height - width))
    {
        p4 = 0;
    } else {
        p4 = cur_temps[i + width];
    }

    new_temps[i] = diff_const * ((p1 + p2 + p3 + p4) / 4 - cur_temps[i]) + cur_temps[i];
}
