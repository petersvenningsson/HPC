## Relevant concepts

- **I/O bottlenecks**.  We have learned in previous assignments that writing and reading from files can be slow and 
  one wants to limit the number of read/writes to a file.
  
- **Data type conversion**.  Conversion from string to int can be costly and one has to choose the most effcient method here.
  
- **Multithreading**.  In this case, we will achieve this using OpenMP.

- **Synchronization**.  Or the lack thereof. We want to avoid explicit synchronization if possible.

- **Memory management**.  We will want to load as much of the input data into the memory as possible without breaking 
   the constraits given in the assignment. This to further minimize the read/write operations needed.

- **Data types**.  We will need to pay close attention to what data types we choose to limit memroy usage and increase 
   performance in calculations.

## Intended program layout

1. **Count the number of coordinates**:
First step is to count the number of coordinates. This will later be used to create arrays and block sizes (how  much we will load into the memory at once) of the correct size. We plan to implement this using a simple for loop over the number of characters in the input file and increment a counter each time we encounter a newline. Initial tests have shown this to be pretty fast for large files.

2. **Calculate block size**:
We will only load a certain number of coordinates into the memory at any given time. We do this to ensure that we won't go over the specified memory limitations given in the assignment. We set an upper limit to the block size, 10 000. With this limit, we will have a  maximum of 10 000 * 3 * 2 * 4 bytes for coordinates stored. This is well within the limit: 1024^3 - 10 000 * 3 * 2 * 4 = 1 073 501 824 bytes, so we will have plenty of memory left for other components. 4 here stands for the data type, int, which we will store each sub-coordinate with. A coordinate will be a 3-tuple of ints, hence the 3. And the 2 comes from that we will load 2 sets of 10 000 coordinates at the same time, more on this on the next step.

3. **Read coordinates to be computed**
We visualize all the possible combinations of coordinates to compute the distance between in a matrix. We only need to compute one half of this matrix and we choose to compute the upper half above the diagonal (not including the diagonal). We then take a block of coordinates on the x axis and a block of coordinates on the y axis of this matrix and read those into memory from the file. We read the coordinates from string and save them as integers, by skipping to read the dot. In effect, they become multiplied by 100 so we will have to keep this in mind for future calculations when we want to present the results.

4. **Calculate distances**
With some coordinates read into memory, we can then compute the distances between these. We will parallize this part using OpenMP. We compute the distance using the square root of the difference in each sub-coordinate.
All the distances are then saved to an array with the maximum length being equal to the maxium distance possible * 100. This means that we can save the frequencies of the distances in sort of a hash map which is very fast.

5. **Print the result**
Finally, we just have to loop through the array and print the results. The array is already sorted, per design. And we dont have to convert anything to floats, we should just be able printf them into the correct format. We will see!
