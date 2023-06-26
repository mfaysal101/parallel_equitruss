# parallel-equitruss

# Software requirements/built and executed in
(1) Minimum CMake VERSION 3.0, the CMake version in our local compute node is 3.10.2

(2) g++ version 7.5.0 or higher

(3) Linux kernel version "5.4.0-147-generic"

(4) OS "Ubuntu 18.04.6 LTS"

# Pulling the code from the GitHub repository
(1) Pull the code from the main branch of the following code repository "https://github.com/mfaysal101/parallel_equitruss"

(2) There is a directory name gapbs-agile which is empty, please navigate to that directory and pull the code from the branch named "equitruss_infrastructure" from the same GitHub repository

(3) Please make sure the directory hierarchy is maintained, inside the gapbs-agile directory, all the source files and directories should be as seen from "https://github.com/mfaysal101/parallel_equitruss/tree/equitruss_infrastructure". There should not be any other directory in between them, e.g., the "src" directory should be gapbs-agile/src and not like gapbs-agile/other_directory/src.


# How to prepare the datasets for the experiments
Input datasets collected from "https://snap.stanford.edu/data/index.html"
(1) Download the undirected networks from SNAP in the table titled "Networks with ground-truth communities"

(2) For testing purposes, smaller networks, e.g., com-Amazon or com-DBLP are suggested

(3) Getting rid of any comments (if exist) in the network file is a must.

(4) If the starting vertex in the network is "1" please renumber the vertices to start from "0". You may use the to_zero_index.cpp script to do that under the "util" directory. After building the file in g++ using the following command "g++ to_zero_index.cpp -o myexec", run the command "./myexec input_network_name output_network_name".

(5) For the network file, generate their k-truss-ness using any k-truss decomposition application out there and write them in a file in the following format "u  v   k" where, (u,v) is an edge and u < v and k is the trussness.

(6) Alternatively, you can use the sample dataset available from the following repository "https://github.com/mfaysal101/sample_equitruss_dataset" and skip steps 1-5

# How to build the codebase
(1) Navigate into the directory where "parallel_equitruss" is stored in a terminal

(2) Make a directory named "build" by issuing the command "mkdir build"

(3) Navigate into the build directory by "cd build"

(4) Issue the following command to create a Make file by CMake "cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=On .."

(5) Then issue the command "make"

(6) This will build and link the project

# How to execute the code
(1) Following is a simple bash script to execute the code with 32 OpenMP threads. Copy and paste the following script in a file with any name (myscript), save it with execution permission, or (chmod +x myscript). Please note the extension of the network file, change the extension of the edge list network
file downloaded from SNAP from .txt to .el.
" !/bin/bash
export OMP_NUM_THREADS=32
./EquiTruss -sf path-2-input-graph/com-amazon.el -t path-2-ktruss-file/trussd_amazon.txt > amazon-output-32-threads.txt "

(2) The output file "amazon-output-32-threads.txt" should have all the timing information for running different parts of the code (SV connected component to create super node, Afforest connected components to create super node, Super Edge creation time, etc.).

(3) Depending on the shared memory platform used, the experiment should finish within 2 minutes for the Amazon network. Actual timing information for different parts of the code can be found in the generated output file mentioned above.

(4) In the output file generated above, there should be a bunch of printed statements at the end of the file describing the time taken for different components of the application.

(5) To interpret those time information, let us assume "totalExecutionTime" has line number 1.

(6) Since there are 3 different phases/versions of our implementation and all of them are put together in the same application, to get the first version called Baseline EquiTruss (in the paper), add line# (5+6+7+8+9+12+15+17+19+20) to get the total execution time.

(7) Similarly, for the version called Optimal EquiTruss (in the paper), add line# (5+6+7+10+13+15+17+19+20), and for the version called Afforest EquiTruss, add line#
(5+6+7+11+13+15+17+19+20) to get the total execution times for the respective versions.

(8) These results are illustrated for different numbers of threads in Figures 7, 9, and 10 for larger networks such as YouTube,
LiveJournal, and Orkut. 


For further clarification or any difficulty in following the steps
mentioned above, please reach out to the following email address
"faysal@unlv.nevada.edu".
