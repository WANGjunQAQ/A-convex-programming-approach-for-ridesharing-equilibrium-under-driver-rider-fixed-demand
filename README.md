# A-convex-programming-approach-for-ridesharing-equilibrium-under-driver-rider-fixed-demand
This repository is created for storing the data and code of an algorithm implemented to solve ridesharing equilibrium under driver/rider fixed demand. In order to verify the effectiveness and efficiency of this approach, several experiments with varied network scale are conducted. 

In the folder named data, all input of the algorithm are included: net, origin_driver_OD, origin_rider_OD. And for each experiment, the information of net structure are listed in the net.csv. In addition, the origin driver and rider demand could be found in the origin_driver_OD.csv and origin_rider_OD.csv respectively.

The resulting folder denoted with result contains the output of the algorithm. The file named with driver_rider_flow distribution shows the flow distribution among drivers and riders, and the pattern are of the following form:

                            solo driver flow    rider1    rider2    rider3
                            
                driver1
                
                driver2
                
                driver3
                

epsilon1.csv and epsilon2.csv, these two files keep tracks of the running process of the algorithm. And they are of very similar form: both of them have just two columns, the first lists the iteration number of the outer loop, and the latter corresponds the value of epsilon associated the outer iteration.

riders_price.csv just shows the resulting prices that riders have to pay to drivers as a compensation for detour time and some other inconvenience.

run_duration.csv keep a record of total running time of whole algorithm. 

tolerance.csv just stores the parameters setting of the algorithm.

The folder JAVA_CODE is the java code to implement the algorithm.
