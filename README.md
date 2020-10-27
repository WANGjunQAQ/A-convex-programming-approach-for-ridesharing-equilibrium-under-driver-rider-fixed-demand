# A-convex-programming-approach-for-ridesharing-equilibrium-under-driver-rider-fixed-demand
This repository is created for storing the data and code of an algorithm implemented to solve ridesharing equilibrium under driver/rider fixed demand. In order to verify the effectiveness and efficiency of this approach, several experiments with varied network scale are conducted. 

In the folder named data, all input of the algorithm are included: net, origin_driver_OD, origin_rider_OD. And for each experiment, the information of net structure are listed in the net.csv. In addition, the origin driver and rider demand could be found in the origin_driver_OD.csv and origin_rider_OD.csv respectively.

The resulting folder denoted with result contains the output of the algorithm. The file named with driver_rider_flow distribution shows the flow distribution among drivers and riders, and the pattern are of the following form:

                            solo driver flow    rider1    rider2    rider3
                            
                driver1
                
                driver2
                
                driver3
                

