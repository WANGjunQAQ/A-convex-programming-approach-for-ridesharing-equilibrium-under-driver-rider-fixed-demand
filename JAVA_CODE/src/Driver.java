
public class Driver {
	int driver_index;
	int start_vertex_index;
	int end_vertex_index;
	double demand;
	double[] flow_distribution; //0 for solo driver, 1 for the rider_index 1, ...
	
	//used for the imporved direction finding
	int shared_rider; //0 denoted for solo driver, 1 for the rider_index 1,...	
	
	public Driver(int driver_index, int start_vertex_index, int end_vertex_index, double demand, int rider_number){
		this.driver_index = driver_index;
		this.start_vertex_index = start_vertex_index;
		this.end_vertex_index = end_vertex_index;
		this.demand = demand;
		this.flow_distribution = new double[rider_number + 1];
		for(int i = 0; i < rider_number + 1; i ++){
			this.flow_distribution[i] = 0;
		}
	}
}
