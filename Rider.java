
public class Rider {
	int rider_index;
	int start_vertex_index;
	int end_vertex_index;
	double demand;
	double satisfied_demand;
	double price;
	
	public Rider(int rider_index, int start_vertex_index, int end_vertex_index, double demand){
		this.rider_index = rider_index;
		this.start_vertex_index = start_vertex_index;
		this.end_vertex_index = end_vertex_index;
		this.demand = demand;
		this.satisfied_demand = 0;
		
		this.price = 55;
	}
}
