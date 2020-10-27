
public class Link {
	int link_index;
	Vertex start_vertex;
	Vertex end_vertex;
	int type;
	int original_link_index;
	double length;
	double capacity;
	double free_flow_time;
	double flow;
	double current_cost;
	
	public Link(int link_index){
		this.link_index = link_index;
		this.flow = 0;
	}
}
