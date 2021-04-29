import java.util.ArrayList;
import java.util.List;

public class Solution {
	int vertex_number;
	int link_number;
	int driver_number;
	int rider_number;
	List<Driver> drivers = new ArrayList<> ();
	List<Rider> riders = new ArrayList<> ();
	Graph graph;
	
	double t0;
	double delta;
	
	//derived value
	double total_cost;
	
	public Solution(int vertex_number, int link_number, int driver_number, int rider_number, IO io){
		this.vertex_number = vertex_number;
		this.link_number = link_number;
		this.driver_number = driver_number;
		this.rider_number = rider_number;
		this.graph = new Graph(this.vertex_number, this.link_number);
		double[][] net = io.get_net_data();
		graph.establish_network(net);
		
		initialize_driver(io);
		initialize_rider(io);
	}
	
	public void initialize_driver(IO io){
		double[][] driver_od = io.get_driver_OD();
		for(int i = 0; i < this.driver_number; i ++){
			int start_vertex_index = (int)driver_od[i][0];
			int end_vertex_index = (int)driver_od[i][1];
			double demand = driver_od[i][2];
			Driver driver = new Driver(i + 1, start_vertex_index, end_vertex_index, demand, this.rider_number);
			this.drivers.add(driver);
			graph.vertexes.get(start_vertex_index - 1).drivers.add(driver); 
		}
	}
	
	public void initialize_rider(IO io){
		double[][] rider_od = io.get_rider_OD();
		for(int i = 0; i < this.rider_number; i ++){
			int start_vertex_index = (int)rider_od[i][0];
			int end_vertex_index = (int)rider_od[i][1];
			double demand = rider_od[i][2];
			Rider rider = new Rider(i + 1, start_vertex_index, end_vertex_index, demand);
			this.riders.add(rider);
			graph.vertexes.get(start_vertex_index - 1).riders_start.add(rider);
			graph.vertexes.get(end_vertex_index - 1).riders_end.add(rider);
		}
	}
	
	//get shortest path between all possible vertex combination
	public Path[] get_all_shortest_path(){
		Path[][] driver_start_rider_start = new Path[this.driver_number][this.rider_number];
		Path[][] rider_end_driver_end = new Path[this.rider_number][this.driver_number];
		Path[] drivers = new Path[this.driver_number];
		Path[] riders = new Path[this.rider_number]; //since link has its own direction, the path for each rider need to be obtained separately
		
		for(int i = 0; i < this.vertex_number; i ++){
			Vertex vertexi = this.graph.vertexes.get(i);
			if(vertexi.drivers.isEmpty() && vertexi.riders_end.isEmpty() && vertexi.riders_start.isEmpty()){
				continue;
			}
			this.graph.find_shortest_path_from_one_vertex(i + 1);
			
			if(!vertexi.drivers.isEmpty()){
				for(Driver driveri : vertexi.drivers){
					int driver_index = driveri.driver_index;
					for(int riderj = 0; riderj < this.rider_number; riderj ++){
						int riderj_start_vertex_index = this.riders.get(riderj).start_vertex_index;
						Vertex riderj_start_vertex = this.graph.vertexes.get(riderj_start_vertex_index - 1);
						
						Path pathij = new Path(i + 1, riderj_start_vertex_index);
						Vertex current_vertex = riderj_start_vertex;
						while(current_vertex.previous_link != null){
							pathij.links.add(current_vertex.previous_link);
							current_vertex = current_vertex.previous_link.start_vertex;
						}
						pathij.path_cost = riderj_start_vertex.current_cost;
						driver_start_rider_start[driver_index - 1][riderj] = pathij;
					}
					
					//calculate the path cost of solo driver
					int driveri_end_vertex_index = driveri.end_vertex_index;
					Vertex driveri_end_vertex = this.graph.vertexes.get(driveri_end_vertex_index - 1);
					Path path_driveri = new Path(i + 1, driveri_end_vertex_index);
					Vertex current_vertex_driveri = driveri_end_vertex;
					while(current_vertex_driveri.previous_link != null){
						path_driveri.links.add(current_vertex_driveri.previous_link);
						current_vertex_driveri = current_vertex_driveri.previous_link.start_vertex;
					}
					path_driveri.path_cost = driveri_end_vertex.current_cost;
					drivers[driver_index - 1] = path_driveri;
				}
			}
			if(!vertexi.riders_end.isEmpty()){
				for(Rider rideri : vertexi.riders_end){
					int rider_index = rideri.rider_index;
					for(int driverj = 0; driverj < this.driver_number; driverj ++){
						int driverj_end_vertex_index = this.drivers.get(driverj).end_vertex_index;
						Vertex driverj_end_vertex = this.graph.vertexes.get(driverj_end_vertex_index - 1);
						
						Path pathij_rd = new Path(i + 1, driverj_end_vertex_index);
						Vertex current_vertex_rd = driverj_end_vertex;
						while(current_vertex_rd.previous_link != null){
							pathij_rd.links.add(current_vertex_rd.previous_link);
							current_vertex_rd = current_vertex_rd.previous_link.start_vertex;
						}
						pathij_rd.path_cost = driverj_end_vertex.current_cost;
						rider_end_driver_end[rider_index - 1][driverj] = pathij_rd;
					}					
				}				
			}
			if(!vertexi.riders_start.isEmpty()){
				for(Rider rideri_s : vertexi.riders_start){
					//calculate the path cost of single rider
					int rider_index = rideri_s.rider_index;
					int rideri_end_vertex_index = rideri_s.end_vertex_index;
					Vertex rideri_end_vertex = this.graph.vertexes.get(rideri_end_vertex_index - 1);
					Path path_rideri = new Path(rideri_end_vertex_index, i + 1);
					Vertex current_vertex_rideri = rideri_end_vertex;
					while(current_vertex_rideri.previous_link != null){
						path_rideri.links.add(current_vertex_rideri.previous_link);
						current_vertex_rideri = current_vertex_rideri.previous_link.start_vertex;
					}
					path_rideri.path_cost = rideri_end_vertex.current_cost;
					riders[rider_index - 1] = path_rideri;
				}
			}
		}
		
		//find a path for each driver
		Path[] driver_pathes = new Path[this.driver_number];
		int index = 0;
		for(Driver driveri : this.drivers){
			Path driveri_path = get_shortest_path_single_driver(driveri, driver_start_rider_start, rider_end_driver_end, drivers, riders);
			driver_pathes[index] = driveri_path;
			index = index + 1;
		}
		return driver_pathes;
	}
	
	//get the shortest path between driver and rider
	public Path get_shortest_path_single_driver(Driver driver, Path[][] driver_start_rider_start, Path[][] rider_end_driver_end, Path[] drivers, Path[] riders){
		int driver_index = driver.driver_index;
		double solo_driver_cost = drivers[driver_index - 1].path_cost;
		double min_path_cost = solo_driver_cost;
		int rider_index = 0;//0 denote solo driver
		for(int rideri = 0; rideri< this.rider_number; rideri ++){
			double path_cost = driver_start_rider_start[driver_index - 1][rideri].path_cost + riders[rideri].path_cost + rider_end_driver_end[rideri][driver_index - 1].path_cost;
			path_cost = path_cost + delta + 2 * t0 - this.riders.get(rideri).price;
			
			if(path_cost < min_path_cost){
				rider_index = rideri + 1;
				min_path_cost = path_cost;
			}
		}
		
		if(rider_index == 0){
			Path path_driver_solo = drivers[driver_index - 1];
			path_driver_solo.solo_rider_indentifier = 0;
			return path_driver_solo;
		}else{
			Path path_driver_rider = new Path(driver.start_vertex_index, driver.end_vertex_index);
			path_driver_rider.solo_rider_indentifier = rider_index;
			for(Link link : driver_start_rider_start[driver_index - 1][rider_index - 1].links){
				path_driver_rider.links.add(link);
			}
			for(Link link : riders[rider_index - 1].links){
				path_driver_rider.links.add(link);
			}
			for(Link link : rider_end_driver_end[rider_index - 1][driver_index - 1].links){
				path_driver_rider.links.add(link);
			}
			path_driver_rider.path_cost = min_path_cost;
			return path_driver_rider;
		}
	}
	
	//before traffic assignment, we should initialize the link flow of the solution
	public void initialize_link_rider_flow(){
		for(Link linki : this.graph.links){
			linki.flow = 0;
		}
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			this.riders.get(rideri).satisfied_demand = 0;
		}
	}
	
	//put all the flow on the given shortest path for each driver
	public void get_improve_direction(Path[] driver_path){
		for(int pi = 0; pi < this.driver_number; pi ++){
			Path pathi = driver_path[pi];
			Driver driveri = this.drivers.get(pi);
			driveri.shared_rider = pathi.solo_rider_indentifier;
			if(driveri.shared_rider != 0){
				Rider shared_rider = this.riders.get(pathi.solo_rider_indentifier - 1);
				shared_rider.satisfied_demand = shared_rider.satisfied_demand + driveri.demand;
			}
			
			for(Link link : pathi.links){
				link.flow = link.flow + driveri.demand;
			}
			
			driveri.current_shortest_path = pathi;
		}
	}
	
	//get the total cost of the current solution
	public void get_total_cost(){
		this.total_cost = 0;
		for(Link link : this.graph.links){
			this.total_cost = this.total_cost + 4 * (link.free_flow_time * link.flow + 0.03 * link.free_flow_time * Math.pow(link.flow, 5)/Math.pow(link.capacity, 4));
		}
		for(Driver driveri : this.drivers){
			this.total_cost = this.total_cost + (2 * this.t0 + this.delta) * (driveri.demand - driveri.flow_distribution[0]);
		}
		for(Rider rideri : this.riders){
			this.total_cost = this.total_cost + rideri.price * (rideri.demand - rideri.satisfied_demand);
		}
	}
	
	public double get_total_cost_new(){
		double inner_total_cost = 0;
		for(Link link : this.graph.links){
			inner_total_cost = inner_total_cost + link.flow * link.current_cost;
		}
		for(Driver driveri : this.drivers){
			inner_total_cost = inner_total_cost + (2 * this.t0 + this.delta) * (driveri.demand - driveri.flow_distribution[0]);
		}
		for(Rider rideri : this.riders){
			inner_total_cost = inner_total_cost - rideri.price * rideri.satisfied_demand;
		}
		return inner_total_cost;
	}
	
	public double get_total_cost_for_improved_direction_only(){
		double inner_total_cost_for_improved_direction = 0;
		for(Driver driveri : this.drivers){
			
			inner_total_cost_for_improved_direction = inner_total_cost_for_improved_direction + driveri.demand * driveri.current_shortest_path.path_cost;
		}
		
		return inner_total_cost_for_improved_direction;
	}
	
	//calculate cost for link under current link flow
	public void calculate_link_cost(){
		for(Link link : this.graph.links){
			link.current_cost = 4 * link.free_flow_time * (1 + 0.15 * Math.pow(link.flow/link.capacity, 4));
		}
	}
}
