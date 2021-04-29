import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;


public class Graph {
	int vertex_number;
	int link_number;
	List<Vertex> vertexes = new ArrayList<>();
	List<Link> links = new ArrayList<>();
	
	PriorityQueue<Vertex> queue = new PriorityQueue<>(idComparator);
	
	public Graph(int node_number, int link_number){
		this.vertex_number = node_number;
		this.link_number = link_number;
		
		for(int i = 0; i < this.vertex_number; i++){
			Vertex nodei = new Vertex(i + 1);
			vertexes.add(nodei);
		}
		
		for(int i = 0; i < this.link_number; i ++){
			Link linki = new Link(i + 1);
			links.add(linki);
		}
	}
	
	//establish the network
	public void establish_network(double[][] net){
		for(int i = 0; i < net.length; i ++){
			int start_node_index = (int)(net[i][0] - 1);
			int end_node_index = (int)(net[i][1] - 1);
			double capacity = net[i][2];
			double length = net[i][3];
			double free_flow_time = net[i][4];
			//int link_index = (int)(net[i][5] - 1);
			Link linki = this.links.get(i);
			Vertex start_node = this.vertexes.get(start_node_index);
			Vertex end_node = this.vertexes.get(end_node_index);
			linki.start_vertex = start_node;
			linki.end_vertex = end_node;
			linki.length = length;
			linki.capacity = capacity;
			//linki.free_flow_time = linki.length;  // this setting may be changed later
			linki.free_flow_time = free_flow_time; //for simple network
			linki.current_cost = linki.free_flow_time;
			start_node.adjacent_links.add(linki);
		}
	}
	
	public static Comparator<Vertex> idComparator= new Comparator<Vertex>(){
		public int compare(Vertex e1,Vertex e2){
			if(e1.current_cost-e2.current_cost>0){
				return 1;
			}else if(e1.current_cost==e2.current_cost){
				return 0;
			}else{
				return -1;
			}
		}
	};
	
	
	//graph operations------------------------------------------------------------
	//initialize the graph for the purpose of path finding
	public void init_graph(){
		for(int i = 0; i< this.vertex_number; i++){
			this.vertexes.get(i).current_cost = Double.MAX_VALUE;
		}
	}
	
	//find the shortest path from one vertex to all other vertexes by Dijkstra algorithm
	public void find_shortest_path_from_one_vertex(int vertex_index){
		init_graph();
		Vertex init_vertex = this.vertexes.get(vertex_index - 1);
		init_vertex.current_cost = 0;
		init_vertex.previous_link = null;
		this.queue.add(init_vertex);
		
		while(!queue.isEmpty()){
			Vertex current_vertex = queue.poll();
			double current_cost = current_vertex.current_cost;
			List<Link> links = current_vertex.adjacent_links;
			for(Link link : links){
				Vertex end_vertex = link.end_vertex;
				if(current_cost + link.current_cost < end_vertex.current_cost){
					end_vertex.current_cost = current_cost + link.current_cost;
					end_vertex.previous_link = link;
					
					if(!queue.contains(end_vertex)){
						queue.add(end_vertex);
					}
				}
			}
		}
	}
}
