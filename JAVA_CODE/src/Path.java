import java.util.ArrayList;
import java.util.List;

public class Path {
	int start_vertex_index;
	int end_vertex_index;
	List<Link> links = new ArrayList<> ();
	double path_cost;
	int solo_rider_indentifier;//0 for solo driver, 1 for ridesharing with rider 1,...
	
	public Path(int start_vertex_index, int end_vertex_index){
		this.start_vertex_index = start_vertex_index;
		this.end_vertex_index = end_vertex_index;
	}
}
