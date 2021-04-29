import java.util.ArrayList;
import java.util.List;

public class Vertex {
	int vertex_index;
	List<Link> adjacent_links = new ArrayList<>();
	
	double current_cost;//used for finding the shortest path
	Link previous_link; //used for path finding
	List<Driver> drivers = new ArrayList<>();//the set of driver who starts from this vertex
	List<Rider> riders_start = new ArrayList<>();//the set of rider who starts from this vertex
	List<Rider> riders_end = new ArrayList<>(); // the set of rider who end with this vertex
	
	public Vertex(int vertex_index){
		this.vertex_index = vertex_index;
	}
}
