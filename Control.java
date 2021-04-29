import java.util.List;

public class Control {
	//-----------------------------------------background-----------------------------------------
	//simple_network : vertex_number = 3, link_number = 3, driver_number = 3, rider_number = 3
	//siousxFall: vertex_number = 24, link_number = 76, driver_number = 20, rider_number = 20
	//chicago: vertex_number = 933, link_number = 2950, driver_number = 100, rider_number = 100
	//Anaheim: vertex_number = 416, link_number = 914, driver_number = 50, rider_number = 50
	String[] case_names = {"simple_network", "siouxFalls", "chicago", "chicago_200", "chicago_500", "Anaheim", "Anaheim1"};
	int[] vertex_numbers = {3, 24, 933, 933, 933, 416, 416};
	int[] link_numbers = {3, 76, 2950, 2950, 2950, 914, 914};
	int[] driver_od_numbers = {3, 20, 100, 200, 500, 50, 50};
	int[] rider_od_numbers = {3, 20, 100, 200, 500, 50, 50};
	
	int case_index = 3; 
	//0:simple_network; 1: siouxall; 2: chicago; 3: chicago_200; 4: chicago_500; 5: Anaheim;
	
	int vertex_number = vertex_numbers[case_index];
	int link_number = link_numbers[case_index];
	int driver_number = driver_od_numbers[case_index];
	int rider_number = rider_od_numbers[case_index];	
	String case_name = case_names[case_index];
	
	//input filename	
	String filename_net = "D:\\document\\equalibrium\\data5\\" + case_name + "\\net.csv";
	String filename_driver = "D:\\document\\equalibrium\\data5\\" + case_name + "\\origin_driver_OD.csv";
	String filename_rider = "D:\\document\\equalibrium\\data5\\" + case_name + "\\origin_rider_OD.csv";
	
	//output filename
	String filename_driver_rider_flow_distribution = "D:\\document\\equalibrium\\result6\\" + case_name + "\\driver_rider_flow_distribution.csv";
	String filename_riders_price = "D:\\document\\equalibrium\\result6\\" + case_name + "\\riders_price.csv";
	String filename_run_duration = "D:\\document\\equalibrium\\result6\\" + case_name + "\\run_duration.csv";
	String filename_epsilon2 = "D:\\document\\equalibrium\\result6\\" + case_name + "\\epsilon2.csv";
	String filename_epsilon1 = "D:\\document\\equalibrium\\result6\\" + case_name + "\\epsilon1.csv";
	String filename_theta_p = "D:\\document\\equalibrium\\result6\\" + case_name + "\\theta_p.csv";
	String filename_tolerance = "D:\\document\\equalibrium\\result6\\" + case_name + "\\tolerance.csv";
	String filename_net_flow = "D:\\document\\equalibrium\\result6\\" + case_name + "\\net_flow.csv";
	String filename_final_solution_converge_value = "D:\\document\\equalibrium\\result6\\" + case_name + "\\converge_value.csv";
	
	//parameter setting
	double t0 = 2;
	double delta = 5;
	
	//parameter control
	double inner_k_scale_up = 2;
	int inner_maximum_iteration_for_single_solution = 30;
	int inner_maximum_line_search_iteration_number = 100;
	double inner_tolerance = 0.0005;
	double alpha_tolerance = 0.001;
	
	double demand_gap_tolerance = 0.005;//associated with epsilon1
	
	double out_alpha_hat = 2;
	double out_k = 2;
	double out_tolerance = 0.0005; //associated with epsilon2
	
	//temporary variable
	double total1 = 0;
	
	//objective variable
	Solution solution1;
	Solution solution2;
	Solution solution3;
	Solution solution_temp;
	IO io;
	
	//result and parameter needs to be printed out
	long run_duration;
	double[][] epsilon2_array = new double[3000][3]; //column names = ['out_iteration', 'absolute value of gap', 'ebusilon2']
	int epsilon2_array_index = 0;
	
	double[][] epsilon1_array = new double[3000][2];
	int epsilon1_array_index = 0;
	
	double[][] theta_p_array = new double[3000][3];
	int theta_p_array_index = 0;
	
	public Control(){		
		//initialize IO
		this.io = new IO(this.vertex_number, this.link_number, this.driver_number, this.rider_number);
		io.filename_net = this.filename_net;
		io.filename_driver = this.filename_driver;
		io.filename_rider = this.filename_rider;
		io.filename_driver_rider_flow_distribution = this.filename_driver_rider_flow_distribution;
		io.filename_riders_price = this.filename_riders_price;
		io.filename_run_duration = this.filename_run_duration;
		io.filename_epsilon2 = this.filename_epsilon2;
		io.filename_epsilon1 = this.filename_epsilon1;
		io.filename_theta_p = this.filename_theta_p;
		io.filename_tolerance = this.filename_tolerance;
		io.filename_net_flow = this.filename_net_flow;
		io.filename_final_solution_converge_value = this.filename_final_solution_converge_value;
		
		io.print_all_tolerance(this.inner_tolerance, this.demand_gap_tolerance, this.out_tolerance);
		
		this.solution1 = new Solution(this.vertex_number, this.link_number, this.driver_number, this.rider_number, io);
		this.solution2 = new Solution(this.vertex_number, this.link_number, this.driver_number, this.rider_number, io);
		this.solution3 = new Solution(this.vertex_number, this.link_number, this.driver_number, this.rider_number, io);
		this.solution_temp = new Solution(this.vertex_number, this.link_number, this.driver_number, this.rider_number, io);
		solution1.t0 = this.t0;
		solution2.t0 = this.t0;
		solution3.t0 = this.t0;
		solution_temp.t0 = this.t0;
		solution1.delta = this.delta;
		solution2.delta = this.delta;
		solution3.delta = this.delta;
		solution_temp.delta = this.delta;
	}
	
	public static void main(String[] args){
		long start_time = System.currentTimeMillis();
		
		Control control = new Control();
		control.get_improved_solution_direction(control.solution1);
		control.get_improved_solution_direction(control.solution2);
		control.get_improved_solution_direction(control.solution3);
		control.get_improved_solution_direction(control.solution_temp);
		
		control.initialize_initial_solution(control.solution1);
		control.initialize_initial_solution(control.solution2);
		control.initialize_initial_solution(control.solution3);
		control.initialize_initial_solution(control.solution_temp);
		
		control.outer_loop(control.solution1, control.solution2, control.solution_temp);
		long end_time = System.currentTimeMillis();
		control.run_duration = end_time - start_time;
		//control.simple_network_check(control.solution1);
		control.solution1.calculate_link_cost();
		control.copy_solution1_to_solution2_price_total_cost(control.solution1, control.solution3);
		control.result_handle(control.io);
	}
	
	public void outer_loop(Solution solution1, Solution solution2, Solution solution_temp){
		int out_iteration = 0;
		int out_iteration1 = 0;
		boolean flag_convergence_out = false;
		boolean flag_lambda = false;
		double old_theta_p = solution1.total_cost;
		double alpha_right = 1;
		double lambda_denominator = 1;
		
		while(!(flag_convergence_out && flag_lambda && out_iteration > 300)){
			out_iteration1 = out_iteration;
			double inner_gap = 0;

			//the following ten lines which optimize the setting of inner tolerance are specified for large scale 
			//instances including chicago, chicago_200, chicago_500. And other smaller instances don't need these
			//setting, therefore can be ignored.
			//-------------------------------setting for large scale instances only: start----------------------------------------
			if(out_iteration < 2000){
				this.inner_tolerance = 0.5;
			}else if(out_iteration < 4300){
				this.inner_tolerance = 0.05;
			}else{
				this.inner_tolerance = 0.0005;
			}
			if(out_iteration == 4300){
				lambda_denominator = lambda_denominator * 10;
			}
			//-------------------------------setting for large scale instances only: end----------------------------------------
			
			
			inner_gap = inner_loop_include_line_search1_with_bounded_right(solution1, solution2, solution_temp, 10, out_iteration, alpha_right);
			construct_ergodic_sequence(this.solution3, solution1, out_iteration + 1);
			double z =calculate_z(this.solution3);
			
			//print the convergence process (theta_p)
			solution1.get_total_cost();
			double theta_p_absolute_value = Math.abs(old_theta_p - solution1.total_cost);
			double theta_p = Math.abs(theta_p_absolute_value/solution1.total_cost);
			old_theta_p = solution1.total_cost;
			theta_p_array[theta_p_array_index][0] = out_iteration;
			theta_p_array[theta_p_array_index][1] = theta_p_absolute_value;
			theta_p_array[theta_p_array_index][2] = theta_p;
			theta_p_array_index = theta_p_array_index + 1;
			if(theta_p_array_index == this.theta_p_array.length){
				this.io.print_theta_p(theta_p_array, theta_p_array_index);
				theta_p_array_index = 0;
			}
			
			//print the convergence process (epsilon2)
			double absolute_value = Math.abs(z - solution1.total_cost);
			double epsilon2 = absolute_value/z;

			epsilon2_array[epsilon2_array_index][0] = out_iteration;
			epsilon2_array[epsilon2_array_index][1] = absolute_value;
			epsilon2_array[epsilon2_array_index][2] = epsilon2;
			epsilon2_array_index = epsilon2_array_index + 1;
			if(epsilon2_array_index == this.epsilon2_array.length){
				this.io.print_epsilon2(epsilon2_array, epsilon2_array_index);
				epsilon2_array_index = 0;
			}
			
			if(epsilon2 < this.out_tolerance){
				flag_convergence_out = true;
			}else{
				flag_convergence_out = false;
			}
			update_lambda_with_adjustable_denominator(solution1, out_iteration1 + 1, lambda_denominator);
			double epsilon1 = check_solution_epsilon1(this.solution3);
			
			//print the convergence process (epsilon1)
			epsilon1_array[epsilon1_array_index][0] = out_iteration;
			epsilon1_array[epsilon1_array_index][1] = epsilon1;
			epsilon1_array_index = epsilon1_array_index + 1;
			if(epsilon1_array_index == this.epsilon1_array.length){
				this.io.print_epsilon1(epsilon1_array, epsilon1_array_index);
				epsilon1_array_index = 0;
			}
			if(epsilon1 < this.demand_gap_tolerance && inner_gap < 0.0005){
				flag_lambda = true;
			}else{
				flag_lambda = false;
			}
			
			double converge_value = 0;
			copy_solution1_to_solution2_link_flow_demand(solution1, solution2);
			copy_solution1_to_solution2_price_total_cost(solution1, solution2);
			copy_solution1_to_solution2_link_flow_demand(solution1, solution_temp);
			copy_solution1_to_solution2_price_total_cost(solution1, solution_temp);
			out_iteration = out_iteration + 1;	
			System.out.println(out_iteration + "   " + epsilon2 + "  " + converge_value  + "  " + inner_gap + "   " + old_theta_p + "  " + lambda_denominator + "   " + alpha_right + "   " + flag_lambda);
			if(out_iteration > 100000 - 1){
				break;
			}
		}
		this.io.print_epsilon2(epsilon2_array, epsilon2_array_index);
		this.io.print_epsilon1(epsilon1_array, epsilon1_array_index);
		this.io.print_theta_p(theta_p_array, theta_p_array_index);
	}
	
	public double check_solution_epsilon1(Solution solution){
		double enumerator = 0;
		double denominator = 0;
		for(Rider rider : solution.riders){
			double gap = rider.demand  - rider.satisfied_demand;
			if(gap < 0){
				gap = 0;
			}
			enumerator = enumerator + Math.abs(gap);
			denominator = denominator + rider.demand;
		}
		return enumerator/denominator;
	}
	
	public void update_lambda_with_adjustable_denominator(Solution solution, double out_iteration, double denominator){
		for(Rider rider : solution.riders){
			double current_price = rider.price;
			double gap = rider.demand - rider.satisfied_demand;
			gap = gap / denominator;
			double updated_step = this.out_alpha_hat * gap /out_iteration;
			if(out_iteration < 50){
				if(updated_step > 1){
					updated_step = 1;
				}else if(updated_step < -1){
					updated_step = -1;
				}
			}
			if(current_price + updated_step < 0){
				rider.price = 0;
			}else{
				rider.price = current_price + updated_step;
			}
		}
	}
	
	public void construct_ergodic_sequence(Solution solution3, Solution solution1, double out_iteration){
		double total2 = this.total1 + Math.pow(out_iteration, this.out_k);
		double total3 = Math.pow(out_iteration, this.out_k);
		
		double factor1 = this.total1/total2;
		double factor2 = total3/total2;
		this.total1 = total2;
		
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			Rider rideri_solution3 = solution3.riders.get(rideri);
			Rider rideri_solution1 = solution1.riders.get(rideri);
			rideri_solution3.satisfied_demand = factor1 * rideri_solution3.satisfied_demand + factor2 * rideri_solution1.satisfied_demand;
		}
		for(int driveri = 0; driveri < this.driver_number; driveri ++){
			Driver driveri_solution3 = solution3.drivers.get(driveri);
			Driver driveri_solution1 = solution1.drivers.get(driveri);
			for(int j = 0; j < this.rider_number + 1; j ++){
				driveri_solution3.flow_distribution[j] = factor1 * driveri_solution3.flow_distribution[j] + factor2 * driveri_solution1.flow_distribution[j];
			}
		}
		for(int linki = 0; linki < this.link_number; linki ++){
			Link linki_solution3 = solution3.graph.links.get(linki);
			Link linki_solution1 = solution1.graph.links.get(linki);
			linki_solution3.flow = factor1 * linki_solution3.flow + factor2 * linki_solution1.flow;
		}
	}
	
	public double calculate_z(Solution solution3){
		double total_cost = 0;
		List<Link> links = solution3.graph.links;
		List<Driver> drivers = solution3.drivers;
		for(Link link : links){
			total_cost = total_cost + 4 * (link.free_flow_time * link.flow + 0.03 * link.free_flow_time * Math.pow(link.flow, 5)/Math.pow(link.capacity, 4));
		}
		for(Driver driveri : drivers){
			total_cost = total_cost + (2 * this.t0 + this.delta) * (driveri.demand - driveri.flow_distribution[0]);
		}
		return total_cost;
	}
	
	public boolean find_minimal_along_given_direction_for_inner_loop_with_bounded_right(Solution original_solution, Solution improved_direction, Solution solution_temp, double alpha_right){
		boolean flag_optimized = false;
		double[] link_flow = new double[this.link_number];
		double[] riders_satisfied_demand = new double[this.rider_number];
		copy_improved_direction_to_array(improved_direction, link_flow, riders_satisfied_demand);
		
		double alpha_left = 0;
		double alpha_gap = alpha_right - alpha_left;
		double alpha_left_rate = (3 - Math.sqrt(5))/2;
		double alpha_right_rate = (Math.sqrt(5) - 1)/2;
		double alpha_left_point = alpha_left + alpha_left_rate * alpha_gap;
		double alpha_right_point = alpha_left + alpha_right_rate * alpha_gap;
		update_solution_from_improved_array_for_factor(original_solution, improved_direction, alpha_left_point, link_flow, riders_satisfied_demand);
		update_solution_from_improved_array_for_factor(original_solution, solution_temp, alpha_right_point, link_flow, riders_satisfied_demand);
		improved_direction.get_total_cost();
		solution_temp.get_total_cost();
		double alpha_middle_left_cost = improved_direction.total_cost;
		double alpha_middle_right_cost = solution_temp.total_cost;
		while(alpha_gap > this.alpha_tolerance){
						
			if(alpha_middle_left_cost < alpha_middle_right_cost){
				alpha_right = alpha_right_point;
				alpha_middle_right_cost = alpha_middle_left_cost;
				alpha_right_point = alpha_left_point;
				alpha_gap = alpha_right - alpha_left;
				alpha_left_point = alpha_left + alpha_left_rate * alpha_gap;
				update_solution_from_improved_array_for_factor(original_solution, improved_direction, alpha_left_point, link_flow, riders_satisfied_demand);
				improved_direction.get_total_cost();
				alpha_middle_left_cost = improved_direction.total_cost;
			}else{
				flag_optimized = true;
				alpha_left = alpha_left_point;
				alpha_middle_left_cost = alpha_middle_right_cost;
				alpha_left_point = alpha_right_point;
				alpha_gap = alpha_right - alpha_left;
				alpha_right_point = alpha_left + alpha_right_rate * alpha_gap;
				update_solution_from_improved_array_for_factor(original_solution, solution_temp, alpha_right_point, link_flow, riders_satisfied_demand);
				solution_temp.get_total_cost();
				alpha_middle_right_cost = solution_temp.total_cost;
			}
		}
		return flag_optimized;
	}
	
	public double inner_loop_include_line_search1_with_bounded_right(Solution original_solution, Solution improved_direction, Solution solution_temp_line, double k, int out_iteration, double alpha_right){
		double converge_value = 0;
		original_solution.calculate_link_cost();
		improved_direction.calculate_link_cost();
		original_solution.get_total_cost();
		
		double inner_line_search_iteration = 0;
		boolean flag_convergence = false;
		
		
		while(!flag_convergence){
			if(out_iteration < 3000){
				if(inner_line_search_iteration > this.inner_maximum_line_search_iteration_number/2){
					break;
				}
			}else if(out_iteration < 5000){
				if(inner_line_search_iteration > this.inner_maximum_line_search_iteration_number * 1.5){
					break;
				}
			}else{
				if(inner_line_search_iteration > this.inner_maximum_line_search_iteration_number * 7){
					break;
				}
			}
			inner_line_search_iteration = inner_line_search_iteration + 1;
			get_improved_solution_direction(improved_direction);
			double total_new_cost_for_improved_direction = improved_direction.get_total_cost_for_improved_direction_only();
			//copy_improved_direction_to_array(improved_direction, link_flow, riders_satisfied_demand);
			
			double total_new_cost = original_solution.get_total_cost_new();
			converge_value = (total_new_cost - total_new_cost_for_improved_direction)/Math.abs(total_new_cost);
			if(converge_value < this.inner_tolerance){
				flag_convergence = true;
			}
			//System.out.println(converge_value);
			
			boolean flag_optimized = find_minimal_along_given_direction_for_inner_loop_with_bounded_right(original_solution, improved_direction, solution_temp_line, alpha_right);
			if(flag_optimized == true){
				Solution solution_temp = improved_direction;
				improved_direction = original_solution;
				original_solution = solution_temp;
				original_solution.calculate_link_cost();
				copy_solution1_to_solution2_links(original_solution, improved_direction);
			}else{
				this.alpha_tolerance = this.alpha_tolerance * 0.1;
			}
		}
		return converge_value;
	}
	
	public void update_solution(Solution solution1, Solution improved_direction, double iteration){
		double factor = 1/iteration;
		
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			improved_direction.riders.get(rideri).satisfied_demand = (1 - factor) * solution1.riders.get(rideri).satisfied_demand + factor * improved_direction.riders.get(rideri).satisfied_demand; 
		}
		for(int driveri = 0; driveri < this.driver_number; driveri ++){
			Driver driver = improved_direction.drivers.get(driveri);
			int driver_shared_rider = driver.shared_rider;
			driver.flow_distribution[driver_shared_rider] = (1 - factor) * solution1.drivers.get(driveri).flow_distribution[driver_shared_rider] + factor * driver.demand;
		}
		for(int linki = 0; linki < this.link_number; linki ++){
			improved_direction.graph.links.get(linki).flow = (1 - factor) * solution1.graph.links.get(linki).flow + factor * improved_direction.graph.links.get(linki).flow;
		}
	}
	
	public void update_solution_from_improved_array(Solution solution1, Solution improved_direction, double iteration, double[] links_flow, double[] riders_satisfied_demand){
		double factor = 1/iteration;
		
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			improved_direction.riders.get(rideri).satisfied_demand = (1 - factor) * solution1.riders.get(rideri).satisfied_demand + factor * riders_satisfied_demand[rideri]; 
		}
		for(int driveri = 0; driveri < this.driver_number; driveri ++){
			Driver driver = improved_direction.drivers.get(driveri);
			for(int riderj = 0; riderj < this.rider_number + 1; riderj ++){
				driver.flow_distribution[riderj] = (1 - factor) * solution1.drivers.get(driveri).flow_distribution[riderj];
			}
			int driver_shared_rider = driver.shared_rider;
			driver.flow_distribution[driver_shared_rider] = driver.flow_distribution[driver_shared_rider] + factor * driver.demand;
		}
		for(int linki = 0; linki < this.link_number; linki ++){
			improved_direction.graph.links.get(linki).flow = (1 - factor) * solution1.graph.links.get(linki).flow + factor * links_flow[linki];
		}
	}
	
	public void update_solution_from_improved_array_for_factor(Solution solution1, Solution improved_direction, double factor1, double[] links_flow, double[] riders_satisfied_demand){
		double factor = factor1;
		
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			improved_direction.riders.get(rideri).satisfied_demand = (1 - factor) * solution1.riders.get(rideri).satisfied_demand + factor * riders_satisfied_demand[rideri]; 
		}
		for(int driveri = 0; driveri < this.driver_number; driveri ++){
			Driver driver = improved_direction.drivers.get(driveri);
			for(int riderj = 0; riderj < this.rider_number + 1; riderj ++){
				driver.flow_distribution[riderj] = (1 - factor) * solution1.drivers.get(driveri).flow_distribution[riderj];
			}
			int driver_shared_rider = driver.shared_rider;
			driver.flow_distribution[driver_shared_rider] = driver.flow_distribution[driver_shared_rider] + factor * driver.demand;
		}
		for(int linki = 0; linki < this.link_number; linki ++){
			improved_direction.graph.links.get(linki).flow = (1 - factor) * solution1.graph.links.get(linki).flow + factor * links_flow[linki];
		}
	}
		
	public void copy_solution1_to_solution2_links(Solution original_solution, Solution derived_solution){
		List<Link> derived_links = derived_solution.graph.links;
		List<Link> original_links = original_solution.graph.links;
		for(int linki = 0; linki < this.link_number; linki++){
			derived_links.get(linki).current_cost = original_links.get(linki).current_cost;
		}
	}
	
	public void copy_solution1_to_solution2_price_total_cost(Solution original_solution, Solution derived_solution){
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			derived_solution.riders.get(rideri).price = original_solution.riders.get(rideri).price;
		}
		derived_solution.total_cost = original_solution.total_cost;
	}
	
	public void copy_solution1_to_solution2_link_flow_demand(Solution original_solution, Solution derived_solution){
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			derived_solution.riders.get(rideri).satisfied_demand = original_solution.riders.get(rideri).satisfied_demand;
		}
		for(int linki = 0; linki < this.link_number; linki ++){
			derived_solution.graph.links.get(linki).flow = original_solution.graph.links.get(linki).flow;
		}
		for(int driveri = 0; driveri < this.driver_number; driveri ++){
			for(int j = 0; j < this.rider_number + 1; j ++){
				derived_solution.drivers.get(driveri).flow_distribution[j] = original_solution.drivers.get(driveri).flow_distribution[j];
			}
		}
	}
	
	public void copy_improved_direction_to_array(Solution solution, double[] links_flow, double[] riders_satisfied_demand){
		for(int linki = 0; linki < this.link_number; linki++){
			links_flow[linki] = solution.graph.links.get(linki).flow;
		}
		for(int rideri = 0; rideri < this.rider_number; rideri ++){
			riders_satisfied_demand[rideri] = solution.riders.get(rideri).satisfied_demand;
		}
	}
	
	public void get_improved_solution_direction(Solution solution){	
		//simple_network_check(solution);
		Path[] drivers_path = solution.get_all_shortest_path();
		solution.initialize_link_rider_flow();
		solution.get_improve_direction(drivers_path);
	}
	
	public void result_handle(IO io){
		io.print_driver_rider_flow_distribution(this.solution3);
		io.print_riders_price(this.solution1);
		io.print_run_duration(this.run_duration);
		io.print_net_flow(this.solution3);
	}
	
	//this method is used for initializing the initial solution
	public void initialize_initial_solution(Solution solution){
		for(Driver driver : solution.drivers){
			driver.flow_distribution[driver.shared_rider] = driver.demand;
		}
		solution.calculate_link_cost();
		solution.get_total_cost();
	}
	
	public void simple_network_check(Solution solution){
		double linka_cost = solution.graph.links.get(0).current_cost;
		double linkb_cost = solution.graph.links.get(1).current_cost;
		double linkc_cost = solution.graph.links.get(2).current_cost;
		
		double rider1_price = solution.riders.get(0).price;
		double rider2_price = solution.riders.get(1).price;
		double rider3_price = solution.riders.get(2).price;
		
		double tra1 = linka_cost + 9 - rider1_price;
		double tra2 = linka_cost + linkb_cost + linkc_cost + 9 - rider2_price;
		double tra3 = linka_cost + linkb_cost + linkc_cost + 9 - rider3_price;
		double tra4 = linka_cost;
		
		double tra5 = linka_cost + linkc_cost + 9 - rider1_price;
		double tra6 = linka_cost + linkc_cost + 9 - rider3_price;
		double tra7 = linka_cost + linkc_cost + 9 - rider2_price;
		double tra8 = linka_cost + linkc_cost;
		
		double tra9 = linkc_cost + 9 - rider3_price;
		double tra10 = linkc_cost;
		
		System.out.println(tra1 + "\t" + tra2 + "\t" + tra3 + "\t" + tra4);
		System.out.println(tra5 + "\t" + tra6 + "\t" + tra7 + "\t" + tra8);
		System.out.println(tra9 + "\t" + tra10);
	}
}
