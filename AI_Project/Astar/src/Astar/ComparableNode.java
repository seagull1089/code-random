package Astar;

import java.util.ArrayList;


public class ComparableNode implements Comparable<ComparableNode>{

	public Node node;
	public double g_cost;
	public double h_cost;
	public double cost;
	public ArrayList<String> path; 
	public ComparableNode(Node n,ArrayList<String> n_path,double g_cost,double h_cost){
		this.node = n;
		this.path  = n_path;
		this.g_cost = g_cost;
		this.h_cost = h_cost;
		this.cost = g_cost+h_cost;
		this.path.add(n.name);
		
	}
	
	public int compareTo(ComparableNode c){
		if(this.cost - c.cost < 0){
			return -1;
		}else if(this.cost-c.cost>0){
			return 1;
		}
		return 0;
	}
	

}
