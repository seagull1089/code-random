package Astar;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.PriorityQueue;

public class AlgoAStar {
	public static Result getMinPath(String s1,String s2){
		
		Result result = new Result();
		ArrayList<String> path = null;
		HashSet<String> closedNodes = new HashSet<String>();
		double start_time = System.currentTimeMillis();
		double pathCost = 0;
		int maxNodesExpanded =0;
		int maxNodesinPQ = 0;
		if(MapGraph.getNode(s1)== null){
			path = new ArrayList<String>();
		}else if(s1.equals(s2)){
			path = new ArrayList<String>();
			path.add(s1);
			pathCost =0;
		}else{
			PriorityQueue<ComparableNode> pq = new PriorityQueue<ComparableNode>();
			Node tmp_node1 = MapGraph.getNode(s1);
			Iterator it = tmp_node1.adjacentNodes.keySet().iterator();
			while(it.hasNext()){
				String adj_node_name =(String) it.next();
				ArrayList<String> curr_path = new ArrayList<String>();
				curr_path.add(tmp_node1.name);
				
				double g_cost = tmp_node1.adjacentNodes.get(adj_node_name).doubleValue();
			//	System.out.println(" Initial cost = " + " " + s1 + " " + adj_node_name + " " + g_cost);
				double h_cost =  Astar.Distance.getDistance(adj_node_name, s2);
				
				Node adjacentNode = MapGraph.getNode(adj_node_name);
				
				if(adjacentNode == null){
					continue;
				}
				
				ComparableNode cmp_node = new ComparableNode(adjacentNode,curr_path,g_cost,h_cost);
			//	System.out.println(" Adding  " + adjacentNode.name + " Cost = "+ cmp_node.cost);
				pq.add(cmp_node);
				//maxNodesExpanded++;
				maxNodesinPQ = (maxNodesinPQ >pq.size()) ? maxNodesinPQ : pq.size();
			}
			closedNodes.add(s1);
			while(!pq.isEmpty()){
				ComparableNode cmp_node_tmp = pq.remove();
				Node n1 = cmp_node_tmp.node;
				//System.out.println("Node = " + n1.name + " Cost ==  "+ cmp_node_tmp.g_cost);
				
				if(n1.name.equals(s2)){
					path = cmp_node_tmp.path;
					pathCost = cmp_node_tmp.g_cost;
			//		System.out.println(" Found ..");
					break;
				}
				
				if(closedNodes.contains(n1.name)){
					continue;
				}
				
			
				
				//maxNodesExpanded++;
				Iterator it1 = n1.adjacentNodes.keySet().iterator();
				while(it1.hasNext()){
					
					    String adj_node = (String)it1.next();
					  
					    
					    if(closedNodes.contains(adj_node)){
							continue;
						}
					
					    
						Node tmp_node = MapGraph.getNode(adj_node);
						ArrayList<String> curr_path = new ArrayList<String>(cmp_node_tmp.path);
						
						/*
						if(tmp_node == null){
							tmp_node = new Node(adj_node);
						}*/
						
						double g_cost  = cmp_node_tmp.g_cost + n1.adjacentNodes.get(adj_node).doubleValue();
						double h_cost = Astar.Distance.getDistance(adj_node, s2);
						for( String s:curr_path){
							//System.out.println(" Inside : " + s);
							
						}
						
						pq.add(new ComparableNode(tmp_node,curr_path,g_cost,h_cost));
						
						maxNodesinPQ = (maxNodesinPQ >pq.size()) ? maxNodesinPQ : pq.size();
				}
						
				closedNodes.add(n1.name);
			}  // end of while loop
			
		}
		
		double end_time = System.currentTimeMillis();
		result.path = (path == null) ? new ArrayList<String>() : path ;
		result.pathCost = pathCost;
		result.maxnodesExpanded = closedNodes.size();
		result.maxnodesInPQ = maxNodesinPQ;
		result.timeTaken = (end_time - start_time)/1000;
		return result;
	}
	
	public static Result getMinPath(int x1,int y1,int x2,int y2){
		String s1 = ""+x1+":"+y1;
		String s2 = ""+x2+":"+y2;
		return getMinPath(s1,s2);
	}

}
