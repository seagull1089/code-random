package Astar;


import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Set;

public class AlgoIDAStar {
	public static double stepFactor =1.5;

	public void setStepFactor(double max){
			stepFactor = max;
	}

	public static Result getMinPath(String s1,String s2){
			
			double mincost = 0,maxcost = 1.2*Astar.Distance.getDistance(s1, s2);
			HashSet<String> closedNodes = new HashSet<String>();
			
			double start_time = System.currentTimeMillis();
			Result rs = new Result();
			
			
			//System.out.println("Initial Cost Estimate " + maxcost);
			
			if(MapGraph.getNode(s1)== null){
				rs.path = new ArrayList<String>();
				
			}else if(s1.equals(s2)){
				rs.path = new ArrayList<String>();
				rs.path.add(s1);
			}else{
				Node tmp_node1 = MapGraph.getNode(s1);
				Iterator it = tmp_node1.adjacentNodes.keySet().iterator();
				PriorityQueue<ComparableNode> pq = new PriorityQueue<ComparableNode>();
				
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
				//	maxNodesExpanded++;
					rs.maxnodesInPQ = (rs.maxnodesInPQ >pq.size()) ? rs.maxnodesInPQ : pq.size();
				}
				closedNodes.add(s1);
				
				
				PriorityQueue<ComparableNode> stepSearch = new PriorityQueue<ComparableNode>(pq);
				PriorityQueue<ComparableNode> residualNodes = new PriorityQueue<ComparableNode>();
				
				while(!searchPQ(stepSearch, s2, mincost, maxcost, rs,residualNodes,closedNodes)){
					
					if(stepSearch.isEmpty() && residualNodes.isEmpty()){
						rs.path = new ArrayList<String>();
						break;
					}
					mincost = maxcost;
					maxcost = maxcost*stepFactor;
					stepSearch = new PriorityQueue<ComparableNode>(residualNodes);
					residualNodes = new PriorityQueue<ComparableNode>();
				}
				
			}	
			rs.maxnodesExpanded = closedNodes.size();
			double end_time = System.currentTimeMillis();
			rs.timeTaken = (end_time - start_time)/1000;
			return rs;
		}
		
	private static boolean searchPQ(PriorityQueue<ComparableNode> pq,String destNode,
			 double mincost,double maxcost,Result rs,
			 PriorityQueue<ComparableNode> residualNodes,Set<String> closedNodes){
		
			boolean found = false;
			
			while(!pq.isEmpty()){
				ComparableNode cmpnode = pq.remove();
				Node node = cmpnode.node;
				double cost  = cmpnode.cost;
			
				if(cost > maxcost){
					residualNodes.add(cmpnode);
					continue;
				}
				
				if(node.name.equals(destNode)){
					
						rs.path = cmpnode.path;
						rs.pathCost  = cmpnode.g_cost;
						found = true;
						break;
				}
				if(closedNodes.contains(node.name)){
					continue;
				}
					
				Iterator it1 = node.adjacentNodes.keySet().iterator();
				while(it1.hasNext()){
					String adj_node = (String)it1.next();
					if(closedNodes.contains(adj_node)){
							continue;
					}
					Node tmp_node = MapGraph.getNode(adj_node);
					ArrayList<String> curr_path = new ArrayList<String>(cmpnode.path);
					if(tmp_node == null){
						continue;
					}
					
					double g_cost=0,h_cost=0;	
					g_cost  = cmpnode.g_cost + node.adjacentNodes.get(adj_node).doubleValue();
					h_cost = Astar.Distance.getDistance(adj_node, destNode);
					
					pq.add(new ComparableNode(tmp_node,curr_path,g_cost,h_cost));
					
					rs.maxnodesInPQ = (rs.maxnodesInPQ >pq.size()) ? rs.maxnodesInPQ : pq.size();
				}
				closedNodes.add(node.name);
				
			}
			
			return found;
		}
		
	public static Result getMinPath(int x1,int y1,int x2,int y2){
			String s1 = ""+x1+":"+y1;
			String s2 = ""+x2+":"+y2;
			return getMinPath(s1,s2);
		}

		
}

	
