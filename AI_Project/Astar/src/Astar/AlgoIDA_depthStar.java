package Astar;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.PriorityQueue;

public class AlgoIDA_depthStar {
	public static int maxDepthStep =3;

	public void setmaxDepthStep(int max){
		maxDepthStep = max;
	}

	public static Result getMinPath(String s1,String s2){
		
		Result rs = new Result();
		ArrayList<String> path;
		int depthTraversed = 0;
		int currmaxdepth = 3;
		if(s1.equals(s2)){
		
			path = new ArrayList<String>();
			path.add(s1);
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
			//	maxNodesinPQ = (maxNodesinPQ >pq.size()) ? maxNodesinPQ : pq.size();
			}

			while(!searchPQ(new PriorityQueue<ComparableNode>(pq), s2, depthTraversed, currmaxdepth, rs)){
				if(depthTraversed > 1300){
					break;
				}
				depthTraversed = currmaxdepth;
				currmaxdepth += maxDepthStep;
				
			}
			
			
		}
		
		return rs;
	}
	
	private static boolean searchPQ(PriorityQueue<ComparableNode> pq,String destNode,int mindepth,int maxdepth,Result rs){
		boolean found = false;
		System.out.println(" destNode " + destNode + " mindepth " + mindepth + " maxdepth " + maxdepth );
		while(!pq.isEmpty()){
			ComparableNode cmpnode = pq.remove();
			Node node = cmpnode.node;
			int depth = cmpnode.path.size();
			System.out.println(" Node opened = "+ node.name);
			if(depth > maxdepth){
				continue;
			}
			
			if(node.name.equals(destNode)){
					System.out.println(" Found " + destNode);
					rs.path = cmpnode.path;
					rs.pathCost  = cmpnode.g_cost;
					found = true;
					break;
			}
					
			Iterator it1 = node.adjacentNodes.keySet().iterator();
			while(it1.hasNext()){
				String adj_node = (String)it1.next();
				Node tmp_node = MapGraph.getNode(adj_node);
				ArrayList<String> curr_path = new ArrayList<String>(cmpnode.path);
				if(tmp_node == null){
					continue;
				}
				
				double g_cost=0,h_cost=0;
				if(depth > mindepth){
					g_cost  = cmpnode.g_cost + node.adjacentNodes.get(adj_node).doubleValue();
					h_cost = Astar.Distance.getDistance(adj_node, destNode);
				}
				
				for( String s:curr_path){
					//System.out.println(" Inside : " + s);
					
				}
				
				pq.add(new ComparableNode(tmp_node,curr_path,g_cost,h_cost));
			//	maxNodesExpanded++;
				//maxNodesinPQ = (maxNodesinPQ >pq.size()) ? maxNodesinPQ : pq.size();
			}
			
		}
		
		return found;
	}
	
	public static Result getMinPath(int x1,int y1,int x2,int y2){
		String s1 = ""+x1+":"+y1;
		String s2 = ""+x2+":"+y2;
		return getMinPath(s1,s2);
	}

	
}
