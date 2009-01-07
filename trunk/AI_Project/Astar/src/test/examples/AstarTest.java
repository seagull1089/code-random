package test.examples;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import Astar.*;

public class AstarTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		MapGraph.LoadData("/Users/ragip002/AI_Project/MapData.txt");
		/*Node n1 = MapGraph.getNode("1699:7706");
		System.out.println("Name "+ n1.name);
		Iterator it =  n1.adjacentNodes.keySet().iterator();
		while(it.hasNext()){
			String key = (String) it.next();
			System.out.println(key + ":" + n1.adjacentNodes.get(key));
		}*/
		//System.out.println(" Distance -- " + Astar.Distance.getDistance(1121, 7568 ,982, 7504));
/*
		for (String name : MapGraph.getAllNodeNames()){
		//	System.out.println("Name of the node = " + name);
			HashSet<String> hasSeen = new HashSet<String>();
			Node tmp1 = MapGraph.getNode(name);
			ArrayList<String> path = new ArrayList<String>();
			while(tmp1.adjacentNodes!=null){
			//	System.out.println(" searching thru " + tmp1.name);
				Iterator it = tmp1.adjacentNodes.keySet().iterator();
				boolean termNotFound = false; 
				while(it.hasNext()){
					String tmpName = (String)it.next();
					if(hasSeen.contains(tmpName)){
						continue;
					}
				//	System.out.println(" Adding " + tmpName);
					path.add(tmpName);
					hasSeen.add(tmpName);
					tmp1 = MapGraph.getNode(tmpName);
					if(tmp1 != null)
						termNotFound = true;
					break;
				}
				
				if(!termNotFound){
					break; 
				}
				
			}
			
			//System.out.println(" path Length = " + path.size() + " Init = " + path.get(0) + " Goal =" + path.get(path.size()-1));
			System.out.println(" " + path.get(0) + " ," + path.get(path.size()-1));
			for(String s: path ){
				System.out.print(" " + s);
			}
			//System.out.println();
		}
	*/	
	//	1999,7740,1740,7601
		
		Result rs = AlgoAStar.getMinPath(1999,7740,1740,7601);
	//	Result rs1 = AlgoAStar.getMinPath(1999,7740,1740,7601);
		System.out.println("maxNodes in PQ ="+ rs.maxnodesInPQ);
		System.out.println("maxNodesExapnded =" + rs.maxnodesExpanded);
		System.out.println("pathCost = "+ rs.pathCost);

		for(String s:rs.path){
			
			System.out.print(" " + s );
		//	System.out.println(" Distance = " +  Astar.Distance.getDistance(tmp_s, s));
			
		}
/*		
		System.out.println("maxNodes in PQ ="+ rs1.maxnodesInPQ);
		System.out.println("maxNodesExapnded =" + rs1.maxnodesExpanded);
		System.out.println("pathCost = "+ rs1.pathCost);
		for(String s:rs1.path){
			
			System.out.print(" " + s );
		//	System.out.println(" Distance = " +  Astar.Distance.getDistance(tmp_s, s));
			
		}
		
	*/	
		//System.out.println(" Distance -- " + Astar.Distance.getDistance("1522:7879", "1087:8996"));
	}

}
