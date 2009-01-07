package Astar;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Set;


public class MapGraph {
	
	static HashMap<String,Node> mapGraph;
	static boolean initialized = false;

	public static void LoadData(String filename){
		if(!initialized){
			try {
				mapGraph = new HashMap<String, Node>();
				
				String str_line;
				BufferedReader in = new BufferedReader(new FileReader(filename));
				while((str_line = in.readLine())!=null){
					String int_str[] = str_line.split(" ");
					int oneway = new Integer(int_str[0]).intValue();
					int x1 = new Integer(int_str[1]).intValue();
					int y1 = new Integer(int_str[2]).intValue();
					int x2 = new Integer(int_str[3]).intValue();
					int y2 = new Integer(int_str[4]).intValue();
					
					//System.out.println(oneway + ":" + x1 + ":" + y1 + ":" + x2 + ":" + y2);
					String point1 = ""+x1+":"+y1;
					String point2 = ""+x2+":"+y2;
					//System.out.println(point1 + " || " + point2);
					if(oneway ==1 || oneway == 2){
						if(!mapGraph.containsKey(point1)){
							mapGraph.put(point1, new Node(x1,y1));
						}
						if(!mapGraph.containsKey(point2)){
							mapGraph.put(point2, new Node(x2,y2));
						}
						Node tmpNode = mapGraph.get(point1);
						tmpNode.addAdjacentNode(x2,y2);
					}
					if(oneway == 2){
						/*if(!mapGraph.containsKey(point2)){
							mapGraph.put(point2, new Node(x2,y2));
						}*/
						Node tmpNode = mapGraph.get(point2);
						tmpNode.addAdjacentNode(x1,y1);
					}
					
				}
				in.close();
				initialized = true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}
	}

	public static Set<String> getAllNodeNames(){
		return mapGraph.keySet();
	}
	
	public static Node getNode(String name){
		return mapGraph.get(name);
	}
	
}
