package Astar;

import java.util.HashMap;

public class Node {
	public String name;
	public int coordX,coordY;
	public HashMap<String,Float> adjacentNodes;

	public Node(int x,int y){
		this.name = ""+x+":"+y;
		this.coordX = x;
		this.coordY = y;
		adjacentNodes = new HashMap<String,Float>();
	}
	
	public Node(String Name){
		String coords[] = Name.split(":");
		this.coordX = new Integer(coords[0]).intValue();
		this.coordY = new Integer(coords[1]).intValue();
		this.name = Name;
		adjacentNodes = new HashMap<String,Float>();
	}
	
	public void addAdjacentNode(int x,int y,Float distance){
		String nodeName = ""+x+":"+y;
		adjacentNodes.put(nodeName, distance);
	}
	
	public void addAdjacentNode(String Name){
		String coords[] = Name.split(":");
		int x = new Integer(coords[0]).intValue();
		int y = new Integer(coords[1]).intValue();
		addAdjacentNode(x, y);
	}
	
	public void addAdjacentNode(int x,int y){
		double d = Astar.Distance.getDistance(coordX,coordY,x,y);
		addAdjacentNode(x, y, new Float(d));
	}
	
}
