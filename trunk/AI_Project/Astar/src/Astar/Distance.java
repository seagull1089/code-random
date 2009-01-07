package Astar;

public class Distance {

	public static double getDistance(int x1 , int y1 , int x2,int y2){
		double d =  Math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
		//System.out.println(" " + x1 + " " + y1 + " " + x2 + " " + y2 + " dis:" + d);
		return d;
		
	}
	
	public static double getDistance(String s1,String s2){
		String point1[] = s1.split(":");
		int x1 = new Integer(point1[0]).intValue();
		int y1 = new Integer(point1[1]).intValue();
		String point2[] = s2.split(":");
		int x2 = new Integer(point2[0]).intValue();
		int y2 = new Integer(point2[1]).intValue();
		return getDistance(x1,y1,x2,y2);
	}
	
}
