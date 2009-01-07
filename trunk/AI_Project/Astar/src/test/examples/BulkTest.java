package test.examples;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.SortedSet;

import Astar.AlgoAStar;
import Astar.AlgoIDAStar;
import Astar.MapGraph;
import Astar.Result;

public class BulkTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		/*
		String file = "/Users/ragip002/AI_Project/stringPath.txt";
		String graphFile = "/Users/ragip002/AI_Project/MapData.txt";
		String outFile = "out.txt";
		*/
	
		if(args.length !=2) {
			System.out.println("Invalid number of arguments, the GraphFile and the searchNodesFile should be given as commandline args");
			System.exit(1);
		}
	
		String graphFile = args[0];
		String searchNodesFile = args[1];
		String outFile = "out.txt";
	


		try{
			double start_time = System.currentTimeMillis();
			MapGraph.LoadData(graphFile);
			String in_str;
		//	HashMap<Integer, ArrayList<Result>> ResultSet = new HashMap<Integer,ArrayList<Result>>();
			BufferedReader in  = new BufferedReader(new FileReader(searchNodesFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
		        System.out.println("Calculating the paths:");
			int numPaths = 0;	
			while((in_str = in.readLine())!=null){
				String nodes[] = in_str.split(",");
				
				Result rsa = AlgoAStar.getMinPath(nodes[0], nodes[1]);
				Result rsida = AlgoIDAStar.getMinPath(nodes[0], nodes[1]);
				/*
				Integer key = new Integer(rsa.path.size());
				
				if(!ResultSet.containsKey(key) ) {
					ResultSet.put(key, new ArrayList<Result>());
				}
				
				ResultSet.get(key).add(rsa);
				ResultSet.get(key).add(rsida);
				*/
				
				out.write("StartNode = " + nodes[0] + " ,GoalNode =" + nodes[1] + ", PathLength =" + rsa.path.size()+"\n");
				out.write(" AStar Cost         = 	" + rsa.pathCost  +
					 		", IDAStar cost          = " + rsida.pathCost +"\n");
				out.write(" AStar maxNodesExp  = 	" + rsa.maxnodesExpanded  + 
							", IDAStar maxnodesExp   = " + rsida.maxnodesExpanded + "\n");
				out.write(" AStar maxNodesInPQ = 	" + rsa.maxnodesInPQ  + 	
							", IDAStar maxnodesInPQ  = " + rsida.maxnodesInPQ+"\n");
				out.write(" AStar timeTaken    = 	" + rsa.timeTaken  + 	
							", IDAStar timeTaken     = " + rsida.timeTaken+"\n");
			    
				
				out.write("AStar  path :");
				for(String s : rsa.path){
					out.write( "->"+ s);
				}
				out.write("\n");
				out.write("IDAStar path :");
				for(String s : rsida.path){
					out.write( "->"+ s);
				}
				out.write("\n");
				out.write("\n");
		
				System.out.print(".");	
				/*
				Object[] keys = ResultSet.keySet().toArray();
				for(Object k : keys){
					System.out.println(" " + (Integer)k );
				}
				System.out.println();
				*/
				numPaths++;
			}
			double end_time = System.currentTimeMillis();
			System.out.println("!!Done!!");
			System.out.println("Total number of paths Calculated 	  = " + numPaths);
			System.out.println("Total time Taken to compute all paths = " + (end_time-start_time)/1000);
			System.out.println("The output is gathered into file : 	  " + outFile);
			in.close();
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}

	}

}
