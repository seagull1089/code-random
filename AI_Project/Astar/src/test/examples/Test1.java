package test.examples;


import java.io.BufferedReader;

import java.io.FileReader;

public class Test1 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String filename = "/Users/ragip002/AI_Project/MapData.txt";
		try{
			String map_line;
			int count = 0;
			BufferedReader in = new 
				BufferedReader(new FileReader(filename));
			while((map_line = in.readLine())!=null){
				System.out.println(map_line);
				++count;
			}
			in.close();
			System.out.println("Final Count = "+count);
		}catch(Exception e){
			e.printStackTrace();
		}

	}
}