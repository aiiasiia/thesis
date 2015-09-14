/**
  * PMS8outputReader.java
  * Prints out l,d,run,time(s) of a multi-run PMS8 printout. One row per run.
  * Note that the reading is dependent on dataset file naming convention.
  *
  * @author Aia Sia
  * @version 1.0 9/09/2015
*/

import java.util.*;
import java.io.*;


public class PMS8outputReader {

	public static void main(String args[]) throws Exception {
		Scanner input = new Scanner (new File ( args[0] ));
		
		while(input.hasNextLine()) {
			String[] line = input.nextLine().split(" ");
			if( line[0].equals("Time") && line[1].equals("on")) {
				String[] filename = line[6].split("/");
				System.out.println(filename[1] + "," + line[4].substring(0, line[4].length()-2) );
			}
		}

		input.close();
	}
}