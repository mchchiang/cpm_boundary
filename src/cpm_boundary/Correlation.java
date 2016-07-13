package cpm_boundary;

import java.io.*;

public class Correlation {
	public static void main (String [] args) throws IOException {
		int x = Integer.parseInt(args[0]);
		int y = Integer.parseInt(args[1]);
		int startRow = Integer.parseInt(args[2]);
		String filename = args[3];
		String outputFile = args[4];

		double [][] data = new double [y][x];

		BufferedReader reader = new BufferedReader(new FileReader(filename));

		int count = 0;
		String line;
		String [] array;
		while (reader.ready() && count < y){
			line = reader.readLine();
			array = line.trim().split("\\s++");
			if (array.length > 0 && !(array.length == 1 && array[0].equals(""))){
				if (array.length != x){
					System.out.println("error in line " + count + " array length: " 
							+ array.length + " x: " + x);
				}
				for (int i = 0; i < x; i++){
					data[count][i] = Double.parseDouble(array[i]);
				}
				count++;
			}
		}
		reader.close();
		

		PrintWriter writer = new PrintWriter(
				new BufferedWriter(new FileWriter(outputFile)));

		double sum;
		int tmax;
		double [] result = new double [y-startRow];
		double avgSq = 0.0;
		for (int i = 0; i < y-startRow; i++){
			avgSq += data[i][1];
		}
		avgSq /= (double) (y-startRow);
		avgSq *= avgSq;
		
		for (int tau = 0; tau < y-startRow; tau++){
			if (tau % 10000 == 0){
				System.out.println(tau);
			}
			sum = 0.0;
			tmax = y - tau - startRow;
			for (int t = 0; t < tmax; t++){
				sum += data[t][1] * data[tau+t][1];
			}
			sum /= (double) tmax;
			result[tau] = sum - avgSq;
		}

		double norm = result[0];
		for (int tau = 0; tau < result.length; tau++){
			writer.printf("%d %.5f %.5f\n", tau, result[tau], result[tau] / norm);	
		}

		writer.close();
	}
}
