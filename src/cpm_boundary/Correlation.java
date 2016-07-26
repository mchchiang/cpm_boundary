package cpm_boundary;

import java.io.*;

public class Correlation {
	public static void main (String [] args) throws IOException {
		int x = Integer.parseInt(args[0]);
		int y = Integer.parseInt(args[1]);
		int startRow = Integer.parseInt(args[2]);
		int maxTau = Integer.parseInt(args[3]);
		String filename = args[4];
		String outputFile = args[5];
		
		if (maxTau > y-startRow) maxTau = y-startRow;

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

		double sum1, sum2;
		int tmax;
		double [][] result = new double [2][maxTau];
		double avgSq1 = 0.0;
		double avgSq2 = 0.0;
		for (int i = 0; i < y-startRow; i++){
			avgSq1 += data[i][1];
			avgSq2 += data[i][2];
		}
		avgSq1 /= (double) (y-startRow);
		avgSq2 /= (double) (y-startRow);
		avgSq1 *= avgSq1;
		avgSq2 *= avgSq2;
		
		for (int tau = 0; tau < maxTau; tau++){
			if (tau % 10000 == 0){
				System.out.println(tau);
			}
			sum1 = 0.0;
			sum2 = 0.0;
			tmax = y - tau - startRow;
			for (int t = 0; t < tmax; t++){
				sum1 += data[t][1] * data[tau+t][1];
				sum2 += data[t][2] * data[tau+t][2];
			}
			sum1 /= (double) tmax;
			sum2 /= (double) tmax;
			result[0][tau] = sum1 - avgSq1;
			result[1][tau] = sum2 - avgSq2;
		}

		double norm1 = result[0][0];
		double norm2 = result[1][0];
		for (int tau = 0; tau < result.length; tau++){
			writer.printf("%d %.5f %.5f %.5f %.5f\n", tau, 
					result[0][tau], result[0][tau] / norm1,
					result[1][tau], result[1][tau] / norm2);	
		}

		writer.close();
	}
}
