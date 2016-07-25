package cpm_boundary;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class PerimeterToAreaRatio {
	public static void main (String [] args) throws IOException {
		String filename = args[0];
		String outputFile = args[1];

		BufferedReader reader = new BufferedReader(new FileReader(filename));

		int count = 0;
		String [] array;
		String line = reader.readLine();
		String [] info = line.trim().split("\\s++");

		int nx = Integer.parseInt(info[0]);//lattice width
		int ny = Integer.parseInt(info[1]);//lattice height
		int q = Integer.parseInt(info[2]);//number of cells (exclude zero)

		int [][] spin = new int [nx][ny];

		while (reader.ready() && count < ny){
			line = reader.readLine();
			array = line.trim().split("\\s++");
			if (array.length > 0 && !(array.length == 1 && array[0].equals(""))){
				if (array.length != nx){
					System.out.println("error in line " + count + " array length: " 
							+ array.length + " x: " + nx);
				}
				for (int i = 0; i < nx; i++){
					spin[count][i] = Integer.parseInt(array[i]);
				}
				count++;
			}
		}
		reader.close();
		
		ArrayList<Integer> area = new ArrayList<Integer>(q+1);
		ArrayList<Integer> perimeter = new ArrayList<Integer>(q+1);
		for (int i = 0; i <= q; i++){
			area.add(0);
			perimeter.add(0);
		}
		
		int s;
		int p;
		
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				s = spin[i][j];
				area.set(s, area.get(s)+1);
				p = 0;
				if (spin[i][j] != spin[iup(i,nx)][j]) p++;
				if (spin[i][j] != spin[idown(i,nx)][j]) p++;
				if (spin[i][j] != spin[i][jup(j,ny)]) p++;
				if (spin[i][j] != spin[i][jdown(j,ny)]) p++;
				perimeter.set(s, perimeter.get(s)+p);
			}
		}
		
		double avg = 0.0;
		double avgSq = 0.0;
		
		PrintWriter writer = new PrintWriter(
				new BufferedWriter(new FileWriter(outputFile)));
		int a;
		double ratio;
		for (int i = 1; i <= q; i++){
			a = area.get(i);
			p = perimeter.get(i);
			ratio = (double) p / Math.sqrt(a);
			writer.printf("%d %d %d %.5f\n", i, p, a, ratio);
			avg += ratio;
			avgSq += ratio * ratio;
		}
		avg /= (double) q;
		avgSq /= (double) q;
		writer.printf("%.5f %.5f\n", avg, Math.sqrt(avgSq - avg * avg));
		writer.close();
	}

	//periodic boundary methods
	private static int iup(int i, int nx){
		if (i == nx-1) return 0;
		return i+1;
	}

	private static int idown(int i, int nx){
		if (i == 0) return nx-1;
		return i-1;
	}

	private static int jup(int j, int ny){
		if (j == ny-1) return 0;
		return j+1;
	}

	private static int jdown(int j, int ny){
		if (j == 0) return ny-1;
		return j-1;
	}
}
