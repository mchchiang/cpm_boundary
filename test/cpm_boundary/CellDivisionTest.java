package cpm_boundary;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

public class CellDivisionTest {
	private final int nx = 4;
	private final int ny = 4;
	private final int q = 4;
	private final double alpha = 2.0;
	private final double beta = 16.0;
	private final double tol = 0.000001;
	private final double temperature = 1;
	private final double lambda = 1;
	private final double motility = 1.0;
	private final int seed = -1;
	private final int numOfSweeps = 0;
	private final int nequil = 0;
	private final double rotateDiff = 1.0;

	//convert primitive double array to arraylist
	public ArrayList<Double> toList(double [] array){
		ArrayList<Double> list = new ArrayList<Double>();
		for (double a : array){
			list.add(a);
		}
		return list;
	}

	@Test
	public void testGyrationTensor1XX(){
		int [][] spin = new int [][]{
				{0,0,0,0,0,0,0,0},
				{0,0,1,1,1,0,0,0},
				{0,1,1,1,1,1,0,0},
				{0,1,1,1,1,1,1,0},
				{0,0,1,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		double [] areaTarget = new double [] {1.0,1.0};
		CellPottsModel model = new CellPottsModel(8, 8, 1, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		model.initSpin(spin);
		model.updateCM();
		
		CellDivision division = new CellDivision();
		assertEquals("return wrong gyration tensor comp xx value",
				2.1376, division.gyrationTensor(model, 1, 0, 0), tol);
	}

	@Test
	public void testGyrationTensor1YY(){
		int [][] spin = new int [][]{
				{0,0,0,0,0,0,0,0},
				{0,0,1,1,1,0,0,0},
				{0,1,1,1,1,1,0,0},
				{0,1,1,1,1,1,1,0},
				{0,0,1,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		double [] areaTarget = new double [] {1.0,1.0};
		CellPottsModel model = new CellPottsModel(8, 8, 1, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		model.initSpin(spin);
		model.updateCM();
		
		CellDivision division = new CellDivision();
		assertEquals("return wrong gyration tensor comp yy value",
				1.9264, division.gyrationTensor(model, 1, 1, 1), tol);
	}

	@Test
	public void testGyrationTensor1XY(){
		int [][] spin = new int [][]{
				{0,0,0,0,0,0,0,0},
				{0,0,1,1,1,0,0,0},
				{0,1,1,1,1,1,0,0},
				{0,1,1,1,1,1,1,0},
				{0,0,1,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		double [] areaTarget = new double [] {1.0,1.0};
		CellPottsModel model = new CellPottsModel(8, 8, 1, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		model.initSpin(spin);
		model.updateCM();
		
		CellDivision division = new CellDivision();
		assertEquals("return wrong gyration tensor comp xy value",
				0.3392, division.gyrationTensor(model, 1, 0, 1), tol);
	}

	@Test
	public void testGyrationTensor1YX(){
		int [][] spin = new int [][]{
				{0,0,0,0,0,0,0,0},
				{0,0,1,1,1,0,0,0},
				{0,1,1,1,1,1,0,0},
				{0,1,1,1,1,1,1,0},
				{0,0,1,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		double [] areaTarget = new double [] {1.0,1.0};
		CellPottsModel model = new CellPottsModel(8, 8, 1, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);

		model.initSpin(spin);
		model.updateCM();
		
		CellDivision division = new CellDivision();
		assertEquals("return wrong gyration tensor comp xy value",
				0.3392, division.gyrationTensor(model, 1, 1, 0), tol);
	}

	@Test
	public void testGetMajorAxis1(){
		int [][] spin = new int [][]{
				{0,0,0,0,0,0,0,0},
				{0,0,1,1,1,0,0,0},
				{0,1,1,1,1,1,0,0},
				{0,1,1,1,1,1,1,0},
				{0,0,1,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		double [] areaTarget = new double [] {1.0,1.0};
		CellPottsModel model = new CellPottsModel(8, 8, 1, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);

		model.initSpin(spin);
		model.updateCM();
		
		CellDivision division = new CellDivision();
		double [] vec = division.getMajorAxis(model, 1);
		assertEquals("return wrong value for major axis",
				0.80537229, vec[0], tol);
		assertEquals("return wrong value for major axis",
				0.592769327, vec[1], tol);
	}

	@Test
	public void testSplitCell1(){
		int [][] spin = new int [][]{
				{0,0,0,0,0,0,0,0},
				{0,0,1,1,1,0,0,0},
				{0,1,1,1,1,1,0,0},
				{0,1,1,1,1,1,1,0},
				{0,0,1,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		double [] areaTarget = new double [] {1.0,1.0};
		CellPottsModel model = new CellPottsModel(8, 8, 1, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		CellDivision division = new CellDivision();
		
		model.addCPMExtension(division);
		model.initSpin(spin);
		model.updateCM();
		
		division.addCellDivisionListener(model);
		
		division.splitCell(model, 1, 1);

		int [][] expectedSpin = new int [][] {
				{0,0,0,0,0,0,0,0},
				{0,0,2,2,2,0,0,0},
				{0,2,2,2,2,2,0,0},
				{0,2,2,2,1,1,1,0},
				{0,0,2,1,1,1,1,0},
				{0,0,1,1,1,1,0,0},
				{0,0,0,1,1,0,0,0},
				{0,0,0,0,0,0,0,0}
		};

		for (int i = 0; i < spin.length; i++){
			for (int j = 0; j < spin[0].length; j++){
				spin[i][j] = model.getSpin(i, j);
			}
		}
		assertArrayEquals("splitted cell incorrectly", expectedSpin, spin);
	}

}
