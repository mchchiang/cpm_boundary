package cpm_boundary;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

public class CellPottsModelTest {
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
	
	private final CellPottsModel defaultModel = new CellPottsModel(
			nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
			seed, numOfSweeps, nequil, false);
	
	//convert primitive double array to arraylist
	public ArrayList<Double> toList(double [] array){
		ArrayList<Double> list = new ArrayList<Double>();
		for (double a : array){
			list.add(a);
		}
		return list;
	}
	
	//Test accessor methods

	@Test
	public void testGetSpin1(){
		defaultModel.initSpin();
		defaultModel.setSpin(2, 3, 1);
		assertEquals("Problem in setting/getting spin",
				1, defaultModel.getSpin(2, 3));
	}
	
	@Test
	public void testSetAlpha1(){
		double newAlpha = 123.2;
		defaultModel.setAlpha(newAlpha);
		assertEquals("Set a wrong value for alpha",
				newAlpha, defaultModel.getAlpha(), tol);
	}
	
	@Test
	public void testGetAlpha1(){
		assertEquals("Returned a wrong value for alpha",
				alpha, defaultModel.getAlpha(), tol);
	}
	
	@Test
	public void testSetBeta1(){
		double newBeta = 435.5;
		defaultModel.setBeta(newBeta);
		assertEquals("Set a wrong value for beta",
				newBeta, defaultModel.getBeta(), tol);
	}
	
	@Test
	public void testGetBeta1(){
		assertEquals("Returned a wrong value for beta",
				beta, defaultModel.getBeta(), tol);
	}
	
	@Test
	public void testSetLambda1(){
		double newLambda = 8937.23;
		defaultModel.setLambda(newLambda);
		assertEquals("Set a wrong value for lambda",
				newLambda, defaultModel.getLambda(), tol);
	}
	
	@Test
	public void testGetLambda1(){
		assertEquals("Returned a wrong value for lambda",
				lambda, defaultModel.getLambda(), tol);
	}
	
	@Test
	public void testGetNumOfRows1(){
		assertEquals("Returned a wrong value for the number of rows in lattice",
				ny, defaultModel.getNumOfRows());
	}
	
	@Test
	public void testGetNumOfColumns1(){
		assertEquals("Returned a wrong value for the number of columns in lattice",
				nx, defaultModel.getNumOfColumns());
	}
	
	@Test
	public void testSetTemp1(){
		double newTemp = 4351.3;
		defaultModel.setTemp(newTemp);
		assertEquals("Set a wrong value for the temperature",
				newTemp, defaultModel.getTemp(), tol);
	}
	
	@Test
	public void testGetTemp1(){
		assertEquals("Returned a wrong value for the temperature",
				temperature, defaultModel.getTemp(), tol);
	}
	
	@Test
	public void testSetRotateDiff1(){
		double newRotateDiff = 4351.3;
		defaultModel.setRotateDiff(newRotateDiff);
		assertEquals("Set a wrong value for the rotational diffusion coefficient",
				newRotateDiff, defaultModel.getRotateDiff(), tol);
	}
	
	@Test
	public void testGetRotateDiff1(){
		assertEquals("Returned a wrong value for the rotational diffusion coefficient",
				rotateDiff, defaultModel.getRotateDiff(), tol);
	}
	
	@Test
	public void testSetNEquil1(){
		int newNequil = 3472;
		defaultModel.setNEquil(newNequil);
		assertEquals("Set a wrong value for the number of steps to reach equilibrium (nequil)",
				newNequil, defaultModel.getNEquil());
	}
	
	@Test
	public void testSetNEquil2(){
		int newNequil = -239047;
		defaultModel.setNEquil(newNequil);
		assertEquals("Set a wrong value for the number of steps to reach equilibrium (nequil)",
				nequil, defaultModel.getNEquil());
	}
	
	@Test
	public void testGetNEquil1(){
		assertEquals("Returned a wrong value for the number of steps to reach equilibrium (nequil)",
				nequil, defaultModel.getNEquil(), tol);
	}
	
	@Test
	public void testSetNumOfSweeps1(){
		int newNumOfSweeps = 45668;
		defaultModel.setNumOfSweeps(newNumOfSweeps);
		assertEquals("Set a wrong value for the number of sweeps",
				newNumOfSweeps, defaultModel.getNumOfSweeps());
	}
	
	@Test
	public void testSetNumOfSweeps2(){
		int newNumOfSweeps = -234789239;
		defaultModel.setNumOfSweeps(newNumOfSweeps);
		assertEquals("Set a wrong value for the number of sweeps",
				numOfSweeps, defaultModel.getNumOfSweeps());
	}
	
	@Test
	public void testGetNumOfSweeps1(){
		assertEquals("Returned a wrong value for the number of sweeps",
				numOfSweeps, defaultModel.getNumOfSweeps());
	}
	
	@Test
	public void testSetMotilityConst1(){
		double newMotiltyConst = 34985.423;
		defaultModel.setMotilityConst(newMotiltyConst);
		assertEquals("Set a wrong value for the motility constant",
				newMotiltyConst, defaultModel.getMotilityConst(), tol);
	}
	
	@Test
	public void testSetMotilityConst2(){
		double newMotiltyConst = -238947.34598;
		defaultModel.setMotilityConst(newMotiltyConst);
		assertEquals("Set a wrong value for the motility constant",
				motility, defaultModel.getMotilityConst(), tol);
	}
	
	@Test
	public void testGetMotilityConst1(){
		assertEquals("Returned a wrong value for the motility constant",
				motility, defaultModel.getMotilityConst(), tol);
	}
	
	@Test
	public void testSetMotility1(){
		double newMotility = 90812.324;
		defaultModel.setMotility(1, newMotility);
		assertEquals("Set a wrong value for the motility constant",
				newMotility, defaultModel.getMotility(1), tol);
	}
	
	@Test
	public void testSetMotility2(){
		double newMotility = 905.232;
		defaultModel.setMotility(3, newMotility);
		assertEquals("Set a wrong value for the motility constant",
				newMotility, defaultModel.getMotility(3), tol);
	}
	
	@Test
	public void testSetMotility3(){
		double newMotility = -8234.34534;
		defaultModel.setMotility(1, newMotility);
		assertEquals("Set a wrong value for the motility constant",
				0.0, defaultModel.getMotility(1), tol);
	}
	
	@Test
	public void testGetMotility1(){
		assertEquals("Returned a wrong value for the motility constant for spin 1",
				0.0, defaultModel.getMotility(1), tol);
		assertEquals("Returned a wrong value for the motility constant for spin 2",
				0.0, defaultModel.getMotility(2), tol);
		assertEquals("Returned a wrong value for the motility constant for spin 3",
				0.0, defaultModel.getMotility(3), tol);
		assertEquals("Returned a wrong value for the motility constant for spin 4",
				0.0, defaultModel.getMotility(4), tol);
	}
	
	@Test
	public void testSetAverageInterval1(){
		int avgInt = 45668;
		defaultModel.setAverageInterval(avgInt);
		assertEquals("Set a wrong value for the number of sweeps",
				avgInt, defaultModel.getAverageInterval());
	}
	
	@Test
	public void testSetAverageInterval2(){
		int avgInt = -234;
		defaultModel.setAverageInterval(avgInt);
		assertEquals("Set a wrong value for the number of sweeps",
				1, defaultModel.getAverageInterval());
	}
	
	@Test
	public void testGetAverageInterval1(){
		assertEquals("Returned a wrong value for the number of sweeps",
				1, defaultModel.getAverageInterval());
	}
	
	@Test
	public void testIsRunning1(){
		assertFalse("Returned true even when the model is not running", 
				defaultModel.isRunning());
	}
	
	@Test
	public void testIsPaused1(){
		assertFalse("Returned true even when the model is not running", 
				defaultModel.isPaused());
	}
	
	@Test
	public void testGetTotalTypesOfSpin1(){
		assertEquals("Returned a wrong value for the total types of spin",
				q+1, defaultModel.getTypesOfSpin());
	}
	
	@Test
	public void testGetTotalEnergy1(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Return a wrong value for the total energy",
				(19+23) * alpha + 2 * lambda, model.getTotalEnergy(), tol);
	}
	
	//Test init methods
	@Test
	public void testInitSpin1(){
		int [][] expectedSpin = new int [][]{
				{1,1,1,1},
				{1,2,2,1},
				{1,2,2,1},
				{1,1,1,1}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(expectedSpin);
		
		int nx = expectedSpin.length;
		int ny = expectedSpin[0].length;
		
		int [][] actualSpin = new int [nx][ny];
		
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				actualSpin[i][j] = model.getSpin(i, j);
			}
		}
		
		assertArrayEquals("Initialise spins incorrectly", 
				expectedSpin, actualSpin);
	}
	
	@Test
	public void testInitSpin2(){
		int [][] expectedSpin = new int [][] {
				{1,2,2,1},
				{1,2,3,4},
				{4,3,3,4},
				{1,2,3,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(expectedSpin);
		
		int nx = expectedSpin.length;
		int ny = expectedSpin[0].length;
		
		int [][] actualSpin = new int [nx][ny];
		
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				actualSpin[i][j] = model.getSpin(i, j);
			}
		}
		
		assertArrayEquals("Initialise spins incorrectly", 
				expectedSpin, actualSpin);
	}
	
	@Test
	public void testInitSpin3(){
		int [][] expectedSpin = new int [][] {
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(expectedSpin);
		
		int nx = expectedSpin.length;
		int ny = expectedSpin[0].length;
		
		int [][] actualSpin = new int [nx][ny];
		
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				actualSpin[i][j] = model.getSpin(i, j);
			}
		}
		
		assertArrayEquals("Initialise spins incorrectly", 
				expectedSpin, actualSpin);
	}
	
	@Test
	public void testInitSpin4(){
		int [][] expectedSpin = new int [][]{
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
		
		model.initSpin(expectedSpin);
		
		int nx = expectedSpin.length;
		int ny = expectedSpin[0].length;
		
		int [][] actualSpin = new int [nx][ny];
		
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				actualSpin[i][j] = model.getSpin(i, j);
			}
		}
		
		assertArrayEquals("Initialise spins incorrectly", 
				expectedSpin, actualSpin);
	}
	
	@Test
	public void testInitMotility1(){
		int numOfMotileCells = 2;
		defaultModel.initMotility(numOfMotileCells);
		double cellMotility;
		int count = 0;
		for (int i = 1; i < defaultModel.getTypesOfSpin(); i++){
			cellMotility = defaultModel.getMotility(i);
			if (cellMotility > 0.0){
				count++;
				assertEquals("Returned wrong motility strength for spin " + i,
					motility, cellMotility, tol);
			} else {
				assertEquals("Returned wrong motility strength for spin " + i,
						0.0, cellMotility, tol);
			}
		}
		assertEquals("Initialise wrong number of motile cells", 
				numOfMotileCells, count);
	}
	
	@Test
	public void testInitMotility2(){
		int numOfMotileCells = (int) (q * 0.5);
		defaultModel.initMotility(0.5);
		double cellMotility;
		int count = 0;
		for (int i = 1; i < defaultModel.getTypesOfSpin(); i++){
			cellMotility = defaultModel.getMotility(i);
			if (cellMotility > 0.0){
				count++;
				assertEquals("Returned wrong motility strength for spin " + i,
					motility, cellMotility, tol);
			} else {
				assertEquals("Returned wrong motility strength for spin " + i,
						0.0, cellMotility, tol);
			}
		}
		assertEquals("Initialise wrong number of motile cells", 
				numOfMotileCells, count);
	}
	
	//Test other methods	
	@Test
	public void testPottsEnergy1(){
		assertEquals("Returned wrong energy value for a pair of same spins",
				defaultModel.pottsEnergy(3,3), 0.0, tol);
	}
	
	@Test
	public void testPottsEnergy2(){
		assertEquals("Returned wrong energy value for a pair that involves q = 0",
				defaultModel.pottsEnergy(0,2), beta, tol);
	}
	
	@Test
	public void testPottsEnergy3(){
		assertEquals("Returned wrong energy value",
				defaultModel.pottsEnergy(1,2), alpha, tol);
	}
	
	@Test
	public void testNegDeltaE1(){
		int [][] spin = new int [][] {
				{1,2,2,1},
				{1,2,3,4},
				{4,3,3,4},
				{1,2,3,4}
		};
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0, 4.0};
		double lambda = 2;
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		//simulate a switch in spin
		model.setSpin(2, 1, 4);
		
		double negDE = model.negDeltaE(2, 1, 4, 4.0, 4.0, 3.0, 5.0);
		assertEquals("Returned wrong negative delta E value",
				-4.0, negDE, tol);
	}
	
	@Test
	public void testNegDeltaE2(){
		int [][] spin = new int [][] {
				{1,2,2,1},
				{1,2,3,4},
				{4,3,3,4},
				{1,2,3,4}
		};
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0, 4.0};
		double lambda = 2;
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		//simulate a switch in spin
		model.setSpin(2, 1, 4);
		
		double negDE = model.negDeltaE(2, 1, 4, 4.0, 4.0, 3.0, 5.0);
		assertEquals("Returned wrong negative delta E value",
				-4.0, negDE, tol);
	}
	
	@Test
	public void testCalculateCM1_1Xa(){
		/* test for the configuration
				{1,1,1,1}
				{1,2,2,1}
				{1,2,2,1}
				{1,1,1,1}
		*/
		
		CellPottsModel model = new CellPottsModel(
				nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		ArrayList<Integer> xPos = new ArrayList<Integer>();
		xPos.add(1);
		xPos.add(1);
		xPos.add(2);
		xPos.add(2);
		assertEquals("Returned wrong xcm value for spin 2",
				2.0, model.calculateCM(xPos, 4), tol);
		
	}
	
	//same test as testCalculateCM1a but with internal array list
	@Test
	public void testCalculateCM1_1Xb(){
		int [][] spin = new int [][]{
				{1,1,1,1},
				{1,2,2,1},
				{1,2,2,1},
				{1,1,1,1}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 2",
				2.0, model.calculateCM(model.getSpinXPos(2), 4), tol);
	}
	
	
	@Test
	public void testCalculateCM2_1X(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 3",
				0.5, model.calculateCM(model.getSpinXPos(3), 4), tol);
	}
	
	@Test
	public void testCalculateCM2_1Y(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 3",
				1.75, model.calculateCM(model.getSpinYPos(3), 4), tol);
	}

	@Test
	public void testCalculateCM2_2X(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 4",
				19.0/6.0, model.calculateCM(model.getSpinXPos(4), 4), tol);
	}
	
	@Test
	public void testCalculateCM2_2Y(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 4",
				1.0/6.0, model.calculateCM(model.getSpinYPos(4), 4), tol);
	}
	
	@Test
	public void testInitCM2(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned wrong xcm value for spin 3",
				0.5, model.getXCM(3), tol);
		assertEquals("Returned wrong ycm value for spin 3",
				1.75, model.getYCM(3), tol);
		assertEquals("Returned wrong xcm value for spin 4",
				19.0/6.0, model.getXCM(4), tol);
		assertEquals("Returned wrong ycm value for spin 4",
				1.0/6.0, model.getYCM(4), tol);
	}
	
	@Test
	public void testCalculateCM3_1X(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 1",
				3.5, model.calculateCM(model.getSpinXPos(1), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_1Y(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 1",
				5.5, model.calculateCM(model.getSpinYPos(1), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_2X(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 2",
				17.0/14.0, model.calculateCM(model.getSpinXPos(2), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_2Y(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 2",
				37.0/14.0, model.calculateCM(model.getSpinYPos(2), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_3X(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 3",
				103.0/18.0, model.calculateCM(model.getSpinXPos(3), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_3Y(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 3",
				13.0/18.0, model.calculateCM(model.getSpinYPos(3), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_4X(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 4",
				5.0/6.0, model.calculateCM(model.getSpinXPos(4), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_4Y(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 4",
				4.5, model.calculateCM(model.getSpinYPos(4), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_5X(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong xcm value for spin 5",
				3.5, model.calculateCM(model.getSpinXPos(5), 6), tol);
	}
	
	@Test
	public void testCalculateCM3_5Y(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		assertEquals("Returned wrong ycm value for spin 5",
				2.5, model.calculateCM(model.getSpinYPos(5), 6), tol);
	}
	
	@Test
	public void testInitCM3(){
		int [][] spin = new int [][]{
				{3,2,2,2,4,3},
				{3,2,2,2,4,4},
				{1,5,5,2,4,1},
				{1,5,5,5,1,1},
				{3,3,5,5,1,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {7.0, 7.0, 7.0, 7.0, 7.0};
		CellPottsModel model = new CellPottsModel(6, 6, 5, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned wrong xcm value for spin 1",
				3.5, model.getXCM(1), tol);
		assertEquals("Returned wrong ycm value for spin 1",
				5.5, model.getYCM(1), tol);
		assertEquals("Returned wrong xcm value for spin 2",
				17.0/14.0, model.getXCM(2), tol);
		assertEquals("Returned wrong ycm value for spin 2",
				37.0/14.0, model.getYCM(2), tol);
		assertEquals("Returned wrong xcm value for spin 3",
				103.0/18.0, model.getXCM(3), tol);
		assertEquals("Returned wrong ycm value for spin 3",
				13.0/18.0, model.getYCM(3), tol);
		assertEquals("Returned wrong xcm value for spin 4",
				5.0/6.0, model.getXCM(4), tol);
		assertEquals("Returned wrong ycm value for spin 4",
				4.5, model.getYCM(4), tol);
		assertEquals("Returned wrong xcm value for spin 5",
				3.5, model.getXCM(5), tol);
		assertEquals("Returned wrong ycm value for spin 5",
				2.5, model.getYCM(5), tol);
	}
	
	@Test
	public void testCalculateCM4_1X(){
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
		assertEquals("Returned wrong xcm value for spin 1",
				3.82, model.calculateCM(model.getSpinXPos(1), 8), tol);
	}
	
	@Test
	public void testCalculateCM4_1Y(){
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
		assertEquals("Returned wrong ycm value for spin 1",
				3.94, model.calculateCM(model.getSpinYPos(1), 8), tol);
	}
	
	@Test
	@Deprecated
	public void testCalculateDeltaCM1a(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertEquals("Returned wrong delta xcm value",
				-1.0/3.0, model.calculateDeltaCM(1, 1, 3, true)[0], tol);
	}
	
	@Test
	@Deprecated
	public void testCalculateDeltaCM1b(){
		int [][] spin = new int [][]{
				{1,3,3,1},
				{1,3,2,1},
				{4,2,2,2},
				{4,3,2,4}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(4, 4, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertEquals("Returned wrong delta ycm value",
				7.0/3.0-2.25, model.calculateDeltaCM(1, 1, 3, true)[1], tol);
	}
	
	@Test
	public void testHasSameNeighbours1(){
		int [][] spin = new int [][]{
				{3,3,4,4,4,3},
				{3,3,4,4,3,3},
				{1,1,2,2,2,1},
				{1,1,2,2,2,1},
				{1,2,2,2,2,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(6, 6, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertTrue("Expect to have same neighbours", model.hasSameNeighbours(0, 0));	
	}
	
	@Test
	public void testHasSameNeighbours2(){
		int [][] spin = new int [][]{
				{3,3,4,4,4,3},
				{3,3,4,4,3,3},
				{1,1,2,2,2,1},
				{1,1,2,2,2,1},
				{1,2,2,2,2,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(6, 6, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertTrue("Expect to have same neighbours", model.hasSameNeighbours(3, 3));
	}
	
	@Test
	public void testHasSameNeighbours3(){
		int [][] spin = new int [][]{
				{3,3,4,4,4,3},
				{3,3,4,4,3,3},
				{1,1,2,2,2,1},
				{1,1,2,2,2,1},
				{1,2,2,2,2,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(6, 6, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertFalse("Expect to have same neighbours", model.hasSameNeighbours(3, 0));
	}
	
	@Test
	public void testHasSameNeighbours4(){
		int [][] spin = new int [][]{
				{3,3,4,4,4,3},
				{3,3,4,4,3,3},
				{1,1,2,2,2,1},
				{1,1,2,2,2,1},
				{1,2,2,2,2,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(6, 6, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertFalse("Expect to have same neighbours", model.hasSameNeighbours(1, 3));
	}
	
	@Test
	public void testHasSameNeighbours5(){
		int [][] spin = new int [][]{
				{3,3,4,4,4,3},
				{3,3,4,4,3,3},
				{1,1,2,2,2,1},
				{1,1,2,2,2,1},
				{1,2,2,2,2,1},
				{3,3,3,4,4,3}
		};
		
		double [] areaTarget = new double [] {4.0, 4.0, 4.0, 4.0};
		CellPottsModel model = new CellPottsModel(6, 6, 4, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		
		assertFalse("Expect to have same neighbours", model.hasSameNeighbours(4, 2));
	}
	
	
	@Test
	public void testXDiff1(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				-2.0 , model.xDiff(2, 4), tol);
	}
	
	@Test
	public void testXDiff2(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				3.0 , model.xDiff(5, 2), tol);
	}
	
	//boundary cases
	@Test
	public void testXDiff3(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				3.0 , model.xDiff(2, 6), tol);
	}
	
	@Test
	public void testXDiff4(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				-2.0 , model.xDiff(6, 1), tol);
	}
	
	@Test
	public void testYDiff1(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				-2.0 , model.yDiff(2, 4), tol);
	}
	
	@Test
	public void testYDiff2(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				3.0 , model.yDiff(5, 2), tol);
	}
	
	//boundary cases
	@Test
	public void testYDiff3(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				3.0 , model.xDiff(2, 6), tol);
	}
	
	@Test
	public void testYDiff4(){
		CellPottsModel model = new CellPottsModel(
				7, 7, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		model.initSpin();
		
		assertEquals("Return wrong xDiff value",
				-2.0 , model.yDiff(6, 1), tol);
	}
	
	@Test
	public void testDot1(){
		CellPottsModel model = new CellPottsModel(
				nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		assertEquals("Return wrong dot product",
				0.0, model.dot(2.0, 0.0, 0.0, 7.0), tol);
	}
	
	@Test
	public void testDot2(){
		CellPottsModel model = new CellPottsModel(
				nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		assertEquals("Return wrong dot product",
				154.0, model.dot(15.0, 2.0, 8.0, 17.0), tol);
	}
	
	@Test
	public void testMag2_1(){
		CellPottsModel model = new CellPottsModel(
				nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		assertEquals("Return wrong magnitude squared value",
				0.0, model.mag2(0.0, 0.0), tol);
	}
	
	@Test
	public void testMag2_2(){
		CellPottsModel model = new CellPottsModel(
				nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		assertEquals("Return wrong magnitude squared value",
				25.0, model.mag2(3.0, 4.0), tol);
	}
	
	@Test
	public void testMag2_3(){
		CellPottsModel model = new CellPottsModel(
				nx, ny, q, temperature, lambda, alpha, beta, motility, rotateDiff,
				seed, numOfSweeps, nequil, false);
		assertEquals("Return wrong magnitude squared value",
				169.0, model.mag2(5.0, 12.0), tol);
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
		model.initCM();
		assertEquals("return wrong gyration tensor comp xx value",
				2.1376, model.gyrationTensor(1, 0, 0), tol);
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
		model.initCM();
		assertEquals("return wrong gyration tensor comp yy value",
				1.9264, model.gyrationTensor(1, 1, 1), tol);
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
		model.initCM();
		assertEquals("return wrong gyration tensor comp xy value",
				0.3392, model.gyrationTensor(1, 0, 1), tol);
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
		model.initCM();
		assertEquals("return wrong gyration tensor comp xy value",
				0.3392, model.gyrationTensor(1, 1, 0), tol);
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
		model.initCM();
		double [] vec = model.getMajorAxis(1);
		assertEquals("return wrong value for major axis",
				0.80537229, vec[0], tol);
		assertEquals("return wrong value for major axis",
				0.592769327, vec[1], tol);
		assertEquals("eigenvectors are not orthonormal",
				0.0, vec[0] * vec[2] + vec[1] * vec[3], tol);
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
		
		model.initSpin(spin);
		model.initCM();
		model.splitCell(1, 1);

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
	
	@Test
	public void testCalculatePositionSum1_1X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(1);
		assertEquals("Returned a wrong value for the x position sum for spin 1",
				5.0, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_1Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(1);
		assertEquals("Returned a wrong value for the y position sum for spin 1",
				3.0, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_2X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(2);
		assertEquals("Returned a wrong value for the x position sum for spin 2",
				31.5, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_2Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(2);
		assertEquals("Returned a wrong value for the y position sum for spin 2",
				76.5, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_3X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(3);
		assertEquals("Returned a wrong value for the x position sum for spin 3",
				8.0, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_3Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(3);
		assertEquals("Returned a wrong value for the y position sum for spin 3",
				16.0, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_4X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(4);
		assertEquals("Returned a wrong value for the x position sum for spin 4",
				38.5, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_4Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(4);
		assertEquals("Returned a wrong value for the y position sum for spin 4",
				52.5, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_6X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(6);
		assertEquals("Returned a wrong value for the x position sum for spin 6",
				28, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_6Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(6);
		assertEquals("Returned a wrong value for the y position sum for spin 6",
				22, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_7X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(7);
		assertEquals("Returned a wrong value for the x position sum for spin 7",
				0.0, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_7Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(7);
		assertEquals("Returned a wrong value for the y position sum for spin 7",
				17.0, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_8X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(8);
		assertEquals("Returned a wrong value for the x position sum for spin 8",
				38.5, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_8Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(8);
		assertEquals("Returned a wrong value for the y position sum for spin 8",
				31.5, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_9X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> xPos = model.getSpinXPos(9);
		assertEquals("Returned a wrong value for the x position sum for spin 9",
				22.0, model.calculatePositionSum(xPos, 8), tol);
	}
	
	@Test
	public void testCalculatePositionSum1_9Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		ArrayList<Integer> yPos = model.getSpinYPos(9);
		assertEquals("Returned a wrong value for the y position sum for spin 9",
			7.0, model.calculatePositionSum(yPos, 8), tol);
	}
	
	@Test
	public void testCalculateDXCM1_1(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 1",
				1.0/11.0, model.calculateDXCM(1, 1, false), tol);
	}
	
	@Test
	public void testCalculateDYCM1_1(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 1",
				-9.0/55.0, model.calculateDYCM(1, 6, false), tol);
	}
	
	@Test
	public void testCalculateDXCM1_2(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 2",
				3.0/22.0, model.calculateDXCM(2, 1, true), tol);
	}
	
	@Test
	public void testCalculateDYCM1_2(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 2",
				1.0/22.0, model.calculateDYCM(2, 6, true), tol);
	}
	
	@Test
	public void testCalculateDXCM1_3(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 3",
				-0.1, model.calculateDXCM(3, 1, false), tol);
	}
	
	@Test
	public void testCalculateDYCM1_3(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 3",
				-0.3, model.calculateDYCM(3, 2, false), tol);
	}
	
	@Test
	public void testCalculateDXCM1_4(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 4",
				-1.0/6.0, model.calculateDXCM(4, 6, true), tol);
	}
	
	@Test
	public void testCalculateDYCM1_4(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 4",
				-1.0/6.0, model.calculateDYCM(4, 0, true), tol);
	}
	
	@Test
	public void testCalculateDXCM1_6(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 6",
				0.0, model.calculateDXCM(6, 3, true), tol);
	}
	
	@Test
	public void testCalculateDYCM1_6(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 6",
				-0.25, model.calculateDYCM(6, 4, true), tol);
	}
	
	@Test
	public void testCalculateDXCM1_7(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 7",
				-0.3, model.calculateDXCM(7, 1, true), tol);
	}
	
	@Test
	public void testCalculateDYCM1_7(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 7",
				1.0/15.0, model.calculateDYCM(7, 2, true), tol);
	}
	
	@Test
	public void testCalculateDXCM1_8(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 8",
				-0.25, model.calculateDXCM(8, 3, false), tol);
	}
	
	@Test
	public void testCalculateDYCM1_8(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 8",
				0.0, model.calculateDYCM(8, 4, false), tol);
	}
	
	@Test
	public void testCalculateDXCM1_9(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dxcm for spin 9",
				0.2, model.calculateDXCM(9, 6, false), tol);
	}
	
	@Test
	public void testCalculateDYCM1_9(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		assertEquals("Returned a wrong value for dycm for spin 9",
				-0.25, model.calculateDYCM(9, 0, false), tol);
	}
	
	@Test
	public void testUpdateCM1_1X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(1, 1, 6, false);
		assertEquals("Didn't update the x cm for spin 1 correctly",
				13.0/22.0, model.getXCM(1), tol);
	}
	
	@Test
	public void testUpdateCM1_1Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(1, 1, 6, false);
		assertEquals("Didn't update the y cm for spin 1 correctly",
				3.0/22.0, model.getYCM(1), tol);
	}
	
	@Test
	public void testUpdateCM1_2X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(2, 1, 6, true);
		assertEquals("Didn't update the x cm for spin 2 correctly",
				3.0, model.getXCM(2), tol);
	}
	
	@Test
	public void testUpdateCM1_2Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(2, 1, 6, true);
		assertEquals("Didn't update the y cm for spin 2 correctly",
				7.0, model.getYCM(2), tol);
	}
	
	@Test
	public void testUpdateCM1_3X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(3, 1, 2, false);
		assertEquals("Didn't update the x cm for spin 3 correctly",
				1.9, model.getXCM(3), tol);
	}
	
	@Test
	public void testUpdateCM1_3Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(3, 1, 2, false);
		assertEquals("Didn't update the y cm for spin 3 correctly",
				3.7, model.getYCM(3), tol);
	}
	
	@Test
	public void testUpdateCM1_4X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(4, 6, 0, true);
		assertEquals("Didn't update the x cm for spin 4 correctly",
				16.0/3.0, model.getXCM(4), tol);
	}
	
	@Test
	public void testUpdateCM1_4Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(4, 6, 0, true);
		assertEquals("Didn't update the y cm for spin 4 correctly",
				22.0/3.0, model.getYCM(4), tol);
	}
	
	@Test
	public void testUpdateCM1_6X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(6, 3, 4, true);
		assertEquals("Didn't update the x cm for spin 6 correctly",
				3.5, model.getXCM(6), tol);
	}
	
	@Test
	public void testUpdateCM1_6Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(6, 3, 4, true);
		assertEquals("Didn't update the y cm for spin 6 correctly",
				2.5, model.getYCM(6), tol);
	}
	
	@Test
	public void testUpdateCM1_7X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(7, 1, 2, true);
		assertEquals("Didn't update the x cm for spin 7 correctly",
				7.7, model.getXCM(7), tol);
	}
	
	@Test
	public void testUpdateCM1_7Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(7, 1, 2, true);
		assertEquals("Didn't update the y cm for spin 7 correctly",
				2.9, model.getYCM(7), tol);
	}
	
	@Test
	public void testUpdateCM1_8X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(8, 3, 4, false);
		assertEquals("Didn't update the x cm for spin 8 correctly",
				5.25, model.getXCM(8), tol);
	}
	
	@Test
	public void testUpdateCM1_8Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(8, 3, 4, false);
		assertEquals("Didn't update the y cm for spin 8 correctly",
				4.5, model.getYCM(8), tol);
	}
	
	@Test
	public void testUpdateCM1_9X(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(9, 6, 0, false);
		assertEquals("Didn't update the x cm for spin 9 correctly",
				5.7, model.getXCM(9), tol);
	}
	
	@Test
	public void testUpdateCM1_9Y(){
		int [][] spin = new int [][] {
				{1,1,7,7,5,5,1,1},
				{1,1,7,3,3,2,2,1},
				{2,6,6,3,3,2,2,2},
				{2,6,6,6,6,2,2,2},
				{2,9,6,6,8,8,4,4},
				{4,9,9,8,8,8,4,4},
				{4,9,7,8,8,5,5,4},
				{1,1,7,7,5,5,5,1}
		};
		
		double [] areaTarget = new double [] {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
		CellPottsModel model = new CellPottsModel(8, 8, 9, toList(areaTarget), 
				temperature, lambda, alpha, beta, motility, seed);
		
		model.initSpin(spin);
		model.initCM();
		model.updateCM(9, 6, 0, false);
		assertEquals("Didn't update the y cm for spin 9 correctly",
				1.5, model.getYCM(9), tol);
	}
}
