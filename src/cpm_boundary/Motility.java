package cpm_boundary;

import java.util.ArrayList;
import java.util.Random;

public class Motility implements CPMEnergyExtension, CellDivisionListener {
	
	private ArrayList<Double> px;
	private ArrayList<Double> py;
	private ArrayList<Double> theta;
	private ArrayList<Double> motility;
	private double rotateDiff = 0.1;
	private Random rand;

	@Override
	public void init(CellPottsModel model) {
		px = new ArrayList<Double>();
		py = new ArrayList<Double>();
		theta = new ArrayList<Double>();
		motility = new ArrayList<Double>();
		rand = new Random();
		
		for (int i = 0; i < model.getNumberOfSpins(); i++){
			addNewCell();
		}
	}

	@Override
	public void addNewCell() {
		double angle = rand.nextDouble() * 2 * Math.PI;
		theta.add(angle);
		px.add(Math.cos(angle));
		py.add(Math.sin(angle));
	}

	@Override
	public void execute(CellPottsModel model, int time) {
		//update polarity
		double angle;
		for (int i = 1; i < model.getNumberOfSpins(); i++){
			angle = theta.get(i) + Math.sqrt(6 * rotateDiff) * (rand.nextDouble()*2-1);
			theta.set(i, angle);
			px.set(i, Math.cos(angle));
			py.set(i, Math.sin(angle));
		}
	}

	@Override
	public double getDeltaE(CellPottsModel model, 
			int iold, int jold, int inew, int jnew, 
			int oldSpin, int newSpin) {
		double energy = 0.0;
		double [] dcmOld = model.calculateDeltaCM(iold, jold, oldSpin, true);
		double [] dcmNew = model.calculateDeltaCM(iold, jold, newSpin, false);
		energy += motility.get(oldSpin) * model.dot(dcmOld[0], dcmOld[1], px.get(oldSpin), py.get(oldSpin));
		energy += motility.get(newSpin) * model.dot(dcmNew[0], dcmNew[1], px.get(newSpin), py.get(newSpin));
		return energy;
	}

	@Override
	public double getEnergy(CellPottsModel model) {
		return 0.0;
	}
	
	@Override
	public boolean executePriorEquil(){
		return false;
	}
	
	//Accessor methods
	public void setMotility(int cellIndex, double mu){
		motility.set(cellIndex, mu);
	}
	
	public double getMotility(int cellIndex){
		return motility.get(cellIndex);
	}
	
}
