package cpm_boundary;

import java.util.ArrayList;
import java.util.LinkedList;

public class AverageDisplacement extends CPMMeasurement 
implements CellDivisionListener {

	private ArrayList<LinkedList<Double>> dxData;
	private ArrayList<LinkedList<Double>> dyData;
	private ArrayList<Double> sumDX;
	private ArrayList<Double> sumDY;
	private ArrayList<Double> avgDX;
	private ArrayList<Double> avgDY;
	private ArrayList<Double> avgD;
	private int avgInt;
	
	@Override
	public void init(CellPottsModel model){
		dxData = new ArrayList<LinkedList<Double>>();
		dyData = new ArrayList<LinkedList<Double>>();
		sumDX = new ArrayList<Double>();
		sumDY = new ArrayList<Double>();
		avgDX = new ArrayList<Double>();
		avgDY = new ArrayList<Double>();
		avgD = new ArrayList<Double>();
		avgInt = 100;
		for (int i = 0; i < model.getNumberOfSpins(); i++){
			addNewCell();
		}
	}

	@Override
	public void measure(CellPottsModel model, int time) {
		LinkedList<Double> dxList;
		LinkedList<Double> dyList;

		double dx, dy, sumX, sumY, avgX, avgY, avgDis;

		for (int i = 1; i < model.getNumberOfSpins(); i++){
			dxList = dxData.get(i);
			dyList = dyData.get(i);

			dx = model.getDX(i);
			dy = model.getDY(i);

			dxList.add(dx);
			dyList.add(dy);

			sumX =  sumDX.get(i) + dx;
			sumY =  sumDY.get(i) + dy;

			sumDX.set(i,sumX);
			sumDY.set(i,sumY);

			if (dxList.size() >= avgInt){
				avgX = sumX / (double) dxList.size();
				avgY = sumY / (double) dyList.size();

				avgDX.set(i, avgX);
				avgDY.set(i, avgY);

				avgDis = Math.sqrt(avgX * avgX + avgY * avgY);
				avgD.set(i, avgDis);

				sumDX.set(i, sumX - dxList.remove());
				sumDY.set(i, sumY - dyList.remove());
			}
		}
		
		notifyCPMMeasurementListener(time);
	}

	public double getAvgDX(int cellIndex){
		return avgDX.get(cellIndex);
	}

	public double getAvgDY(int cellIndex){
		return avgDY.get(cellIndex);
	}
	
	public double getAvgD(int cellIndex){
		return avgD.get(cellIndex);
	}

	@Override
	public int getDataSize() {
		return avgDX.size();
	}

	@Override
	public void addNewCell() {
		dxData.add(new LinkedList<Double>());
		dyData.add(new LinkedList<Double>());
		sumDX.add(0.0);
		sumDY.add(0.0);
		avgDX.add(0.0);
		avgDY.add(0.0);
		avgD.add(0.0);
	}

	@Override
	public boolean measurePriorEquil() {
		return false;
	}
}
