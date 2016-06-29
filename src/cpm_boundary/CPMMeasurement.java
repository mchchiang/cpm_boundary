package cpm_boundary;

import java.util.ArrayList;

public abstract class CPMMeasurement {
	
	private ArrayList<CPMMeasurementListener> listeners;
	
	public abstract int getDataSize();
	
	public abstract void init(CellPottsModel model);
	public abstract boolean measurePriorEquil();
	protected abstract void measure(CellPottsModel model, int time);
	
	public CPMMeasurement(){
		listeners = new ArrayList<CPMMeasurementListener>();
	}
	
	public void makeMeasurement(CellPottsModel model, int time) throws Exception {
		int dataSize = getDataSize();
		int numOfSpins = model.getNumberOfSpins();
		if (dataSize != numOfSpins){
			throw new Exception("Data size not equal to number of cells in model:"
					+ " Data size: " + dataSize + " Number of cells:" + numOfSpins);
		}
		measure(model,time);
	}
	
	public void addCPMMeasurementListener(CPMMeasurementListener l){
		if (!listeners.contains(l)){
			listeners.add(l);
		}
	}
	
	public void removeCPMMeasurementListener(CPMMeasurementListener l){
		if (listeners != null){
			listeners.remove(l);
		}
	}
	
	public void notifyCPMMeasurementListener(int time){
		for (CPMMeasurementListener l : listeners){
			l.update(this, time);
		}
	}
}
