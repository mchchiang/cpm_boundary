package cpm_boundary;

import java.util.ArrayList;
import java.util.Random;

public class CellDivision implements CPMExtension {

	private int dormantPeriod = 100;
	private ArrayList<Integer> lastDivisionTime;
	private double divisionConst = 5E7;
	private ArrayList<CellDivisionListener> listeners;
	private Object listenersLock = new Object();
	private Random rand;
	
	public CellDivision(){
		listeners = new ArrayList<CellDivisionListener>();
	}
	
	@Override
	public void init(CellPottsModel model) {
		rand = new Random();
		lastDivisionTime = new ArrayList<Integer>();
		listeners = new ArrayList<CellDivisionListener>();
		for (int i = 0; i < model.getNumberOfSpins(); i++){
			lastDivisionTime.add(model.getNEquil() + 
					(int) (model.getNEquil() * rand.nextDouble()));
		}
	}

	@Override
	public void execute(CellPottsModel model, int time) {
		for (int i = 1; i < model.getNumberOfSpins(); i++){
			if (shouldDivide(model, time, i, rand.nextDouble())){
				splitCell(model, time, i);
			}
		}
	}
	
	@Override
	public boolean executePriorEquil(){
		return false;
	}
	
	public void addNewCell() {
		lastDivisionTime.add(0);
	}
	
	
	public void addCellDivisionListener(CellDivisionListener l){
		//synchronized(listenersLock){
			if (!listeners.contains(l)){
				listeners.add(l);
			}
		//}
	}
	
	public void removeCellDivisionListener(CellDivisionListener l){
		//synchronized(listenersLock){
			listeners.remove(l);
		//}
	}
	
	public void notifyCellDivisionListener(){
		//synchronized(listenersLock){
			for (CellDivisionListener l : listeners){
				System.out.println(l);
				l.addNewCell();
			}
		//}
	}
	
	public boolean shouldDivide(CellPottsModel model, int time, int cellIndex, double r){
		int dt = time - lastDivisionTime.get(cellIndex);
		
		if (dt < dormantPeriod || model.getArea(cellIndex) < 1.0 * model.getCellArea()){
			return false;
		} else {
			double sq = model.getBeta() / model.getAlpha() / (double) dt;
			sq *= sq; 
			double prob = 1.0 / (divisionConst * sq + 1);
			if (prob > r){
				return true;
			} else {
				return false;
			}
		}
	}
	
	public void splitCell(CellPottsModel model, int time, int cellIndex){
		System.out.println("split cell");
		addNewCell();
		notifyCellDivisionListener();
		
		int q = model.getNumberOfSpins() - 1;
		double cellArea = model.getCellArea();
		
		double [] w = getMajorAxis(model, cellIndex);

		ArrayList<Vector2D> cellPos = model.getSpinPos(cellIndex);
		
		int x, y;
		double xDiff, yDiff, value;
		double x0 = model.getXCM(cellIndex);
		double y0 = model.getYCM(cellIndex);

		for (Vector2D pt : cellPos){
			x = pt.getX();
			y = pt.getY();
			xDiff = model.xDiff(x + 0.5, x0);
			yDiff = model.yDiff(y + 0.5, y0);
			value = model.dot(xDiff, yDiff, w[0], w[1]);

			if (value < 0){
				model.setSpin(x, y, q);
			}
		}
		
		model.setAreaTarget(cellIndex, cellArea);
		model.setAreaTarget(q, cellArea);
		lastDivisionTime.set(cellIndex, time);
		lastDivisionTime.set(q, time);
	}

	/*
	 * store the eigenvectors of the gyration tensor of a cell as 
	 * rows in a matrix
	 */
	public double [] getMajorAxis(CellPottsModel model, int cellIndex){
		//compute eigenvalues
		double sxx = gyrationTensor(model, cellIndex, 0, 0);
		double syy = gyrationTensor(model, cellIndex, 1, 1);
		double sxy = gyrationTensor(model, cellIndex, 0, 1);
		double sum = sxx + syy;
		double sqrt = Math.sqrt(sum * sum - 4.0 * (sxx * syy - sxy * sxy));

		/*
		 * compute the largest eigenvalue as that corresponds to the eigenvector
		 * that describes the major axis
		 */
		double eigenval = (sum + sqrt) / 2.0;
		double [] eigenvec = new double [2];
		double x, y, len;
		x = 1.0;
		y = (eigenval - sxx) / sxy;
		len = Math.sqrt(model.mag2(x,y));
		x /= len;
		y /= len;
		eigenvec[0] = x;
		eigenvec[1] = y;

		return eigenvec;
	}

	//return the (i,j) element of the gyration tensor
	public double gyrationTensor(CellPottsModel model, int cellIndex, int i, int j){
		double avg = 0.0;
		ArrayList<Vector2D> pos = model.getSpinPos(cellIndex);
		int n = pos.size();
		if (i == 0 && j == 0){
			double x0 = model.getXCM(cellIndex);
			double dx;
			for (int k = 0; k < n; k++){
				dx = model.xDiff(pos.get(k).getX() + 0.5, x0);
				avg += dx * dx;
			}
			return avg / (double) n;

		} else if (i == 1 && j == 1){
			double y0 = model.getYCM(cellIndex);
			double dy;
			for (int k = 0; k < n; k++){
				dy = model.yDiff(pos.get(k).getY() + 0.5, y0);
				avg += dy * dy;
			}
			return avg / (double) n;

		} else if (i == 0 && j == 1 || i == 1 && j == 0){
			double x0 = model.getXCM(cellIndex);
			double y0 = model.getYCM(cellIndex);
			double dx, dy;
			Vector2D pt;
			for (int k = 0; k < n; k++){
				pt = pos.get(k);
				dx = model.xDiff(pt.getX() + 0.5, x0);
				dy = model.yDiff(pt.getY() + 0.5, y0);
				avg += dx * dy;
			}
			return avg / (double) n;
		}
		return 0.0;
	}
}
