package cpm_boundary;

public class Haptotaxis implements CPMEnergyExtension {

	//variables for haptotaxis
	private double [][] fibrin;
	private double kh = 40;
	private double kn = Math.log(100.0) / 2000.0;
	private double ki = 2 * kn;


	@Override
	public void init(CellPottsModel model) {
		//init fibrin concentration
		int nx = model.getNumOfColumns();
		int ny = model.getNumOfRows();
		fibrin = new double [nx][ny];
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				fibrin[i][j] = 1.0;
			}
		}
	}

	@Override
	public void execute(CellPottsModel model, int time) {
		int nx = model.getNumOfColumns();
		int ny = model.getNumOfRows();
		
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				if (model.getSpin(i, j) != 0){
					fibrin[i][j] *= Math.exp(-ki);
				} else if (model.getSpin(i, model.jup(j)) != 0 ||
						model.getSpin(i, model.jdown(j)) != 0 ||
						model.getSpin(model.iup(i), j) != 0 ||
						model.getSpin(model.idown(i), j) != 0 ||
						model.getSpin(model.iup(i), model.jup(j)) != 0 ||
						model.getSpin(model.iup(i), model.jdown(j)) != 0 ||
						model.getSpin(model.idown(i), model.jup(j)) != 0 ||
						model.getSpin(model.idown(i), model.jdown(j)) != 0){
					fibrin[i][j] *= Math.exp(-kn);
				} 
			}
		}
	}
	
	@Override
	public boolean executePriorEquil(){
		return false;
	}

	@Override
	public double getDeltaE(CellPottsModel model, 
			int iold, int jold, int inew, int jnew, 
			int oldSpin, int newSpin) {
		return kh * (fibrin[inew][jnew] - fibrin[iold][jold]);
	}

	@Override
	public double getEnergy(CellPottsModel model) {
		// TODO Auto-generated method stub
		return 0;
	}

}
