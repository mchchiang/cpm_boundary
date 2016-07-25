package cpm_boundary;

/**
 * Write the lattice configuration after the final MCS to file
 * 
 * @author Michael Chiang
 *
 */
public class SpinWriter extends DataWriter {
	
	private int numOfSweeps;
	private int numOfCells;
	
	/**
	 * Initialise the spin writer
	 * @param sweeps number of sweeps or MCS for the simulation
	 */
	public SpinWriter(int sweeps, int q){
		numOfSweeps = sweeps;
		numOfCells = q;
	}
	
	@Override
	public void writeData(CellPottsModel model, int time) {
		if (numOfSweeps == time){
			int nx = model.getNumOfColumns();
			int ny = model.getNumOfRows();
			writer.println(nx + " " + ny + " " + numOfCells);
			for (int i = 0; i < ny; i++){
				for (int j = 0; j < nx; j++){
					writer.print(model.getSpin(j, i) + " ");
				}
				writer.println();
			}
		}
	}

}
