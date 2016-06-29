package cpm_boundary;

public interface CPMEnergyExtension extends CPMExtension {
	
	public double getDeltaE(CellPottsModel model, int iold, int jold, 
			int inew, int jnew, int oldSpin, int newSpin);
	
	public double getEnergy(CellPottsModel model);
}
