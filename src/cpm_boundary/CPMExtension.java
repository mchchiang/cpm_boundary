package cpm_boundary;

public interface CPMExtension {
	
	public void init(CellPottsModel model);
	
	public void execute(CellPottsModel model, int time);
	
	public boolean executePriorEquil();
	
}
