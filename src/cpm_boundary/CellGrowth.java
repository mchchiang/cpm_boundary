package cpm_boundary;

public class CellGrowth implements CPMExtension {

	@Override
	public void init(CellPottsModel model) {}

	@Override
	public void execute(CellPottsModel model, int time) {
		int s = model.getSpin(model.getNumOfColumns()/2, model.getNumOfRows()/2);
		model.setAreaTarget(s, model.getAreaTarget(s) * 1.01);
	}

	@Override
	public boolean executePriorEquil() {
		return false;
	}

}
