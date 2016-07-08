package cpm_boundary;

public class RoughnessWriter extends DataWriter {

	@Override
	public void writeData(CellPottsModel model, int time) {
		writer.printf("%d %.5f\n", time, model.calculateRoughness());
	}

}
