package cpm_boundary;

public class RoughnessWriter extends DataWriter {

	@Override
	public void writeData(CellPottsModel model, int time) {
		double [] roughness = model.calculateRoughness();
		writer.printf("%d %.5f %.5f\n", time, roughness[0], roughness[1]);
	}

}
