package cpm_boundary;

public class PerimeterToAreaRatioWriter extends DataWriter {
	@Override
	public void openWriter(String filename){
		super.openWriter(filename);
		writer.println("# time, ratio, ratio_stdev, ratio_stderror");
	}
	
	@Override
	public void writeData(CellPottsModel model, int time) {
		double [] ratioData = model.getPerimeterToAreaRatio();
		writer.printf("%d %.5f %.5f %.5f\n", 
				time, ratioData[0], ratioData[1], ratioData[2]);
	}

}
