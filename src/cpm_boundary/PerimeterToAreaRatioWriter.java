package cpm_boundary;

public class PerimeterToAreaRatioWriter extends DataWriter {
	@Override
	public void openWriter(String filename){
		super.openWriter(filename);
		writer.println("# time, ratio for cell 490\n");
	}
	
	@Override
	public void writeData(CellPottsModel model, int time) {
		//double [] ratioData = model.getPerimeterToAreaRatio();
		writer.printf("%d %.5f\n", 
				time, model.getPerimeterToAreaRatio(490));
	}

}
