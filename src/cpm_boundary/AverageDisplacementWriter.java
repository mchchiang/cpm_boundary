package cpm_boundary;

public class AverageDisplacementWriter extends DataWriter {
	
	private int measureTime;
	
	public AverageDisplacementWriter(int measureTime){
		this.measureTime = measureTime;
	}

	@Override
	public void writeData(CellPottsModel model, int time) {
		if (time == measureTime){
			System.out.println("writing data");
			//find growing cell
			int spin = model.getSpin(model.getNumOfColumns()/2, model.getNumOfRows()/2);
			double r, x, y, dx, dy;
			double x0 = model.getXCM(spin);
			double y0 = model.getYCM(spin);
			for (int i = 1; i < model.getTypesOfSpin(); i++){
				if (i != spin){
					x = model.getXCM(i);
					y = model.getYCM(i);
					dx = model.xDiff(x, x0);
					dy = model.yDiff(y, y0);
					r = Math.sqrt(model.mag2(dx, dy));
					writer.printf("%.5f %.5f %.5f %.5f\n", r, 
							model.getAvgDX(i), 
							model.getAvgDY(i), 
							model.getAvgD(i));
				}
			}
		}
	}
}
