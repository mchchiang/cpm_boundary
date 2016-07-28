package cpm_boundary;

/**
 * Write the mean square displacement of the cells (from their initial
 * positions) at each MCS to file
 * 
 * @author Michael Chiang
 *
 */
public class R2Writer extends DataWriter {
	
	private int intervalIndex;
	
	public R2Writer(){
		this(0);
	}
	
	public R2Writer(int intervalIndex){
		this.intervalIndex = intervalIndex;
	}
	
	@Override
	public void openWriter(String filename){
		super.openWriter(filename);
		writer.println("# time, r2, r2_error\n");
	}
	
	@Override
	public void writeData(CellPottsModel model, int time) {
		if (time >= model.getR2StartPoint(intervalIndex)){
			double [] r2Data = model.calculateR2(intervalIndex);
			writer.printf("%d %.8f %.8f\n", time, r2Data[0], r2Data[1]);	
		}
	}
}
