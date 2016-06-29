package cpm_boundary;

/**
 * CMWriter.java
 * 
 * Write the center of mass of each cell at each MCS to file
 * 
 * @author Michael Chiang
 *
 */
public class CMWriter extends DataWriter {

	@Override
	public void writeData(CellPottsModel model, int time) {
		writer.println();
		int q = model.getNumberOfSpins();
		writer.printf("%d ", time);
		for (int i = 0; i < q; i++){
			writer.printf("%.8f %.8f ", model.getXCM(i), model.getYCM(i));
		}		
	}

}
