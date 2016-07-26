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
		int q = model.getTypesOfSpin();
		writer.printf("%d ", time);
		for (int i = 1; i < q; i++){
			writer.printf("%.8f %.8f ", model.getXCM(i), model.getYCM(i));
		}	
		writer.println();
	}

}
