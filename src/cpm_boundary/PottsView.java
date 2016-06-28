package cpm_boundary;

import java.awt.*;

import javax.swing.*;

/**
 * PottsView.java
 * 
 * Handle the window operations for the GUI. Main class of the program
 * for simulations with visualisation.
 * 
 * @author Michael Chiang
 *
 */

@SuppressWarnings("serial")
public class PottsView extends JFrame {
	private PottsViewPanel viewPanel;
	private DisplacementPanel disPanel;
	private PottsControlPanel controlPanel;
	private EnergyPanel energyPanel;
	
	//constructor
	public PottsView(){
		this.setSize(1700, 500);
		this.setTitle("Cellular Potts Model");
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		this.setResizable(true);
		
		viewPanel = new PottsViewPanel(new NullModel());
		controlPanel = new PottsControlPanel(this);
		disPanel = new DisplacementPanel();
		energyPanel = new EnergyPanel();
		
		JPanel displayPanel = new JPanel();
		displayPanel.setLayout(new GridLayout(0,3));
		displayPanel.add(viewPanel);
		displayPanel.add(disPanel);
		displayPanel.add(energyPanel);
		
		Container content = this.getContentPane();
		content.add(displayPanel, BorderLayout.CENTER);
		content.add(controlPanel, BorderLayout.SOUTH);
		
		this.setVisible(true);
	}
	
	public static void main (String [] args){
		SwingUtilities.invokeLater(new Runnable() {
		    @Override
		    public void run() {
		    	new PottsView(); 
		    }
		});
	}
	
	/**
	 * Set the model to be display on screen
	 * @param model
	 */
	public void setModel(SpinModel model){
		viewPanel.setModel(model);
		disPanel.setModel((CellPottsModel) model); 
		energyPanel.setModel((CellPottsModel) model);
	}
	
	/**
	 * Start drawing the model configuration to screen
	 */
	public void initImage(){
		viewPanel.initImage();
		disPanel.startDrawingImage();
	}
	
	/**
	 * Stop drawing the model configuration to screen
	 */
	public void stopDrawingImage(){
		viewPanel.stopDrawingImage();
		disPanel.stopDrawingImage();
	}
}
