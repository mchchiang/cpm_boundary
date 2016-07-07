package cpm_boundary;

import javax.swing.*;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Observable;
import java.util.Observer;
import java.util.Timer;
import java.util.TimerTask;

/**
 * PottsViewPanel.java
 * 
 * Paint the lattice configuration on screen
 * 
 * @author Michael Chiang
 *
 */

@SuppressWarnings("serial")
public class PottsViewPanel extends JPanel implements Observer {

	private SpinModel model;
	private ArrayList<Color> colours;
	private BufferedImage fg = null;
	private Timer timer = null;


	/**
	 * Initialise the view panel
	 * @param model the spin model to be displayed on screen
	 */
	public PottsViewPanel(SpinModel model){
		this.model = model;
		this.model.addObserver(this);
		setColours();
	}

	/**
	 * Set the spin model to be displayed on screen
	 * @param model spin model
	 */
	public void setModel(SpinModel model){
		this.model.deleteObserver(this);
		this.model = model;
		this.model.addObserver(this);
		setColours();
		repaint();
	}

	/**
	 * Set the colour for each spin
	 */
	public void setColours(){
		int typesOfSpin = model.getTypesOfSpin();
		colours = new ArrayList<Color>();
		
		if (model instanceof CellPottsModel){
			CellPottsModel cpm = (CellPottsModel) model;
			for (int i = 0; i < typesOfSpin; i++){
				if (i == 0){
					colours.add(Color.WHITE);
				} else if (cpm.getMotility(i) > 0.0) {
					colours.add(Color.BLACK);
				} else {
					colours.add(generateColour());
				}
			}
		} else {
			for (int i = 0; i < typesOfSpin; i++){
				if (i == 0){
					colours.add(Color.WHITE);
				} else {
					colours.add(generateColour());
				}
			}
		}
	}

	/**
	 * Begin painting the model configuration on screen with fixed 
	 * update rate
	 */
	public void initImage(){
		int width = model.getNumOfColumns();
		int height = model.getNumOfRows();

		fg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

		//draw an initial image of the cell
		for (int i = 0; i < width; i++){
			for (int j = 0; j < height; j++){
				fg.setRGB(i, j, colours.get(model.getSpin(i, j)).getRGB());
			}
		}

		final JPanel panel = this;
		panel.getGraphics().drawImage(fg, 0, 
				panel.getInsets().top, panel.getWidth(), 
				panel.getHeight() - panel.getInsets().top, null);

		timer = new Timer();
		timer.scheduleAtFixedRate(new TimerTask() {
			public void run() {
				panel.getGraphics().drawImage(
						fg, 0, panel.getInsets().top, 
						panel.getWidth(), 
						panel.getHeight() - panel.getInsets().top, null);
			}
		}, 0, 33);
	}

	/**
	 * Stop painting the model configuration on screen
	 */
	public void stopDrawingImage(){
		timer.cancel();
		timer = null;
	}

	@Override
	public void paint(Graphics g){
		super.paint(g);
		g.drawImage(
				fg, 0, this.getInsets().top, 
				this.getWidth(), 
				this.getHeight() - this.getInsets().top, null);
	}

	@Override
	public void update(Observable o, Object spin) {
		Object [] spinIndices = (Object []) spin;
		drawSpin((Integer) spinIndices[0], (Integer) spinIndices[1]);
	}	

	/**
	 * Update the spin value of a particular lattice site on screen
	 * @param i column index of the site
	 * @param j row index of the site
	 */
	public void drawSpin(int i, int j){
		updateColourMap();
		fg.setRGB(i, j, colours.get(model.getSpin(i, j)).getRGB());
	}
	
	//check if new cells are created and need new colours
	public void updateColourMap(){
		int diff = model.getTypesOfSpin() - colours.size();
		if (diff > 0){
			for (int i = 0; i < diff; i++){
				colours.add(generateColour());
			}
		}
	}
	
	public Color generateColour(){
		float s, b, h;
		h = (float) Math.random();
		b = (float) (Math.random() * 0.5 + 0.5);
		s = (float) (Math.random() * 0.5 + 0.5);
		return new Color(Color.HSBtoRGB(h, s, b));
	}
}
