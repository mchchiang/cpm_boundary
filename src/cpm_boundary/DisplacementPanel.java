package cpm_boundary;

import javax.swing.*;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Timer;
import java.util.TimerTask;

@SuppressWarnings("serial")
public class DisplacementPanel extends JPanel implements DataListener {

	private CellPottsModel model = null;
	private BufferedImage fg = null;
	private Graphics2D fgGraphics = null;
	private Object lock = new Object();			
	private Timer timer = null;

	private int arrowSize = 2;
	private int width;
	private int height;

	private final double qPI = 0.25*Math.PI;

	private ArrayList<Integer> avgDColour;
	private double maxD;


	public void setModel(CellPottsModel model){
		if (this.model != null){
			this.model.removeDataListener(this);
		}
		this.model = model;
		this.model.addDataListener(this);

		width = model.getNumOfColumns() * arrowSize;
		height = model.getNumOfRows() * arrowSize;

		//for the case of null model with zero width and height
		if (width == 0) width = 1;
		if (height == 0) height = 1;

		avgDColour = new ArrayList<Integer>();
		maxD = 0.1;

		for (int i = 0; i < model.getTypesOfSpin(); i++){
			addNewCell();
		}

		fg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

		fgGraphics = fg.createGraphics();
		fgGraphics.setBackground(Color.WHITE);
		fgGraphics.setColor(Color.BLACK);
		
		repaint();
	}

	/**
	 * Begin painting the model configuration on screen with fixed 
	 * update rate
	 */
	public void startDrawingImage(){
		if (model != null){
			final JPanel panel = this;
			timer = new Timer();
			timer.scheduleAtFixedRate(new TimerTask() {
				public void run() {
					synchronized(lock){
						drawImage(panel.getGraphics());
					}
				}
			}, 0, 33);
		}
	}

	/**
	 * Stop painting the model configuration on screen
	 */
	public void stopDrawingImage(){
		timer.cancel();
		timer = null;
		repaint();
	}

	public void drawImage(Graphics g){
		g.drawImage(
				fg, 0, this.getInsets().top,
				this.getWidth(),
				this.getHeight() - this.getInsets().top, null);
	}

	@Override
	public void paint(Graphics g){
		super.paint(g);
		g.drawImage(
				fg, 0, this.getInsets().top, 
				this.getWidth(), 
				this.getHeight() - this.getInsets().top, null);
	}

	public void drawArrow(int x1, int y1, double l1, double a) {
		int x2 = x1 + (int) Math.round(l1*Math.cos(a));
		int y2 = y1 + (int) Math.round(l1*Math.sin(a));
		fgGraphics.drawLine(x1, y1, x2, y2);// tail
		double l2 = 0.25*l1;
		double am = a - qPI;
		double ap = a + qPI;
		// tip:
		fgGraphics.drawLine(x2, y2, x2 - (int) Math.round(l2*Math.cos(am)), 
				y2 - (int) Math.round(l2*Math.sin(am)));
		fgGraphics.drawLine(x2, y2, x2 - (int) Math.round(l2*Math.cos(ap)), 
				y2 - (int) Math.round(l2*Math.sin(ap)));
	}

	@Override
	public void update(CellPottsModel model, int time) {

		double diff = model.getTypesOfSpin() - avgDColour.size();

		while (diff > 0){
			addNewCell();
			diff--;
		}
		
		float h = 0.3f, s = 0.0f, b = 1.0f;
		for (int i = 1; i < model.getTypesOfSpin(); i++){
			s = (float) (model.getAvgD(i) / maxD);
			avgDColour.set(i, Color.getHSBColor(h, s, b).getRGB());
		}

		synchronized(lock){
			//draw the magnitude of displacement
			for (int i = 0; i < width; i++){
				for (int j = 0; j < height; j++){
					fg.setRGB(i, j, avgDColour.get(model.getSpin(i/arrowSize, j/arrowSize)));
				}
			}

			//draw the displacement arrows
			double x, y, a, dx, dy;

			for (int i = 1; i < model.getTypesOfSpin(); i++){
				x = model.getXCM(i);
				y = model.getYCM(i);
				dx = model.getAvgDX(i);
				dy = model.getAvgDY(i);
				
				if (dx != 0.0 && dy != 0.0){
					a = Math.atan2(dy, dx);
					//System.out.println(a);
					drawArrow((int) (x * arrowSize), (int) (y * arrowSize), 
							arrowSize * 5, a);
				}
			}
		}
	}

	public void addNewCell(){
		avgDColour.add(Color.WHITE.getRGB());
	}
}
