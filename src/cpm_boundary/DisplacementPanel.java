package cpm_boundary;

import javax.swing.*;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;

@SuppressWarnings("serial")
public class DisplacementPanel extends JPanel implements DataListener {

	private CellPottsModel model = null;
	private BufferedImage fg = null;
	private Graphics2D fgGraphics = null;
	private Object lock = new Object();			
	private Timer timer = null;

	private int arrowSize = 10;
	private int width;
	private int height;

	private final double qPI = 0.25*Math.PI;

	public DisplacementPanel(){

	}

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

		fg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

		fgGraphics = fg.createGraphics();
		fgGraphics.setBackground(Color.WHITE);
		fgGraphics.setColor(Color.BLACK);

		updateArrow();
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

	public void updateArrow(){ 
		synchronized(lock){
			//clear current arrow
			fgGraphics.clearRect(0, 0, width, height);

			//draw new arrows
			double x, y, dx, dy, a;
			int nx = model.getNumOfColumns();
			int ny = model.getNumOfRows();

			//don't draw arrows at boundary
			for (int i = 1; i < model.getTypesOfSpin(); i++){
				dx = model.getDX(i);
				dy = model.getDY(i);
				x = model.getXCM(i);
				y = model.getYCM(i);
				a = Math.atan2(dy,dx);
				//if (x + dx < nx && x + dx > 0 &&
				//		y + dy < ny && y + dy > 0){
					drawArrow((int) (x * arrowSize), (int) (y * arrowSize), 
							arrowSize*10, a);
				//}
			}
		}
	}

	@Override
	public void update(CellPottsModel model, int time) {
		updateArrow();
	}
}
