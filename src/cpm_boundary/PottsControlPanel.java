package cpm_boundary;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 * PottsControlPanel.java
 * 
 * The parameter control panel on the GUI. It processes the parameter
 * inputs by the user and initialise the simulation when the "run"
 * button is clicked.
 * 
 * @author Michael Chiang
 *
 */
@SuppressWarnings("serial")
public class PottsControlPanel extends JPanel implements ActionListener {
	private JPanel simParamsPanel;
	private JPanel modelParamsPanel;

	private JLabel lblWidth;
	private JTextField txtWidth;

	private JLabel lblHeight;
	private JTextField txtHeight;

	private JLabel lblQ;
	private JTextField txtQ;

	private JLabel lblAlpha;
	private JTextField txtAlpha;

	private JLabel lblBeta;
	private JTextField txtBeta;

	private JLabel lblMotility;
	private JTextField txtMotility;

	private JLabel lblRotateDiff;
	private JTextField txtRotateDiff;

	private JLabel lblLambda;
	private JTextField txtLambda;

	private JLabel lblTemp;
	private JTextField txtTemp;

	private JLabel lblGrowthRate;
	private JTextField txtGrowthRate;

	private JLabel lblFracOccupied;
	private JTextField txtFracOccupied;

	private JLabel lblNumOfSteps;
	private JTextField txtNumOfSteps;

	private JLabel lblNEquil;
	private JTextField txtNEquil;

	private JButton btnRun;
	private JButton btnStop;
	private JButton btnPause;

	private PottsView view;

	private CellPottsModel model;
	private MyThread runThread;

	public PottsControlPanel(PottsView view){
		this.view = view;

		lblWidth = new JLabel("Width: ");
		txtWidth = new JTextField(3);

		lblHeight = new JLabel("Height: ");
		txtHeight = new JTextField(3);

		lblQ = new JLabel("Number of Cells: ");
		txtQ = new JTextField(3);

		txtLambda = new JTextField(3);
		lblLambda = new JLabel("Lambda: ");

		txtTemp = new JTextField(3);
		lblTemp = new JLabel("Temp: ");

		txtAlpha = new JTextField(3);
		lblAlpha = new JLabel("Alpha: ");

		txtBeta = new JTextField(3);
		lblBeta = new JLabel("Beta: ");

		txtMotility = new JTextField(3);
		lblMotility = new JLabel("Motility: ");

		txtRotateDiff = new JTextField(3);
		lblRotateDiff = new JLabel("Rotate Diff: ");

		txtGrowthRate = new JTextField(3);
		lblGrowthRate = new JLabel("Growth Rate: ");

		txtNEquil = new JTextField(3);
		lblNEquil = new JLabel("nequil: ");

		txtFracOccupied = new JTextField(3);
		lblFracOccupied = new JLabel("Fraction occupied: ");

		txtNumOfSteps = new JTextField(3);
		lblNumOfSteps = new JLabel("Steps: ");

		btnRun = new JButton("Run");
		btnRun.addActionListener(this);

		btnStop = new JButton("Stop");
		btnStop.addActionListener(this);
		btnStop.setEnabled(false);

		btnPause = new JButton("Pause");
		btnPause.addActionListener(this);
		btnPause.setEnabled(false);

		modelParamsPanel = new JPanel();
		modelParamsPanel.add(lblAlpha);
		modelParamsPanel.add(txtAlpha);
		modelParamsPanel.add(lblBeta);
		modelParamsPanel.add(txtBeta);
		modelParamsPanel.add(lblLambda);
		modelParamsPanel.add(txtLambda);
		modelParamsPanel.add(lblTemp);
		modelParamsPanel.add(txtTemp);
		modelParamsPanel.add(lblMotility);
		modelParamsPanel.add(txtMotility);
		modelParamsPanel.add(lblRotateDiff);
		modelParamsPanel.add(txtRotateDiff);
		modelParamsPanel.add(lblGrowthRate);
		modelParamsPanel.add(txtGrowthRate);

		simParamsPanel = new JPanel();
		simParamsPanel.add(lblWidth);
		simParamsPanel.add(txtWidth);
		simParamsPanel.add(lblHeight);
		simParamsPanel.add(txtHeight);
		simParamsPanel.add(lblQ);
		simParamsPanel.add(txtQ);
		simParamsPanel.add(lblFracOccupied);
		simParamsPanel.add(txtFracOccupied);
		simParamsPanel.add(lblNumOfSteps);
		simParamsPanel.add(txtNumOfSteps);
		simParamsPanel.add(lblNEquil);
		simParamsPanel.add(txtNEquil);
		simParamsPanel.add(btnRun);
		simParamsPanel.add(btnStop);
		simParamsPanel.add(btnPause);

		setLayout(new BorderLayout());
		add(modelParamsPanel, BorderLayout.NORTH);
		add(simParamsPanel, BorderLayout.SOUTH);
	}


	@Override
	public void actionPerformed(ActionEvent e) {
		/*
		 * run the simulation with the parameters entered by the user 
		 * when the "run" button is clicked
		 */
		
		if (e.getSource() == btnRun){

			int nx = Integer.parseInt(txtWidth.getText());
			int ny = Integer.parseInt(txtHeight.getText());
			int q = Integer.parseInt(txtQ.getText());
			double temp = Double.parseDouble(txtTemp.getText());
			double lambda = Double.parseDouble(txtLambda.getText());
			double alpha = Double.parseDouble(txtAlpha.getText());
			double beta = Double.parseDouble(txtBeta.getText());
			double motility = Double.parseDouble(txtMotility.getText());
			double rotateDiff = Double.parseDouble(txtRotateDiff.getText());
			double growthRate = Double.parseDouble(txtGrowthRate.getText());
			double fracOccupied = Double.parseDouble(txtFracOccupied.getText());
			int seed = -1;
			int numOfSweeps = Integer.parseInt(txtNumOfSteps.getText());
			int nequil = Integer.parseInt(txtNEquil.getText());

			model = new CellPottsModel(
					nx, ny, q, temp, lambda, alpha, beta, motility, rotateDiff,
					growthRate, fracOccupied, seed, numOfSweeps, nequil, new DataWriter [] {}, true);

			runThread = new MyThread();
			runThread.start();

		} else if (e.getSource() == btnStop){
			btnPause.setEnabled(false);
			btnPause.setText("Pause");
			model.stop();
		} else if (e.getSource() == btnPause){
			if (model.isPaused()){
				model.resume();
				btnPause.setText("Pause");
			} else {
				btnPause.setText("Resume");
				model.pause();
			}
		}
	}
	
	class MyThread extends Thread {
		@Override
		public void run(){
			btnRun.setEnabled(false);
			btnStop.setEnabled(true);
			btnPause.setEnabled(true);
			
			model.initSpin();
			model.initPolarity();

			view.setModel(model);
			view.initImage();

			model.run();
			view.stopDrawingImage();

			btnRun.setEnabled(true);
			btnStop.setEnabled(false);
			btnPause.setEnabled(false);
		}
	}
}
