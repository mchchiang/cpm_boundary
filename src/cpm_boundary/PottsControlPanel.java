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
	
	private JLabel lblNumOfSteps;
	private JTextField txtNumOfSteps;
	
	private JButton btnRun;
	
	private PottsView view;
	
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
		
		txtNumOfSteps = new JTextField(3);
		lblNumOfSteps = new JLabel("Steps: ");
		
		btnRun = new JButton("Run");
		btnRun.addActionListener(this);
		
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
		
		simParamsPanel = new JPanel();
		simParamsPanel.add(lblWidth);
		simParamsPanel.add(txtWidth);
		simParamsPanel.add(lblHeight);
		simParamsPanel.add(txtHeight);
		simParamsPanel.add(lblQ);
		simParamsPanel.add(txtQ);
		simParamsPanel.add(lblNumOfSteps);
		simParamsPanel.add(txtNumOfSteps);
		simParamsPanel.add(btnRun);
		
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
			Thread runthread = new Thread(){
				@Override
				public void run(){					
					int nx = Integer.parseInt(txtWidth.getText());
					int ny = Integer.parseInt(txtHeight.getText());
					int q = Integer.parseInt(txtQ.getText());
					double temp = Double.parseDouble(txtTemp.getText());
					double lambda = Double.parseDouble(txtLambda.getText());
					double alpha = Double.parseDouble(txtAlpha.getText());
					double beta = Double.parseDouble(txtBeta.getText());
					double motility = Double.parseDouble(txtMotility.getText());
					double rotateDiff = Double.parseDouble(txtRotateDiff.getText());
					int seed = -1;
					int numOfSweeps = Integer.parseInt(txtNumOfSteps.getText());
					int nequil = 0;
					
					btnRun.setEnabled(false);
					//SpinReader reader = new SpinReader();
					//reader.openReader("init_spin.dat");
					CellPottsModel model = new CellPottsModel(
							nx, ny, q, temp, lambda, alpha, beta, motility, rotateDiff,
							seed, numOfSweeps, nequil, new DataWriter [] {}, true);
					model.initSpin();
					//model.initSpin(reader.readSpins());
					model.initPolarity();
					
					view.setModel(model);
					view.initImage();
					
					model.run();
					
					view.stopDrawingImage();
					
					btnRun.setEnabled(true);
				}
			};
			runthread.start();
		}
	}
}
