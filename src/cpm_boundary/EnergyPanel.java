package cpm_boundary;

import java.awt.*;

import javax.swing.*;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.*;

@SuppressWarnings("serial")
public class EnergyPanel extends JPanel implements DataListener {
	
	private CellPottsModel model = null;
	private XYSeries series = null;
	private XYSeriesCollection data = null;
	private JFreeChart chart = null;
	private ChartPanel chartPanel = null;
	private JLabel lblAcceptRate;
	private JLabel lblPerimeterToAreaRatio;
	private JPanel dataPanel;
	
	public void setModel(CellPottsModel model){
		if (this.model != null){
			this.model.removeDataListener(this);
		}
		this.model = model;
		this.model.addDataListener(this);
		init();
	}
	
	public void init(){
		if (chartPanel != null){
			this.remove(chartPanel);
			this.remove(lblAcceptRate);
		}
		series = new XYSeries("Energy data");
		data = new XYSeriesCollection();
		data.addSeries(series);
		chart = ChartFactory.createXYLineChart(
				"Energy data",
				"Time (MCS)",
				"Energy (kT)",
				data,
				PlotOrientation.VERTICAL,
				true,
				true,
				false);
		chartPanel = new ChartPanel(chart);
		lblAcceptRate = new JLabel("Accept Rate: ");
		lblPerimeterToAreaRatio = new JLabel("Perimeter to Area Ratio: ");
		dataPanel = new JPanel();
		dataPanel.add(lblAcceptRate);
		dataPanel.add(lblPerimeterToAreaRatio);
		this.setLayout(new BorderLayout());
		this.add(chartPanel, BorderLayout.CENTER);
		this.add(dataPanel, BorderLayout.SOUTH);
		this.validate();
	}

	@Override
	public void update(CellPottsModel model, int time) {
		series.add(time, model.calculateR2()[0]);
		if (time % 10 == 0){
			lblAcceptRate.setText(String.format("Accept Rate: %.5f",model.getAcceptRate()));
			lblPerimeterToAreaRatio.setText(
					String.format("Perimeter to Area Ratio: %.5f", model.getPerimeterToAreaRatio()[0]));
			this.validate();
		}
	}
}
