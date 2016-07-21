package cpm_boundary;

import java.awt.*;

import javax.swing.*;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.*;

@SuppressWarnings("serial")
public class EnergyPanel extends JPanel implements DataListener {
	
	private CellPottsModel model = null;
	private XYSeries series1 = null;
	private XYSeries series2 = null;
	private XYSeriesCollection data = null;
	private JFreeChart chart = null;
	private ChartPanel chartPanel = null;
	
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
		}
		series1 = new XYSeries("Roughness data stdv");
		series2 = new XYSeries("Roughness data diff");
		data = new XYSeriesCollection();
		data.addSeries(series1);
		data.addSeries(series2);
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
		this.setLayout(new BorderLayout());
		this.add(chartPanel, BorderLayout.CENTER);
		this.validate();
	}

	@Override
	public void update(CellPottsModel model, int time) {
		double [] roughness = model.calculateRoughness();
		series1.add(time, roughness[0]);
		series2.add(time, roughness[1]);
	}
}
